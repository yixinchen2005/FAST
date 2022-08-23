#---------------------------pytorch----------------------------------
import torch
import torch.multiprocessing as mp
from torch.utils.tensorboard import SummaryWriter
from torch.nn import CrossEntropyLoss, LogSoftmax, NLLLoss
#--------------------------------------------------------------------
#---------------------------pytorch geometric------------------------
from torch_geometric.data import Data, Batch
#--------------------------------------------------------------------
#------------------------others-------------------------
import argparse
import os
import random
from datetime import datetime
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.metrics import precision_recall_curve, average_precision_score
#-------------------------------------------------------
from lib.dataloader import Graphs
from lib.utils import load_dataset
from lib.models import MLP


parser = argparse.ArgumentParser()
parser.set_defaults(
                    min_size = 3,
                    max_size = 29,
                    input_dim = 84,
                    hidden_dim = 64,
                    output_dim = 64,
                    num_layers = 4,
                    lr = 1e-4,
                    dropout = 0.0,
                    batch_size = 10,
                    num_batches = 500000,
                    num_workers = 4,
                    dataset = '',
                    eval_interval= 1000,
                    layer_type = 'GIN',
                    model_path = 'saved_model/balanced_',
                    val_size= 1000,
                    test = False,
                    skip = 'all'
                    )
args = parser.parse_args()

device = torch.device('cuda') if torch.cuda.is_available() \
                              else torch.device("cpu")
args.model_path += 'model_' + args.layer_type + '_' + str(args.input_dim) +\
                   '_' + str(args.hidden_dim) + '_' + str(args.output_dim) +\
                   '_' + args.dataset + '.pt'


def batchmaker(target, pos, neg):
    target = Batch.from_data_list(target)
    pos = Batch.from_data_list(pos)
    neg = Batch.from_data_list(neg)
    return target, pos, neg
    
def validation(model, testset, logger, num_batches, epoch):
    model.eval()
    labels = torch.tensor(([1]*args.batch_size + [0] * args.batch_size) *\
                                 len(testset)).to(device)
    #print(all_labels)
    #print(len(testset))
    all_preds = []
    all_raw_preds = []
    all_labels = []
    
    for target, pos, neg in testset:
        with torch.no_grad():
            target_emb, pos_emb, neg_emb = \
                            (model.emb_model(target.to(device)),
                             model.emb_model(pos.to(device)), 
                             model.emb_model(neg.to(device)))
                                          
            emb_t = torch.cat((target_emb, target_emb),
                               dim=0)
            emb_q = torch.cat((pos_emb, neg_emb), dim=0)
            pred = model(emb_t, emb_q)

            raw_pred = pred[:,1]
            pred = pred.argmax(dim=-1)
            
        #-----------------------------------------------------------------------
        
        all_raw_preds.append(raw_pred)
        all_preds.append(pred)
        
        #-----------------------------------------------------------------------
    pred = torch.cat(all_preds, dim = -1).to(device)
    #labels = torch.cat(all_labels, dim = -1).to(device)
    raw_pred = torch.cat(all_raw_preds, dim=-1).to(device)
    
    acc = torch.mean((pred == labels).type(torch.float))
    prec = (torch.sum(pred * labels).item() / torch.sum(pred).item() if
        torch.sum(pred) > 0 else float("NaN"))
    recall = (torch.sum(pred * labels).item() /
        torch.sum(labels).item() if torch.sum(labels) > 0 else
        float("NaN"))
    labels = labels.detach().cpu().numpy()
    raw_pred = raw_pred.detach().cpu().numpy()
    pred = pred.detach().cpu().numpy()
    auroc = roc_auc_score(labels, raw_pred)
    avg_prec = average_precision_score(labels, raw_pred)
    tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
  
    print('-------------------------------------------------------------------')    
    print("\n{}".format(str(datetime.now())))
    print("Validation. Epoch {}. Acc: {:.4f}. "
        "P: {:.4f}. R: {:.4f}. AUROC: {:.4f}. AP: {:.4f}.\n     "
        "TN: {}. FP: {}. FN: {}. TP: {}".format(epoch,
            acc, prec, recall, auroc, avg_prec,
            tn, fp, fn, tp))
                 
    if not args.test:
        logger.add_scalar("Accuracy/test", acc, num_batches)
        logger.add_scalar("Precision/test", prec, num_batches)
        logger.add_scalar("Recall/test", recall, num_batches)
        logger.add_scalar("AUROC/test", auroc, num_batches)
        logger.add_scalar("AvgPrec/test", avg_prec, num_batches)
        logger.add_scalar("TP/test", tp, num_batches)
        logger.add_scalar("TN/test", tn, num_batches)
        logger.add_scalar("FP/test", fp, num_batches)
        logger.add_scalar("FN/test", fn, num_batches)
        print("Saving {}".format(args.model_path))
        torch.save(model.state_dict(), args.model_path)
    
    
def train(model, train_set, in_queue, out_queue):
    done = False
    trainset = Graphs(train_set, args.min_size, args.max_size,
                            args.eval_interval)
    optimizer = torch.optim.Adam(model.parameters(), lr = args.lr)
    
    labels = torch.tensor([1]*args.batch_size + [0] * args.batch_size).to(device)
    while not done:
        msg, _ = in_queue.get()
        if msg == "done":
            done = True
            break
        target, pos, neg = [], [], []
        for i in range(args.batch_size):
            t, p, n = trainset.__getitem__(1)
            target.append(t)
            pos.append(p)
            neg.append(n)
        target, pos, neg = batchmaker(target, pos, neg)
        model.train()
        model.zero_grad()
            
        target_emb, pos_emb, neg_emb = \
                        (model.emb_model(target.to(device)),
                         model.emb_model(pos.to(device)), 
                         model.emb_model(neg.to(device)))
                                   
        emb_t = torch.cat((target_emb, target_emb), dim = 0)
        emb_q = torch.cat((pos_emb, neg_emb), dim = 0)
            
        pred = model(emb_t, emb_q)
        criterion = NLLLoss()
        loss = criterion(pred, labels)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        pred = pred.argmax(dim=-1)
        acc = torch.mean((pred == labels).type(torch.float))
        out_queue.put(("step",
                     (loss.item(), acc.item())))

def train_multiprocess():
    if not os.path.exists(os.path.dirname(args.model_path)):
        os.makedirs(os.path.dirname(args.model_path))
        
    in_queue, out_queue = mp.Queue(), mp.Queue()
    
    trainset, test = load_dataset()
    
    print("Using dataset {}".format(args.dataset))
    
    model = MLP(args.input_dim, args.hidden_dim, args)
    model.share_memory()
    model.to(device)

    testset = Graphs(test, args.min_size, args.max_size,
                            args.val_size)
    tests = []
    for i in range(args.val_size//args.batch_size):
        print(i)
        target, pos, neg = [], [], []
        for j in range(args.batch_size):
            t, p, n = testset.__getitem__(1)
            target.append(t)
            pos.append(p)
            neg.append(n)
        target, pos, neg = batchmaker(target, pos, neg)
        
        target, pos, neg = (target.to(torch.device('cpu')),
                            pos.to(torch.device('cpu')),
                            neg.to(torch.device('cpu')))
        tests.append((target, pos, neg))
    print('Validation generated.')

    
    record_keys = ['layer_type','num_layers','input_dim',
                   'hidden_dim','output_dim','dataset', 'method_type']
    comment = ".".join(["{}={}".format(k, v) \
              for k, v in sorted(vars(args).items()) if k in record_keys])
    logger = SummaryWriter(comment = comment)

    workers = []
    for i in range(args.num_workers):
        worker = mp.Process(target = train, args = (model, trainset, \
                            in_queue, out_queue))
        worker.start()
        workers.append(worker)
        
    if args.test:
        validation(model, tests, logger, 0, 0, True)
    else:
        num_batches = 0
        for epoch in range(args.num_batches // args.eval_interval):
            for i in range(args.eval_interval):
                in_queue.put(('step', None))
            for i in range(args.eval_interval):
                msg, params = out_queue.get()
                train_loss, train_acc = params
                print("Batch {}. Loss: {:.4f}. Training acc: {:.4f}.".format(
                    num_batches, train_loss, train_acc),
                    end="              \r")
                logger.add_scalar("Loss/train", train_loss, num_batches)
                logger.add_scalar("Accuracy/train", train_acc, num_batches)
                num_batches += 1
            validation(model, tests, logger, num_batches, epoch)

    for i in range(args.num_workers):
        in_queue.put(("done", None))
    for worker in workers:
        worker.join()

def main():
    mp.set_start_method("spawn", force=True)
    train_multiprocess()
    
if __name__ == '__main__':
    main()
