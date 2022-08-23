import argparse
import torch
import os
import networkx as nx
import pickle

from torch_geometric.data import Batch
from torch_geometric.utils import to_networkx
#-------------------------------------------------------------------------------
from lib.models import MLP
from lib.utils import load_splited, load_query
from lib.dataloader import gen_k_hop_subgraph


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
                    dataset = 'splited',
                    eval_interval= 1000,
                    layer_type = 'GIN',
                    model_path = 'saved_model/balanced_',
                    val_size= 1000,
                    skip = 'all'
                    )
args = parser.parse_args()

device = torch.device('cuda') if torch.cuda.is_available() \
                              else torch.device("cpu")

args.model_path += 'model_' + args.layer_type + '_' + str(args.input_dim) +\
                   '_' + str(args.hidden_dim) + '_' + str(args.output_dim) +\
                   '_' + args.dataset + '.pt'

#-------------------------------------------------------------------------------

def main():
    if not os.path.exists("./candidates/"):
        os.makedirs("./candidates/")

    model = MLP(args.input_dim, args.hidden_dim, args)
    model.load_state_dict(torch.load(args.model_path,
            map_location=device))
    model.eval()
    print('Pre-trained model loaded. ')
    
    name = 'path'
    candidates = set()
    for k in range(10):
        if k != 1:
            num = k
            dic, data = load_pruned(name, num)

            query = load_query(name, str(num))
            candidate = set()
            embs = []
            for i in range(query.num_nodes):   
                print('comparing query node : ', i)
                center, q = gen_k_hop_subgraph(i, query)
                q = Batch.from_data_list([q]).to(device)
                emb_q = model.emb_model(q)
                all_preds = []
                for j in range(data.num_nodes):
                    #print(j)
                    if len(embs) < data.num_nodes:
                        center, target = gen_k_hop_subgraph(j, data)
                        target = Batch.from_data_list([target]).to(device)
                        emb_t = model.emb_model(target)
                        #embs.append(emb_t)
                    else:
                        emb_t = embs[j].to(device)
                        #print(emb_t)
                    pred = model(emb_t, emb_q)
                        
                    raw_pred = model.predict(pred)
                    raw_pred = raw_pred[:,1]
                    all_preds.append(raw_pred.unsqueeze(1).item())

                    pred = pred.argmax(dim=-1)
                    if j % 10000 == 0:
                        print(j)
                    if pred.item() == 1:
                        #print(j, ' : ', data.y[j].item())
                        tmp = data.y[j].item()
                        candidate.update([tmp])
                rpath = './candidates/' + name + '_' + str(num) + '_node_' + str(i)
                if not os.path.exists(rpath):
                    f = open(rpath, 'x')
                    f.close()
            
                f = open(rpath, 'wb')
                pickle.dump(candidate, f)
                candidates.update(candidate)
        rpath = './candidates/' + name + '_' + str(num)
        if not os.path.exists(rpath):
            f = open(rpath, 'x')
            f.close()
            
            f = open(rpath, 'wb')
        pickle.dump(candidates, f)

    
if __name__ == '__main__':
    main()
