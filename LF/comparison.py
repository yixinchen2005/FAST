import pickle
from lib.utils import load_splited
from torch_geometric.utils import subgraph

def load_gt(gtype, num):
    path = './candidate/' + gtype + '_' + str(num) + '_gt'
    with open(path, 'rb') as f:
        gt = pickle.load(f)
    return gt

def load_node_set(gtype, num, nnum):
    path = './candidates/' + gtype + '_' + str(num) + '_node_' + str(nnum)
    with open(path, 'rb') as f:
        nodeset = pickle.load(f)
    return nodeset
    
def load_candidates(gtype, num):
    path = './candidates/' + gtype + '_' + str(num)
    with open(path, 'rb') as f:
        nodeset = pickle.load(f)
    return nodeset

name = 'path'
avg_acc = 0.0
avg_nodes = 0.0
avg_edges = 0.0
avg_after_nodes = 0.0
avg_after_edges = 0.0
n_g = 0.0
for num in range(10):
    print('***********************************************************************')
    print('on query ', name, ' ', str(num))
    gt = load_gt(name, num)
    n_g += 1.0
    candidates = load_candidates(name, num)
    #print(candidates)
    candidate = []
    for i in range(8):
        tmp = load_node_set(name, num, i)
        #print(len(tmp))
        candidate.append(tmp)
    union = candidate[7]
    for i in range(7):
        union = candidate[i] | union
    dic, data = load_splited(name, num)
    all_neg = data.num_nodes - len(gt)
    all_pos = len(gt)
    # confusion matrix union
    # tp pos sample and pos pred, both in union and gt
    # tn neg sample and neg pred, neither in union or gt
    # fp neg sample and pos pred, in union not in gt
    # fn pos sample and neg pred, in gt but not in union
    tp = len(union & gt)
    tn = data.num_nodes - len(gt|union)
    fp = len(union - gt)
    fn = len(gt - union)
    acc = (tp + tn)/data.num_nodes
    #---------------------------
    print('num candidates (union) : ', len(union))
    print('num gt : ', len(gt))
    print('---------------------------------------------------')
    print('false neg (union) : ', fn)
    print('true neg (union) : ', tn)
    print('false pos (union) : ', fp)
    print('true pos (union) : ', tp)
    #print('jaccard sim (union) : ', tp / len(gt | union))
    print('Acc. (union) : ',  acc)
    avg_acc += acc

    print(data)
    nodes = []
    union = list(union)
    for i in range(len(union)):
        #print(dic)
        nodes.append(dic[union[i]])

    print('num node after pruning (union): ', len(nodes))
    #print(nodes[0:10])   
    edge_index, _ = subgraph(nodes, data.edge_index)
    print('num edge after pruning (union): ', len(edge_index[0]))
    
    avg_nodes += data.num_nodes
    avg_edges += len(data.edge_index[0])
    avg_after_nodes += len(nodes)
    avg_after_edges += len(edge_index[0])
    print('***********************************************************************')
print('On queries : ', name)
print('avg. acc : ', avg_acc / n_g)
print('avg. num nodes : ', avg_nodes / n_g)
print('avg. num edges : ', avg_edges / n_g)
print('avg. num nodes after pruning : ', avg_after_nodes / n_g)
print('avg. num edges after pruning : ', avg_after_edges / n_g)
print(n_g) 
