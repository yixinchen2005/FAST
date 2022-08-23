#---------------------------pytorch----------------------------------
import torch
from torch.utils.data import Dataset, DataLoader
#--------------------------------------------------------------------
#---------------------------pytorch geometric------------------------
from torch_geometric.data import Batch, Data
from torch_geometric.utils import to_networkx, k_hop_subgraph, subgraph, \
                                  dropout_adj, contains_isolated_nodes, \
                                  sort_edge_index
#--------------------------------------------------------------------

#---------------------------home made lib----------------------------
from lib.pruned_utils import load_pruned
#--------------------------------------------------------------------
#------------------------others--------------------------------------
import random
import numpy as np
import scipy.stats as stats
import networkx as nx
#--------------------------------------------------------------------
            
def discrete_sampling(dataset):
    probs = np.array([data.num_nodes for data in dataset], dtype=np.float)
    probs /= np.sum(probs)
    dist = stats.rv_discrete(values = (np.arange(len(dataset)), probs))
    return dist

def gen_k_hop_subgraph(idx, data, relabel = True):
    #print('idx : ', idx)
    subset, edge_index,_, _ = k_hop_subgraph(idx, 3, data.edge_index, 
                                             relabel_nodes = False)
    #print('idx label : ', torch.argmax(data.x[idx]).item() + 1)
    if relabel: 
        x = None
        y = None
        for i in range(len(subset)):
            if subset[i] == idx :
                center = i
            if x == None and data.x != None:
                x = data.x[subset[i]].unsqueeze(0)
            elif x != None:
                x = torch.cat((x, data.x[subset[i]].unsqueeze(0)), 0)
            if y == None and data.y != None:
                y = data.y[subset[i]].unsqueeze(0)
            elif y != None:
                y = torch.cat((y, data.y[subset[i]].unsqueeze(0)), 0)
        edge_index, _ = subgraph(subset, data.edge_index, 
                                 relabel_nodes = True)
    else:
        x = data.x
        y = data.y
        center = idx
    #print('idx label : ', torch.argmax(data.x[idx]).item() + 1)
    if y == None:
        return center, Data(x = x, edge_index = edge_index)
    else :
        return center, Data(x = x, edge_index = edge_index, y = y)

def gen_subgraph(subset, data):
    x = None
    y = None
    for i in range(len(subset)):
        if x == None and data.x != None:
            x = data.x[subset[i]].unsqueeze(0)
        elif x != None:
            x = torch.cat((x, data.x[subset[i]].unsqueeze(0)), 0)
        if y == None and data.y != None:
            y = data.y[subset[i]].unsqueeze(0)
        elif y != None:
            y = torch.cat((y, data.y[subset[i]].unsqueeze(0)), 0)
    edge_index, _ = subgraph(subset, data.edge_index, 
                             relabel_nodes = True)
    return Data(x = x, edge_index = edge_index, y = y)
                  
def sample_nodes(dataset, prob, size, fix_id = None):
    while True:
        if fix_id == None:
            idx = prob.rvs()
        else:
            idx = fix_id
        graph = to_networkx(dataset[idx])
        start_node = random.choice(list(graph.nodes))
        neigh = [start_node]
        graph = nx.ego_graph(graph, start_node, radius = 3)
        frontier = list(set(graph.neighbors(start_node)) - set(neigh))
        visited = set([start_node])
        while len(neigh) < size and frontier:
            new_node = random.choice(list(frontier))
            #new_node = max(sorted(frontier))
            assert new_node not in neigh
            neigh.append(new_node)
            visited.add(new_node)
            frontier += list(graph.neighbors(new_node))
            frontier = [x for x in frontier if x not in visited]
        if len(neigh) == size:
            #    pos = gen_subgraph(neigh, dataset[idx])
            return idx, start_node, neigh

class Graphs(Dataset):
    def __init__(self, dataset, min_size, max_size, batch_size):
        self.max_size = max_size
        self.min_size = min_size
        self.dataset = dataset
        self.batch_size = batch_size
        #print(self.dataset)
        self.sampling_prob = discrete_sampling(self.dataset)
        
    
    def __getitem__(self, idx):
        size = random.randint(self.min_size, self.max_size)
        idx, center, pos_nodes = sample_nodes(self.dataset, 
                                              self.sampling_prob,
                                              size)
        size = random.randint(self.min_size, 
                              self.max_size)
        center, target = gen_k_hop_subgraph(center, self.dataset[idx])
        pos = gen_subgraph(pos_nodes, self.dataset[idx])
        while True:
            idneg, neg_center, neg_nodes = sample_nodes(self.dataset, 
                                                    self.sampling_prob,
                                                    size)
            if torch.argmax(self.dataset[idneg].x[neg_center]).item() != \
                torch.argmax(target.x[center]).item():
                #print('pass')
                neg = gen_subgraph(neg_nodes, self.dataset[idneg])
                
                return target, pos, neg
            
    def __len__(self):
        #return len(self.dataset)
        return self.batch_size
                                              
