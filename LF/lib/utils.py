import torch
import torch.nn.functional as F
import random
import os
import pickle

from torch_geometric.data import Data
from torch_geometric.utils import remove_self_loops
from torch_sparse import coalesce

def load_dataset():
    dataset = []
    for i in range(10):
        if i != 1:
            _, data = load_splited('cycle', i)
            dataset.append(data)
            print('----------------------')
            print(data.num_nodes)
        _, data = load_splited('path', i)
        print(data.num_nodes)
        dataset.append(data)
        _, data = load_splited('tree', i)
        dataset.append(data)
        print(data.num_nodes)
        _, data = load_splited('other', i)
        dataset.append(data)
        print(data.num_nodes)
        print('----------------------')
        print(data.num_nodes)
    #print(len(dataset))
    #random.shuffle(dataset)
    idx = int(len(dataset) * 0.9)
    return dataset[0:idx], dataset[idx:]

def load_splited(gtype, graph_num):
    path = './graphs/HR_' + gtype + '_' + str(graph_num) + '/pruned_graph'
    f = open(path, "r")
    nodes = []
    edges = []
    node_labels = []
    edge_labels = []
    dic = {}
    num = 0
    din, dout = [],[]
    while f:
        line = f.readline()
        if len(line) <= 0:
            break
        splited = line.split()
        if splited[0] == 'v':
            nodes.append(int(splited[1]))
            node_labels.append(int(splited[2]))
            dic[int(splited[1])] = num
            num += 1
        elif splited[0] == 'e':
            #edges.append((int(splited[1]),int(splited[2])))
            #edges.append((int(splited[2]),int(splited[1])))
            din.append(int(splited[1]))
            dout.append(int(splited[2]))
    f.close()
    #print(dout[-3:])
    for i in range(len(din)):
        edges.append((dic[din[i]], dic[dout[i]]))
        edges.append((dic[dout[i]], dic[din[i]]))
    edges = torch.tensor(edges, dtype=torch.long).squeeze().t()
    labels = torch.tensor(node_labels, dtype=torch.long)
    if labels.dim() == 1:
        labels = labels.unsqueeze(-1)
    labels = labels - labels.min(dim = 0)[0]
    labels = labels.unbind(dim = -1)
    '''
        x size
    '''
    labels = [F.one_hot(x, num_classes=84) for x in labels]
    x = torch.cat(labels, dim = -1).to(torch.float)
    num_nodes = edges.max().item() + 1 if x is None else x.size(0)
    edge_index, _ = remove_self_loops(edges)
    edge_index, _ = coalesce(edge_index, None,  num_nodes, num_nodes)
    y = torch.tensor(nodes, dtype = torch.int)
    data = Data(x = x,
                edge_index = edge_index,
                y = y)
    return dic, data


def load_query(gtype, name):
    path = './examples/'
    if gtype == 'cycle':
        path += 'cycle_none/cycle_8_' + name + '.graph'
    elif gtype == 'path':
        path += 'path_none/path_8_' + name + '.graph'
    elif gtype == 'tree':
        path += 'tree_none/tree_8_' + name + '.graph'
    else:
        path += 'other_none/other_8_' + name + '.graph'
    f = open(path, "r")
    nodes = []
    edges = []
    node_labels = []
    edge_labels = []
    line = f.readline()
    while f:
        line = f.readline()
        if len(line) <= 0:
            break
        splited = line.split()
        if splited[0] == 'v':
            nodes.append(int(splited[1]))
            node_labels.append(int(splited[2]))
        elif splited[0] == 'e':
            edges.append([int(splited[1]),int(splited[2])])
            edges.append([int(splited[2]),int(splited[1])])
    f.close()
    edges = torch.tensor(edges, dtype=torch.long).squeeze().t()
    labels = torch.tensor(node_labels, dtype=torch.long)
    if labels.dim() == 1:
        labels = labels.unsqueeze(-1)
    labels = labels - labels.min(dim = 0)[0]
    labels = labels.unbind(dim = -1)
    '''
        x size
    '''
    labels = [F.one_hot(x, num_classes=84) for x in labels]
    x = torch.cat(labels, dim = -1).to(torch.float)
    num_nodes = edges.max().item() + 1 if x is None else x.size(0)
    edge_index, _ = remove_self_loops(edges)
    edge_index, _ = coalesce(edge_index, None,  num_nodes, num_nodes)
    data = Data(x = x,
                edge_index = edge_index)
    return data