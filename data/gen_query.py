import argparse
from ast import walk
from email.errors import HeaderParseError
from tracemalloc import start 
import networkx as nx
import ipdb
import random
import os
import math
import ipdb

parser = argparse.ArgumentParser()
parser.add_argument('--edge_list')
parser.add_argument('--label_map')
parser.add_argument('--nodes')
parser.add_argument('--num')
parser.add_argument('--output_dir')
parser.add_argument('--type', help='opttion: [path, tree ,cycle, other]')

args = parser.parse_args()

'''
constraint generation part
'''

def search(G, path, result, result_indx):
    if path[-1] == path[0] and len(path)>3:
        tmp = path[:-1]
        tmp.sort()
        if tmp not in result_indx:
            result_indx.append(tmp)
            result.append(path[:-1])
        return
    for n in G.neighbors(path[-1]):
        if n not in path[1:]:
            path_ = path[:]
            path_.append(n)
            search(G, path_, result, result_indx)
    
def find_all_cycle(G, start):
    result = []
    search(G, [start], result, [])
    return result

def get_CC(G, start, rand=False):
    cc = find_all_cycle(G, start)
    if len(cc) == 0:
        return None
    r = []
    if rand == True:
        return random.choice(cc)
    for c in cc:
        label_map = nx.get_node_attributes(G,'label')
        labels = set([label_map[node] for node in c])
        r.append((c, float(len(labels))/float(len(c))))
    r = sorted(r, key=lambda x:x[1], reverse=True)
    return r[0][0]

def get_PC(G, start, rand=False):
    target = []
    label_map = nx.get_node_attributes(G,'label')
    for node in G.nodes():
        if label_map[node] == label_map[start] and node != start:
            target.append(node)
    if len(target) == 0:
        return None
    paths = []
    for n in target:
        for p in nx.all_shortest_paths(G, source=start, target=n):
            if len(p)>2:
                paths.append(p)
        # paths += [p for p in nx.all_shortest_paths(G, source=start, target=n)]
    if len(paths) == 0:
        return None
    r = []
    if rand==True:
        return random.choice(paths)
    for path in paths:
        labels = set([label_map[node] for node in path])
        r.append((path, float(len(labels))/float(len(path))))
    r = sorted(r, key=lambda x:x[1], reverse=True)
    return r[0][0]

def generat_NLCC(G, label_map):
    label_dist = {}
    for n in label_map:
        label = label_map[n]
        if label not in label_dist:
            label_dist[label] = 0
        label_dist[label] += 1
    degree_list = []
    for node in G.nodes():
        degree_list.append((node, label_dist[G.nodes[node]['label']]))
    degree_list = sorted(degree_list, key=lambda x:x[1], reverse=False)
    selected_start_vertex = [ele[0] for ele in degree_list[:len(degree_list)//3]]
    NLCC = []
    for node in selected_start_vertex:
        cc = get_CC(G, node)
        if cc is not None:
            NLCC.append(("cc", cc))
        else:
            pc = get_PC(G, node)
            if pc is not None:
                NLCC.append(("pc", pc))
    return NLCC
    
'''
query generation part
'''

def random_walk(adj, start_id, label_map):
    walk = [start_id]
    id_map = {start_id:0}
    G = nx.Graph()
    G.add_node(0, label=label_map[start_id])
    G_edges = []
    for i in range(int(args.nodes)-1):
        can = adj[walk[i]]-set(walk)
        if len(can) == 0:
            return None, None
        else:
            next_node = random.choice(list(can))
            id_map[next_node] = len(id_map)
            G.add_node(id_map[next_node], label=label_map[next_node])
            conn_edge = adj[next_node] & set(walk)
            for e in conn_edge:
                G_edges.append([next_node, e])
            walk.append(next_node)
    for e in G_edges:
        G.add_edge(id_map[e[0]], id_map[e[1]])
    return walk, G

def sample(adj, label_map):
    # start_id = random.randint(0, len(label_map)-1)
    start_id = random.choice(list(label_map.keys()))
    path, G  = random_walk(adj, start_id, label_map)
    if path == None:
        return None
    if args.type == 'path':
        g = nx.Graph()
        for i in range(len(path)):
            g.add_node(i, label=label_map[path[i]])
        for i in range(len(path)-1):
            g.add_edge(i, i+1)
        return g
    elif args.type == 'tree':
        
        tree = nx.tree.minimum_spanning_edges(G, algorithm='prim', data=True)
        tree = nx.Graph(tree)
        for n in tree.nodes():
            tree.nodes[n]['label'] = label_map[path[n]]
        return tree
        # return g
    elif args.type == 'cycle':
        if G.has_edge(0, len(path)-1):
            g = nx.Graph()
            for i in range(len(path)):
                g.add_node(i, label=label_map[path[i]])
            for i in range(len(path)-1):
                g.add_edge(i, i+1)
            g.add_edge(len(path)-1, 0)
            return g
        else:
            return None
    else:
        tree = nx.tree.minimum_spanning_edges(G, algorithm='prim', data=True)
        tree = nx.Graph(tree)
        non_tree_edge = set(G.edges())-set(tree.edges())
        non_tree_edge = list(non_tree_edge)
        if len(non_tree_edge) == 0:
            return None
        selected_non_tree_edge = random.sample(non_tree_edge, min(len(non_tree_edge), len(path)))
        for e in selected_non_tree_edge:
            tree.add_edge(e[0],e[1])
        for n in tree.nodes():
            tree.nodes[n]['label'] = label_map[path[n]]
        return tree

def generate(edge_list, label_map):
    adj = {}
    label = {}
    with open(edge_list) as f:
        for line in f.readlines():
            ele = [int(e) for e in line.strip('\n').split(" ")]
            if ele[0] not in adj:
                adj[ele[0]] = set([])
            if ele[1] not in adj:
                adj[ele[1]] = set([])
            adj[ele[0]].add(ele[1])
            adj[ele[1]].add(ele[0])
    if label_map == None:
        for n in adj:
            label[n] = math.ceil((math.log(len(adj[n])+1, 2)))
    else:
        with open(label_map) as f:
            for line in f.readlines():
                ele = [int(e) for e in line.strip('\n').split(" ")] 
                label[ele[0]] = ele[1]

    count = 0
    print("done loading")
    while count<int(args.num):
        g = sample(adj, label)
        print("sampling")
        if g != None:
            # write the query
            if not os.path.exists(args.output_dir):
                os.system('mkdir -p {0}'.format(args.output_dir))
            filename = os.path.join(args.output_dir, "{0}_{1}_{2}.graph".format(args.type, args.nodes, count))
            NLCC = generat_NLCC(g, label)
            if(len(NLCC)) == 0:
                continue
            vertices = [n for n in g.nodes()]
            vertices.sort()
            edges = [(e[0], e[1]) for e in g.edges()]
            edges.sort()
            with open(filename, 'w') as f:
                f.write("t # 0\n")
                for n in vertices:
                    f.write("v {0} {1}\n".format(n, g.nodes[n]["label"]))
                for e in edges:
                    f.write("e {0} {1}\n".format(e[0], e[1]))
                # for nlcc in NLCC:
                #     if nlcc[0] == "cc":
                #         cc = " ".join([str(i) for i in nlcc[1]+[nlcc[1][0]]])
                #         f.write("c {0} : {0} : {0} : 1 : 0 : 1\n".format(cc))
                #     elif nlcc[0] == "pc":
                #         f.write("c {0} : {0} : {0} : 0 : 0 : 1\n".format(" ".join([str(i) for i in nlcc[1]])))
                # f.write("c diameter : {0}\n".format(nx.diameter(g)))
                f.write("t # -1\n")
            count += 1


def generate_sync(graph_file):
    adj = {}
    label = {}
    with open(graph_file) as f:
        for line in f.readlines():
            ele = line.strip('\n').split(" ")
            if ele[0] == 'v':
                label[int(ele[1])] = int(ele[2])
            elif ele[0] == 'e':
                if int(ele[1]) not in adj:
                    adj[int(ele[1])] = set([])
                if int(ele[2]) not in adj:
                    adj[int(ele[2])] = set([])
                adj[int(ele[1])].add(int(ele[2]))
                adj[int(ele[2])].add(int(ele[1]))

    count = 0
    print("done loading")

    start_idx = max([int(f.split("_")[2].rstrip(".graph")) for f in os.listdir(args.output_dir)])+1 if os.path.exists(args.output_dir) and len([f.split("_")[2] for f in os.listdir(args.output_dir)])>0 else 0
    while count<int(args.num):
        g = sample(adj, label)
        print("sampling")
        if g != None:
            # write the query
            if not os.path.exists(args.output_dir):
                os.system('mkdir -p {0}'.format(args.output_dir))
            filename = os.path.join(args.output_dir, "{0}_{1}_{2}.graph".format(args.type, args.nodes, start_idx+count))
            vertices = [n for n in g.nodes()]
            vertices.sort()
            edges = [(e[0], e[1]) for e in g.edges()]
            edges.sort()
            with open(filename, 'w') as f:
                f.write("t # 0\n")
                for n in vertices:
                    f.write("v {0} {1}\n".format(n, g.nodes[n]["label"]))
                for e in edges:
                    f.write("e {0} {1}\n".format(e[0], e[1]))
                f.write("t # -1\n")
            count += 1

generate(args.edge_list, args.label_map)

# this is the generator for synthetic graphs and will automatically assign the labels following the way in PruneJuice
# generate_sync(args.edge_list)
