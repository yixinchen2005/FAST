import argparse
from cProfile import label
from logging import root 
import networkx as nx
import ipdb
import random
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--input_label_map')
parser.add_argument('--output')
parser.add_argument('--mode', help='options: [random, heuristic, none]')
parser.add_argument('--para', help='potion of the selected vertices')

args = parser.parse_args()

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
        paths += [p for p in nx.all_shortest_paths(G, source=start, target=n)]
    r = []
    if rand==True:
        return random.choice(paths)
    for path in paths:
        labels = set([label_map[node] for node in path])
        r.append((path, float(len(labels))/float(len(path))))
    r = sorted(r, key=lambda x:x[1], reverse=True)
    return r[0][0]

def print_nlcc(nlcc):
    l = []
    tail = ""
    if nlcc[0] == "cc":
        l = [i for i in nlcc[1]+[nlcc[1][0]]]
        tail = "1 : 0 : 1"
    elif nlcc[0] == "pc":
        l = [i for i in nlcc[1]]
        tail = "0 : 0 : 1"
    scan_map = {}
    offset = 0
    for ele in l:
        if ele not in scan_map:
            scan_map[ele] = offset
        offset += 1
    l2 = [scan_map[i] for i in l]
    return "c {0} : {1} : {2} : {3}".format(" ".join([str(i) for i in l]), " ".join([str(i) for i in l2]), " ".join(["0" for i in range(len(l))]), tail)
        

def generate_heuristic(input_file, output_file, data_label):
    # analyse the data graph
    label_dist = {}
    with open(data_label) as f:
        for line in f.readlines():
            ele = line.strip("\n").split()
            if ele[1] not in label_dist:
                label_dist[ele[1]] = 0
            label_dist[ele[1]] += 1

    # read the graph
    G = nx.Graph()
    vertices = []
    edges = []
    with open(input_file) as f:
        for line in f.readlines():
            ele = line.strip('\n').split()
            if line.startswith('v'):
                G.add_node(int(ele[1]), label = ele[2])
                vertices.append(line)
            elif line.startswith('e'):
                G.add_edge(int(ele[1]), int(ele[2]))
                edges.append(line)

    # degree_list = []
    # for node in G.nodes():
    #     degree_list.append((node, G.degree(node)))
    # degree_list = sorted(degree_list, key=lambda x:x[1], reverse=True)
    degree_list = []
    for node in G.nodes():
        degree_list.append((node, label_dist[G.nodes[node]['label']]))
    degree_list = sorted(degree_list, key=lambda x:x[1], reverse=False)
    print(degree_list)

    selected_start_vertex = [ele[0] for ele in degree_list[:len(degree_list)//int(args.para)]]
    NLCC = []
    for node in selected_start_vertex:
        cc = get_CC(G, node)
        if cc is not None:
            NLCC.append(("cc", cc))
        else:
            pc = get_PC(G, node)
            if pc is not None:
                NLCC.append(("pc", pc))

    with open(output_file, 'w') as f:
        f.write("t # 0\n")
        vertices.sort()
        edges.sort()
        for v in vertices:
            f.write(v)
        for e in edges:
            f.write(e)
        for nlcc in NLCC:
            f.write(print_nlcc(nlcc)+"\n")
            # if nlcc[0] == "cc":
            #     cc = " ".join([str(i) for i in nlcc[1]+[nlcc[1][0]]])
            #     f.write("c {0} : {0} : {0} : 1 : 0 : 1\n".format(cc))
            # elif nlcc[0] == "pc":
            #     f.write("c {0} : {0} : {0} : 0 : 0 : 1\n".format(" ".join([str(i) for i in nlcc[1]])))
        f.write("c diameter : {0}\n".format(nx.diameter(G)))
        f.write("t # -1\n")

def generate_random(input_file, output_file):
    G = nx.Graph()
    vertices = []
    edges = []
    with open(input_file) as f:
        for line in f.readlines():
            ele = line.strip('\n').split()
            if line.startswith('v'):
                G.add_node(int(ele[1]), label = ele[2])
                vertices.append(line)
            elif line.startswith('e'):
                G.add_edge(int(ele[1]), int(ele[2]))
                edges.append(line)

    start_vertices = random.sample([k for k in G.nodes()], len(G.nodes())//int(args.para))
    print(start_vertices)
    NLCC = []
    for node in start_vertices:
        ch = random.choice(["CC", "PC"])
        if ch=="CC":
            cc = get_CC(G, node, True)
            if cc is not None:
                NLCC.append(("cc", cc))
            else:
                pc = get_PC(G, node, True)
                if pc is not None:
                    NLCC.append(("pc", pc))
        else:
            pc = get_PC(G, node, True)
            if pc is not None:
                NLCC.append(("pc", pc))
            else:
                cc = get_CC(G, node, True)
                if cc is not None:
                    NLCC.append(("cc", cc))
    vertices.sort()
    edges.sort()
    with open(output_file, 'w') as f:
        f.write("t # 0\n")
        for v in vertices:
            f.write(v)
        for e in edges:
            f.write(e)
        for nlcc in NLCC:
            f.write(print_nlcc(nlcc)+"\n")
            # if nlcc[0] == "cc":
            #     cc = " ".join([str(i) for i in nlcc[1]+[nlcc[1][0]]])
            #     f.write("c {0} : {0} : {0} : 1 : 0 : 1\n".format(cc))
            # elif nlcc[0] == "pc":
            #     f.write("c {0} : {0} : {0} : 0 : 0 : 1\n".format(" ".join([str(i) for i in nlcc[1]])))
        f.write("c diameter : {0}\n".format(nx.diameter(G)))
        f.write("t # -1\n")

def generate_none(input_file, output_file):
    G = nx.Graph()
    vertices = []
    edges = []
    with open(input_file) as f:
        for line in f.readlines():
            ele = line.strip('\n').split()
            if line.startswith('v'):
                G.add_node(int(ele[1]), label = ele[2])
                vertices.append(line)
            elif line.startswith('e'):
                G.add_edge(int(ele[1]), int(ele[2]))
                edges.append(line)
    vertices = [n for n in G.nodes()]
    vertices.sort()
    edges = [(e[0], e[1]) for e in G.edges()]
    edges.sort()
    vertices.sort()
    edges.sort()
    with open(output_file, 'w') as f:
        f.write("t # 0\n")
        for v in vertices:
            f.write(v)
        for e in edges:
            f.write(e)
        f.write("c diameter : {0}\n".format(nx.diameter(G)))
        f.write("t # -1\n")


if not os.path.exists(args.output):
    os.mkdir(args.output)
for file in os.listdir(args.input):
    if not file.endswith(".graph"):
        continue
    filename = os.path.join(args.input, file)
    output_filename = os.path.join(args.output, file)
    print(output_filename)
    if args.mode == "heuristic":
        generate_heuristic(filename, output_filename, args.input_label_map)
    elif args.mode == "random":
        generate_random(filename, output_filename)
    elif args.mode == "none":
        generate_none(filename, output_filename)
    


# root_path = "./dataset/queries/twitch-game"
# data_graph = "./dataset/twitch-game/twitch_labels"

# generate_heuristic("./dataset/queries/gemsec/HR/cycle/cycle_8_0.graph", "test", data_graph)
# mode = "heuristic"
# for t in ["cycle", "other", "path", "tree"]:
#     if not os.path.exists(os.path.join(root_path, t+"_"+mode)):
#         os.mkdir(os.path.join(root_path, t+"_"+mode))
#     for q in os.listdir(os.path.join(root_path, t)):
#         if q.endswith(".graph"):
#             query_file = os.path.join(os.path.join(root_path, t), q)
#             output_file = os.path.join(os.path.join(root_path, t+"_"+mode), q)
#             print(query_file)
#             print(output_file)
#             if args.mode == "heuristic":
#                 generate_heuristic(query_file, output_file, data_graph)
#             elif args.mode == "random":
#                 generate_random(query_file, output_file)
        

# if args.mode == "heuristic":
#     generate_heuristic()
# else:
#     generate_random()
