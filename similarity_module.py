#! /usr/bin/env python3
from igraph import Graph
import math
from files import stat_assos


def read_interactome(interactomefile, weights, directed):
    """
    read a file into an igraph object using igraph read_Ncol method
    :param interactomefile: a file path of a "ncol" file contains interactome
    :param weights: same as "weights" in igraph object construction method:
    True, False, "auto", "if_present"
    :param directed: True or False
    :return: an igraph object
    """
    with open(interactomefile, mode='r') as f:
        g = Graph.Read_Ncol(f, names=True, weights=weights, directed=directed)
        return g


def clear_interactome(interactomefile, cleared_file, sep='\t'):
    """
    remove self loops and parallel edges in the interactomefile,
    save the cleared data in a new file
    :param interactomefile: a file path of a "ncol" file contains interactome
    :param cleared_file: file path of the cleared data to write
    :param sep: delimiter
    :return: None
    """
    node_pairs = {}
    with open(interactomefile, mode='r') as f:
        for line in f:
            words = line.split(sep)
            node1 = words[0].strip()
            node2 = words[1].strip()
            if node1 != node2:
                if node1 in node_pairs.keys() and node2 in node_pairs.keys():
                    if node2 not in node_pairs[node1] and node1 not in node_pairs[node2]:
                        node_pairs[node1].add(node2)
                elif node1 not in node_pairs.keys() and node2 not in node_pairs.keys():
                    node_pairs[node1] = set()
                    node_pairs[node1].add(node2)
                elif node1 in node_pairs.keys():
                    node_pairs[node1].add(node2)
                else:
                    node_pairs[node2].add(node1)
            else:
                print(line, end="")
    with open(cleared_file, mode='w') as wf:
        for n in node_pairs.keys():
            for m in node_pairs[n]:
                wf.write(n+sep+m+"\n")


def density(dgassos, graph):
    """
    given a dict which contains a list of node groups, and an igraph object,
    return a dict which gives density of these node groups based on this igraph object.
    :param dgassos: a dict object which keys are module names and values are module nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are their density scores
    """
    ddensityscore = {}
    for name in dgassos.keys():
        nodes = graph.vs.select(name_in=dgassos[name])
        ddensityscore[name] = graph.subgraph(nodes).density()
    return ddensityscore


def similarity_cal_module_1(dgassos, graph, gncutoff=0):
    """
    module cal method 1
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gncutoff: gene number cut off, only diseases whose number of associated
    genes in graph is no less than gncutoff will be calculated
    :return: a dict, (key-value: string-dict<string-float>)
    """
    gvs = set(graph.vs['name'])

    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = gvs.intersection(dgassos[d])
        if len(dgleft) > gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    sims = {}
    for i in range(0, len(diseases)):
        sims[diseases[i]] = {}
        for j in range(i, len(diseases)):
            gsi = dgassos_new[diseases[i]]
            gsj = dgassos_new[diseases[j]]
            gsintersect = gsi.intersection(gsj)
            gsid = gsi.difference(gsj)
            gsjd = gsj.difference(gsi)
            conncount = 0
            if len(gsintersect) != 0:
                for g in gsid:
                    conncount += connected_count(graph, g, gsintersect)
                for g in gsjd:
                    conncount += connected_count(graph, g, gsintersect)
            for g in gsid:
                conncount += connected_count(graph, g, gsjd)
            sims[diseases[i]][diseases[j]] = (len(gsintersect)**2+conncount)/(len(gsi)*len(gsj))
        print(i, "done..")
    return sims


def connected_count(graph, node, nodes):
    count = 0
    for n in nodes:
        if graph.are_connected(node, n):
            count += 1
    return count


def similarity_cal_module_2(dgassos, graph, gncutoff=0):
    """
    module cal method 2
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gncutoff: gene number cut off, only diseases whose number of associated
    genes in graph is no less than gncutoff will be calculated
    :return: a dict, (key-value: string-dict<string-float>)
    """
    gvs = set(graph.vs['name'])

    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = gvs.intersection(dgassos[d])
        if len(dgleft) > gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))
    genetrandis = {}
    sims = {}
    gene_pairscount = 0
    for i in range(0, len(diseases)):
        sims[diseases[i]] = {}
        print(i, diseases[i], len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            gsi = dgassos_new[diseases[i]]
            gsj = dgassos_new[diseases[j]]
            sumtdis = 0.0
            for gi in gsi:
                for gj in gsj:
                    gene_pairscount += 1
                    if gi in genetrandis.keys() and gj in genetrandis[gi].keys():
                        sumtdis += genetrandis[gi][gj]
                    elif gj in genetrandis.keys() and gi in genetrandis[gj].keys():
                        sumtdis += genetrandis[gj][gi]
                    else:
                        if gi not in genetrandis.keys():
                            genetrandis[gi] = {}
                        tdistemp = transformed_distance(
                            graph.shortest_paths(source=gi, target=gj, weights=None)[0][0])
                        sumtdis += tdistemp
                        genetrandis[gi][gj] = tdistemp
            sims[diseases[i]][diseases[j]] = sumtdis/(len(gsi)*len(gsj))
        print(" done..")
    stat_assos(genetrandis)
    print("gene pairs calculated:", gene_pairscount)
    return sims


def transformed_distance(shortestpathdis=0, a=1, b=1):
    return a*math.exp(-b*shortestpathdis)


def diameter(dgassos, graph):
    """
    given a dict which contains a list of node groups, and an igraph object,
    return a dict which gives diameters of each node group based
    on this igraph object.
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are diameters of modules
    """
    ddiameter = {}
    for d in dgassos.keys():
        nodes = graph.vs.select(name_in=dgassos[d])
        sps = graph.shortest_paths(source=nodes, target=nodes, weights=None)
        maxsps = []
        for line in sps:
            maxsps.append(max(line))
        ddiameter[d] = max(maxsps)
    return ddiameter


def avg_degree(dgassos, graph):
    """
    given a dict which contains a list of node groups, and an igraph object,
    return a dict which gives average degrees of these nodes in node groups based
    on this igraph object.
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are average degrees of nodes
    in modules
    """
    davgdegree = {}
    for d in dgassos.keys():
        nodes = graph.vs.select(name_in=dgassos[d])
        indegree = sum(graph.subgraph(nodes).degree())
        davgdegree[d] = indegree/len(nodes)
    return davgdegree


def avg_shortestpath(dgassos, graph):
    """
    given a dict which contains a list of node groups, and an igraph object,
    return a dict which gives average shortest path length of these nodes in each
    node group based on this igraph object.
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are their shortest
    path length values
    """
    dsp = {}
    for d in dgassos.keys():
        nodes = graph.vs.select(name_in=dgassos[d])
        sps = graph.shortest_paths(source=nodes, target=nodes, weights=None)
        sumsps = 0
        for ls in sps:
            sumsps += sum(ls)
        if len(nodes) > 1:
            dsp[d] = sumsps/(len(nodes)*(len(nodes)-1))
        else:
            dsp[d] = 0
    return dsp


def lambda_module(dgassos, graph):
    """
    get lambda value from the "lambda-module"
    :param dgassos: a dict object which keys are module names and values are module nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are their lambda values
    """
    dlambda = {}
    for name in dgassos.keys():
        nodes = graph.vs.select(name_in=dgassos[name])
        din = sum(graph.subgraph(nodes).degree())
        dout = sum(graph.degree(vertices=nodes))
        dout = dout-din
        if dout == 0:
            dlambda[name] = float('inf')
        else:
            dlambda[name] = din/dout
    return dlambda
