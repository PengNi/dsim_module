#! /usr/bin/env python3
import math
import time
import numpy
from scipy import spatial
from igraph import Graph
from files import stat_assos
from files import invert_dict
from common_use import hypergeometric_test


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


def similarity_cal_module_1(dgassos, graph, gncutoff=1):
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
        if len(dgleft) >= gncutoff:
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


def similarity_cal_module_2(dgassos, graph, gncutoff=1):
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
        if len(dgleft) >= gncutoff:
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


def similarity_cal_module_3(dgassos, graph, gncutoff=1):
    """
    module cal method 3
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
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    sims = {}
    for i in range(0, len(diseases)):
        sims[diseases[i]] = {}
        for j in range(i, len(diseases)):
            gsi_ori = dgassos[diseases[i]]
            gsj_ori = dgassos[diseases[j]]
            gsintersect_ori = gsi_ori.intersection(gsj_ori)
            conncount = 0
            if diseases[i] in dgassos_new.keys() and diseases[j] in dgassos_new.keys():
                gsi = dgassos_new[diseases[i]]
                gsj = dgassos_new[diseases[j]]
                gsintersect = gsi.intersection(gsj)
                gsid = gsi.difference(gsj)
                gsjd = gsj.difference(gsi)
                if len(gsintersect) != 0:
                    for g in gsid:
                        conncount += connected_count(graph, g, gsintersect)
                    for g in gsjd:
                        conncount += connected_count(graph, g, gsintersect)
                for g in gsid:
                    conncount += connected_count(graph, g, gsjd)
            sims[diseases[i]][diseases[j]] = (len(gsintersect_ori) ** 2 + conncount) / (len(gsi_ori) * len(gsj_ori))
        print(i, "done..")
    return sims


def similarity_cal_module_4(dgassos, graph, gncutoff=1):
    """
    module cal method 4
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
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    allgenes = set()
    for d in diseases:
        allgenes |= dgassos_new[d]
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
            avgic = (ic(len(gsi), len(allgenes)) + ic(len(gsj), len(allgenes)))/2
            sims[diseases[i]][diseases[j]] = (len(gsintersect)**2+conncount)/(len(gsi)*len(gsj))*avgic
        print(i, "done..")
    return sims


def ic(lowern, uppern):
    return -math.log(lowern/uppern, 2)


def similarity_cal_module_5(dgassos, graph, gncutoff=1):
    """
    module cal method 5
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
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    sims = {}
    for i in range(0, len(diseases)-1):
        sims[diseases[i]] = {}
        for j in range(i+1, len(diseases)):
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
            sim = (len(gsintersect)**2+conncount)/(len(gsi)*len(gsj))
            sims[diseases[i]][diseases[j]] = sim
            if diseases[j] not in sims.keys():
                sims[diseases[j]] = {}
            sims[diseases[j]][diseases[i]] = sim
        print(i, "done..")
    normsims = {}
    maxsims = {}
    for d in diseases:
        maxsims[d] = max([i for i in sims[d].values()])
    for i in range(0, len(diseases)-1):
        normsims[diseases[i]] = {}
        for j in range(i+1, len(diseases)):
            avgmaxsim = (maxsims[diseases[i]] + maxsims[diseases[j]])/2
            if avgmaxsim == 0:
                normsims[diseases[i]][diseases[j]] = 0
            else:
                normsims[diseases[i]][diseases[j]] = sims[diseases[i]][diseases[j]]/avgmaxsim
        print(i, "done...")
    return normsims


def similarity_cal_module_6(dgassos, graph, gncutoff=1):
    """
    module cal method 6
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gncutoff: gene number cut off, only diseases whose number of associated
    genes in graph is no less than gncutoff will be calculated
    :return: a dict, (key-value: string-dict<string-float>)
    """
    # only use lcc-----
    glcc = graph.clusters(mode='WEAK').giant()
    gvs = set(glcc.vs['name'])

    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = gvs.intersection(dgassos[d])
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    glccnodes = list(gvs)
    gn2loc = {}
    for i in range(0, len(glccnodes)):
        gn2loc[glccnodes[i]] = i
    adjmat = numpy.zeros(shape=(len(glccnodes), len(glccnodes)), dtype=numpy.int8)
    for i in range(0, len(glccnodes) - 1):
        for j in range(i + 1, len(glccnodes)):
            if glcc.are_connected(glccnodes[i], glccnodes[j]):
                adjmat[i, j] = 1
                adjmat[j, i] = 1
    adjmat2nd = numpy.full(shape=(len(glccnodes), len(glccnodes)), fill_value=-1, dtype=numpy.int8)
    print('sim cal beginning...')
    sims = {}
    for i in range(0, len(diseases)):
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        sims[diseases[i]] = {}
        for j in range(i, len(diseases)):
            gsi = dgassos_new[diseases[i]]
            gsj = dgassos_new[diseases[j]]
            gsilen = len(gsi)
            gsjlen = len(gsj)
            gsid = gsi.difference(gsj)
            gsjd = gsj.difference(gsi)
            sim = 2 * len(gsi.intersection(gsj))
            for g in gsid:
                cc_one, cc_two = 0, 0
                gloc = gn2loc[g]
                for otherg in gsj:
                    othergloc = gn2loc[otherg]
                    cc_one += adjmat[gloc, othergloc]
                    if adjmat2nd[gloc, othergloc] == -1:
                        adjmat2nd[gloc, othergloc] = numpy.dot(adjmat[gloc], adjmat[othergloc])
                        adjmat2nd[othergloc, gloc] = adjmat2nd[gloc, othergloc]
                    if adjmat2nd[gloc, othergloc] > 0:
                        cc_two += 1
                sim += (cc_one * 2 + cc_two) / (gsjlen * 4)
            for g in gsjd:
                cc_one, cc_two = 0, 0
                gloc = gn2loc[g]
                for otherg in gsi:
                    othergloc = gn2loc[otherg]
                    cc_one += adjmat[gloc, othergloc]
                    if adjmat2nd[gloc, othergloc] == -1:
                        adjmat2nd[gloc, othergloc] = numpy.dot(adjmat[gloc], adjmat[othergloc])
                        adjmat2nd[othergloc, gloc] = adjmat2nd[gloc, othergloc]
                    if adjmat2nd[gloc, othergloc] > 0:
                        cc_two += 1
                sim += (cc_one * 2 + cc_two) / (gsilen * 4)
            sims[diseases[i]][diseases[j]] = sim / (len(gsi) + len(gsj))
        print("---------------------------------------cost time:", str(time.time() - now))
    return sims


def similarity_cal_module_7(dgassos, graph, gncutoff=1):
    """
    module cal method 7
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
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    nodenames = list(gvs)
    sps = graph.shortest_paths(source=nodenames, target=nodenames, weights=None, mode=3)
    gn2loc = {}
    for i in range(0, len(nodenames)):
        gn2loc[nodenames[i]] = i


def similarity_cal_spavgn(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """
    use average shortest path and normalize to cal
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gfilter: True or False, use graph to filter disease gene assos or not
    :param gncutoff: gene number cut off, when filter is True, only diseases whose
    number of associated genes in graph is no less than gncutoff will be calculated
    :param transformdistance: True or False, use transform distance or divide
    :return: a dict, (key-value: string-dict<string-float>)
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")
    dselfavg = {}
    for d in diseases:
        avgavg = 0.0
        for g in dgassos_new[d]:
            avgavg += sim_geneset2gene_avg(g, dgassos_new[d], sim_gene2gene)
        dselfavg[d] = avgavg/len(dgassos_new[d])

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[i]], sim_gene2gene)
            osim = (simsum / (len(dgassos_new[diseases[i]]) + len(dgassos_new[diseases[j]])))
            navg = (dselfavg[diseases[i]] + dselfavg[diseases[j]]) / 2

            result[diseases[i]][diseases[j]] = osim/navg
        print("---------------------------------------cost time:", str(time.time()-now))
    return result


def similarity_cal_spavgn_circle(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """
    use average shortest path and normalize to cal
    :param dgassos: a dict object which keys are module names and values are module
    nodes sets
    :param graph: an igraph object
    :param gfilter: True or False, use graph to filter disease gene assos or not
    :param gncutoff: gene number cut off, when filter is True, only diseases whose
    number of associated genes in graph is no less than gncutoff will be calculated
    :param transformdistance: True or False, use transform distance or divide
    :return: a dict, (key-value: string-dict<string-float>)
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")
    dselfavg = {}
    for d in diseases:
        avgavg = 0.0
        for g in dgassos_new[d]:
            avgavg += sim_geneset2gene_avg(g, dgassos_new[d], sim_gene2gene)
        dselfavg[d] = avgavg/len(dgassos_new[d])

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_avg(g, dgassos_new[diseases[i]], sim_gene2gene)
            osim = (simsum / (len(dgassos_new[diseases[i]]) + len(dgassos_new[diseases[j]])))
            navg = (dselfavg[diseases[i]] + dselfavg[diseases[j]]) / 2
            navg = math.sqrt(1-pow(navg-1, 2))

            result[diseases[i]][diseases[j]] = osim/navg
        print("---------------------------------------cost time:", str(time.time()-now))
    return result


def similarity_cal_spmaxn(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """

    :param dgassos:
    :param graph:
    :param gfilter:
    :param gncutoff:
    :param transformdistance:
    :return:
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")
    dselfavg = {}
    for d in diseases:
        avgavg = 0.0
        for g in dgassos_new[d]:
            avgavg += sim_geneset2gene_avg(g, dgassos_new[d], sim_gene2gene)
        dselfavg[d] = avgavg / len(dgassos_new[d])

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[i]], sim_gene2gene)
            osim = (simsum / (len(dgassos_new[diseases[i]]) + len(dgassos_new[diseases[j]])))
            navg = (dselfavg[diseases[i]] + dselfavg[diseases[j]]) / 2

            result[diseases[i]][diseases[j]] = osim / navg
        print("---------------------------------------cost time:", str(time.time() - now))
    return result


def similarity_cal_spmax(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """

    :param dgassos:
    :param graph:
    :param gfilter:
    :param gncutoff:
    :param transformdistance:
    :return:
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[i]], sim_gene2gene)
            result[diseases[i]][diseases[j]] = (simsum /
                                                (len(dgassos_new[diseases[i]]) +
                                                 len(dgassos_new[diseases[j]])))
        print("---------------------------------------cost time:", str(time.time() - now))
    return result


def similarity_cal_spmaxn_circle(dgassos, graph, gfilter=False, gncutoff=1, transformdistance=True):
    """

    :param dgassos:
    :param graph:
    :param gfilter:
    :param gncutoff:
    :param transformdistance:
    :return:
    """
    gvs = set(graph.vs['name'])

    if gfilter:
        dgassos_new = {}
        for d in dgassos.keys():
            dgleft = gvs.intersection(dgassos[d])
            if len(dgleft) >= gncutoff:
                dgassos_new[d] = dgleft
    else:
        dgassos_new = dgassos
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    dgs = set()
    for d in diseases:
        dgs |= set(dgassos_new[d])
    print("disease genes num:", len(dgs))
    sim_gene2gene = sim_gene2gene_shortestpath(dgs, graph, transformdistance)
    print("gene2gene sim cal done..")
    dselfavg = {}
    for d in diseases:
        avgavg = 0.0
        for g in dgassos_new[d]:
            avgavg += sim_geneset2gene_avg(g, dgassos_new[d], sim_gene2gene)
        dselfavg[d] = avgavg / len(dgassos_new[d])

    result = {}
    for i in range(0, len(diseases)):
        result[diseases[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(dgassos_new[diseases[i]]))
        for j in range(i, len(diseases)):
            simsum = 0.0
            for g in dgassos_new[diseases[i]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[j]], sim_gene2gene)
            for g in dgassos_new[diseases[j]]:
                simsum += sim_geneset2gene_max(g, dgassos_new[diseases[i]], sim_gene2gene)
            osim = (simsum / (len(dgassos_new[diseases[i]]) + len(dgassos_new[diseases[j]])))
            navg = (dselfavg[diseases[i]] + dselfavg[diseases[j]]) / 2
            navg = math.sqrt(1 - pow(navg - 1, 2))
            result[diseases[i]][diseases[j]] = osim / navg
        print("---------------------------------------cost time:", str(time.time() - now))
    return result


def similarity_cal_pathwayvector(dgassos, pgassos, graph, simway='SpMin', vectorway='Cosine', gncutoff=1):
    """

    :param dgassos:
    :param pgassos:
    :param graph:
    :param simway:
    :param vectorway:
    :param gncutoff:
    :return:
    """
    print('nodes:', len(graph.vs()), 'edges:', len(graph.es()))
    # ---only use the lcc-----------
    glcc = graph.clusters(mode='WEAK').giant()
    gnodes = set(glcc.vs['name'])
    print('lcc nodes:', len(glcc.vs()), 'edges:', len(glcc.es()))
    stat_assos(dgassos)
    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = gnodes.intersection(dgassos[d])
        if len(dgleft) >= gncutoff:
            dgassos_new[d] = dgleft
    print("disease gene assos left: ", end='')
    stat_assos(dgassos_new)
    stat_assos(pgassos)
    pgassos_new = {}
    for p in pgassos.keys():
        pgleft = gnodes.intersection(pgassos[p])
        if len(pgleft) >= 1:
            pgassos_new[p] = pgleft
    print("pathway gene assos left: ", end='')
    stat_assos(pgassos_new)

    dpvectors = disease2pahtwayvectors(dgassos_new, pgassos_new, glcc, simway)
    return sim_dpvectors(dpvectors, vectorway)


def disease2pahtwayvectors(dgassos, pgassos, graph, simway='SpMin'):
    # node 2 node similarity---------------------------
    nodenames = graph.vs['name']
    sps = graph.shortest_paths(source=nodenames, target=nodenames, weights=None, mode=3)
    g2loc = {}
    for i in range(0, len(nodenames)):
        g2loc[nodenames[i]] = i
    # -------------------------------------------------
    pathways = list(pgassos.keys())
    dpvector = {}
    counter = 0
    if simway == 'SpMin':
        for d in dgassos.keys():
            dpvector[d] = []
            dgenes = dgassos[d]
            for p in pathways:
                pgenes = pgassos[p]
                simsum = 0.0
                for g in dgenes:
                    simsum += minsp_geneset2gene(g, pgenes, g2loc, sps)
                for g in pgenes:
                    simsum += minsp_geneset2gene(g, dgenes, g2loc, sps)
                dpvector[d].append(simsum/(len(dgenes) + len(pgenes)))
            counter += 1
            print("disease2pathwayvectors:", counter)
    elif simway == 'SpTransMax':
        spstrans = []
        for i in range(0, len(nodenames)):
            spstrans.append([])
            for j in range(0, len(nodenames)):
                spstrans[i].append(transformed_distance(sps[i][j]))

        for d in dgassos.keys():
            dpvector[d] = []
            dgenes = dgassos[d]
            for p in pathways:
                pgenes = pgassos[p]
                simsum = 0.0
                for g in dgenes:
                    simsum += maxtrans_geneset2gene(g, pgenes, g2loc, spstrans)
                for g in pgenes:
                    simsum += maxtrans_geneset2gene(g, dgenes, g2loc, spstrans)
                dpvector[d].append(simsum/(len(dgenes) + len(pgenes)))
            counter += 1
            print("disease2pathwayvectors:", counter)
    elif simway == 'SpTransMax_tfidf':
        spstrans = []
        for i in range(0, len(nodenames)):
            spstrans.append([])
            for j in range(0, len(nodenames)):
                spstrans[i].append(transformed_distance(sps[i][j]))
        # ---get disease2pathway assos---
        dpassos = hypergeometric_test(dgassos, pgassos)
        print('disease pathway assos: ', end='')
        stat_assos(dpassos)
        pdassos = invert_dict(dpassos)
        dsnum = len(dpassos.keys())
        pidf = {}
        for p in pdassos.keys():
            pidf[p] = math.log2(dsnum/len(pdassos[p]))
        pathways = list(pdassos.keys())
        p2loc = {}
        for i in range(0, len(pathways)):
            p2loc[pathways[i]] = i
        # -------------------------------
        for d in dpassos.keys():
            dpvector[d] = [0.0] * len(pathways)
            dgenes = dgassos[d]
            for p in dpassos[d]:
                pgenes = pgassos[p]
                simsum = 0.0
                for g in dgenes:
                    simsum += maxtrans_geneset2gene(g, pgenes, g2loc, spstrans)
                for g in pgenes:
                    simsum += maxtrans_geneset2gene(g, dgenes, g2loc, spstrans)
                weightdp = simsum / (len(dgenes) + len(pgenes))
                dpvector[d][p2loc[p]] = pidf[p] * weightdp
            counter += 1
            print("disease2pathwayvectors:", counter)
    elif simway == 'SpAvg':
        pass
    else:
        print('no matched strategy!')
        return None
    return dpvector


def sim_dpvectors(dpvectors, vectorway='Cosine'):
    sims = {}
    ds = list(dpvectors.keys())
    if vectorway == 'Cosine':
        for i in range(0, len(ds)):
            sims[ds[i]] = {}
            for j in range(i, len(ds)):
                sims[ds[i]][ds[j]] = 1 - spatial.distance.cosine(dpvectors[ds[i]], dpvectors[ds[j]])
            print('sim_dpvectors:', i)
    elif vectorway == 'mutualinformation':
        pass
    elif vectorway == 'euclidean':
        for i in range(0, len(ds)):
            sims[ds[i]] = {}
            for j in range(i, len(ds)):
                sims[ds[i]][ds[j]] = spatial.distance.euclidean(dpvectors[ds[i]], dpvectors[ds[j]])
            print('sim_dpvectors:', i)
    elif vectorway == 'euclidean_exp':
        for i in range(0, len(ds)):
            sims[ds[i]] = {}
            for j in range(i, len(ds)):
                sims[ds[i]][ds[j]] = \
                    transformed_distance(spatial.distance.euclidean(dpvectors[ds[i]], dpvectors[ds[j]]))
            print('sim_dpvectors:', i)
    else:
        print('no vectorway specified!')
        return None
    return sims


def minsp_geneset2gene(g1, g2set, gene2loc, spmatrix):
    result = float('inf')
    for g2 in g2set:
        sptemp = spmatrix[gene2loc[g1]][gene2loc[g2]]
        if sptemp < result:
            result = sptemp
    return result


def maxtrans_geneset2gene(g1, g2set, gene2loc, spmatrix):
    result = float('-inf')
    for g2 in g2set:
        sptemp = spmatrix[gene2loc[g1]][gene2loc[g2]]
        if sptemp > result:
            result = sptemp
    return result


def sim_gene2gene_shortestpath(dgs, graph, transformdistance=True):
    nodenames = graph.vs['name']
    sps = graph.shortest_paths(source=nodenames, target=nodenames, weights=None, mode=3)
    result = {}
    for n in nodenames:
        result[n] = {}
    if transformdistance is True:
        for i in range(0, len(nodenames)):
            for j in range(i, len(nodenames)):
                result[nodenames[i]][nodenames[j]] = transformed_distance(sps[i][j])
                result[nodenames[j]][nodenames[i]] = result[nodenames[i]][nodenames[j]]
    else:
        for i in range(0, len(nodenames)):
            for j in range(i, len(nodenames)):
                result[nodenames[i]][nodenames[j]] = 1 / 2**(sps[i][j])
                result[nodenames[j]][nodenames[i]] = result[nodenames[i]][nodenames[j]]
    # supply
    dglefts = set(dgs).difference(set(nodenames))
    for dg in dglefts:
        result[dg] = {}
        for eg in dgs:
            result[dg][eg] = 0.0
        result[dg][dg] = 1.0
    dginters = set(dgs).intersection(set(nodenames))
    for dg in dginters:
        for eg in dglefts:
            result[dg][eg] = 0.0
    return result


def sim_geneset2gene_avg(g, gset, sim_gene2gene):
    sims = []
    for i in gset:
        sims.append(sim_gene2gene[g][i])
    if len(sims) > 0:
        return sum(sims)/len(sims)
    else:
        return 0.0


def sim_geneset2gene_max(g, gset, sim_gene2gene):
    result = 0.0
    for i in gset:
        if sim_gene2gene[g][i] > result:
            result = sim_gene2gene[g][i]
    return result


def katzwalklength(g, k=5, fileprefix='data/interactome_science/interactome_adjmat'):
    """

    :param g: igraph graph object
    :param k: cutoff of walk length
    :param fileprefix:
    :return:
    """
    gvs = set(g.vs['name'])
    gnodes = list(gvs)
    # node2loc = {}
    # for i in range(0, len(gnodes)):
    #     node2loc[gnodes[i]] = i
    adjmat = numpy.zeros(shape=(len(gnodes), len(gnodes)), dtype=numpy.int8)
    for i in range(0, len(gnodes)):
        for j in range(i, len(gnodes)):
            if g.are_connected(gnodes[i], gnodes[j]):
                adjmat[i, j] = 1
                adjmat[j, i] = 1
    print('similarity_module: katzwalklength:', 'preprocess done..')
    kadjmat = 1
    for m in range(1, k+1):
        kadjmat = numpy.dot(kadjmat, adjmat)
        with open(fileprefix+str(m)+'.tsv', mode='w') as wf:
            for g in gnodes:
                wf.write('\t'+g)
            wf.write('\n')
            for i in range(0, len(gnodes)):
                wf.write(gnodes[i])
                for j in range(0, len(gnodes)):
                    wf.write('\t' + str(kadjmat[i, j]))
                wf.write('\n')
        print('similarity_module: katzwalklength:', m)


def similarity_cal_katz(dgassos, gene2loc, adjmats, beta=None):
    """

    :param dgassos:
    :param gene2loc:
    :param adjmats:
    :param beta:
    :return:
    """
    if beta is None:
        beta = numpy.linalg.norm(adjmats[0], 2)
        beta = 1 / beta
        print(beta)
        # beta = 10**(-6)
    genes = set(gene2loc.keys())
    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = genes.intersection(dgassos[d])
        if len(dgleft) >= 1:
            dgassos_new[d] = dgleft
    diseases = list(dgassos_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    katzmat = numpy.zeros(shape=(len(genes), len(genes)), dtype=numpy.float32)
    for i in range(0, len(genes)):
        for j in range(i, len(genes)):
            for c in range(0, len(adjmats)):
                katzmat[i, j] += beta**(c+1) * adjmats[c][i, j]
            katzmat[j, i] = katzmat[i, j]

    print("disease sim cal beginning..")
    sims = {}
    for i in range(0, len(diseases)):
        sims[diseases[i]] = {}
        gsi = dgassos_new[diseases[i]]
        loci = [gene2loc[g] for g in gsi]
        for j in range(i, len(diseases)):
            gsj = dgassos_new[diseases[j]]
            locj = [gene2loc[g] for g in gsj]
            simtemp = 0.0
            for loc in loci:
                simtemp += max(katzmat[loc, locj])
            for loc in locj:
                simtemp += max(katzmat[loc, loci])
            sims[diseases[i]][diseases[j]] = simtemp / (len(gsi) + len(gsj))
        print("similarity_cal_katz:", i)
    return sims


def read_adjmatrix(filepath):
    """

    :param filepath:
    :return:
    """
    adjlist = []
    with open(filepath, mode='r') as rf:
        next(rf)
        for line in rf:
            words = line.strip().split('\t')
            adjlist.append([int(i) for i in words[1:len(words)]])
    return numpy.array(adjlist)


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


def get_lcc(graphfile):
    """
    given a graph ncol file path, return the largest connected component
    :param graphfile: path of a graph ncol file
    :return: a igraph.Graph object, the lcc
    """
    with open(graphfile, 'r') as f:
        g = Graph.Read_Ncol(f, names=True, weights=False, directed=False)
    print('nodes:', len(g.vs))
    print('edges:', len(g.es))
    glcc = g.clusters(mode='WEAK').giant()
    print('lcc nodes:', len(glcc.vs))
    print('lcc edges:', len(glcc.es))
    return glcc


def write_lcc():
    glcc = get_lcc('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop.tab')
    # glccnodes = glcc.vs['name']
    # glccnodes2loc = {}
    # for i in range(0, len(glccnodes)):
    #     glccnodes2loc[glccnodes[i]] = i
    # with open('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop_maxcc_nodes_numbered.tab', 'w') as f:
    #     for gln in glccnodes:
    #         f.write(gln+'\t'+str(glccnodes2loc[gln])+'\n')
    #
    # with open('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop_maxcc_edges_numbered.tab', 'w') as f:
    #     gles = glcc.es
    #     f.write(str(len(glcc.vs))+' '+str(len(glcc.es))+'\n')
    #     for e in gles:
    #         ns = glcc.vs[e.source]['name']
    #         nt = glcc.vs[e.target]['name']
    #         f.write(str(glccnodes2loc[ns])+' '+str(glccnodes2loc[nt])+'\n')
    pass
