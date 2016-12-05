#! /usr/bin/env python3
from random import randint
import numpy as np


def transformed_distance(shortestpathdis=0, a=1, b=1):
    return a*np.exp(-b*shortestpathdis)


def sim_geneset2gene_avg(g, gset, sim_gene2gene):
    sims = []
    for i in gset:
        sims.append(sim_gene2gene[g][i])
    return sum(sims)/len(sims)


def similarity_cal(d2g, d2t, tissuegraphs):
    tissue_genesims = {}
    tissuegeneset = {}
    for t in tissuegraphs.keys():
        nodenames = tissuegraphs[t].vs['name']
        tissuegeneset[t] = set(nodenames)
        sps = tissuegraphs[t].shortest_paths(source=nodenames, target=nodenames, weights=None, mode=3)
        result = {}
        for n in nodenames:
            result[n] = {}
        for i in range(0, len(nodenames)):
            for j in range(i, len(nodenames)):
                result[nodenames[i]][nodenames[j]] = transformed_distance(sps[i][j])
                result[nodenames[j]][nodenames[i]] = result[nodenames[i]][nodenames[j]]
        tissue_genesims[t] = result
        print(t)
    print("tissue genesim done..")
    d2tg = {}
    dselfsim = {}
    for d in d2t.keys():
        d2tg[d] = {}
        dselfsim[d] = {}
        for t in d2t[d]:
            d2tg[d][t] = d2g[d].intersection(tissuegeneset[t])
            avgsim = 0.0
            for g in d2tg[d][t]:
                avgsim += sim_geneset2gene_avg(g, d2tg[d][t], tissue_genesims[t])
            dselfsim[d][t] = avgsim / len(d2tg[d][t])
    print("disease selfsim done..")
    sims = {}
    ds = list(d2tg.keys())
    for i in range(0, len(ds)):
        sims[ds[i]] = {}
        for j in range(i, len(ds)):
            commont = d2t[ds[i]].intersection(d2t[ds[j]])
            if len(commont) == 0:
                sims[ds[i]][ds[j]] = 0.0
            else:
                simij = []
                for t in commont:
                    gsi = d2tg[ds[i]][t]
                    gsj = d2tg[ds[j]][t]
                    sim = 0.0
                    for g in gsi:
                        sim += sim_geneset2gene_avg(g, gsj, tissue_genesims[t])
                    for g in gsj:
                        sim += sim_geneset2gene_avg(g, gsi, tissue_genesims[t])
                    sim /= (len(gsi) + len(gsj))
                    simij.append(sim * 2 / (dselfsim[ds[i]][t] + dselfsim[ds[j]][t]))
                sims[ds[i]][ds[j]] = max(simij)
        print('sim_cal:', i)
    return sims


def find_diseasetissueassos(tissuegraphs, d2g, gcutoff=5, zcutoff=1.6):
    """

    :param tissuegraphs:
    :param d2g:
    :param gcutoff:
    :param zcutoff:
    :return:
    """
    d2t = {}
    for d in d2g.keys():
        d2t[d] = set()
    for t in tissuegraphs.keys():
        tgraph = tissuegraphs[t]
        tvs = tgraph.vs['name']
        for d in d2g.keys():
            dgs_t = d2g[d].intersection(set(tvs))
            lendgs_t = len(dgs_t)
            if lendgs_t >= gcutoff:
                len_dlcc_t = tgraph.subgraph(tgraph.vs.select(name_in=dgs_t)).clusters(mode='WEAK').giant().vcount()
                len_randlcc_t = []
                for i in range(0, 1000):
                    randloc = set()
                    while len(randloc) < lendgs_t:
                        randloc.add(randint(0, len(tvs)-1))
                    randnodes = [tvs[i] for i in randloc]
                    len_randlcc_t.append(tgraph.subgraph(tgraph.vs.select(name_in=randnodes)).clusters(mode='WEAK')
                                         .giant().vcount())
                zscore = (len_dlcc_t - np.mean(len_randlcc_t)) / np.std(len_randlcc_t)
                if zscore >= zcutoff:
                    d2t[d].add(t)
        print(t, len(tvs))
    ds = list(d2t.keys())
    for d in ds:
        if len(d2t[d]) == 0:
            del d2t[d]
    return d2t


def find_tissuesubgraphs(tissue2gene, graph):
    """

    :param tissue2gene:
    :param graph:
    :return: dict, tissue-igraph object
    """
    tissuegraphs = {}
    for tissue in tissue2gene.keys():
        nodes = graph.vs.select(name_in=tissue2gene[tissue])
        tissuegraphs[tissue] = graph.subgraph(nodes)
    return tissuegraphs


def find_tissuegenes(tissuegeneexp, cutoff=1.0):
    """

    :param tissuegeneexp: filename of tissue specific gene expression data,
    e.g. 'data/tissuespec_srep/srep35241-s4.txt'
    :param cutoff: gene exp z-score cutoff
    :return: dict, tissue-set of genes, string-set<string>
    """
    gene_exp = {}
    with open(tissuegeneexp, mode='r') as f:
        fstwords = next(f).strip().split('\t')
        tissues = [fstwords[i] for i in range(1, len(fstwords))]
        for line in f:
            words = line.strip().split('\t')
            valvector = []
            for i in range(1, len(words)):
                valvector.append(float(words[i].strip()))
            gene_exp[words[0].strip()] = valvector
    tissue2gene = {}
    for tissue in tissues:
        tissue2gene[tissue] = set()
    for gene in gene_exp.keys():
        geneexpvalue = gene_exp[gene]
        for i in range(0, len(geneexpvalue)):
            if geneexpvalue[i] >= cutoff:
                tissue2gene[tissues[i]].add(gene)
    return tissue2gene
