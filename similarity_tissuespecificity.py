#! /usr/bin/env python3
from random import randint
import numpy as np


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
