#! /usr/bin/env python3
"""calculate similarity of gene or disease pairs based on GO"""


def similarity_cal(dgassos, g2gsimmat):
    diseases = list(dgassos.keys())
    sims = {}
    for i in range(0, len(diseases)):
        sims[diseases[i]] = {}
        gsi = dgassos[diseases[i]]
        for j in range(i, len(diseases)):
            gsj = dgassos[diseases[j]]

            gsij = gsi.intersection(gsj)
            gsid = gsi.difference(gsj)
            gsjd = gsj.difference(gsi)
            simtemp = 2 * len(gsij)
            for g in gsid:
                simtemp += sim_geneset2gene_max(g, gsj, g2gsimmat)
            for g in gsjd:
                simtemp += sim_geneset2gene_max(g, gsi, g2gsimmat)
            sims[diseases[i]][diseases[j]] = simtemp/(len(gsi) + len(gsj))
        print('similarity_cal:', i)
    return sims


def simmatrix(genes, g2gsim):
    print('simmatrix begin..')
    sims = {}
    genes = list(genes)
    for g in genes:
        sims[g] = {}
    for i in range(0, len(genes)):
        for j in range(i, len(genes)):
            if genes[i] in g2gsim.keys() and genes[j] in g2gsim[genes[i]].keys():
                sims[genes[i]][genes[j]] = g2gsim[genes[i]][genes[j]]
            elif genes[j] in g2gsim.keys() and genes[i] in g2gsim[genes[j]].keys():
                sims[genes[i]][genes[j]] = g2gsim[genes[j]][genes[i]]
            else:
                if genes[i] == genes[j]:
                    sims[genes[i]][genes[j]] = 1.0
                else:
                    sims[genes[i]][genes[j]] = 0.0
            sims[genes[j]][genes[i]] = sims[genes[i]][genes[j]]
    print('simmatrix end..')
    return sims


def simmatrix_combine(genes, g2gsim1, g2gsim2, alpha=0.5):
    print('simmatrix begin..')
    sims = {}
    genes = list(genes)
    for g in genes:
        sims[g] = {}
    for i in range(0, len(genes)):
        for j in range(i, len(genes)):
            sims[genes[i]][genes[j]] = alpha * findsimvalue(genes[i], genes[j], g2gsim1) + \
                                       (1 - alpha) * findsimvalue(genes[i], genes[j], g2gsim2)
            sims[genes[j]][genes[i]] = sims[genes[i]][genes[j]]
    print('simmatrix end..')
    return sims


def sim_geneset2gene_max(g, gset, sim_gene2gene):
    result = 0.0
    for i in gset:
        if sim_gene2gene[g][i] > result:
            result = sim_gene2gene[g][i]
    return result


def findsimvalue(g1, g2, gsims):
    if g1 in gsims.keys() and g2 in gsims[g1].keys():
        ggsim = gsims[g1][g2]
    elif g2 in gsims.keys() and g1 in gsims[g2].keys():
        ggsim = gsims[g2][g1]
    else:
        if g1 == g2:
            ggsim = 1.0
        else:
            ggsim = 0.0
    return ggsim
