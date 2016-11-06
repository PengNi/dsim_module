#! /usr/bin/env python3
import math
from files import invert_dict


def sim_geneset2gene(dgassos):
    """

    :param dgassos:
    :return:
    """
    gdassos = invert_dict(dgassos)
    ds = list(dgassos.keys())
    gs = set()
    for d in ds:
        gs.update(dgassos[d])

    geneprob = {}
    for g in gs:
        geneprob[g] = len(gdassos[g]) / len(ds)

    genepairprob = {}
    gs2gscore = {}
    counter = 0
    for d in ds:
        gs2gscore[d] = {}
        dgs = dgassos[d]
        for g in gs:
            dgmi = 0.0
            for dg in dgs:
                if (dg, g) in genepairprob.keys():
                    pdgg = genepairprob[(dg, g)]
                elif (g, dg) in genepairprob.keys():
                    pdgg = genepairprob[(g, dg)]
                else:
                    internum = len(gdassos[dg].intersection(gdassos[g]))
                    pdgg = internum / len(ds)
                    genepairprob[(g, dg)] = pdgg
                if pdgg != 0.0:
                    dgmi += pdgg * math.log2(pdgg / (geneprob[g] * geneprob[dg]))
            gs2gscore[d][g] = dgmi
        counter += 1
        print('sim_geneset2gene:', counter)
    return gs2gscore


def similarity_mutiinfomation(dgassos):
    """

    :param dgassos:
    :return: simdict
    """
    gs2gscore = sim_geneset2gene(dgassos)
    ds = list(dgassos.keys())
    sims = {}
    for i in range(0, len(ds)):
        sims[ds[i]] = {}
        for j in range(i, len(ds)):
            sims[ds[i]][ds[j]] = 0
        print("similarity_mutulinformation:", i)
    return sims
