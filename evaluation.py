#! /usr/bin/env python3
"""Evaluation for disease-disease association computing methods"""
from files import read_sims
from files import stat_assos
from copy import deepcopy
from random import randint


def eva_70benchmarkpairs(methodfilepaths, benchmarkpairs, times=1):
    """
    use 70 disease pairs as benchmark set to evaluation disease similarity methods
    reference: Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring disease
    similarity by integrating semantic and gene functional association[J]. PloS one,
    2014, 9(6): e99415.
    :param methodfilepaths: a list of file paths of methods' results
    :param benchmarkpairs: list of tuples which each tuple contains two diseases,
    can be obtained from file "data/ground_truth_68_disease_pairs_umlsid.tsv"
    :param times: times to replicate the evaluation, one time produces one result which
    can be used to draw ROC
    :return: a list of dicts, one dict contains one result,
    dict: {disease1: {disease2:{method1: simvalue, method2: simvalue, ..., label: 0/1}, }, }
    """
    sims = {}
    for fp in methodfilepaths:
        sims[fp] = read_sims(fp)
    print("reading sim files completed!")

    dpairs = {}
    for dk in sims[methodfilepaths[0]].keys():
        dpairs[dk] = set(sims[methodfilepaths[0]][dk].keys())
    for i in range(1, len(methodfilepaths)):
        simtemp = sims[methodfilepaths[i]]
        for dpk in dpairs.keys():
            scopy = deepcopy(dpairs[dpk])
            for s in scopy:
                if not ((dpk in simtemp.keys() and s in simtemp[dpk].keys()) or
                        (s in simtemp.keys() and dpk in simtemp[s].keys())):
                    dpairs[dpk].remove(s)
    ds = list(dpairs.keys())
    for d in ds:
        dpairs[d].discard(d)
        if len(dpairs[d]) == 0:
            del dpairs[d]
    print("number of disease pairs which all methods can cal: ", end='')
    stat_assos(dpairs)
    for d in dpairs.keys():
        dpairs[d] = list(dpairs[d])

    bencallist = []
    for d1, d2 in benchmarkpairs:
        if ((d1 in dpairs.keys() and d2 in dpairs[d1]) or
                (d2 in dpairs.keys() and d1 in dpairs[d2])):
            bencallist.append((d1, d2))
    print("number of benchmark disease pairs which all methods can cal:", len(bencallist))

    ress = []
    for i in range(0, times):
        ress.append(eva_70benchmarkpairs_onetime(sims, dpairs, bencallist))
    return ress


def findsimvalue(d1, d2, dsims):
    """
    find d1 and d2's simvalue in dsims
    :param d1: string
    :param d2: string
    :param dsims: dict, key-value: string-dict<string-value> ({entity1: {entity2: sim, }, }
    :return: float
    """
    if d1 in dsims.keys() and d2 in dsims[d1].keys():
        return dsims[d1][d2]
    else:
        return dsims[d2][d1]


def eva_70benchmarkpairs_onetime(allsims, alldpairs, bencallist):
    """
    use 70benchmarkpairs to evaluate for one time
    :param allsims: dict, {sim1: { d1: {d2: simvalue, }, }, sim2: {}, }
    :param alldpairs: dict, {d1: [d2, d3], d2: [], }
    :param bencallist: list of tuples
    :return: dict, dict patterm is
    {disease1: {disease2:{method1: simvalue, method2: simvalue, ..., label: 0/1}, }, }
    """
    d1s = list(alldpairs.keys())
    randomlist = []
    randomsetlen = 10 * len(bencallist)
    while len(randomlist) < randomsetlen:
        d1 = d1s[randint(0, len(d1s)-1)]
        d2 = alldpairs[d1][randint(0, len(alldpairs[d1])-1)]
        if (((d1, d2) not in bencallist) and ((d2, d1) not in bencallist) and
                ((d1, d2) not in randomlist) and ((d2, d1) not in randomlist)):
            randomlist.append((d1, d2))

    res = {}
    for d1, d2 in bencallist:
        if d1 not in res.keys():
            res[d1] = {}
        res[d1][d2] = {}
        for sim in allsims.keys():
            res[d1][d2][sim] = findsimvalue(d1, d2, allsims[sim])
        res[d1][d2]['label'] = 1
    for d1, d2 in randomlist:
        if d1 not in res.keys():
            res[d1] = {}
        res[d1][d2] = {}
        for sim in allsims.keys():
            res[d1][d2][sim] = findsimvalue(d1, d2, allsims[sim])
        res[d1][d2]['label'] = 0
    return res