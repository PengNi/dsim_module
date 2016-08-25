#! /usr/bin/env python3
"""Evaluation for disease-disease association computing methods"""
from files import read_sims
from files import stat_assos
from copy import deepcopy
from random import randint
from operator import itemgetter


def eva_readsims(methodfilepaths):
    """
    read sims from all files in methodfilepaths
    :param methodfilepaths: a list of file paths (strings)
    :return: dict, keys are file paths, values are results from method
    read_sims()
    """
    sims = {}
    for fp in methodfilepaths:
        sims[fp] = read_sims(fp)
    print("reading sim files completed!")
    return sims


def eva_groundtruth(sims, groundtruthfilepath, topn=5000):
    """
    evaluate disease sim methods from methodfilepaths based on the disease sim
    scores in groundtruthfilepath
    :param sims: got from method eva_readsims(), contains sims of methods need
    to be evaluated
    :param groundtruthfilepath: file path of ground truth disease sim scores
    :param topn: get topn pairs have the highest scores for evaluating
    :return: to be decided
    """
    groundtruth_sim = read_sims(groundtruthfilepath)
    print("reading groundtrthfile completed!")

    dpairs = {}
    for dk in groundtruth_sim.keys():
        dpairs[dk] = set(groundtruth_sim[dk].keys())
    methodfilepaths = list(sims.keys())
    for i in range(0, len(methodfilepaths)):
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

    topnpairs = {}
    for fp in methodfilepaths:
        simtemp = []
        for d in dpairs.keys():
            for p in dpairs[d]:
                simtemp.append((d, p, findsimvalue(d, p, sims[fp])))
        simtemp = sorted(simtemp, key=itemgetter(2), reverse=True)
        topnpairs[fp] = []
        if len(simtemp) < topn:
            topn = len(simtemp)
        for i in range(0, topn):
            topnpairs[fp].append((simtemp[i][0], simtemp[i][1], simtemp[i][2],
                                  findsimvalue(simtemp[i][0], simtemp[i][1], groundtruth_sim)))
    return topnpairs


def eva_tprfprs(scoredicts):
    """
    transform scores and labels to tprs and fprs
    :param scoredicts: a list of dicts, one dict contains one result,
    dict: {disease1: {disease2:{method1: simvalue, method2: simvalue, ..., label: 0/1}, }, }
    (from eva_70benchmarkpairs())
    :return: a list of dicts, each dict contains tpr and tpr for one or more methods
    """
    ress = []
    for scoredict in scoredicts:
        ress.append(eva_tprfpr(eva_ranking(scoredict)))
    return ress


def eva_ranking(scoredict):
    """
    based on scores and label in scoredict, ranking disease pairs in the purpose of
    calculating tpr and fpr
    :param scoredict: dict: {disease1: {disease2:{method1: simvalue, method2: simvalue, ...,
    label: 0/1}, }, } (in this case, scoredict is from eva_70benchmarkpairs_onetime())
    :return: a dict which each value is a sorted list based on a sort strategy,
    {method1: [(disease1, disease2, simscore, label), ], method2: [], }
    """
    scorenlabel = {}
    d1 = list(scoredict.keys())[0]
    d2 = list(scoredict[d1].keys())[0]
    methodnames = set(scoredict[d1][d2].keys())
    methodnames.discard('label')
    for m in methodnames:
        scorenlabel[m] = []

    for d1 in scoredict.keys():
        for d2 in scoredict[d1].keys():
            for m in methodnames:
                scorenlabel[m].append((d1, d2, scoredict[d1][d2][m], scoredict[d1][d2]['label']))
    for m in methodnames:
        # sort strategy---------------------------------------------------------
        scorenlabel[m] = sorted(scorenlabel[m], key=itemgetter(3))
        scorenlabel[m] = sorted(scorenlabel[m], key=itemgetter(2), reverse=True)
        # ----------------------------------------------------------------------
    return scorenlabel


def eva_tprfpr(scorenlabel):
    """
    based on all kind of scores and label in scoredict, calculate the tpr and fpr
    of each kind of score
    :param scorenlabel: a dict which each value is a sorted list based on a sort strategy,
    {method1: [(disease1, disease2, simscore, label), ], method2: [], }
    :return: dict, {method1: [(float1, float2), ], method2: [], } (float1 is fpr, float2 is tpr)
    """
    tpfp = {}
    methodnames = set(scorenlabel.keys())
    for m in methodnames:
        tpfp[m] = []
        sltemp = scorenlabel[m]
        slable = []
        for i in range(0, len(sltemp)):
            slable.append(sltemp[i][3])
        conditionp = sum(slable)
        conditionf = len(slable) - conditionp
        tp = 0
        tpfp[m].append((0.0, 0.0))
        for i in range(0, len(slable)):
            tp += slable[i]
            fp = i+1-tp
            tpfp[m].append((fp/conditionf, tp/conditionp))
    return tpfp


def eva_70benchmarkpairs(sims, benchmarkpairs, times=1):
    """
    use 70 disease pairs as benchmark set to evaluation disease similarity methods
    reference: Cheng L, Li J, Ju P, et al. SemFunSim: a new method for measuring disease
    similarity by integrating semantic and gene functional association[J]. PloS one,
    2014, 9(6): e99415.
    :param sims: got from method eva_readsims(), contains sims of methods need
    to be evaluated
    :param benchmarkpairs: list of tuples which each tuple contains two diseases,
    can be obtained from file "data/ground_truth_68_disease_pairs_umlsid.tsv"
    :param times: times to replicate the evaluation, one time produces one result which
    can be used to draw ROC
    :return: a list of dicts, one dict contains one result,
    dict: {disease1: {disease2:{method1: simvalue, method2: simvalue, ..., label: 0/1}, }, }
    """
    methodfilepaths = list(sims.keys())
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
    :return: dict, dict pattern is
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
