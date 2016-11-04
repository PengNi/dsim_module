#! /usr/bin/env python3
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
from files import stat_assos


def normalize_simdict(simdict):
    """
    normalize a simdict, xc = (x-minx)/(max-minx)
    :param simdict: a dict contains sim values, (key-value: string-dict<string-float>)
    :return: a dict, (key-value: string-dict<string-float>)
    """
    minv = float('inf')
    maxv = float('-inf')
    for d in simdict.keys():
        for k in simdict[d].keys():
            if minv > simdict[d][k]:
                minv = simdict[d][k]
            if maxv < simdict[d][k]:
                maxv = simdict[d][k]
    nsimdict = {}
    print(maxv, minv)
    delta = maxv - minv
    for d in simdict.keys():
        nsimdict[d] = {}
        for k in simdict[d].keys():
            nsimdict[d][k] = (simdict[d][k] - minv)/delta
    return nsimdict


def hypergeometric_test(acassos, bcassos, mtmethod='fdr_bh', pvcutoff=0.05):
    """
    based on A-C and B-C associations, using hypergeometric test to get
    A-B associations
    :param acassos: dict, string-set<string>,
    e.g. A-C (e.g. disease-gene) associations
    :param bcassos: dict, string-set<string>,
    e.g. B-C (e.g. pathway-gene) associations
    :param mtmethod: multiple tests p-value correction method,
    `bonferroni` : one-step correction
    `sidak` : one-step correction
    `holm-sidak` : step down method using Sidak adjustments
    `holm` : step-down method using Bonferroni adjustments
    `simes-hochberg` : step-up method  (independent)
    `hommel` : closed method based on Simes tests (non-negative)
    `fdr_bh` : Benjamini/Hochberg  (non-negative)
    `fdr_by` : Benjamini/Yekutieli (negative)
    `fdr_tsbh` : two stage fdr correction (non-negative)
    `fdr_tsbky` : two stage fdr correction (non-negative)
    :param pvcutoff, p-value cutoff to tell if a test is significant or not
    :return: dict, string-set<string>, A-B (e.g. disease-pathway) associations
    """
    print("A-C assos: ", end='')
    stat_assos(acassos)
    print("B-C assos: ", end='')
    stat_assos(bcassos)
    bcs = set()
    for b in bcassos.keys():
        bcs |= bcassos[b]
    # -----------------------------------------------
    a_s = list(acassos.keys())
    acassos_new = {}
    for a in a_s:
        acassos_new[a] = bcs.intersection(acassos[a])
        if len(acassos_new[a]) == 0:
            del acassos_new[a]
    print("filtered A-C: ", end='')
    stat_assos(acassos_new)
    # -----------------------------------------------
    bcslen = len(bcs)
    a_s = list(acassos_new.keys())
    bs = list(bcassos.keys())
    abassos = {}
    for a in a_s:
        aclen = len(acassos_new[a])
        pvalue = []
        bs_cal = []
        for b in bs:
            abclen = len(acassos_new[a].intersection(bcassos[b]))
            if abclen > 0:
                bs_cal.append(b)
                pglen = len(bcassos[b])
                pvalue.append(1 - hypergeom.cdf(abclen-1, bcslen, pglen, aclen))
        if len(bs_cal) > 0:
            qvalue_result = multipletests(pvalue, alpha=pvcutoff, method=mtmethod,
                                          is_sorted=False, returnsorted=False)
            abassos[a] = set()
            for i in range(0, len(bs_cal)):
                if qvalue_result[0][i]:
                    abassos[a].add(bs_cal[i])
            if len(abassos[a]) == 0:
                del abassos[a]
    print("common_use.hypergeometric_test(): done..")
    return abassos
