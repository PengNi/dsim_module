#! /usr/bin/env python3
"""calculate similarity of disease pairs based on comorbidity data.
Reference:Hidalgo C A, Blumm N, Barabási A L, et al. A dynamic network approach
for the study of human phenotypes[J]. PLoS Comput Biol, 2009, 5(4): e1000353.
"""


# read "AllNet3.net" file
def read_rr_allnet3(filepath, rrvalue_col):
    """
    read "AllNet3.net" file to get each disease pair's relative risk value
    :param filepath: "AllNet3.net"
    :param rrvalue_col: the col where the rr values are at
    :return: a dict, key-value is {disease1:{disease2: <float>, }, }
    """
    rrs = {}
    rrvalue_col -= 1
    with open(filepath, mode='r') as f:
        for line in f:
            words = line.split("\t")
            disease1 = words[0].strip()
            disease2 = words[1].strip()
            rr = float(words[rrvalue_col].strip())
            if disease1 not in rrs.keys():
                rrs[disease1] = {}
            rrs[disease1][disease2] = rr
    return rrs


# read "AllNet3.net" file
def read_phicor_ttestval_allnet3(filepath, phivalue_col, ttestvalue_col):
    """
    read "AllNet3.net" file to get each disease pair's phi-correlation and its corresponding
    t-test value
    :param filepath: "AllNet3.net"
    :param phivalue_col: the col where the phi-correlation values are at
    :param ttestvalue_col: the col where the t-test values are at
    :return: a dict, key-value is {disease1:{disease2:{phi: <float>, ttest: <float>}, }, }
    """
    phis = {}
    phivalue_col -= 1
    ttestvalue_col -= 1
    with open(filepath, mode='r') as f:
        for line in f:
            words = line.split("\t")
            disease1 = words[0].strip()
            disease2 = words[1].strip()
            phi = float(words[phivalue_col].strip())
            ttest = float(words[ttestvalue_col].strip())
            if disease1 not in phis.keys():
                phis[disease1] = {}
            phis[disease1][disease2] = {'phi': phi, 'ttest': ttest}
    return phis
