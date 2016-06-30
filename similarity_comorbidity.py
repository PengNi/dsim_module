#! /usr/bin/env python3
"""calculate similarity of disease pairs based on comorbidity data.
Reference:Hidalgo C A, Blumm N, Barabási A L, et al. A dynamic network approach
for the study of human phenotypes[J]. PLoS Comput Biol, 2009, 5(4): e1000353.
"""


# read "AllNet3.net" file
def read_phicor_ttestval_allnet3(filepath):
    """
    read "AllNet3.net" file to get each disease pair's phi-correlation and its corresponding
    ttest value
    :param filepath: "AllNet3.net"
    :return: a dict, key-value pattern is {disease1:{disease2:{phi: <float>, ttest: <float>}, }, }
    """
    phis = {}
    with open(filepath, mode='r') as f:
        for line in f:
            words = line.split("\t")
            disease1 = words[0].strip()
            disease2 = words[1].strip()
            phi = float(words[8].strip())
            ttest = float(words[9].strip())
            if disease1 not in phis.keys():
                phis[disease1] = {}
            phis[disease1][disease2] = {'phi': phi, 'ttest': ttest}
    return phis
