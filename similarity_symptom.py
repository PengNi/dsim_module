#! /usr/bin/env python3
"""calculate similarity of disease pairs based on disease-symptom associations.
Reference: Zhou X Z, Menche J, Barabási A L, et al. Human symptoms–disease network[J].
Nature communications, 2014, 5.
"""


# read "ncomms5212-s5.txt" file
def read_terms_similarity(filepath):
    """
    read "ncomms5212-s5.txt" file to get symptom similarity between mesh terms
    :param filepath: "ncomms5212-s5.txt"
    :return: a dict (key-value: string-dict<string-float>)
    """
    sim = {}
    with open(filepath, mode='r') as f:
        next(f)
        for line in f:
            words = line.split("\t")
            disease1 = words[0].strip()
            disease2 = words[1].strip()
            simvalue = float(words[2].strip())
            if disease1 not in sim.keys():
                sim[disease1] = {}
            sim[disease1][disease2] = simvalue
    return sim
