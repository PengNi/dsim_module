#! /usr/bin/env python3
"""calculate similarity of disease pairs based on disease-symptom associations.
Reference: Zhou X Z, Menche J, Barabási A L, et al. Human symptoms–disease network[J].
Nature communications, 2014, 5.
"""


def termnamesim2termidsim(termname_sim, termname2termid):
    """
    the "ncomms5212-s5.txt" file has symptom similairty between meah term names,
    this method converts termnames to ids based on termname-termid mapping relationships.
    :param termname_sim: dict obtained from method read_terms_similarity()
    :param termname2termid: dict, has termname2termid mappings, obtained from
    "MeshTreeHierarchy.csv"
    :return: dcit (key-value: string-dict<string-float>), which contains similarities
    between mesh term ids.
    """
    sim = {}
    for n1 in termname_sim.keys():
        for n2 in termname_sim[n1].keys:
            if n1 in termname2termid.keys() and n2 in termname2termid.keys():
                if termname2termid[n1] not in sim.keys():
                    sim[termname2termid[n1]] = {}
                sim[termname2termid[n1]][termname2termid[n2]] = termname_sim[n1][n2]
    return sim


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
