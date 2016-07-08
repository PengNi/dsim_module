#! /usr/bin/env python3
"""calculate similarity of disease pairs based on disease-symptom associations.
Reference: Zhou X Z, Menche J, Barabási A L, et al. Human symptoms–disease network[J].
Nature communications, 2014, 5.
"""
from copy import deepcopy


def termnamemap2umlsid(meshnames, umlsdiseases, meshname2meshid):
    """
    convert mesh term names to umls ids
    :param meshnames: set contains mesh names, may be contains all
    mesh names from "ncomms5212-s5.txt" file
    :param umlsdiseases: an UmlsDiseases object
    :param meshname2meshid: dict from method read_mappings("MeshTreeHierarchy.csv")
    :return: a dict, key-value: string(mesh term name)-string(umls id)
    """
    meshname2umlsid = {}
    umlsname2umlsid = {}
    umlsname2umlsid_nor = {}
    meshid2umlsid = {}
    umlsdis = umlsdiseases.getumlsdiseases()
    for u in umlsdis.keys():
        umlsname2umlsid_nor[umlsdis[u].getname()] = umlsdis[u].getid()
        umlsname2umlsid[umlsdis[u].getname().lower()] = umlsdis[u].getid()
        for m in umlsdis[u].getmeshids():
            if m not in meshid2umlsid.keys():
                meshid2umlsid[m] = set()
            meshid2umlsid[m].add(umlsdis[u].getid())
    # print("termnamemap2umlsid(): umlsname2umlsid:", end="")
    # stat_maps(umlsname2umlsid)
    # print("termnamemap2umlsid(): meshid2umlsid:", end="")
    # stat_assos(meshid2umlsid)

    perfect = 0
    umlsnamemap = 0
    meshidmap = 0
    good = 0
    bad = 0
    badbad = 0
    notgood = 0
    notgoodeither = 0
    for mn in meshnames:
        umlsidtemp = ""
        umlsidtemps = set()
        if mn.lower() in umlsname2umlsid.keys():
            umlsidtemp = umlsname2umlsid[mn.lower()]

        meshname2meshid_lowerkeys = {}
        for m2m in meshname2meshid.keys():
            meshname2meshid_lowerkeys[m2m.lower()] = deepcopy(meshname2meshid[m2m])

        if mn.lower() in meshname2meshid_lowerkeys.keys():
            meshidtemp = "mesh:"+meshname2meshid_lowerkeys[mn.lower()]
            if meshidtemp in meshid2umlsid.keys():
                umlsidtemps = meshid2umlsid[meshidtemp]

        if umlsidtemp == "":
            if len(umlsidtemps) == 0:
                badbad += 1
            elif len(umlsidtemps) == 1:
                meshidmap += 1
                meshname2umlsid[mn] = set()
                meshname2umlsid[mn].add(list(umlsidtemps)[0])
            else:
                bad += 1
                meshname2umlsid[mn] = deepcopy(umlsidtemps)
        else:
            # if mn not in meshname2meshid and mn.lower() in meshname2meshid_lowerkeys:
            #     print("meshname:", mn, ";", mn.lower())
            # if mn not in umlsname2umlsid_nor and mn.lower() in umlsname2umlsid:
            #     print("umlsname:", mn, ";", mn.lower())
            if len(umlsidtemps) == 0:
                umlsnamemap += 1
                meshname2umlsid[mn] = set()
                meshname2umlsid[mn].add(umlsidtemp)
            elif len(umlsidtemps) == 1:
                if umlsidtemp == list(umlsidtemps)[0]:
                    perfect += 1
                    meshname2umlsid[mn] = set()
                    meshname2umlsid[mn].add(umlsidtemp)
                else:
                    notgood += 1
                    meshname2umlsid[mn] = set()
                    meshname2umlsid[mn].add(umlsidtemp)
                    # print("11 but not the same:", mn, umlsidtemp, umlsidtemps)
            else:
                if umlsidtemp in umlsidtemps:
                    good += 1
                    meshname2umlsid[mn] = set()
                    meshname2umlsid[mn].add(umlsidtemp)
                else:
                    notgoodeither += 1
    print("using umlsnames to map got 0 result, using mesh ids to map got"
          " 0 result:", badbad)
    print("using umlsnames to map got 0 result, using mesh ids to map got"
          " 1 result:", meshidmap)
    print("using umlsnames to map got 0 result, using mesh ids to map got"
          " multi results:", bad)
    print("using umlsnames to map got 1 result, using mesh ids to map got"
          " 0 result:", umlsnamemap)
    print("using umlsnames to map got 1 result, using mesh ids to map got"
          " 1 result, and these 2 results are the same:", perfect)
    print("using umlsnames to map got 1 result, using mesh ids to map got"
          " 1 result, and these 2 results are different:", notgood)
    print("using umlsnames to map got 1 result, using mesh ids to map got"
          " multi results, and 1 result by umlsnames is in multi results:", good)
    print("using umlsnames to map got 1 result, using mesh ids to map got"
          " multi results, but 1 result by umlsnames is not in multi results:", notgoodeither)
    return meshname2umlsid


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
