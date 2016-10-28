#! /usr/bin/env python3
"""functions about read from and write to files"""
import copy


def read_one_col(filepath, col, header=False, sep="\t", encoding='utf-8'):
    """
    read a column from a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param encoding: encoding
    :param sep: delimiter, a string, default "\t"
    :param col: the number of col needs to be extracted
    :return: a list contains elements from the extracted column
    """
    eles = []
    col -= 1
    with open(filepath, encoding=encoding, mode='r') as f:
        if header:
            next(f)
        for line in f:
            words = line.strip().split(sep=sep)
            if len(words) > col:
                eles.append(words[col].strip())
    return eles


def read_assos(filepath, header=False, sep="\t", xcol=1, ycol=2):
    """
    read associations (e.g. disease-gene associations) into key-value format from
    a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param xcol: number of the column contains the dict's key, default 1
    :param ycol: number of the column contains the dict's value, default 2
    :return: a dict object which keys are entities in the first row, values are
    set of entities have relationships with the key.
    """
    assos = {}
    xcol -= 1
    ycol -= 1
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        for line in f:
            words = line.strip().split(sep=sep)
            if len(words) > xcol and len(words) > ycol:
                k = words[xcol].strip()
                if k not in assos:
                    assos[k] = set()
                # # ---------------------
                # if words[ycol].strip() in assos[words[xcol].strip()]:
                #     print("read_assos() duplicated associations:", words[ycol].strip(), words[xcol].strip())
                # # ---------------------
                assos[k].add(words[ycol].strip())
    return assos


def read_mappings(filepath, header=False, sep="\t", xcol=1, ycol=2):
    """
    read mappings (e.g. gene symbol-gene entrez id mapping) into key-value format
    from a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param xcol: number of the column contains the dict's key, default 1
    :param ycol: number of the column contains the dict's value, default 2
    :return: a dict which each key is a string and the corresponding value is a string
    """
    mapping = {}
    xcol -= 1
    ycol -= 1
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        for line in f:
            words = line.strip().split(sep=sep)
            if len(words) > xcol and len(words) > ycol:
                k = words[xcol].strip()
                if k not in mapping:
                    mapping[k] = words[ycol].strip()
                else:
                    # print("read_mappings() the key is duplicated:", words[xcol].strip())
                    if mapping[k] != words[ycol].strip():
                        print("read_mappings() fetal error, the key has different mapping values:",
                              "key:", k, "value:", mapping[k], words[ycol].strip())
    return mapping


def read_sims(filepath, header=False, sep="\t", xcol=1, ycol=2, vcol=3):
    """
    read triplet similarity into key-value format from a table file
    :param filepath: a path of a table file which contains at least 3 columns,
    2 cols are names and the other col is values
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param xcol: number of the column contains entities, default 1
    :param ycol: number of another column contains entities, default 2
    :param vcol: number of the column contains the similarity values, default 3
    :return: dict, key-value: string-dict<string-value> ({entity1: {entity2: sim, }, }
    """
    sim = {}
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        xcol -= 1
        ycol -= 1
        vcol -= 1
        for line in f:
            words = line.split(sep)
            entity1 = words[xcol].strip()
            entity2 = words[ycol].strip()
            simvalue = float(words[vcol].strip())
            if entity1 not in sim.keys():
                sim[entity1] = {}
            sim[entity1][entity2] = simvalue
            if entity1 != entity2 and entity2 in sim.keys() and entity1 in sim[entity2].keys():
                print(entity1, entity2, "are duplicate in the file.")
                del sim[entity2][entity1]
        ks = list(sim.keys())
        for k in ks:
            if len(sim[k]) == 0:
                del sim[k]
    return sim


def read_simmatrix(filepath, rowheader=True, colheader=True, sep="\t"):
    """
    read similarity matrix into key-value format from a table file, the table
    file must contains rownames (rowheader=True) or colnames (colheader=True)
    or both, if there is only one which is True, this method considers that
    the matrix defaults to a symmetric matrix
    :param filepath: a path of a table file which contains a sim matrix
    :param rowheader: a boolean variable indicates if the first col is a head or not,
    default True
    :param colheader: a boolean variable indicates if the first row is a head or not,
    default True
    :param sep: delimiter, a string, default "\t"
    :return: dict, key-value: string-dict<string-value> ({entity1: {entity2: sim, }, }
    """
    sim = {}
    if not rowheader | colheader:
        print("the matrix must have at least one header for row and col names")
        return None
    with open(filepath, mode='r') as f:
        rownames = []
        colnames = []
        starter = 0
        if rowheader:
            starter = 1
            rownames = read_one_col(filepath, 1, colheader, sep)
        if colheader:
            rns = next(f)
            colnames = rns.strip().split(sep)
        if not colheader & rowheader:
            if rowheader:
                colnames = rownames
            else:
                rownames = colnames
        rcount = 0
        for line in f:
            words = line.strip().split(sep)
            if rownames[rcount] not in sim.keys():
                sim[rownames[rcount]] = {}
            for i in range(0, len(colnames)):
                vtemp = float(words[i+starter])
                if colnames[i] in sim.keys() and rownames[rcount] in sim[colnames[i]].keys():
                    if vtemp != sim[colnames[i]][rownames[rcount]]:
                        print(rownames[rcount], colnames[i], "have two sim values.")
                else:
                    sim[rownames[rcount]][colnames[i]] = vtemp
            rcount += 1
    return sim


def write_mappings(dictionary, filepath, header=False, sep="\t"):
    """
    write dict where each key only has one value (int or str or else, not list) to a file
    :param dictionary: dict object
    :param filepath: file path of a file
    :param header: boolean, need a head or not, if True, the head will be "V1'sep'V2",
    default False
    :param sep: delimiter, default "\t"
    :return: none
    """
    with open(filepath, mode='w') as f:
        if header:
            f.write("V1"+sep+"V2\n")
        for k in dictionary.keys():
            f.write(str(k)+sep+str(dictionary[k])+"\n")
    print("write_mappings: writing finished.")


def write_assos(dictionary, filepath, header=False, sep="\t"):
    """
    write dict where each key's value is a set contains entities associated with the key
    to a file
    :param dictionary: dict object
    :param filepath: file path of a file
    :param header: boolean, need a head or not, if True, the head will be "V1'sep'V2",
    default False
    :param sep: delimiter, default "\t"
    :return: none
    """
    with open(filepath, mode='w') as f:
        if header:
            f.write("V1" + sep + "V2\n")
        for k in dictionary.keys():
            for v in dictionary[k]:
                f.write(str(k) + sep + str(v) + "\n")
    print("write_assos: writing finished.")


def write_sims(simdict, filepath, header=False, sep='\t'):
    """
    write dict contains similarities between each two entities (e.g. diseases)
    to a file
    :param simdict: dict, key-value (string-dict<string-float>):
    {entity1:{entity2: sim1, entity3: sim2,},}
    :param filepath: file path to a file
    :param header: boolean, need a head or not, if True, the head will be
    "V1sepV2sepsim",default False
    :param sep: delimiter, default '\t'
    :return: None
    """
    with open(filepath, mode='w') as f:
        if header:
            f.write("V1" + sep + "V2" + sep + "sim\n")
        for k1 in simdict.keys():
            for k2 in simdict[k1].keys():
                f.write(str(k1) + sep + str(k2) + sep + str(simdict[k1][k2]) + "\n")
    print("write_sims: writing finished.")


def write_simmatrix(simdict, filepath, asorder=False, dorder=None, sep='\t'):
    """
    write dict contains similarities between each two entities (e.g. diseases)
    to a file in matrix format
    :param simdict: dict, key-value (string-dict<string-float>):
    {entity1:{entity2: sim1, entity3: sim2,},}
    :param filepath: file path to a file
    "V1sepV2sepsim",default False
    :param asorder: if row name and col name of the output matrix as order
    :param dorder: list of ordered matrix dimnames, if asorder is True, this will work
    :param sep: delimiter, default '\t'
    :return: None
    """
    def findsimvalue(d1, d2, dsims):
        if d1 in dsims.keys() and d2 in dsims[d1].keys():
            return dsims[d1][d2]
        elif d2 in dsims.keys() and d1 in dsims[d2].keys():
            return dsims[d2][d1]
        else:
            if d1 == d2:
                return 1.0
            return 0.0

    dnames = set(simdict.keys())
    for d in simdict.keys():
        dnames.update(set(simdict[d].keys()))
    if asorder:
        dnames = dorder
    else:
        dnames = list(dnames)
    with open(filepath, mode='w') as f:
        for d in dnames:
            f.write(sep + d)
        f.write("\n")
        for da in dnames:
            f.write(da)
            for db in dnames:
                f.write(sep + str(findsimvalue(da, db, simdict)))
            f.write("\n")
        print("write_simmatrix: writing finished.")


def write_slist(slist, filepath, header=False):
    """
    write a list which each element is a string to a file.
    :param slist: a list (set is ok, may be) of strings
    :param filepath: file path to a file
    :param header: boolean, need a head or not, if True, the head will be "V1",
    default False
    :return: None
    """
    with open(filepath, mode='w') as f:
        if header:
            f.write("V1\n")
        for s in slist:
            f.write(s+"\n")
    print("write_slist: writing finished.")


def write_sorteddict(dictlist, filepath, header=False, sep="\t"):
    """
    write a list (converted from dict by method sorted()) which each element is a key-value tuple
    to a file.
    :param dictlist: sorted dict list
    :param filepath: file path to a file
    :param header: boolean, need a head or not, if True, the head will be "V1'sep'V2",
    default False
    :param sep: delimiter, default "\t"
    :return: None
    """
    with open(filepath, mode='w') as f:
        if header:
            f.write("V1"+sep+"V2\n")
        for t in dictlist:
            f.write(str(t[0])+sep+str(t[1])+"\n")
    print("write_sorteddict: writing finished.")


def stat_assos(assos):
    """
    print how many keys/values/associations in the dict assos
    :param assos: a dict (key-value: string-set<string>)
    :return: None
    """
    values = set()
    counta = 0
    for k in assos.keys():
        for v in assos[k]:
            values.add(v)
            counta += 1
    print("stat_assos: keys:", len(assos.keys()), "values:", len(values), "associations:", counta)


def stat_maps(maps):
    """
    print how many keys/values in the dict maps
    :param maps: a dcit (key-value: string-string)
    :return: None
    """
    values = set()
    for k in maps.keys():
        values.add(maps[k])
    print("stat_maps: keys:", len(maps.keys()), "values:", len(values))


def stat_sims(sims):
    """
    print some stats of sims
    :param sims: sims from method read_sims() or read_simmatrix()
    :return: None
    """
    nkeys = len(sims)
    entities = set()
    assos = 0
    for p in sims.keys():
        for q in sims[p].keys():
            entities.add(q)
            assos += 1
    print("stat_sims: keys:", nkeys, "entities:", len(entities), "assos:", assos)


def stat_network(assos):
    """
    print stats of node assos of network
    :param assos: a dict (key-value: string-set<string>), a key means a node
    and the key's value is a set of its associated nodes
    :return: None
    """
    nodes = set()
    edgenum = 0
    for p in assos.keys():
        nodes.add(p)
        for q in assos[p]:
            nodes.add(q)
            edgenum += 1
    print("stat_network: nodes:", len(nodes), "edges:", edgenum)


def invert_dict(xy_dict):
    """
    construct a new dict upon the inputed dict, the new dict's key are the inputed
    dict's values, values are sets of inputed dict's keys
    :param xy_dict: a dict which keys are strings, values are sets of entities associated
    with the key
    :return: a dict (yx_dict)
    """
    yx_dict = {}
    for x in xy_dict.keys():
        for y in xy_dict[x]:
            if y not in yx_dict.keys():
                yx_dict[y] = set()
            yx_dict[y].add(x)
    return yx_dict


def combine_two_assos(dict1, dict2):
    """
    combine two dicts, each dict's key-value is string-set<string>
    :param dict1: a dict
    :param dict2: another dict
    :return:
    """
    dict_all = copy.deepcopy(dict1)
    for k in dict2.keys():
        if k not in dict_all.keys():
            dict_all[k] = set()
        dict_all[k] |= dict2[k]
    return dict_all
