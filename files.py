#! /usr/bin/env python3
"""functions about read from and write to files"""
import copy


def read_one_col(filepath, col, header=False, sep="\t"):
    """
    read a column from a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param col: the number of col needs to be extracted
    :return: a list contains elements from the extracted column
    """
    eles = []
    col -= 1
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        for line in f:
            eles.append(line.split(sep=sep)[col].strip())
    return eles


def read_assos(filepath, header=False, sep="\t", xcol=1, ycol=2):
    """
    read associations (e.g. disease-gene associations) into key-value format from
    a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param xcol: the column number contains the dict's key, default 1
    :param ycol: the column number contains the dict's value, default 2
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
            words = line.split(sep=sep)
            if words[xcol].strip() not in assos:
                assos[words[xcol].strip()] = set()
            # # ---------------------
            # if words[ycol].strip() in assos[words[xcol].strip()]:
            #     print("read_assos() duplicated associations:", words[ycol].strip(), words[xcol].strip())
            # # ---------------------
            assos[words[xcol].strip()].add(words[ycol].strip())
    return assos


def read_mappings(filepath, header=False, sep="\t", xcol=1, ycol=2):
    """
    read mappings (e.g. gene symbol-gene entrez id mapping) into key-value format
    from a table file
    :param filepath: a path of a table file which contains at least two columns
    :param header: a boolean variable indicates if the first row is a head or not,
    default False
    :param sep: delimiter, a string, default "\t"
    :param xcol: the column number contains the dict's key, default 1
    :param ycol: the column number contains the dict's value, default 2
    :return: a dict which each key is a string and the corresponding value is a string
    """
    mapping = {}
    xcol -= 1
    ycol -= 1
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        for line in f:
            words = line.split(sep=sep)
            if words[xcol].strip() not in mapping:
                mapping[words[xcol].strip()] = words[ycol].strip()
            else:
                # print("read_mappings() the key is duplicated:", words[xcol].strip())
                if mapping[words[xcol].strip()] != words[ycol].strip():
                    print("read_mappings() fetal error, the key has different mapping values:",
                          "key:", words[xcol].strip(), "value:", mapping[words[xcol].strip()],
                          words[ycol].strip())
    return mapping


def write_mappings(dictionary, filepath, header=False, sep="\t"):
    """
    write dict where each key only has one value (int or str or else, not list) to a file.
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
    to a file.
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


def write_slist(slist, filepath, header=False):
    """
    write a list  which each element is a string to a file.
    :param slist: a list of strings
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
