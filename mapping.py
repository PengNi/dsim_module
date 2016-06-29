#! /usr/bin/env python3
"""functions about ID mapping"""
import mygene


# entrezid to symbol is defaultly considered as one to one relationship
def entrezid2symbol(entrezid_list, species="human"):
    """
    convert entrezid to gene symbol
    :param entrezid_list: a list of entrezids
    :param species: species, default "huamn"
    :return: a dict which keys are entrezids (as string type) and values are gene symbols
    """
    mg = mygene.MyGeneInfo()
    result = mg.querymany(entrezid_list, scope="entrezgene",
                          fields="entrezgene,symbol", species=species)
    mapping = {}
    for r in result:
        if 'notfound' not in r.keys():
            if str(r['entrezgene']) in mapping.keys():
                print(str(r['entrezgene']), ": the key has duplicate values")
            mapping[str(r['entrezgene'])] = r['symbol']
    return mapping


# probesetid may have multiple entrezids
def probesetid2entrezid(probesetid_list, species="human"):
    """
    convert probesetid to entrezid
    :param probesetid_list: a list of probesetids
    :param species: species, default "human"
    :return: a dict which keys are probesetids (as string type) and values are sets of entrezids
    """
    mg = mygene.MyGeneInfo()
    result = mg.querymany(probesetid_list, scope="reporter",
                          fields="entrezgene", species=species)
    mapping = {}
    for r in result:
        if 'notfound' not in r.keys():
            if r['query'] not in mapping.keys():
                mapping[r['query']] = set()
            mapping[r['query']].add(r['entrezgene'])
    return mapping


def convert_dict_values(dictionary, mapping):
    """
    convert dict value ids upon a mapping dict
    :param dictionary: dict whose key-value is string-set<string>
    :param mapping: dict whose key-value is string-string
    :return:
    """
    newdict = {}
    for k in dictionary.keys():
        newdict[k] = set()
        for v in dictionary[k]:
            if v in mapping.keys():
                newdict[k].add(mapping[v])
    return newdict


def covert_dict_keys(dictionary, mapping):
    """
    convert dict value ids upon a mapping dict
    :param dictionary: dict whose key-value is string-whatever
    :param mapping: dict whose key-value is string-string
    :return:
    """
    newdict = {}
    for k in dictionary.keys():
        if k in mapping.keys():
            newdict[mapping[k]] = dictionary[k]
    return newdict
