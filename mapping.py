#! /usr/bin/env python3
"""functions about ID mapping"""
import mygene
from copy import deepcopy

fields = ['entrezgene', 'symbol', 'uniprot', 'reporter']


def doid2xref(regex, complete, do):
    """
    get doid2xref mapping
    :param regex: "ICD9CM", "ICD10CM", "MSH", "UMLS_CUI" "OMIM" etc.
    :param complete: True or False, get the complete id or not,
    e.g. "OMIM:613376" or "613376"
    :param do: DiseaseOntology
    :return: dict, string-set<string>
    """
    doterms = do.getterms()
    assos = {}
    for d in doterms.keys():
        xrefs = doterms[d].getxrefs()
        for x in xrefs:
            if str(x).startswith(regex):
                if d not in assos.keys():
                    assos[d] = set()
                assos[d].add(str(x))
    if complete:
        return assos
    else:
        assosn = {}
        for d in assos.keys():
            assosn[d] = set()
            for x in assos[d]:
                assosn[d].add(str(x).split(':')[1].strip())
        return assosn


def doid2umlsid(umlsdiseases):
    """
    get mapping between doids and umls ids
    :param umlsdiseases: an UmlsDiseases object
    :return: a dict, key-value: string(doid)-set<string>(umls id)
    """
    do2umls = {}
    umlsdis = umlsdiseases.getumlsdiseases()
    for u in umlsdis.keys():
        for d in umlsdis[u].getdoids():
            if d not in do2umls.keys():
                do2umls[d] = set()
            do2umls[d].add(u)
    return do2umls


def icd9cmid2umlsid_alldigit(umlsdiseases):
    """
    get mapping between icd9cm ids (3 digit ids, 4 digit ids,
    5 digit ids) and umls ids
    :param umlsdiseases: an UmlsDiseases object
    :return: a dict, key-value: string(icd9cmid)-set<string>(umls id)
    """
    icd92umls = {}
    umlsdis = umlsdiseases.getumlsdiseases()
    for u in umlsdis.keys():
        for i in umlsdis[u].geticd9cmids():
            icd9id = i[7:]
            if len(icd9id) < 7:
                if icd9id not in icd92umls.keys():
                    icd92umls[icd9id] = set()
                icd92umls[icd9id].add(umlsdis[u].getid())
    return icd92umls


def icd9cmid2umlsid_3digit(umlsdiseases):
    """
    get mapping between icd9cm 3 digit ids and umls ids
    :param umlsdiseases: an UmlsDiseases object
    :return: a dict, key-value: string(icd9cmid_3digit)-set<string>(umls id)
    """
    icd92umls = {}
    umlsdis = umlsdiseases.getumlsdiseases()
    for u in umlsdis.keys():
        for i in umlsdis[u].geticd9cmids():
            icd9id_3digit = i[7:10]
            if icd9id_3digit not in icd92umls.keys():
                icd92umls[icd9id_3digit] = set()
            icd92umls[icd9id_3digit].add(umlsdis[u].getid())
    return icd92umls


def termname2umlsid(meshnames, umlsdiseases, meshname2meshid):
    """
    get mapping between mesh term names and umls ids
    :param meshnames: set contains mesh names, may be contains all
    mesh names from "ncomms5212-s5.txt" file
    :param umlsdiseases: an UmlsDiseases object
    :param meshname2meshid: dict from method read_mappings("MeshTreeHierarchy.csv")
    :return: a dict, key-value: string(mesh term name)-set<string(umls id), >
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
                    meshname2umlsid[mn] = set()
                    meshname2umlsid[mn].add(umlsidtemp)
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


def symbol2uniprot(symbolid_list, species='human'):
    mg = mygene.MyGeneInfo()
    result = mg.querymany(symbolid_list, scope='symbol',
                          fields='uniprot', species=species)
    mapping = {}
    for r in result:
        if 'notfound' not in r.keys():
            if r['query'] not in mapping.keys():
                mapping[r['query']] = set()
            if 'uniprot' in r.keys() and 'Swiss-Prot' in r['uniprot'].keys():
                if isinstance(r['uniprot']['Swiss-Prot'], str):
                    mapping[r['query']].add(r['uniprot']['Swiss-Prot'])
                else:
                    mapping[r['query']].update(set(r['uniprot']['Swiss-Prot']))
    return mapping


def geneida2geneidb(geneida, geneidb, geneida_list, species='human'):
    """
    convert geneida to geneidb
    :param geneida_list: a list of geneidas
    :param geneida: type of geneida, string
    :param geneidb: type of geneidb, string
    :param species: species, default 'human'
    :return: a dict which keys are geneidas (string) and values are sets of geneidbs
    """
    mg = mygene.MyGeneInfo()
    result = mg.querymany(geneida_list, scope=geneida,
                          fields=geneidb, species=species)
    mapping = {}
    for r in result:
        if 'notfound' not in r.keys() and geneidb in r.keys():
            if r['query'] not in mapping.keys():
                mapping[r['query']] = set()
            mapping[r['query']].add(r[geneidb])
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
    :param mapping: dict whose key-value is string-string/string-set<string>
    :return:
    """
    newdict = {}
    for k in dictionary.keys():
        newdict[k] = set()
        for v in dictionary[k]:
            if v in mapping.keys():
                if isinstance(mapping[v], str):
                    newdict[k].add(mapping[v])
                else:
                    newdict[k].update(mapping[v])
    ks = list(newdict.keys())
    for k in ks:
        if len(newdict[k]) == 0:
            del newdict[k]
    return newdict


def convert_dict_keys(dictionary, mapping):
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
