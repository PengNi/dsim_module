#! /usr/bin/env python3
"""class for umls diseases"""
import turtle_parser


def read_umls_disease_info(dgasso_score_cutoff):
    """generate umls disease information, inclding each disease's umls id, disease name,
    mesh classes, mesh ids, mim numbers, doids, icd9cm ids, hpo ids and associated genes
    (in two form: entrezid, gene symbol)
    @:param dgasso_score_cutoff: the cutoff of the disease-gene association score, only associations
    whose score no less than dgasso_score_cutoff will be read
    """
    umlsdiseases = UmlsDiseases()

    umlsdiseaseid2name = turtle_parser.ttl_parser_disease2name("data/disease.ttl")
    umlsdiseaseid2categories = turtle_parser.ttl_parser_disease2class("data/disease.ttl")
    umlsdiseaseid2meshid = turtle_parser.ttl_parser_umls2mesh_exact("data/ls-umls2mesh-exactMatch.ttl")
    umlsdiseaseid2omimid = turtle_parser.ttl_parser_umls2omim_exact("data/ls-umls2omim-exactMatch.ttl")
    umlsdiseaseid2doid = turtle_parser.ttl_parser_umls2do_exact("data/ls-umls2do-exactMatch.ttl")
    umlsdiseaseid2icd9cm = turtle_parser.ttl_parser_umls2icd9cm_exact("data/ls-umls2icd9cm-exactMatch.ttl")
    umlsdiseaseid2hpoid = turtle_parser.ttl_parser_umls2hpo_exact("data/ls-umls2hpo-exactMatch.ttl")

    umlsdiseaseid2entrezid = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
                                                                dgasso_score_cutoff)
    umlsdiseaseid2genesymbol = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
                                                                  dgasso_score_cutoff, True, False)

    for did in umlsdiseaseid2name.keys():
        if str(did).startswith("umls"):
            umlsdiseases.getumlsdiseases()[str(did)] = UmlsDisease()
            umlsdiseases.getumlsdiseases()[str(did)].setid(str(did))
            umlsdiseases.getumlsdiseases()[str(did)].setname(umlsdiseaseid2name[did])
    # print("number of umls disease ids who have names:", len(umlsdiseases.getumlsdiseases()))

    for did in umlsdiseaseid2categories.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2categories but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setcategories(umlsdiseaseid2categories[did])

    for did in umlsdiseaseid2meshid.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2meshid but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setmeshids(umlsdiseaseid2meshid[did])

    for did in umlsdiseaseid2omimid.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2omimid but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setomimids(umlsdiseaseid2omimid[did])

    for did in umlsdiseaseid2doid.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2doid but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setdoids(umlsdiseaseid2doid[did])

    for did in umlsdiseaseid2icd9cm.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2icd9cm but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].seticd9cmids(umlsdiseaseid2icd9cm[did])

    for did in umlsdiseaseid2hpoid.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2hpoid but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].sethpoids(umlsdiseaseid2hpoid[did])

    for did in umlsdiseaseid2entrezid.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2entrezid but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setentrezgenes(umlsdiseaseid2entrezid[did])

    for did in umlsdiseaseid2genesymbol.keys():
        if did not in umlsdiseases.getumlsdiseases().keys():
            print(did, "is in umlsdiseaseid2genesymbol but not in umlsdiseaseid2name")
        else:
            umlsdiseases.getumlsdiseases()[did].setgenesymbols(umlsdiseaseid2genesymbol[did])

    return umlsdiseases


# read "all_gene_disease_associations.tsv" file
def read_all_gene_disease_associations(filepath, score_cutoff, header=True, entrezgene=True):
    """
    read "all_gene_disease_associations.tsv" file to get disease-gene associations
    :param filepath: "all_gene_disease_associations.tsv"
    :param score_cutoff: the cutoff of the association score, only associations whose score no less
    than score_cutoff will be read
    :param header: a boolean variable indicates if the first row is a head or not,
    default True
    :param entrezgene: True or False, True means the associations between disease ids and entrezgene
    ids will be read, False means the associations between disease ids and gene symbols will be read
    :return: a dict contains disease-gene associations, key-value is string-set<string>
    """
    assos = {}
    genecol = 0
    if not entrezgene:
        genecol = 1
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        for line in f:
            words = line.split(sep="\t")
            disease = words[3].strip()
            if float(words[5].strip()) >= score_cutoff:
                if disease not in assos.keys():
                    assos[disease] = set()
                assos[disease].add(words[genecol].strip())
    return assos


class UmlsDiseases:
    def __init__(self):
        """
        umlsdiseases is a dict, key-value is umlsid-UmlsDisease object
        """
        self._umlsdiseases = {}

    def getumlsdiseases(self):
        return self._umlsdiseases


class UmlsDisease:
    def __init__(self):
        self._id = ""
        self._name = ""
        self._categories = set()
        self._meshids = set()
        self._omimids = set()
        self._doids = set()
        self._icd9cmids = set()
        self._hpoids = set()
        self._gene_entrezids = set()
        self._gene_symbols = set()

    def getid(self):
        return self._id

    def setid(self, idnumber):
        self._id = idnumber

    def getname(self):
        return self._name

    def setname(self, name):
        self._name = name

    def getcategories(self):
        return self._categories

    def setcategories(self, categories):
        self._categories = categories

    def getmeshids(self):
        return self._meshids

    def setmeshids(self, meshids):
        self._meshids = meshids

    def getomimids(self):
        return self._omimids

    def setomimids(self, omimids):
        self._omimids = omimids

    def getdoids(self):
        return self._doids

    def setdoids(self, doids):
        self._doids = doids

    def geticd9cmids(self):
        return self._icd9cmids

    def seticd9cmids(self, icd9cmids):
        self._icd9cmids = icd9cmids

    def gethpoids(self):
        return self._hpoids

    def sethpoids(self, hpoids):
        self._hpoids = hpoids

    def getentrezgenes(self):
        return self._gene_entrezids

    def setentrezgenes(self, entrezgenes):
        self._gene_entrezids = entrezgenes

    def getgenesymbols(self):
        return self._gene_symbols

    def setgenesymbols(self, genesymbols):
        self._gene_symbols = genesymbols

    def print(self):
        print("[umls disease]")
        print("id:", self._id)
        print("name:", self._name)
        print("mesh class:", self._categories)
        print("mesh ids:", self._meshids)
        print("omim ids:", self._omimids)
        print("doids:", self._doids)
        print("icd9cm ids:", self._icd9cmids)
        print("hpo ids: ", self._hpoids)
        print("entrezgenes:", self._gene_entrezids)
        print("gene symbols:", self._gene_symbols)
