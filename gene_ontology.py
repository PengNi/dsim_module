#! /usr/bin/env python3
"""defining a class for Gene Ontology"""
import queue

namespaces = ("biological_process", "molecular_function", "cellular_component")
go_annotation_evidence_codes = ("EXP", "IDA", "IMP", "IGI", "IEP", "ISS", "ISA", "ISM", "ISO")
goid_old2new = {"GO:0097481": "GO:0014069", }


def read_go_annotation_file(filepath, hcec, qualifier):
    """
    read go annotation file
    :param filepath: a path of go annotation file
    :param hcec: boolean virable, True means only annotations with high confidence
    evidence code are used, False means annotations with all evidence code are used
    :param qualifier: boolean virable, True means only annotations with empty
    "qualifier" column are used, False means all annotations are used
    :return: a dict object which keys are GO terms, values are set of gene symbols
    annotated the key.
    """
    assos = {}
    with open(filepath, mode='r') as f:
        for line in f:
            if not line.startswith("!"):
                words = line.split("\t")
                if hcec and words[6].strip() not in go_annotation_evidence_codes:
                    continue
                if qualifier and words[3].strip() != "":
                    continue
                if words[4].strip() not in assos:
                    assos[words[4].strip()] = set()
                assos[words[4].strip()].add(words[2].strip())
        # -------------------------------------------change some alt_ids
        for goid_old in goid_old2new.keys():
            if goid_old in assos.keys():
                if goid_old2new[goid_old] in assos.keys():
                    assos[goid_old2new[goid_old]] |= assos[goid_old]
                else:
                    assos[goid_old2new[goid_old]] = assos[goid_old].copy()
                assos.pop(goid_old, None)
        # -------------------------------------------
    return assos


def add_implicit_annotations(goannotationdict, geneontology):
    """
    add all implicit general annotations by up-propagating the given annotations
    along the full GO tree
    :param goannotationdict: dict where key are gene symbols and values are sets of
    GO terms associated with the key, get from go annotation file
    :param geneontology: a GeneOntology object
    :return: a dict
    """
    assos = {}
    for k in goannotationdict.keys():
        assos[k] = set()
        q = queue.Queue()
        for go in goannotationdict[k]:
            q.put(go)
        while not q.empty():
            temp = q.get()
            assos[k].add(temp)
            if temp in geneontology.getterms():
                for t in geneontology.getterms()[temp].getparents():
                    q.put(t)
            else:
                print(temp, " is not in gene ontology for some reasons")
    return assos


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


class GeneOntology:
    """Gene Ontology"""
    def __init__(self):
        self._terms = {}

    def getterms(self):
        return self._terms

    def readobofile(self, filepath):
        """
        read go.obo file
        :param filepath: path of go.obo file
        """
        with open(filepath, mode='r') as f:
            temp = f.read().splitlines()
            i = 0
            templen = len(temp)
            while i < templen:
                if temp[i] == "[Term]":
                    i += 1
                    goid = ""
                    name = ""
                    namespace = ""
                    parents = set()
                    synonyms = set()
                    while temp[i] != "":
                        info = temp[i].strip().split(":")
                        key = info[0].strip()
                        if key == "id":
                            goid = "GO:"+info[2].strip()
                            i += 1
                            continue
                        if key == "name":
                            name = info[1].strip()
                            i += 1
                            continue
                        if key == "namespace":
                            namespace = info[1].strip()
                            i += 1
                            continue
                        if key == "is_a":
                            parentid = info[2].strip().split("!")
                            parents.add("GO:"+parentid[0].strip())
                            i += 1
                            continue
                        if key == "synonym":
                            synonyms.add(info[1].strip())
                            i += 1
                            continue
                        i += 1
                    if goid in self._terms:
                        self._terms[goid].setname(name)
                        self._terms[goid].setnamespace(namespace)
                        self._terms[goid].setparents(parents)
                        self._terms[goid].setsynonyms(synonyms)
                    else:
                        termtemp = GOTerm()
                        termtemp.setgoid(goid)
                        termtemp.setname(name)
                        termtemp.setnamespace(namespace)
                        termtemp.setparents(parents)
                        termtemp.setsynonyms(synonyms)
                        self._terms[goid] = termtemp
                    for p in parents:
                        if p not in self._terms:
                            termtemp = GOTerm()
                            termtemp.setgoid(p)
                            self._terms[p] = termtemp
                        self._terms[p].addchild(goid)
                i += 1
        print("initialize terms size: ", len(self._terms))


class GOTerm:
    """GO term"""
    def __init__(self):
        self._goid = ""
        self._name = ""
        self._namespace = ""
        self._synonyms = set()
        self._parents = set()
        self._children = set()

    def setgoid(self, goid):
        self._goid = goid

    def getgoid(self):
        return self._goid

    def setname(self, name):
        self._name = name

    def getname(self):
        return self._name

    def setnamespace(self, namespace):
        if namespace in namespaces:
            self._namespace = namespace
        else:
            print("namespace is not right!")

    def getnamespace(self):
        return self._namespace

    def getsynonyms(self):
        return self._synonyms

    def setsynonyms(self, synonyms):
        self._synonyms = synonyms

    def addsynonym(self, synonym):
        if not isinstance(synonym, str):
            raise ValueError("Invalid synonym '{}'".format(synonym))
        self._synonyms.add(synonym)

    def getchildren(self):
        return self._children

    def setchildren(self, children):
        self._children = children

    def addchild(self, child):
        if not isinstance(child, str):
            raise ValueError("Invalid child '{}'".format(child))
        self._children.add(child)

    def getparents(self):
        return self._parents

    def setparents(self, parents):
        self._parents = parents

    def addparent(self, parent):
        if not isinstance(parent, str):
            raise ValueError("Invalid parent '{}'".format(parent))
        self._parents.add(parent)

    def print(self):
        print("[Term]")
        print("id:", self.getgoid())
        print("name:", self.getname())
        print("namespace:", self.getnamespace())
        print("synonyms:", self.getsynonyms())
        print("parents:", self.getparents())
        print("children:", self.getchildren())
