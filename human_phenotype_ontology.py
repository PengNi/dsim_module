#! /usr/bin/env python3
"""defining a class for Human Phenotype Ontology"""
from queue import Queue
import math

hpo_annotation_evidence_codes = ("EXP", "IDA", "IMP", "IGI", "IEP", "ISS", "ISA", "ISM", "ISO")


def terms_information_content(hpo, expp2anno):
    """

    :param hpo:
    :param expp2anno:
    :return:
    """
    termsic = {}
    hpoterms = hpo.getterms()
    allanno = set()
    for p in expp2anno.keys():
        allanno.update(expp2anno[p])
    allannonum = len(allanno)
    for t in hpoterms.keys():
        if t in expp2anno.keys():
            termsic[t] = -math.log2(len(expp2anno[t])/allannonum)
        else:
            termsic[t] = 0.0
    return termsic


def hpo_expp2anno(p2anno, hpo):
    """

    :param p2anno:
    :param hpo:
    :return:
    """
    hpoterms = hpo.getterms()
    expp2anno = {}
    for p in hpoterms.keys():
        poffspring = hpo.get_termsoffspring(p)
        poffspring.add(p)
        expp2anno[p] = set()
        for ospring in poffspring:
            if ospring in p2anno.keys():
                expp2anno[p].update(p2anno[ospring])
    expps = list(expp2anno.keys())
    for expp in expps:
        if len(expp2anno[expp]) == 0:
            del expp2anno[expp]
    return expp2anno


def read_hpo_annotation_file(filepath, db='OMIM', hcec=False):
    """
    "phenotype_annotation.tab"
    :param filepath: file path
    :param db: annotation db name, in the first column of the file
    :param hcec: boolean virable, True means only annotations with high confidence
    evidence code are used, False means annotations with all evidence code are used,
    7th column of the file
    :return:
    """
    p2anno = {}
    with open(filepath, mode='r', encoding='utf-8') as rf:
        for line in rf:
            words = line.strip().split('\t')
            if words[0] != db:
                continue
            if hcec and words[6].strip() not in hpo_annotation_evidence_codes:
                continue
            if words[4].strip() not in p2anno.keys():
                p2anno[words[4].strip()] = set()
            p2anno[words[4].strip()].add(words[1].strip())
    return p2anno


def get_terms_at_layern(n, hpo):
    """
    get all termids of terms at layer n
    :param hpo: HumanPhenotypeOntology
    :param n: integer
    :return: set()
    """
    terms = hpo.getterms()
    result = set()
    for t in terms.keys():
        layerstemp = terms[t].getlayers()
        if len(layerstemp) != 0 and min(layerstemp) == n:
            result.add(t)
    return result


def get_terms2offsprings(terms, hpo):
    """
    get terms-offsprings dict, offsprings include the term itself
    :param terms: set()
    :param hpo: HumanPhenotypeOntology
    :return: dict
    """
    result = {}
    for t in terms:
        result[t] = set()
        result[t].add(t)
        result[t].update(hpo.get_termsoffspring(t))
    return result


def get_terms2ancestors(terms, hpo):
    """

    :param terms:
    :param hpo:
    :return:
    """
    result = {}
    for t in terms:
        result[t] = set()
        result[t].add(t)
        result[t].update(hpo.get_termsancestors(t))
    return result


class HumanPhenotypeOntology:
    """Human Phenotype Ontology"""
    def __init__(self):
        self._terms = {}

    def getterms(self):
        return self._terms

    def readobofile(self, filepath):
        """
        read hp.obo file
        :param filepath: path of hp.obo file
        """
        with open(filepath, mode='r') as f:
            line = next(f)
            while line is not None:
                line = line.strip()
                if line == '[Term]':
                    hpoid = ""
                    name = ""
                    altids = set()
                    synonyms = set()
                    xrefs = set()
                    parents = set()
                    line = next(f).strip()
                    while line != '':
                        words = line.split(':')
                        key = words[0].strip()
                        if key == 'id':
                            hpoid = 'HP:' + words[2].strip()
                        elif key == 'name':
                            name = words[1].strip()
                        elif key == 'alt_id':
                            altids.add('HP:' + words[2].strip())
                        elif key == 'synonym':
                            synonyms.add(words[1].strip().split('"')[1].strip())
                        elif key == 'xref':
                            xrefs.add(words[1].strip() + ':' + words[2].strip().split("\"")[0].strip())
                        elif key == 'is_a':
                            parents.add('HP:' + words[2].strip().split('!')[0].strip())
                        line = next(f).strip()
                    if hpoid in self._terms.keys():
                        self._terms[hpoid].sethpoid(hpoid)
                        self._terms[hpoid].setname(name)
                        self._terms[hpoid].setaltids(altids)
                        self._terms[hpoid].setsynonyms(synonyms)
                        self._terms[hpoid].setxrefs(xrefs)
                        self._terms[hpoid].setparents(parents)
                    else:
                        hpoterm = HPOTerm()
                        hpoterm.sethpoid(hpoid)
                        hpoterm.setname(name)
                        hpoterm.setaltids(altids)
                        hpoterm.setsynonyms(synonyms)
                        hpoterm.setxrefs(xrefs)
                        hpoterm.setparents(parents)
                        self._terms[hpoid] = hpoterm
                    for p in parents:
                        if p not in self._terms.keys():
                            hpoterm = HPOTerm()
                            hpoterm.sethpoid(p)
                            self._terms[p] = hpoterm
                        self._terms[p].addchild(hpoid)
                try:
                    line = next(f)
                except StopIteration:
                    line = None
        print("initialize terms size:", len(self._terms))

    def settermslayers(self):
        for t in self._terms.keys():
            self._terms[t].setlayers(set())
        self._terms['HP:0000001'].addlayer(1)
        q = Queue()
        q.put('HP:0000001')
        while not q.empty():
            parentid = q.get()
            parentlayers = self._terms[parentid].getlayers()
            for c in self._terms[parentid].getchildren():
                ctemp = self._terms[c]
                for pl in parentlayers:
                    ctemp.addlayer(pl+1)
                q.put(c)

    def get_termsoffspring(self, termid):
        offspring = set()
        q = Queue()
        for c in self._terms[termid].getchildren():
            q.put(c)
        while not q.empty():
            nextc = q.get()
            for ncc in self._terms[nextc].getchildren():
                q.put(ncc)
            offspring.add(nextc)
        return offspring

    def get_termsancestors(self, termid):
        ancestors = set()
        q = Queue()
        for c in self._terms[termid].getparents():
            q.put(c)
        while not q.empty():
            nextc = q.get()
            for ncc in self._terms[nextc].getparents():
                q.put(ncc)
            ancestors.add(nextc)
        return ancestors


class HPOTerm:
    """HPO term"""
    def __init__(self):
        self._hpoid = ""
        self._name = ""
        self._altids = set()
        self._synonyms = set()
        self._xrefs = set()
        self._parents = set()
        self._children = set()
        self._layers = set()

    def sethpoid(self, hpoid):
        self._hpoid = hpoid

    def gethpoid(self):
        return self._hpoid

    def setname(self, name):
        self._name = name

    def getname(self):
        return self._name

    def getaltids(self):
        return self._altids

    def setaltids(self, altids):
        self._altids = altids

    def addaltid(self, altid):
        if not isinstance(altid, str):
            raise ValueError("Invalid synonym '{}'".format(altid))
        self._altids.add(altid)

    def getsynonyms(self):
        return self._synonyms

    def setsynonyms(self, synonyms):
        self._synonyms = synonyms

    def addsynonym(self, synonym):
        if not isinstance(synonym, str):
            raise ValueError("Invalid synonym '{}'".format(synonym))
        self._synonyms.add(synonym)

    def getxrefs(self):
        return self._xrefs

    def setxrefs(self, xrefs):
        self._xrefs = xrefs

    def addxref(self, xref):
        if not isinstance(xref, str):
            raise ValueError("Invalid synonym '{}'".format(xref))
        self._xrefs.add(xref)

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

    def getlayers(self):
        return self._layers

    def setlayers(self, layers):
        self._layers = layers

    def addlayer(self, layer):
        self._layers.add(layer)

    def print(self):
        print("[Term]")
        print("id:", self.gethpoid())
        print("name:", self.getname())
        print("layers:", self.getlayers())
        print("alt_ids:", self.getaltids())
        print("synonyms:", self.getsynonyms())
        print("xref:", self.getxrefs())
        print("parents:", self.getparents())
        print("children:", self.getchildren())
