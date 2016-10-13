#! /usr/bin/env python3
"""defining a class for Disease Ontology"""
from queue import Queue


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
    get terms-offsprings dict
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
