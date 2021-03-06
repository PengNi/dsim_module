#! /usr/bin/env python3
"""defining a class for Disease Ontology"""
from queue import Queue


def get_terms_at_layern(n, do):
    """
    get all termids of terms at layer n
    :param do: DiseaseOntology
    :param n: integer
    :return: set()
    """
    terms = do.getterms()
    result = set()
    for t in terms.keys():
        layerstemp = terms[t].getlayers()
        if len(layerstemp) != 0 and min(layerstemp) == n:
            result.add(t)
    return result


def get_terms2offsprings(terms, do):
    """
    get terms-offsprings dict
    :param terms: set()
    :param do: DiseaseOntology
    :return: dict
    """
    result = {}
    for t in terms:
        result[t] = set()
        result[t].add(t)
        result[t].update(do.get_termsoffspring(t))
    return result


class DiseaseOntology:
    """Disease Ontology"""
    def __init__(self):
        self._terms = {}

    def getterms(self):
        return self._terms

    def readobofile(self, filepath):
        """
        read HumanDO.obo file
        :param filepath: path of HumanDO.obo file
        """
        with open(filepath, mode='r') as f:
            line = next(f)
            while line is not None:
                line = line.strip()
                if line == '[Term]':
                    doid = ""
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
                            doid = 'DOID:' + words[2].strip()
                        elif key == 'name':
                            name = words[1].strip()
                        elif key == 'alt_id':
                            altids.add('DOID:' + words[2].strip())
                        elif key == 'synonym':
                            synonyms.add(words[1].strip().split('"')[1].strip())
                        elif key == 'xref':
                            if len(words) < 3:
                                print("readobofile exception line:", line)
                                xrefs.add('OMIM:' + words[1].strip())  # for line = 'xref: 611644' exception
                            else:
                                xrefs.add(words[1].strip() + ':' + words[2].strip())
                        elif key == 'is_a':
                            parents.add('DOID:' + words[2].strip().split('!')[0].strip())
                        line = next(f).strip()
                    if doid in self._terms.keys():
                        self._terms[doid].setdoid(doid)
                        self._terms[doid].setname(name)
                        self._terms[doid].setaltids(altids)
                        self._terms[doid].setsynonyms(synonyms)
                        self._terms[doid].setxrefs(xrefs)
                        self._terms[doid].setparents(parents)
                    else:
                        doterm = DOTerm()
                        doterm.setdoid(doid)
                        doterm.setname(name)
                        doterm.setaltids(altids)
                        doterm.setsynonyms(synonyms)
                        doterm.setxrefs(xrefs)
                        doterm.setparents(parents)
                        self._terms[doid] = doterm
                    for p in parents:
                        if p not in self._terms.keys():
                            doterm = DOTerm()
                            doterm.setdoid(p)
                            self._terms[p] = doterm
                        self._terms[p].addchild(doid)
                try:
                    line = next(f)
                except StopIteration:
                    line = None
        print("initialize terms size:", len(self._terms))

    def settermslayers(self):
        for t in self._terms.keys():
            self._terms[t].setlayers(set())
        self._terms['DOID:4'].addlayer(1)
        q = Queue()
        q.put('DOID:4')
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


class DOTerm:
    """DO term"""
    def __init__(self):
        self._doid = ""
        self._name = ""
        self._altids = set()
        self._synonyms = set()
        self._xrefs = set()
        self._parents = set()
        self._children = set()
        self._layers = set()

    def setdoid(self, doid):
        self._doid = doid

    def getdoid(self):
        return self._doid

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
        print("id:", self.getdoid())
        print("name:", self.getname())
        print("layers:", self.getlayers())
        print("alt_ids:", self.getaltids())
        print("synonyms:", self.getsynonyms())
        print("xref:", self.getxrefs())
        print("parents:", self.getparents())
        print("children:", self.getchildren())
