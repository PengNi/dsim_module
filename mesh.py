#! /usr/bin/env python3
"""class for MeSH terms"""
from queue import Queue


def mesh_ancestorscontri(calterms, meshobj, delta=0.5):
    meshterms = meshobj.getterms()
    terms_cancal = set(calterms).intersection(set(meshterms.keys()))

    ancestorscontribution = {}
    for term in terms_cancal:
        ancestorscontribution[term] = {}
        ancestorscontribution[term][term] = 1.0
        q = Queue()
        for c in meshterms[term].getparents():
            q.put(c)
        while not q.empty():
            nextcal = q.get()
            for c in meshterms[nextcal].getparents():
                q.put(c)
            nextcal_contri = 0.0
            for nextcal_child in meshterms[nextcal].getchildren():
                if nextcal_child in ancestorscontribution[term].keys():
                    if nextcal_contri < ancestorscontribution[term][nextcal_child]:
                        nextcal_contri = ancestorscontribution[term][nextcal_child]
            nextcal_contri *= delta
            if nextcal not in ancestorscontribution[term].keys():
                ancestorscontribution[term][nextcal] = nextcal_contri
            elif nextcal_contri > ancestorscontribution[term][nextcal]:
                ancestorscontribution[term][nextcal] = nextcal_contri
    return ancestorscontribution


def dvas(calterms, ancestorscontribution, meshobj):
    meshterms = meshobj.getterms()
    terms_cancal = set(calterms).intersection(set(meshterms.keys()))
    dva = {}
    for term in terms_cancal:
        dvaterm = 0.0
        for t in ancestorscontribution[term].keys():
            dvaterm += ancestorscontribution[term][t]
        dva[term] = dvaterm
    return dva


def misim_cal(calterms, meshobj, delta=0.5):
    """

    :param calterms:
    :param meshobj:
    :param delta:
    :return:
    """
    meshterms = meshobj.getterms()
    terms_cancal = set(calterms).intersection(set(meshterms.keys()))
    print('there are', len(terms_cancal), 'diseases can be calculated.')

    ancestorscontribution = {}
    for term in terms_cancal:
        ancestorscontribution[term] = {}
        ancestorscontribution[term][term] = 1.0
        q = Queue()
        for c in meshterms[term].getparents():
            q.put(c)
        while not q.empty():
            nextcal = q.get()
            for c in meshterms[nextcal].getparents():
                q.put(c)
            nextcal_contri = 0.0
            for nextcal_child in meshterms[nextcal].getchildren():
                if nextcal_child in ancestorscontribution[term].keys():
                    if nextcal_contri < ancestorscontribution[term][nextcal_child]:
                        nextcal_contri = ancestorscontribution[term][nextcal_child]
            nextcal_contri *= delta
            if nextcal not in ancestorscontribution[term].keys():
                ancestorscontribution[term][nextcal] = nextcal_contri
            elif nextcal_contri > ancestorscontribution[term][nextcal]:
                ancestorscontribution[term][nextcal] = nextcal_contri

    dva = {}
    for term in terms_cancal:
        dvaterm = 0.0
        for t in ancestorscontribution[term].keys():
            dvaterm += ancestorscontribution[term][t]
        dva[term] = dvaterm

    sims = {}
    termlist = list(terms_cancal)
    for i in range(0, len(termlist)):
        sims[termlist[i]] = {}
        for j in range(i, len(termlist)):
            common_ancestors = \
                set(ancestorscontribution[termlist[i]].keys()).intersection(set(ancestorscontribution[termlist[j]].keys()))
            numerator = 0.0
            for ancestor in common_ancestors:
                numerator += (ancestorscontribution[termlist[i]][ancestor] + ancestorscontribution[termlist[j]][ancestor])
            sims[termlist[i]][termlist[j]] = numerator / (dva[termlist[i]] + dva[termlist[j]])
    print("misim: finished")
    return sims


class MeSH:
    def __init__(self):
        self._terms = {}

    def getterms(self):
        return self._terms

    def read_meshtreehierarchy_csv(self, filepath, header=True):
        """
        read "MeshTreeHierarchy.csv" file to get each mesh term's heading,
        tree number and unique id
        :param filepath: "MeshTreeHierarchy.csv"
        :param header: if there is a head or not, default True
        :return: None
        """
        with open(filepath, mode='r', encoding="utf-8") as f:
            if header:
                next(f)
            for line in f:
                words = line.split("\t")
                if len(words) < 2:
                    print(words)
                uniqueid = words[1].strip()
                if uniqueid not in self._terms.keys():
                    tempterm = MeSHTerms()
                    tempterm.setheading(words[2].strip())
                    tempterm.setuniqueid(uniqueid)
                    tempterm.addtreenumber(words[0].strip())
                    self._terms[uniqueid] = tempterm
                else:
                    if self._terms[uniqueid].getheading() != words[2].strip():
                        print(" read_meshtreehierarchy_csv() fetal error: the id has different headings:",
                              uniqueid, self._terms[uniqueid].getheading(), words[2].strip())
                    else:
                        self._terms[uniqueid].addtreenumber(words[0].strip())
        print("MeSH read_meshtreehierarchy_csv finished!")

    def read_bin_file(self, filepath):
        """
        read binary file of mesh descriptor, e.g. "d2016.bin"
        :param filepath:
        :return:
        """
        pass

    def setchildrenandparents(self):
        treenumber2id = {}
        for m in self._terms.keys():
            for tnum in self._terms[m].gettreenumbers():
                if tnum in treenumber2id.keys():
                    print(tnum, 'wrong')
                treenumber2id[tnum] = m
        for m in self._terms.keys():
            for tnum in self._terms[m].gettreenumbers():
                tnumtmp = str(tnum).split('.')
                if len(tnumtmp) > 1:
                    parenttmp = tnumtmp[0]
                    for i in range(1, len(tnumtmp) - 1):
                        parenttmp += '.' + tnumtmp[i]
                    self._terms[m].addparent(treenumber2id[parenttmp])
                    self._terms[treenumber2id[parenttmp]].addchild(m)

    def getancetorsbyid(self, termid):
        ancestors = set()
        if termid in self.getterms().keys():
            ancestors.add(termid)
            q = Queue()
            for c in self._terms[termid].getparents():
                q.put(c)
            while not q.empty():
                nextc = q.get()
                for ncc in self._terms[nextc].getparents():
                    q.put(ncc)
                ancestors.add(nextc)
        return ancestors

    def getancestors(self, termids):
        """
        termids is set() or list()
        """
        t2a = {}
        for t in termids:
            tas = self.getancetorsbyid(t)
            if len(tas) > 0:
                t2a[t] = tas
        return t2a


class MeSHTerms:
    def __init__(self):
        self._heading = ""
        self._tree_numbers = set()
        self._unique_id = ""
        self._children = set()
        self._parents = set()

    def getheading(self):
        return self._heading

    def setheading(self, heading):
        self._heading = heading

    def gettreenumbers(self):
        return self._tree_numbers

    def settreenumbers(self, treenumbers):
        self._tree_numbers = treenumbers

    def addtreenumber(self, treenumber):
        self._tree_numbers.add(treenumber)

    def getuniqueid(self):
        return self._unique_id

    def setuniqueid(self, uniqueid):
        self._unique_id = uniqueid

    def getchildren(self):
        return self._children

    def setchildren(self, children):
        self._children = children

    def addchild(self, child):
        self._children.add(child)

    def getparents(self):
        return self._parents

    def setparents(self, parents):
        self._parents = parents

    def addparent(self, parent):
        self._parents.add(parent)

    def print(self):
        print("[mesh term]")
        print("heading:", self._heading)
        print("unique id:", self._unique_id)
        print("tree numbers:", self._tree_numbers)
        print("parents:", self._parents)
        print("children:", self._children)
