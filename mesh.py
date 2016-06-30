#! /usr/bin/env python3
"""class for MeSH terms"""


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


class MeSHTerms:
    def __init__(self):
        self._heading = ""
        self._tree_numbers = set()
        self._unique_id = ""

    def getheading(self):
        return self._heading

    def setheading(self, heading):
        self._heading = heading

    def gettreenumbers(self):
        return self._tree_numbers

    def settreenumber(self, treenumbers):
        self._tree_numbers = treenumbers

    def addtreenumber(self, treenumber):
        self._tree_numbers.add(treenumber)

    def getuniqueid(self):
        return self._unique_id

    def setuniqueid(self, uniqueid):
        self._unique_id = uniqueid

    def print(self):
        print("[mesh term]")
        print("heading:", self._heading)
        print("unique id:", self._unique_id)
        print("tree numbers:", self._tree_numbers)
