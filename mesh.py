#! /usr/bin/env python3
"""class for MeSH terms"""


class MeSH:
    def __init__(self):
        self._terms = {}

    def getterms(self):
        return self._terms

    def read_meshtreehierarchy_csv(self, filepath, header=True):
        """
        read MeshTreeHierarchy.csv file to get each mesh term's heading,
        tree number and unique id
        :param filepath: "MeshTreeHierarchy.csv"
        :param header: if there is a head or not, default True
        :return: None
        """
        with open(filepath, mode='r', encoding="utf-8") as f:
            if header:
                next(f)
            for line in f:
                tempterm = MeSHTerms()
                words = line.split("\t")
                tempterm.setheading(words[2].strip())
                tempterm.setuniqueid(words[1].strip())
                tempterm.settreenumber(words[0].strip())
                self._terms[words[1].strip()] = tempterm
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
        self._tree_number = ""
        self._unique_id = ""

    def getheading(self):
        return self._heading

    def setheading(self, heading):
        self._heading = heading

    def gettreenumber(self):
        return self._tree_number

    def settreenumber(self, treenumber):
        self._tree_number = treenumber

    def getuniqueid(self):
        return self._unique_id

    def setuniqueid(self, uniqueid):
        self._unique_id = uniqueid

    def print(self):
        print("[mesh term]")
        print("heading:", self._heading)
        print("unique id:", self._unique_id)
        print("tree number:", self._tree_number)
