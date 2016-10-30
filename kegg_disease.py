#! /usr/bin/env python3

from openpyxl import load_workbook
import re


class KeggDiseases:
    def __init__(self):
        self._kdterms = {}

    def getterms(self):
        return self._kdterms

    def read_xls(self, filepath):
        wb = load_workbook(filepath)
        ws = wb.get_sheet_by_name(wb.get_sheet_names()[0])
        wsstr = []
        for row in ws:
            rowstr = []
            for cell in row:
                strtemp = str(cell.value).replace(u'\xa0', '\t')
                rowstr.append(strtemp)
            wsstr.append(rowstr)

        pathwayp = re.compile(r'hsa\d\d\d\d\d')
        oldname = ['ICD-10:', 'MeSH:', 'OMIM:', 'ISD-10:', 'MedlinePlus:']
        covdbname = {'ICD-10:': 'icd10cm:',
                     'MeSH:': 'mesh:',
                     'OMIM:': 'omim:',
                     'ISD-10:': 'icd10cm:',
                     'MedlinePlus:': 'MedlinePlus:', }

        def finddbname(text):
            for oname in oldname:
                if str(text).endswith(oname):
                    return covdbname[oname]
            return None

        for row in wsstr:
            kd = KeggDisease()
            for i in range(0, len(row) - 1, 2):
                if row[i] == 'Entry':
                    kd.setname(str(row[i + 1]).split('\t')[0])
                elif row[i] == 'Pathway':
                    pathways = pathwayp.findall(str(row[i + 1]))
                    kd.setpathways(set(pathways))
                elif row[i] == 'Other DBs':
                    words = str(row[i + 1]).strip().split('\t')
                    xrefs = set()
                    if len(words) == 1:
                        print(words)
                    else:
                        for c in range(1, len(words)):
                            prefix = finddbname(words[c - 1].strip())
                            if prefix is not None:
                                strtemp = words[c].strip()
                                for name in oldname:
                                    strtemp = strtemp.split(name)[0]
                                tids = strtemp.strip().split(' ')
                                for tid in tids:
                                    xrefs.add(prefix + tid)
                            else:
                                print('-----------', words)
                    kd.setxrefs(xrefs)
            self._kdterms[kd.getname()] = kd


class KeggDisease:

    def __init__(self):
        self._id = ""
        self._name = ""
        self._category = ""
        self._pathways = set()
        self._xrefs = set()

    def getid(self):
        return self._id

    def setid(self, kid):
        self._id = kid

    def getname(self):
        return self._name

    def setname(self, name):
        self._name = name

    def getcategory(self):
        return self._category

    def setcategory(self, category):
        self._category = category

    def getpathways(self):
        return self._pathways

    def setpathways(self, pathways):
        self._pathways = pathways

    def addpathway(self, pathway):
        self._pathways.add(pathway)

    def getxrefs(self):
        return self._xrefs

    def setxrefs(self, xrefs):
        self._xrefs = xrefs

    def addxref(self, xref):
        self._xrefs.add(xref)


def getmesh2pathways(keggdis):
    if not isinstance(keggdis, KeggDiseases):
        raise ValueError("Invalid parent '{}'".format(keggdis))
    kdterms = keggdis.getterms()
    mesh2pathways = {}
    for kd in kdterms.keys():
        kdpathways = kdterms[kd].getpathways()
        kdxrefs = kdterms[kd].getxrefs()
        if len(kdpathways) != 0:
            for xref in kdxrefs:
                if str(xref).startswith('mesh:'):
                    if xref not in mesh2pathways.keys():
                        mesh2pathways[xref] = set()
                    mesh2pathways[xref].update(kdpathways)
    return mesh2pathways


# if __name__ == '__main__':
#     kds = KeggDiseases()
#     kds.read_xls("data/kegg_disease/aaa.xlsx")
#     print(len(kds.getterms()))
#     mesh2pathway = getmesh2pathways(kds)
