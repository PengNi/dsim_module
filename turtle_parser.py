#! /usr/bin/env python3


# parse "ls-umls2hpo-exactMatch.ttl" flie
def ttl_parser_umls2hpo_exact(filepath):
    """
        read "ls-umls2hpo-exactMatch.ttl" file to get the mapping of disease umls id
        to hpo id
        :param filepath: "ls-umls2hpo-exactMatch.ttl"
        :return: a dict, key-value is string-set<string>, not string-string because
        a umls id may be mapped to many hpo ids
        """
    umls2hpo = {}
    with open(filepath, mode='r') as f:
        umls = ""
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("<http://"):
                umls = "umls:" + stripline.split("/")[6].strip().split(">")[0].strip()
                umls2hpo[umls] = set()
            if stripline.startswith("skos:exactMatch hpo"):
                hpos = stripline.split(",")
                umls2hpo[umls].add(hpos[0].strip().split(" ")[1].strip())
                length = len(hpos)
                if length >= 2:
                    for i in range(1, length - 1):
                        umls2hpo[umls].add(hpos[i].strip())
                    umls2hpo[umls].add(hpos[length - 1].strip().split(".")[0].strip())
            line = f.readline()
    return umls2hpo


# parse "ls-umls2omim-exactMatch.ttl" flie
def ttl_parser_umls2omim_exact(filepath):
    """
    read "ls-umls2omim-exactMatch.ttl" file to get the mapping of disease umls id
    to mim number
    :param filepath: "ls-umls2omim-exactMatch.ttl"
    :return: a dict, key-value is string-set<string>, not string-string because
    a umls id may be mapped to many mim numbers
    """
    umls2omim = {}
    with open(filepath, mode='r') as f:
        umls = ""
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("<http://"):
                umls = "umls:" + stripline.split("/")[6].strip().split(">")[0].strip()
                umls2omim[umls] = set()
            if stripline.startswith("skos:exactMatch "):
                mimnumbers = stripline.split(",")
                for mimnumber in mimnumbers:
                    umls2omim[umls].add("omim:" + mimnumber.strip().split("/")[4].strip().split(">")[0].strip())
            line = f.readline()
    return umls2omim


# parse "ls-umls2icd9cm-exactMatch.ttl" flie
def ttl_parser_umls2icd9cm_exact(filepath):
    """
    read "ls-umls2icd9cm-exactMatch.ttl" file to get the mapping of disease umls id
    to icd9cm id
    :param filepath: "ls-umls2icd9cm-exactMatch.ttl"
    :return: a dict, key-value is string-set<string>, not string-string because
    a umls id may be mapped to many icd9cm ids
    """
    umls2icd9cm = {}
    with open(filepath, mode='r') as f:
        umls = ""
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("<http://"):
                umls = "umls:" + stripline.split("/")[6].strip().split(">")[0].strip()
                umls2icd9cm[umls] = set()
            if stripline.startswith("skos:exactMatch "):
                icd9cmids = stripline.split(",")
                for icd9cmid in icd9cmids:
                    umls2icd9cm[umls].add("icd9cm:"+icd9cmid.strip().split("/")[5].strip().split(">")[0].strip())
            line = f.readline()
    return umls2icd9cm


# parse "ls-umls2do-exactMatch.ttl" flie
def ttl_parser_umls2do_exact(filepath):
    """
    read "ls-umls2do-exactMatch.ttl" file to get the mapping of disease umls id to doid
    :param filepath: "ls-umls2do-exactMatch.ttl"
    :return: a dict, key-value is string-set<string>, not string-string because
    a umls id may be mapped to many do ids
    """
    umls2do = {}
    with open(filepath, mode='r') as f:
        umls = ""
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("<http://"):
                umls = "umls:" + stripline.split("/")[6].strip().split(">")[0].strip()
                umls2do[umls] = set()
            if stripline.startswith("skos:exactMatch "):
                doids = stripline.split(",")
                for doid in doids:
                    umls2do[umls].add(doid.strip().split("/")[4].strip().split(">")[0].strip().lower())
            line = f.readline()
    return umls2do


# parse "ls-umls2mesh-exactMatch.ttl" file
def ttl_parser_umls2mesh_exact(filepath):
    """
    read "ls-umls2mesh-exactMatch.ttl" file to get the mapping of disease umls id to mesh id
    :param filepath: "ls-umls2mesh-exactMatch.ttl"
    :return: a dict, key-value is string-set<string>, not string-string because
    a umls id may be mapped to many mesh ids
    """
    umls2mesh = {}
    with open(filepath, mode='r') as f:
        umls = ""
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("<http://"):
                umls = "umls:" + stripline.split("/")[6].strip().split(">")[0].strip()
                umls2mesh[umls] = set()
            if stripline.startswith("skos:exactMatch mesh"):
                meshes = stripline.split(",")
                umls2mesh[umls].add(meshes[0].strip().split(" ")[1].strip())
                length = len(meshes)
                if length >= 2:
                    for i in range(1, length-1):
                        umls2mesh[umls].add(meshes[i].strip())
                    umls2mesh[umls].add(meshes[length-1].strip().split(".")[0].strip())
            line = f.readline()
    return umls2mesh


# parse "disease.ttl" file
def ttl_parser_disease2name(filepath):
    """
    read "disease.ttl" file to get each disease id and its
    corresponding name/title
    :param filepath: "disease.ttl"
    :return: a dict, key-value is string-string
    """
    diseaseid2name = {}
    disease = ""
    with open(filepath, mode='r') as f:
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("dcterms:identifier \""):
                disease = stripline.split("\"")[1].strip()
                # if disease in diseaseid2name.keys():
                #     print(disease)
            if stripline.startswith("dcterms:title \""):
                diseaseid2name[disease] = stripline.split("\"")[1].strip()
            line = f.readline()
    return diseaseid2name


# parse "disease.ttl" file
def ttl_parser_disease2class(filepath):
    """
    read "disease.ttl" file to get each disease id and classes it belongs to
    :param filepath: "disease.ttl"
    :return: a dict, key-value is string-set<string>
    """
    disease_class = {}
    disease = ""
    with open(filepath, mode='r') as f:
        line = f.readline()
        while line:
            stripline = line.strip()
            if stripline.startswith("dcterms:identifier \""):
                disease = stripline.split("\"")[1].strip()
                if disease.startswith("umls"):
                    disease_class[disease] = set()
            if stripline.startswith("sio:SIO_000095 ") and disease.startswith("umls"):
                categories = stripline.split(",")
                for c in categories:
                    category = c.strip().split("#")[1].strip().split(">")[0]
                    disease_class[disease].add(category)
            line = f.readline()
    return disease_class
