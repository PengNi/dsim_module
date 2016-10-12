#! /usr/bin/env python3
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
from files import stat_assos
from copy import deepcopy


def diseases_similarity_pathway_jaccard(diseases, disease2pathway):
    """
    calculate simialrity of disease pairs based on disease-pathway
    associations, use jaccard index formula.
    :param diseases:
    :param disease2pathway:
    :return:
    """
    diseases = list(set(diseases).intersection(set(disease2pathway.keys())))
    sims = {}
    n = len(diseases)
    for i in range(0, n):
        sims[diseases[i]] = {}
        pathwaysi = disease2pathway[diseases[i]]
        for j in range(i, n):
            pathwaysj = disease2pathway[diseases[j]]
            sims[diseases[i]][diseases[j]] = (len(pathwaysi.intersection(pathwaysj)) /
                                              len(pathwaysi.union(pathwaysj)))
    return sims


def combine_pathway_data(orifilepaths):
    """
    paths of pathway file get from gsea
    :param orifilepaths: paths of pathway file get from gsea
    :return a dict, string-set<string>
    """
    pgassos = {}
    for path in orifilepaths:
        with open(path, mode='r') as rf:
            for line in rf:
                words = line.strip().split('\t')
                paname = words[0].strip()
                genes = set([x.strip() for x in words[2:]])

                alin = False
                for pa in pgassos.keys():
                    if genes == pgassos[pa]:
                        print(pa, pgassos[pa])
                        print(paname, genes)
                        alin = True
                        break

                if alin is False:
                    pgassos[paname] = genes
        print(path, "done..")
    return pgassos
