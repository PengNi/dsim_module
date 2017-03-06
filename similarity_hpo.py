#! /usr/bin/env python3
"""calculate similarity of diseases based hpo"""
import human_phenotype_ontology


def termsimilarity_resnik(t1, t2, p2ancestors, pinforcontent):
    """

    :param t1:
    :param t2:
    :param p2ancestors:
    :param pinforcontent: must contain all terms in hpo
    :return:
    """
    if not (t1 in p2ancestors.keys() and t2 in p2ancestors.keys()):
        return 0.0
    else:
        comancestors = p2ancestors[t1].intersection(p2ancestors[t2])
        res = 0.0
        for c in comancestors:
            if pinforcontent[c] > res:
                res = pinforcontent[c]
        return res


def termssimilarity_resnik(termlist, hpo, p2anno):
    """

    :param termlist:
    :param hpo:
    :param p2anno:
    :return:
    """
    termlist = list(termlist)
    p2ancestors = human_phenotype_ontology.get_terms2ancestors(termlist, hpo)
    expp2anno = human_phenotype_ontology.hpo_expp2anno(p2anno, hpo)
    pinforcontent = human_phenotype_ontology.terms_information_content(hpo, expp2anno)
    sim = {}
    for i in range(0, len(termlist)):
        sim[termlist[i]] = {}
        for j in range(i, len(termlist)):
            sim[termlist[i]][termlist[j]] = termsimilarity_resnik(termlist[i], termlist[j],
                                                                  p2ancestors, pinforcontent)
    return sim


def termsimilarity_lin(t1, t2, p2ancestors, pinforcontent):
    """

    :param t1:
    :param t2:
    :param p2ancestors:
    :param pinforcontent: must contain all terms in hpo
    :return:
    """
    if not (t1 in p2ancestors.keys() and t2 in p2ancestors.keys()):
        return 0.0
    else:
        comancestors = p2ancestors[t1].intersection(p2ancestors[t2])
        res = 0.0
        for c in comancestors:
            if pinforcontent[c] > res:
                res = pinforcontent[c]
        if (pinforcontent[t1] + pinforcontent[t2]) > 0:
            res = 2 * res / (pinforcontent[t1] + pinforcontent[t2])
        else:
            res = 1
        return res


def termssimilarity_lin(termlist, hpo, p2anno):
    """

    :param termlist:
    :param hpo:
    :param p2anno:
    :return:
    """
    termlist = list(termlist)
    p2ancestors = human_phenotype_ontology.get_terms2ancestors(termlist, hpo)
    expp2anno = human_phenotype_ontology.hpo_expp2anno(p2anno, hpo)
    pinforcontent = human_phenotype_ontology.terms_information_content(hpo, expp2anno)
    sim = {}
    for i in range(0, len(termlist)):
        sim[termlist[i]] = {}
        for j in range(i, len(termlist)):
            sim[termlist[i]][termlist[j]] = termsimilarity_lin(termlist[i], termlist[j],
                                                               p2ancestors, pinforcontent)
    return sim


def diseasesimilarity_le(d2p, psim):
    """
    Le D H, Dang V T. Ontology-based disease similarity network for disease gene prediction[J].
    Vietnam Journal of Computer Science, 2016, 3(3): 197-205.
    :param d2p:
    :param psim:
    :return:
    """
    def findsimvalue(p1, p2, dsims):
        if p1 in dsims.keys() and p2 in dsims[p1].keys():
            return dsims[p1][p2]
        elif p2 in dsims.keys() and p1 in dsims[p2].keys():
            return dsims[p2][p1]
        else:
            if p1 == p2:
                return 1.0
            return 0.0
    sim = {}
    selfsim = {}
    diseases = list(d2p.keys())
    for i in range(0, len(diseases)):
        sim[diseases[i]] = {}
        for j in range(i, len(diseases)):
            simvalue = 0.0
            for pi in d2p[diseases[i]]:
                for pj in d2p[diseases[j]]:
                    simtemp = findsimvalue(pi, pj, psim)
                    if simtemp > simvalue:
                        simvalue = simtemp
            sim[diseases[i]][diseases[j]] = simvalue
    for d in diseases:
        selfsim[d] = sim[d][d]
    for d1 in sim.keys():
        for d2 in sim[d1].keys():
            sim[d1][d2] = 2 * sim[d1][d2] / (selfsim[d1] + selfsim[d2])
    return sim


def diseasesimilarity_bma(d2p, psim):
    """

    :param d2p:
    :param psim:
    :return:
    """
    def findsimvalue(p1, p2, dsims):
        if p1 in dsims.keys() and p2 in dsims[p1].keys():
            return dsims[p1][p2]
        elif p2 in dsims.keys() and p1 in dsims[p2].keys():
            return dsims[p2][p1]
        else:
            if p1 == p2:
                return 1.0
            return 0.0
    sim = {}
    diseases = list(d2p.keys())
    for i in range(0, len(diseases)):
        sim[diseases[i]] = {}
        for j in range(i, len(diseases)):
            simvalue = 0.0
            for pi in d2p[diseases[i]]:
                simptemp = 0.0
                for pj in d2p[diseases[j]]:
                    simtemp = findsimvalue(pi, pj, psim)
                    if simtemp > simptemp:
                        simptemp = simtemp
                simvalue += simptemp
            for pi in d2p[diseases[j]]:
                simptemp = 0.0
                for pj in d2p[diseases[i]]:
                    simtemp = findsimvalue(pi, pj, psim)
                    if simtemp > simptemp:
                        simptemp = simtemp
                simvalue += simptemp
            sim[diseases[i]][diseases[j]] = simvalue / (len(d2p[diseases[i]]) + len(d2p[diseases[j]]))
    return sim


def diseasesimilarity_funsimavg(d2p, psim):
    """

    :param d2p:
    :param psim:
    :return:
    """
    def findsimvalue(p1, p2, dsims):
        if p1 in dsims.keys() and p2 in dsims[p1].keys():
            return dsims[p1][p2]
        elif p2 in dsims.keys() and p1 in dsims[p2].keys():
            return dsims[p2][p1]
        else:
            if p1 == p2:
                return 1.0
            return 0.0
    sim = {}
    diseases = list(d2p.keys())
    for i in range(0, len(diseases)):
        sim[diseases[i]] = {}
        for j in range(i, len(diseases)):
            simvalue1, simvalue2 = 0.0, 0.0
            for pi in d2p[diseases[i]]:
                simptemp = 0.0
                for pj in d2p[diseases[j]]:
                    simtemp = findsimvalue(pi, pj, psim)
                    if simtemp > simptemp:
                        simptemp = simtemp
                simvalue1 += simptemp
            simvalue1 /= len(d2p[diseases[i]])
            for pi in d2p[diseases[j]]:
                simptemp = 0.0
                for pj in d2p[diseases[i]]:
                    simtemp = findsimvalue(pi, pj, psim)
                    if simtemp > simptemp:
                        simptemp = simtemp
                simvalue2 += simptemp
            simvalue2 /= len(d2p[diseases[j]])
            sim[diseases[i]][diseases[j]] = (simvalue1 + simvalue2) / 2
    return sim
