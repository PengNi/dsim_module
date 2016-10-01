#! /usr/bin/env python3
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
from files import stat_assos


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


def get_disease_pathway_assos(dgassos, pgassos, mtmethod='fdr_bh', pvcutoff=0.05):
    """
    based on disease-gene and pathway-gene associations, using hypergeometric test to get
    disease-pathway associations
    :param dgassos: dict, string-set<string>, disease-gene associations
    :param pgassos: dict, string-set<string>, pathway-gene associations
    :param mtmethod: multiple tests p-value correlation method,
    `bonferroni` : one-step correction
    `sidak` : one-step correction
    `holm-sidak` : step down method using Sidak adjustments
    `holm` : step-down method using Bonferroni adjustments
    `simes-hochberg` : step-up method  (independent)
    `hommel` : closed method based on Simes tests (non-negative)
    `fdr_bh` : Benjamini/Hochberg  (non-negative)
    `fdr_by` : Benjamini/Yekutieli (negative)
    `fdr_tsbh` : two stage fdr correction (non-negative)
    `fdr_tsbky` : two stage fdr correction (non-negative)
    :param pvcutoff, p-value cutoff to tell if a test is significant or not
    :return: dict, string-set<string>, disease-pathway associations
    """
    print("disease-gene assos: ", end='')
    stat_assos(dgassos)
    print("pathway-gene assos: ", end='')
    stat_assos(pgassos)
    pgenes = set()
    for p in pgassos.keys():
        pgenes |= pgassos[p]
    # -----------------------------------------------
    ds = list(dgassos.keys())
    for d in ds:
        dgassos[d] = pgenes.intersection(dgassos[d])
        if len(dgassos[d]) == 0:
            del dgassos[d]
    print("filtered disease-gene assos: ", end='')
    stat_assos(dgassos)
    # -----------------------------------------------
    pgslen = len(pgenes)
    ds = list(dgassos.keys())
    ps = list(pgassos.keys())
    dpassos = {}
    for d in ds:
        dglen = len(dgassos[d])
        pvalue = []
        ps_cal = []
        for p in ps:
            dpglen = len(dgassos[d].intersection(pgassos[p]))
            if dpglen > 0:
                ps_cal.append(p)
                pglen = len(pgassos[p])
                pvalue.append(1 - hypergeom.cdf(dpglen-1, pgslen, pglen, dglen))
        if len(ps_cal) > 0:
            qvalue_result = multipletests(pvalue, alpha=pvcutoff, method=mtmethod,
                                          is_sorted=False, returnsorted=False)
            dpassos[d] = set()
            for i in range(0, len(ps_cal)):
                if qvalue_result[0][i]:
                    dpassos[d].add(ps_cal[i])
            if len(dpassos[d]) == 0:
                del dpassos[d]
    print("similarity_pathway.get_disease_pathway_assos(): done..")
    return dpassos
