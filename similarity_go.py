#! /usr/bin/env python3
"""calculate similarity of gene or disease pairs based on GO"""


def diseases_similarity_go(diseases, disease2gene, gene2go, go2gene):
    """
    calculate simialrity of disease pairs based on GO.
    :param diseases:
    :param disease2gene:
    :param gene2go:
    :param go2gene:
    :return: a dict (key-value: string-dict<string-float>)
    """
    diseases = list(set(diseases).intersection(set(disease2gene.keys())))
    n = len(diseases)
    sim_result = {}
    for i in range(0, n):
        sim_result[diseases[i]] = {}
        for j in range(i, n):
            simvalue = disease_pair_similarity_go(diseases[i], diseases[j],
                                                  disease2gene, gene2go, go2gene)
            sim_result[diseases[i]][diseases[j]] = simvalue
    return sim_result


def gene_pair_similarity_go(genea, geneb, gene2go, go2gene):
    """
    calculate simialrity of a gene pair based on GO.
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param genea:
    :param geneb:
    :param gene2go:
    :param go2gene:
    :return: a float number or zero
    """
    sharedgos = gene2go[genea].intersection(gene2go[geneb])
    if len(sharedgos) == 0:
        return 0
    num = []
    for go in sharedgos:
        num.append(len(go2gene[go]))
    return 2/min(num)


def disease_pair_similarity_go(diseasea, diseaseb, disease2gene, gene2go, go2gene):
    """
    calculate simialrity of a disease pair based on GO.
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param diseasea:
    :param diseaseb:
    :param disease2gene:
    :param gene2go:
    :param go2gene:
    :return: a float number or zero
    """
    genes_a = disease2gene[diseasea].intersection(list(gene2go.keys()))
    genes_b = disease2gene[diseaseb].intersection(list(gene2go.keys()))

    if len(genes_a) == 0 or len(genes_b) == 0:
        return 0
    else:
        count = 0
        sim_sum = 0
        for a in genes_a:
            for b in genes_b:
                if a != b:
                    count += 1
                    sim_sum += gene_pair_similarity_go(a, b, gene2go, go2gene)
        if count != 0:
            return sim_sum/count
        return 0
