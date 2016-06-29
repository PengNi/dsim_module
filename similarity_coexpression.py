#! /usr/bin/env python3
"""calculate similarity of gene or disease pairs based on gene coexpression"""
from scipy import stats


def gene_pair_similarity_coexp(genea, geneb, gene2expression):
    """
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param genea:
    :param geneb:
    :param gene2expression:
    :return:
    """
    if genea in gene2expression.keys() and geneb in gene2expression.keys():
        return stats.stats.spearmanr(gene2expression[genea], gene2expression[geneb])
    else:
        return 0, 1


def disease_pair_similarity_coexp(diseasea, diseaseb, disease2gene, gene2expression):
    """
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param diseasea:
    :param diseaseb:
    :param disease2gene:
    :param gene2expression:
    :return:
    """
    pass
