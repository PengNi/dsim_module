#! /usr/bin/env python3
"""calculate similarity of gene or disease pairs based on gene coexpression"""
from scipy import stats
from numpy import mean
import time


def read_probeid_expfile(filepath, header=True, sep='\t'):
    """
    read probeid expression file. in this case means "U133AGNF1B.gcrma.avg.cleared.tab"
    Reference: Su A I, Wiltshire T, Batalov S, et al. A gene atlas of the mouse and human
    protein-encoding transcriptomes[J]. Proceedings of the National Academy of Sciences of
    the United States of America, 2004, 101(16): 6062-6067.
    :param filepath: "U133AGNF1B.gcrma.avg.cleared.tab"
    :param header: True or False
    :param sep: delimiter
    :return: a dict (key-value: string-list<float>)
    """
    with open(filepath, mode='r') as f:
        if header:
            next(f)
        probe_exp = {}
        for line in f:
            words = line.split(sep)
            valvector = []
            for i in range(1, len(words)):
                valvector.append(float(words[i].strip()))
            probe_exp[words[0].strip()] = valvector
    return probe_exp


def probeexp2geneexp(probe_exp, gene2probe):
    """
    get genes' expression value vectors based on probes' expression value vectors and
    gene-probe mapping relationships.
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param probe_exp: dict (key-velue: string-list<float>)
    :param gene2probe: dict (key-value: string-set<string>)
    :return: dict (key-velue: string-list<float>)
    """
    gene_exp = {}
    for g in gene2probe.keys():
        max_mean = 0
        chosedp = ''
        for p in gene2probe[g]:
            if p in probe_exp.keys():
                tempmean = mean(probe_exp[p])
                if max_mean < tempmean:
                    max_mean = tempmean
                    chosedp = p
        if chosedp != '':
            gene_exp[g] = probe_exp[chosedp]
    return gene_exp


def diseases_similarity_coexp(diseases, disease2gene, gene2expression):
    """
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param diseases:
    :param disease2gene:
    :param gene2expression:
    :return: a dict (key-value: string-dict<string-float>)
    """
    diseases = list(set(diseases).intersection(set(disease2gene.keys())))
    disease2gene_new = {}
    for d in diseases:
        genestemp = disease2gene[d].intersection(gene2expression.keys())
        if len(genestemp) > 0:
            disease2gene_new[d] = genestemp
    diseases = list(disease2gene_new.keys())
    print("there are {} diseases can be calculated.".format(len(diseases)))

    n = len(diseases)
    sim_result = {}
    print("cal started at:", time.strftime('%X %x %Z'))
    t0 = time.time()
    for i in range(0, n):
        print(i, diseases[i], "gene number:", len(disease2gene_new[diseases[i]]))
        sim_result[diseases[i]] = {}
        for j in range(i, n):
            sim_result[diseases[i]][diseases[j]] = disease_pair_similarity_coexp(diseases[i],
                                                                                 diseases[j],
                                                                                 disease2gene_new,
                                                                                 gene2expression)
        print("------------------------------- cost time:", str(time.time() - t0))
    return sim_result


def gene_pair_similarity_coexp(genea, geneb, gene2expression):
    """
    Reference: Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through
    the incomplete interactome[J]. Science, 2015, 347(6224): 1257601.
    :param genea:
    :param geneb:
    :param gene2expression:
    :return:
    """
    return stats.stats.spearmanr(gene2expression[genea], gene2expression[geneb])


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
    genea = disease2gene[diseasea]
    geneb = disease2gene[diseaseb]

    count = 0
    simsum = 0.0
    for a in genea:
        for b in geneb:
            if a != b:
                gsim = gene_pair_similarity_coexp(a, b, gene2expression)
                simsum += abs(gsim[0])
                count += 1
    if count != 0:
        return simsum/count
    return 0
