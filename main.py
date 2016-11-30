#! /usr/bin/env python3
from copy import deepcopy
from pprint import pprint
import similarity_module
import mapping
import common_use
import experiments
import similarity_genesim
from files import read_one_col
from files import read_mappings
from files import read_assos
from files import read_sims
from files import read_simmatrix
from files import stat_assos
from files import stat_maps
from files import stat_sims
from files import stat_network
from files import write_assos
from files import write_mappings
from files import write_sims
from files import write_simmatrix
from files import write_slist
from files import invert_dict
from umls_disease import read_all_gene_disease_associations
from umls_disease import read_umls_disease_info
from similarity_coexpression import read_probeid_expfile
from similarity_coexpression import probeexp2geneexp
from similarity_coexpression import diseases_similarity_coexp
from similarity_go import diseases_similarity_go
from similarity_pathway import combine_pathway_data
from similarity_pathway import diseases_similarity_pathway_jaccard
from gene_ontology import read_go_annotation_file
from gene_ontology import GeneOntology
from gene_ontology import add_implicit_annotations
from evaluation import eva_readsims
from evaluation import eva_70benchmarkpairs
from evaluation import eva_ranking
from evaluation import eva_rankings
from evaluation import eva_cal_ranks
from evaluation import eva_tprfprs
from evaluation import eva_auc
from evaluation import eva_aucs
from evaluation import eva_groundtruth
from evaluation import eva_roc_getvalidationpairs
from evaluation import eva_get_validationpairsnlabels_comorbidity
from evaluation import eva_test_rankstats_multi
from evaluation import eva_test_pair_rankinfos
from evaluation import eva_test_pair_rankinfo_ranking
from disease_ontology import DiseaseOntology
from disease_ontology import get_terms_at_layern
from disease_ontology import get_terms2offsprings
from human_phenotype_ontology import HumanPhenotypeOntology
from plots import plot_roc
from plots import plot_simshist

namespaces = ("biological_process", "molecular_function", "cellular_component")

evaluation_simfilepaths1 = ['outputs/similarity_funsim_rwrsidd.tsv',
                            'outputs/similarity_hamaneh_rwrsidd_hppinwsl.tsv',
                            'outputs/similarity_experiments_rwr_rwrsidd_hppinwsl.tsv',
                            'outputs/similarity_suntopo_rwrsidd_hppinwosl_triplet.tsv',
                            # 'outputs/similarity_icod_rwrsidd_hppinwsl_triplet.tsv',
                            # 'outputs/similarity_bognew_rwrsidd_triplet.tsv',
                            # 'outputs/similarity_module5_rwrsidd_hppinwsl.tsv',
                            'outputs/similarity_spavgn_trans_rwrsidd_hppinwsl.tsv',
                            # 'outputs/similarity_spmaxexp_rwrsidd_hppinwsl.tsv',
                            # 'outputs/similarity_spmaxn_trans_rwrsidd_hppinwsl.tsv',
                            # 'outputs/similarity_dpathwayEuclideanSpmin_allpathway_rwrsidd_hppinwsl.tsv',
                            # 'outputs/similarity_dpathwayCosSpmaxidf_bkrpathway_rwrsidd_hppinwsl.tsv',
                            # 'outputs/similarity_pathway_jaccard_cp.bkr.v5.1.symbols_rwrsidd_bh005.tsv',
                            # 'outputs/similarity_genefun_rwrsidd_wangbmabp.tsv',
                            # 'outputs/similarity_spavgngenefun_rwrsidd_hppinwsl_wangbmpbp.tsv',
                            ]

shortnames1 = {'outputs/similarity_funsim_rwrsidd.tsv': 'FunSim',
               'outputs/similarity_hamaneh_rwrsidd_hppinwsl.tsv': 'Hamaneh',
               'outputs/similarity_experiments_rwr_rwrsidd_hppinwsl.tsv': 'NetSim',
               'outputs/similarity_suntopo_rwrsidd_hppinwosl_triplet.tsv': 'Sun_topo',
               'outputs/similarity_icod_rwrsidd_hppinwsl_triplet.tsv': 'ICod',
               'outputs/similarity_bognew_rwrsidd_triplet.tsv': 'BOG',
               'outputs/similarity_spavgn_trans_rwrsidd_hppinwsl.tsv': 'ModuleSim',
               'outputs/similarity_module5_rwrsidd_hppinwsl.tsv': 'module5',
               'outputs/similarity_spmaxexp_rwrsidd_hppinwsl.tsv': 'spmax',
               'outputs/similarity_spmaxn_trans_rwrsidd_hppinwsl.tsv': 'spmaxn',
               'outputs/similarity_pathway_jaccard_cp.bkr.v5.1.symbols_rwrsidd_bh005.tsv': 'pathwayj',
               'outputs/similarity_dpathwayCosineSpmax_allpathway_rwrsidd_hppinwsl.tsv': 'pCosSpmaxall',
               'outputs/similarity_dpathwayCosineSpmin_allpathway_rwrsidd_hppinwsl.tsv': 'pCosSpminall',
               'outputs/similarity_dpathwayEuclideanSpmin_allpathway_rwrsidd_hppinwsl.tsv': 'pEucSpminall',
               'outputs/similarity_genefun_rwrsidd_wangbmabp.tsv': 'genefun',
               'outputs/similarity_spavgngenefun_rwrsidd_hppinwsl_wangbmpbp.tsv': 'spavgngenefun',
               }

gtpathlist1 = ['outputs/similarity_pathway_jaccard_cp.bkr.v5.1.symbols_rwrsidd_bh005.tsv',
               'outputs/similarity_domain_jaccard_di2do_sidd.tsv',
               'outputs/similarity_pathway_jaccard_doid2keggpathway.tsv',
               ]

evaluation_simfilepaths2 = [  # str('outputs/similarity_experiments_shortestpath_transformed_less_disgenet_dgcut' +
                              #  'off006_interactome.tsv'),
                            # 'outputs/similarity_suntopo_disgenet_dgcutoff006_interactomemaxcc_triplet.tsv',
                            # 'outputs/similarity_funsim_disgenet_dgcutoff006.tsv',
                            # 'outputs/similarity_bog_disgenet_dgcutoff006_triplet.tsv',
                            # 'outputs/similarity_icod_disgenet_dgcutoff006_interactome_triplet.tsv',
                            # 'outputs/similarity_hamaneh_interactomenumls_dgcuff006.tsv',
                            'outputs/similarity_experiments_rwr_disgenet_dgcutoff006_interactome.tsv',
                            # 'outputs/similarity_module5_disgenet_dgcutoff006_interactome.tsv',
                            'outputs/similarity_spavgn_trans_disgenet_dgcutoff006_interactome.tsv',
                            # 'outputs/similarity_katz5_disgenet_dgcutoff006_interactome_beta05.tsv',
                            # 'outputs/similarity_katz5_disgenet_dgcutoff006_interactome_betadef.tsv',
                            # 'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_beta05.tsv',
                            # 'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_beta025.tsv',
                            'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_betadef.tsv',
                            # 'outputs/similarity_spmaxn_trans_disgenet_dgcutoff006_interactome.tsv',
                            # 'outputs/similarity_genefun_disgenet006_wangbmabp.tsv',
                            # 'outputs/similarity_pathway_jaccard_cp.bkr.v5.1.entrez_disgenet_dgcutoff006_bh005.tsv',
                            ]

shortnames2 = {'outputs/similarity_suntopo_disgenet_dgcutoff006_interactomemaxcc_triplet.tsv': 'Sun_topo',
               'outputs/similarity_funsim_disgenet_dgcutoff006.tsv': 'FunSim',
               'outputs/similarity_bog_disgenet_dgcutoff006_triplet.tsv': 'BOG',
               'outputs/similarity_icod_disgenet_dgcutoff006_interactome_triplet.tsv': 'ICod',
               'outputs/similarity_hamaneh_interactomenumls_dgcuff006.tsv': 'Hamaneh',
               'outputs/similarity_experiments_rwr_disgenet_dgcutoff006_interactome.tsv': 'NetSim',
               'outputs/similarity_experiments_shortestpath_transformed_less_disgenet_'
               'dgcutoff006_interactome.tsv': 'spmax',
               'outputs/similarity_module5_disgenet_dgcutoff006_interactome.tsv': 'module',
               'outputs/similarity_spavgn_trans_disgenet_dgcutoff006_interactome.tsv': 'ModuleSim',
               'outputs/similarity_pathway_jaccard_cp.bkr.v5.1.entrez_disgenet_dgcutoff006_bh005.tsv': 'pathwayj',
               'outputs/similarity_spmaxn_trans_disgenet_dgcutoff006_interactome.tsv': 'spmaxn',
               'outputs/similarity_genefun_disgenet006_wangbmabp.tsv': 'genefun',
               'outputs/similarity_katz5_disgenet_dgcutoff006_interactome_beta05.tsv': 'katz5_05',
               'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_beta05.tsv': 'katz4_05',
               'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_beta025.tsv': 'katz4_025',
               'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_betadef.tsv': 'katz4_def',
               }

gtpathlist2 = ['outputs/similarity_go_bp_umls_dgcutoff006.tsv',
               'outputs/similarity_go_cc_umls_dgcutoff006.tsv',
               'outputs/similarity_go_mf_umls_dgcutoff006.tsv',
               'outputs/similarity_coexp_umls_dgcutoff006.tsv',
               'outputs/similarity_symptom_1445umlsid.tsv',
               'outputs/similarity_symptom_1736umlsid.tsv',
               'outputs/similarity_pathway_jaccard_cp.bkr.v5.1.entrez_disgenet_dgcutoff006_bh005.tsv',
               'outputs/similarity_domain_jaccard_di2do_disgenet_dgcutoff006.tsv',
               ]

evapath_omim = ['D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\similarity_rwr_dgomim_hprd_matrix.tsv',
                'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\similarity_spavgn_dgomim_hprd_matrix.tsv',
                'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\disease_disease_similarity.tsv', ]

shnames_omim = {'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\similarity_rwr_dgomim_hprd_matrix.tsv': 'NetSim',
                'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\similarity_spavgn_dgomim_hprd_matrix.tsv': 'spavgn',
                'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\disease_disease_similarity.tsv': 'MimMiner', }

evaluation_simfilepaths5 = ['outputs/similarity_funsim_rwrsidd.tsv',
                            # 'outputs/similarity_hamaneh_rwrsidd_hppinwsl.tsv',
                            'outputs/similarity_rwr_rwrsidd_interactome.tsv',
                            'outputs/similarity_suntopo_rwrsidd_interactomemaxcc_triplet.tsv',
                            'outputs/similarity_spavgn_rwrsidd_interactome.tsv',
                            ]

shortnames5 = {'outputs/similarity_funsim_rwrsidd.tsv': 'FunSim',
               # 'outputs/similarity_hamaneh_rwrsidd_hppinwsl.tsv',
               'outputs/similarity_rwr_rwrsidd_interactome.tsv': 'NetSim',
               'outputs/similarity_suntopo_rwrsidd_interactomemaxcc_triplet.tsv': 'Sun_topo',
               'outputs/similarity_spavgn_rwrsidd_interactome.tsv': 'ModuleSim',
               }

evaluation_simfilepaths6 = ['outputs/similarity_spavgn_disgenet_dgcutoff000_interactome.tsv',
                            ]

shortnames6 = {'outputs/similarity_spavgn_disgenet_dgcutoff000_interactome.tsv': 'ModuleSim',
               }

evaluation_simfilepaths7 = ['outputs/similarity_funsim_disgenet_curated.tsv',
                            'outputs/similarity_spavgn_disgenet_curated_interactome.tsv',
                            'outputs/similarity_suntopo_disgenet_curated_interactomemaxcc_triplet.tsv',
                            ]

shortnames7 = {'outputs/similarity_funsim_disgenet_curated.tsv': 'FunSim',
               'outputs/similarity_spavgn_disgenet_curated_interactome.tsv': 'ModuleSim',
               'outputs/similarity_suntopo_disgenet_curated_interactomemaxcc_triplet.tsv': 'Sun_topo',
               }


# ---evaluation-----------------------------------------
def evaluation_validationpairs(simpathlist, shortnames, times=1):
    # vpairsnlabels = eva_get_validationpairsnlabels_comorbidity("data/comorbidity/AllNet3_umlsid.net")
    # vpairsnlabels = eva_get_validationpairsnlabels_comorbidity("data/comorbidity/AllNet3_doid.net")
    vpairsnlabels = eva_get_validationpairsnlabels_comorbidity("data/comorbidity/AllNet5.net")

    msims = eva_readsims(simpathlist)
    print("------------scores and lable---------------------------------------------------")
    vpairsinfos = eva_roc_getvalidationpairs(msims, vpairsnlabels, times)
    t = 0
    # for vpairsinfo in vpairsinfos:
    #     if t < 1:
    #         d1 = list(vpairsinfo.keys())[0]
    #         d2 = list(vpairsinfo[d1].keys())[0]
    #         methodnames = list(vpairsinfo[d1][d2].keys())
    #         methodnames.remove('label')
    #         print("d1\td2\tlabel", end='')
    #         for m in methodnames:
    #             print("\t" + shortnames[m], end='')
    #         print()
    #         for d1 in vpairsinfo.keys():
    #             for d2 in vpairsinfo[d1].keys():
    #                 print(d1 + '\t' + d2 + '\t' + str(vpairsinfo[d1][d2]['label']), end='')
    #                 for m in methodnames:
    #                     print("\t" + str(vpairsinfo[d1][d2][m]), end='')
    #                 print()
    #         t += 1
    print("------------pair ranks---------------------------------------------------------")
    vpairs_ranks = eva_rankings(vpairsinfos)
    # t = 0
    # for vpairs_rank in vpairs_ranks:
    #     if t < 1:
    #         ms = list(vpairs_rank.keys())
    #         for m in ms:
    #             print(str(shortnames[m]) + "\t\t\t", end="")
    #         print()
    #         llen = len(vpairs_rank[ms[0]])
    #         for i in range(0, llen):
    #             for m in ms:
    #                 for j in range(0, 4):
    #                     print(str(vpairs_rank[m][i][j]) + "\t", end="")
    #             print()
    #         t += 1
    print("------------tpr and fpr--------------------------------------------------------")
    vpairs_tprfprs = eva_tprfprs(vpairsinfos)
    t = 0
    for vpairs_tprfpr in vpairs_tprfprs:
        if t < 1:
            methodnames = list(vpairs_tprfpr.keys())
            llens = []
            for m in methodnames:
                llens.append(len(vpairs_tprfpr[m]))
                print(str(shortnames[m]) + "_fpr\t" + m + "_tpr\t", end='')
            print()
            maxlen = max(llens)
            for i in range(0, maxlen):
                for m in methodnames:
                    if i < len(vpairs_tprfpr[m]):
                        for j in range(0, 2):
                            print(str(vpairs_tprfpr[m][i][j]) + "\t", end="")
                    else:
                        print('\t\t', end='')
                print()
            auctemp = eva_auc(vpairs_tprfpr)
            plot_roc(vpairs_tprfpr, shortnames, auctemp)
            t += 1
    print("------------auc value----------------------------------------------------------")
    vpairs_aucs = eva_aucs(vpairs_tprfprs)
    auclen = len(vpairs_aucs)
    print("replicated times:", auclen)
    mns = list(vpairs_aucs[0].keys())
    avgaucs = []
    for mn in mns:
        avgauc = 0.0
        for auc in vpairs_aucs:
            avgauc += auc[mn]
        avgaucs.append((str(mn), avgauc / auclen))
    avgaucs = sorted(avgaucs, key=lambda a: a[1], reverse=True)
    for x in avgaucs:
        print(str(shortnames[x[0]]) + "\t" + str(x[1]))


def evaluation_groundtruth(simpathlist, shortnames, gtpaths):

    msims = eva_readsims(simpathlist)

    for g in range(0, len(gtpaths)):
        evagtres = eva_groundtruth(msims, gtpaths[g])
        print(gtpaths[g]+"----------------------------")
        mnames = list(evagtres.keys())
        sumtemp = {}
        for i in range(0, len(mnames)):
            if i == 0:
                print("rank"+"\t"+shortnames[mnames[i]], end='')
            else:
                print("\t"+shortnames[mnames[i]], end='')
            sumtemp[mnames[i]] = 0.0
        print()

        section = []
        for i in range(0, 10):
            section.append(10**i)
        printpoints = []
        for i in range(0, len(section)-1):
            if section[i] < 100:
                istep = section[i]
            else:
                istep = 100
            for j in range(section[i], section[i+1], istep):
                printpoints.append(j)

        index = 0
        for i in range(0, len(evagtres[mnames[0]])):
            for mn in mnames:
                sumtemp[mn] += evagtres[mn][i][3]
            if i+1 == printpoints[index]:
                index += 1
                print(i+1, end='')
                for mn in mnames:
                    print("\t" + str(sumtemp[mn]), end='')
                print()
            elif i == len(evagtres[mnames[0]])-1:
                print(i+1, end='')
                for mn in mnames:
                    print("\t" + str(sumtemp[mn]), end='')
                print()


def evaluation_70benchmarkset(simpathlist, shortnames, mysimloc=0, times=1,
                              bmkpfile="data/benchmarkset_funsim/ground_truth_70_disease_pairs_umlsid.tsv"):

    benchmarkpairs = read_assos(bmkpfile, header=False)
    stat_assos(benchmarkpairs)
    bmptuple = []
    for p in benchmarkpairs.keys():
        for q in benchmarkpairs[p]:
            bmptuple.append((p, q))
    msims = eva_readsims(simpathlist)
    evaress = eva_70benchmarkpairs(msims, bmptuple, times)
    print("---scores and lable---------------------------------------------------------------")
    t = 0
    for er in evaress:
        if t < 1:
            t += 1
            print("-----------------time------------------------------------------------------")
            d1 = list(er.keys())[0]
            d2 = list(er[d1].keys())[0]
            methodnames = list(er[d1][d2].keys())
            methodnames.remove('label')
            print("d1\td2\tlabel", end='')
            for m in methodnames:
                print("\t"+shortnames[m], end='')
            print()
            for d1 in er.keys():
                for d2 in er[d1].keys():
                    print(d1+'\t'+d2+'\t' + str(er[d1][d2]['label']), end='')
                    for m in methodnames:
                        print("\t"+str(er[d1][d2][m]), end='')
                    print()
    print("---ranking disease pairs---------------------------------------------------------")
    t = 0
    for er in evaress:
        if t < 1:
            t += 1
            print("---------------time--------------------------------------------------------")
            ranktemp = eva_cal_ranks(eva_ranking(er))
            ms = list(ranktemp.keys())
            for m in ms:
                print(str(shortnames[m]) + "\t\t\t\t\t", end="")
            print()
            llen = len(ranktemp[ms[0]])
            for i in range(0, llen):
                for m in ms:
                    for j in range(0, 5):
                        print(str(ranktemp[m][i][j])+"\t", end="")
                print()
    print("---tpr fpr----------------------------------------------------------------------")
    tpfprs = eva_tprfprs(evaress)
    t = 0
    for tpfpr in tpfprs:
        if t < 1:
            t += 1
            print("-----------------time------------------------------------------------------")
            methodnames = list(tpfpr.keys())
            llens = []
            for name in methodnames:
                llens.append(len(tpfpr[name]))
                print(str(shortnames[name])+"_fpr\t"+str(shortnames[name])+"_tpr\t", end='')
            print()
            mllen = max(llens)
            for i in range(0, mllen):
                for name in methodnames:
                    if i < len(tpfpr[name]):
                        for j in range(0, 2):
                            print(str(tpfpr[name][i][j])+"\t", end="")
                    else:
                        print("\t\t", end='')
                print()
            auctemp = eva_auc(tpfpr)
            plot_roc(tpfpr, shortnames, auctemp)
    print("---avg auc values------------------------------------------------------------------")
    aucs = eva_aucs(tpfprs)
    auclen = len(aucs)
    mns = list(aucs[0].keys())
    print("--auc at each trial--")
    print('trial', end='')
    for mn in mns:
        print('\t' + shortnames[mn], end='')
    print()
    for i in range(0, len(aucs)):
        print(i+1, end='')
        for mn in mns:
            print('\t' + str(aucs[i][mn]), end='')
        print()
    avgaucs = []
    for mn in mns:
        avgauc = 0.0
        for auc in aucs:
            avgauc += auc[mn]
        avgaucs.append((str(mn), avgauc / auclen))
    avgaucs = sorted(avgaucs, key=lambda a: a[1])
    print("--average auc--")
    for x in avgaucs:
        print(str(shortnames[x[0]])+"\t"+str(x[1]))
    print("---rank stats---------------------------------------------------------------------")
    evarankings = eva_rankings(evaress)
    if bmkpfile == "data/benchmarkset_funsim/ground_truth_70_disease_pairs_umlsid.tsv":
        # ----umls mesh classes------------------------------------
        umlsdiseases = read_umls_disease_info(0.06).getumlsdiseases()
        term2group = {}
        for u in umlsdiseases.keys():
            term2group[u] = umlsdiseases[u].getcategories()
            # ---------------------------------------------------------
    else:
        # ------disease ontology----------------------------------
        do = DiseaseOntology()
        do.readobofile('data/do/HumanDO.obo')
        do.settermslayers()
        layer = 3
        termgroups = get_terms_at_layern(layer, do)
        term2group = invert_dict(get_terms2offsprings(termgroups, do))
        stat_assos(term2group)
        # ---------------------------------------------------------
    # ---same categories stats-----------------
    sameclass = 0
    for p in benchmarkpairs.keys():
        for q in benchmarkpairs[p]:
            if ((p in term2group.keys() and q in term2group.keys()) and
                    (len(term2group[p].intersection(term2group[q])) > 0)):
                sameclass += 1
    print(sameclass, "benchamark pairs are in the same categories.")
    # -----------------------------------------
    print("replicated times:", times)
    for topn in [25, 50, 100, 200, 300, 400, 500]:
        print('------ top', str(topn), '--------------------------')
        rankstats = eva_test_rankstats_multi(evarankings, topn, term2group)
        print("method\ttrue min\ttrue max\ttrue avg\tsame min\tsame max\tsame avg")
        for mn in rankstats.keys():
            truedict = rankstats[mn]['true_labels_num']
            samedict = rankstats[mn]['in_same_class_num']
            print(shortnames[mn], truedict['min'], truedict['max'], truedict['avg'],
                  samedict['min'], samedict['max'], samedict['avg'], sep='\t')
        print('------------------------------------------------')
    print("---rank stats 2---------------------------------------------------------------")
    print("replicated times:", times)
    term2group = {}
    rdtemp = evarankings[0]
    plen = len(rdtemp[list(rdtemp.keys())[0]])
    ms = evarankings[0].keys()
    print('rank', end='')
    for m in ms:
        print('\t' + shortnames[m], end='')
    print()
    for topn in range(0, plen+1):
        rankstats = eva_test_rankstats_multi(evarankings, topn, term2group)
        print(topn, end='')
        for m in ms:
            print('\t' + str(rankstats[m]['true_labels_num']['avg']), end='')
        print()
    print("---rank info------------------------------------------------------------------")
    print("replicated times:", times)
    rankinfos = eva_test_pair_rankinfos(evarankings)
    compsims = []
    for p in simpathlist:
        if p != simpathlist[mysimloc]:
            compsims.append(p)
    riprint = eva_test_pair_rankinfo_ranking(rankinfos, simpathlist[mysimloc], compsims)
    print('disease1', 'disease2', sep='\t', end='')
    for i in range(2, len(riprint['tuplenames'])):
        print('\t' + shortnames[riprint['tuplenames'][i]], end='')
    print()
    tuplelist = riprint['tuplelist']
    for t in tuplelist:
        for i in range(0, len(t) - 1):
            print(str(t[i]) + '\t', end='')
        print(str(t[len(t) - 1]))


def evaluation_sim_histgram(simpathlist, shortnames):
    for simpath in simpathlist:
        sims = read_sims(simpath)
        simsn = common_use.normalize_simdict(sims)
        plot_simshist(simsn, 'hist_' + shortnames[simpath], 50)


def evaluation_sharedgene(simpathlist, shortnames):
    import matplotlib.pyplot as plt
    import numpy as np

    dgassos = read_assos("D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\geneOmimId_diseaseOmimId_in_ppi.txt",
                         False, '\t', 2, 1)

    for simpath in simpathlist:
        sims = read_simmatrix(simpath, True, False)
        sims = common_use.normalize_simdict(sims)
        bardata = {'xtitle': ['0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6',
                              '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0'],
                   'pairssum': [0] * 10,
                   'sgsum': [0] * 10,
                   'ratio': [0.0] * 10}
        pairsclass = {}
        for i in range(0, 10):
            pairsclass[i] = []
        for d1 in sims.keys():
            for d2 in sims[d1].keys():
                value = sims[d1][d2]
                modv = value // 0.1
                if modv < 10.0 and d1 in dgassos.keys() and d2 in dgassos.keys():
                    pairsclass[int(modv)].append((d1, d2))
        for k in pairsclass.keys():
            bardata['pairssum'][k] = len(pairsclass[k])
            for d1, d2 in pairsclass[k]:
                if len(dgassos[d1].intersection(dgassos[d2])) > 0:
                    bardata['sgsum'][k] += 1
        for i in range(0, 10):
            if bardata['pairssum'][i] != 0:
                bardata['ratio'][i] = bardata['sgsum'][i] / bardata['pairssum'][i]

        y_pos = np.arange(len(bardata['xtitle']))
        plt.bar(y_pos, bardata['ratio'], align='center', alpha=0.5, color='green')
        plt.title(shortnames[simpath])
        plt.xticks(y_pos, bardata['xtitle'], rotation=45)
        plt.xlabel("bins")
        plt.ylabel("shared domain ratio")
        plt.savefig('shareddomain_'+shortnames[simpath]+'.png', format='png')
        plt.clf()
# --------------------------------------------------------------------


# ---similarity calculation-------------------------------
def similarity_cal_go():
    disease2gene_symbol = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
                                                             0.06, True, False)
    print("disease gene assos: ", end='')
    stat_assos(disease2gene_symbol)

    geneonto = GeneOntology()
    geneonto.readobofile("go.obo")
    go2gene = read_go_annotation_file("data/gene_association.goa_human")
    print("go2gene assos original: ", end='')
    stat_assos(go2gene)
    go2gene_expand = add_implicit_annotations(go2gene, geneonto)
    print("go2gene assos expanded: ", end='')
    stat_assos(go2gene_expand)
    go2gene_expand_bp = {}
    go2gene_expand_cc = {}
    go2gene_expand_mf = {}
    for go in go2gene_expand.keys():
        if geneonto.getterms()[go].getnamespace() == "biological_process":
            go2gene_expand_bp[go] = deepcopy(go2gene_expand[go])
        if geneonto.getterms()[go].getnamespace() == "cellular_component":
            go2gene_expand_cc[go] = deepcopy(go2gene_expand[go])
        if geneonto.getterms()[go].getnamespace() == "molecular_function":
            go2gene_expand_mf[go] = deepcopy(go2gene_expand[go])
    print("go2gene assos expanded bp: ", end='')
    stat_assos(go2gene_expand_bp)
    print("go2gene assos expanded cc: ", end='')
    stat_assos(go2gene_expand_cc)
    print("go2gene assos expanded mf: ", end='')
    stat_assos(go2gene_expand_mf)

    gene2go_expand_bp = invert_dict(go2gene_expand_bp)
    dsim_go_bp = diseases_similarity_go(list(disease2gene_symbol.keys()), disease2gene_symbol,
                                        gene2go_expand_bp, go2gene_expand_bp)
    write_sims(dsim_go_bp, "similarity_go_bp_umls_dcutoff006.tsv")

    gene2go_expand_cc = invert_dict(go2gene_expand_cc)
    dsim_go_cc = diseases_similarity_go(list(disease2gene_symbol.keys()), disease2gene_symbol,
                                        gene2go_expand_cc, go2gene_expand_cc)
    write_sims(dsim_go_cc, "similarity_go_cc_umls_dcutoff006.tsv")

    gene2go_expand_mf = invert_dict(go2gene_expand_mf)
    dsim_go_mf = diseases_similarity_go(list(disease2gene_symbol.keys()), disease2gene_symbol,
                                        gene2go_expand_mf, go2gene_expand_mf)
    write_sims(dsim_go_mf, "similarity_go_mf_umls_dcutoff006.tsv")


def go2gene_stats():
    geneonto = GeneOntology()
    geneonto.readobofile("go.obo")
    bpcount = 0
    cccount = 0
    mfcount = 0
    for go in geneonto.getterms().keys():
        if geneonto.getterms()[go].getnamespace() == "biological_process":
            bpcount += 1
        if geneonto.getterms()[go].getnamespace() == "cellular_component":
            cccount += 1
        if geneonto.getterms()[go].getnamespace() == "molecular_function":
            mfcount += 1
    print("bp:", bpcount, "cc:", cccount, "mf:", mfcount)

    go2gene = read_go_annotation_file("data/gene_association.goa_human")
    print("go2gene assos original: ", end='')
    stat_assos(go2gene)
    go2gene_expand_bp = {}
    go2gene_expand_cc = {}
    go2gene_expand_mf = {}
    for go in go2gene.keys():
        if geneonto.getterms()[go].getnamespace() == "biological_process":
            go2gene_expand_bp[go] = deepcopy(go2gene[go])
        if geneonto.getterms()[go].getnamespace() == "cellular_component":
            go2gene_expand_cc[go] = deepcopy(go2gene[go])
        if geneonto.getterms()[go].getnamespace() == "molecular_function":
            go2gene_expand_mf[go] = deepcopy(go2gene[go])
    print("go2gene assos original bp: ", end='')
    stat_assos(go2gene_expand_bp)
    print("go2gene assos original cc: ", end='')
    stat_assos(go2gene_expand_cc)
    print("go2gene assos original mf: ", end='')
    stat_assos(go2gene_expand_mf)

    go2gene_expand = add_implicit_annotations(go2gene, geneonto)
    print("go2gene assos expanded: ", end='')
    stat_assos(go2gene_expand)
    go2gene_expand_bp = {}
    go2gene_expand_cc = {}
    go2gene_expand_mf = {}
    for go in go2gene_expand.keys():
        if geneonto.getterms()[go].getnamespace() == "biological_process":
            go2gene_expand_bp[go] = deepcopy(go2gene_expand[go])
        if geneonto.getterms()[go].getnamespace() == "cellular_component":
            go2gene_expand_cc[go] = deepcopy(go2gene_expand[go])
        if geneonto.getterms()[go].getnamespace() == "molecular_function":
            go2gene_expand_mf[go] = deepcopy(go2gene_expand[go])
    print("go2gene assos expanded bp: ", end='')
    stat_assos(go2gene_expand_bp)
    print("go2gene assos expanded cc: ", end='')
    stat_assos(go2gene_expand_cc)
    print("go2gene assos expanded mf: ", end='')
    stat_assos(go2gene_expand_mf)


def similarity_cal_coexpression():
    disease2gene_entrez = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
                                                             0.06, True, True)
    print("disease gene assos: ", end='')
    stat_assos(disease2gene_entrez)
    genecoexp = read_probeid_expfile("data/U133AGNF1B.gcrma.avg.cleared.entrezid.tsv", False)
    print("number of all genes have exp value:s", len(genecoexp))

    alldgs = set()
    for gs in disease2gene_entrez.values():
        for g in gs:
            alldgs.add(g)
    print("numbers of disease genes having exp values:", len(alldgs.intersection(set(genecoexp.keys()))))

    dsim_coexp = diseases_similarity_coexp(set(disease2gene_entrez.keys()), disease2gene_entrez, genecoexp)
    write_sims(dsim_coexp, "similarity_coexp_umls_dgcutoff006.tsv")


def geneid_convert_coexpression():
    gene2probe = read_assos("data/probeid2entrezid_gcrma_anno.tsv", False, "\t", 2, 1)
    stat_assos(gene2probe)
    probeexp = read_probeid_expfile("data/U133AGNF1B.gcrma.avg.cleared.tab", True, '\t')
    print(len(probeexp))
    geneexp = probeexp2geneexp(probeexp, gene2probe)
    print(len(geneexp))
    with open('data/U133AGNF1B.gcrma.avg.cleared.entrezid.tsv', mode='w') as wf:
        for g in geneexp.keys():
            wf.write(g)
            for v in geneexp[g]:
                wf.write('\t'+str(v))
            wf.write('\n')


def similarity_cal_module():
    # disease2gene = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
    #                                                   0.0, True, True)
    # disease2gene = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd_entrezid.tab")
    disease2gene = read_assos('data/disgenet/curated_gene_disease_associations.tsv', True, '\t', 1, 2)
    print("disease gene assos: ", end='')
    stat_assos(disease2gene)

    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    # g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
    #                                        False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    sims = similarity_module.similarity_cal_spavgn(disease2gene, g)
    write_sims(sims, 'outputs/similarity_spavgn_disgenet_curated_interactome.tsv')


def combine_pathway_data_gsea():
    paths = [  # 'data/pathway_gsea/c2.cp.biocarta.v5.1.symbols.gmt',
             # 'data/pathway_gsea/c2.cp.kegg.v5.1.symbols.gmt',
             # 'data/pathway_gsea/c2.cp.reactome.v5.1.symbols.gmt',
             'data/pathway_gsea/c2.all.v5.1.symbols.gmt', ]
    pgassos = combine_pathway_data(paths)
    stat_assos(pgassos)
    # write_assos(pgassos, 'data/pathway_gsea/pgassos.c2.all.v5.1.symbols.tsv')


def get_disease_pathway_assos():
    # disease2gene = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
    #                                                   0.06, True, True)
    # pathway2gene = read_assos("data/pathway_gsea/pgassos.c2.cp.v5.1.entrez.tsv")

    disease2gene = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    pathway2gene = read_assos("data/pathway_gsea/pgassos.c2.cp.bkr.v5.1.symbols.tsv")

    stat_assos(disease2gene)
    stat_assos(pathway2gene)
    dpassos = common_use.hypergeometric_test(disease2gene, pathway2gene, 'fdr_bh', 0.10)
    stat_assos(dpassos)
    pdassos = invert_dict(dpassos)
    for p in pdassos.keys():
        print(p, len(pdassos[p]), sep='\t')
    # write_assos(dpassos,
    #             "data/pathway_gsea/dpassos.c2.cp.v5.1.entrez.disgenet.dgcutoff006.bh005.tsv")
    # write_assos(dpassos, "data/pathway_gsea/dpassos.c2.all.v5.1.symbols.rwrsidd.bh005.tsv")


def similarity_cal_pathway():
    d2p = read_assos("data/kegg_disease/doid2keggpathway.tsv")
    print("disease-pathway assos: ", end='')
    stat_assos(d2p)
    sims = diseases_similarity_pathway_jaccard(list(d2p.keys()), d2p)
    stat_sims(sims)
    write_sims(sims, "outputs/similarity_pathway_jaccard_doid2keggpathway.tsv")


def pathway_in_ppi():
    pathway2gene = read_assos("data/pathway_gsea/pgassos.c2.all.v5.1.symbols.tsv")
    stat_assos(pathway2gene)

    # --read graph---------------------------------------------------------------------------
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
                                           False, False)
    # g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    # ---------------------------------------------------------------------------------------
    print(len(g.vs()), len(g.es()))
    pks = list(pathway2gene.keys())
    for p in pks:
        geneset = pathway2gene[p]
        nodes = g.vs.select(name_in=geneset)
        plcc = g.subgraph(nodes).clusters(mode='WEAK').giant()
        print(p, len(geneset), len(nodes), len(plcc.vs()),
              len(nodes)/len(geneset), len(plcc.vs())/len(nodes), sep='\t')
        if not (len(nodes)/len(geneset) >= 0.8 and len(plcc.vs())/len(nodes) >= 0.4):
            del pathway2gene[p]
    stat_assos(pathway2gene)
    # write_assos(pathway2gene, 'data/pathway_gsea/pgassos.c2.all.v5.1.symbols.0804.hppinwsl.tsv')


def similarity_cal_pathwayvector():
    disease2gene = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    pathway2gene = read_assos('data/pathway_gsea/pgassos.c2.cp.bkr.v5.1.symbols.tsv')
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
                                           False, False)

    sims = similarity_module.similarity_cal_pathwayvector(disease2gene, pathway2gene,
                                                          g, 'SpTransMax_tfidf', 'Cosine')
    write_sims(sims, 'outputs/similarity_dpathwayCosSpmaxidf_bkrpathway_rwrsidd_hppinwsl.tsv')


def get_disease_domain_assos():
    # disease2gene = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    disease2gene = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                      0.06, True, True)
    domain2gene = read_assos("data/org_Hs_eg_db/entrez2symbol2pfam.tsv", True, '\t', 3, 1)

    ddoassos = common_use.hypergeometric_test(disease2gene, domain2gene)
    stat_assos(ddoassos)
    write_assos(ddoassos, "data/org_Hs_eg_db/didoassos_disgenet_dgcutoff006.tsv")


def similarity_cal_domain():
    di2do = read_assos("data/org_Hs_eg_db/didoassos_sidd.tsv", False, "\t", 1, 2)
    print("disease-domain assos: ", end='')
    stat_assos(di2do)
    sims = diseases_similarity_pathway_jaccard(list(di2do.keys()), di2do)
    stat_sims(sims)
    # write_sims(sims, "outputs/similarity_domain_jaccard_di2do_sidd.tsv")


def simlarity_cal_geneseim():
    gfunsim = read_simmatrix('data/test/genesim_disgenet006_wang_bma_gobp.tsv')
    # dgassos = read_assos('data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab')
    dgassos = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                 0.06, True, True)
    stat_assos(dgassos)
    dgs = set()
    for d in dgassos.keys():
        dgs.update(dgassos[d])
    dg2dgsim = similarity_genesim.simmatrix(list(dgs), gfunsim)
    dsims = similarity_genesim.similarity_cal(dgassos, dg2dgsim)
    write_sims(dsims, 'outputs/similarity_genefun_disgenet006_wangbmabp.tsv')


def similarity_cal_combinesim():
    dsim1 = read_sims('outputs/similarity_genefun_rwrsidd_wangbmabp.tsv')
    dsim2 = read_sims('outputs/similarity_spavgn_trans_rwrsidd_hppinwsl.tsv')

    dsim = {}
    for d1 in dsim1.keys():
        dsim[d1] = {}
        for d2 in dsim1[d1].keys():
            dsim[d1][d2] = dsim1[d1][d2] * similarity_genesim.findsimvalue(d1, d2, dsim2)
        print(d1)
    write_sims(dsim, 'outputs/similarity_spavgngenefun_rwrsidd_hppinwsl_wangbmpbp.tsv')


def similarity_cal_katz():
    # g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    # print(len(g.vs()), len(g.es()))
    # similarity_module.katzwalklength(g, 5)
    # ---prepare the data---
    dgassos = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                 0.06, True, True)
    stat_assos(dgassos)
    with open('data/interactome_science/interactome_adjmat1.tsv', mode='r') as rf:
        genes = next(rf).strip().split('\t')
    print('lens of genes:', len(genes))
    gene2loc = {}
    for i in range(0, len(genes)):
        gene2loc[genes[i]] = i
    adjmats = []
    for i in range(1, 6):
        adjmats.append(similarity_module.read_adjmatrix('data/interactome_science/interactome_adjmat'
                                                        + str(i) + '.tsv'))
    # -----------
    sims = similarity_module.similarity_cal_katz(dgassos, gene2loc, adjmats[0:4], None)
    write_sims(sims, 'outputs/similarity_katz4_disgenet_dgcutoff006_interactome_betadef.tsv')
    pass
# -----------------------------------------------------------


# ---mapping-------------------------------------------------
def do2umls_mapping():
    uds = read_umls_disease_info(0)
    do2umls = mapping.doid2umlsid(uds)
    stat_assos(do2umls)
    doid1 = read_one_col("data/benchmarkset_funsim/ground_truth_70_disease_pairs_doid.tsv", 1, True)
    doid2 = read_one_col("data/benchmarkset_funsim/ground_truth_70_disease_pairs_doid.tsv", 2, True)
    doids = set(doid1).union(doid2)
    print("doids:", len(doids))
    do2umls_map = {}
    for d in doids:
        if d in do2umls.keys():
            if len(do2umls[d]) == 1:
                do2umls_map[d] = list(do2umls[d])[0]
            elif d == "DOID:0050700":
                do2umls_map[d] = "umls:C0878544"
            elif d == "DOID:3312":
                do2umls_map[d] = "umls:C0005586"
            elif d == "DOID:6132":
                do2umls_map[d] = "umls:C0006277"
            else:
                print("do2umls_mapping():", d, "has mutiple mappings!")
        elif d == "DOID:83":
            do2umls_map[d] = "umls:C0029531"
    write_mappings(do2umls_map, "data/benchmarkset_funsim/ground_truth_doid2umlsid_47diseases.tsv")
    with open("data/benchmarkset_funsim/ground_truth_70_disease_pairs_doid.tsv", mode='r') as rf:
        with open("data/benchmarkset_funsim/ground_truth_70_disease_pairs_umlsid.tsv", mode='w') as wf:
            for line in rf:
                ids = line.strip().split('\t')
                id1 = ids[0].strip()
                id2 = ids[1].strip()
                if id1 in do2umls_map.keys() and id2 in do2umls_map.keys():
                    wf.write(do2umls_map[id1]+'\t'+do2umls_map[id2]+'\n')


def probeid2entrezid_compare():
    probeid_need2conv = set(read_one_col("data/U133AGNF1B.gcrma.avg.cleared.tab",
                                         1, True))
    print("number of probeids need to be converted:", len(probeid_need2conv))

    p2e_mnd = read_assos("data/probeid2entrezid_mygenendavid.txt")
    print("probeid2entrezid mapping result from mygene and david:\t", end='')
    stat_assos(p2e_mnd)
    p2e_ann = read_assos("data/U133A_annotation.tsv")
    print("probeid2entrezid mapping result from the paper:\t", end='')
    stat_assos(p2e_ann)

    intersects = set(p2e_mnd.keys()).intersection(set(p2e_ann.keys()))
    print("用文章里的mapping和第三方工具mapping都有结果的probeid有", len(intersects), "个")
    print("其中两种方式有不同结果的是下面这些：")
    print("probe\tpaper mapping\tthird tool mapping")
    for p in intersects:
        if p2e_ann[p] != p2e_mnd[p]:
            print(p, p2e_ann[p], p2e_mnd[p], sep='\t')


def get_icd9id2gene_umlsdisgenet():
    # get icd9id 2 gene from disgenet and umls mrconso.rrf
    uds = read_umls_disease_info(0.0)
    icd92umls = mapping.icd9cmid2umlsid_alldigit(uds)
    stat_assos(icd92umls)

    allids_mrconso = allids_mrconsorrf("E:\\UMLS\\2016AA\\META\\MRCONSO.RRF")
    print('mrconso id:', len(allids_mrconso))

    for allid in allids_mrconso:
        if 'icd9cm' in allid.keys():
            for icd9 in allid['icd9cm']:
                icd9id = icd9.split(':')[1]
                if len(icd9id) < 7:
                    if icd9id not in icd92umls.keys():
                        icd92umls[icd9id] = set()
                    icd92umls[icd9id].update(allid['umls'])
    stat_assos(icd92umls)

    umlsterms = uds.getumlsdiseases()
    icd92gene = {}
    for icd9id in icd92umls.keys():
        geneset = set()
        for umlsid in icd92umls[icd9id]:
            if umlsid in umlsterms.keys():
                geneset.update(umlsterms[umlsid].getentrezgenes())
        if len(geneset) > 0:
            icd92gene[icd9id] = geneset
    stat_assos(icd92gene)
    write_assos(icd92gene, "data/comorbidity/icd9cm2entrezgene_disgenet_dgcutoff000.tsv")


def diseaseidmapping_hsdn():
    meshnames1 = read_one_col("ncomms5212-s5.txt", 1, True)
    meshnames2 = read_one_col("ncomms5212-s5.txt", 2, True)

    meshnames = set(meshnames1).union(meshnames2)
    print("meshnames needed to be map:", len(meshnames))
    umlsdiseases = read_umls_disease_info(0)
    meshname2meshid = read_mappings("data/MeshTreeHierarchy.csv", True, "\t", 3, 2)

    meshname2umlsid = mapping.termname2umlsid(meshnames, umlsdiseases, meshname2meshid)
    print("mapped successfully:", len(meshname2umlsid))

    meshname2umlsid_sgl = {}
    for m in meshname2umlsid.keys():
        if len(meshname2umlsid[m]) == 1:
            meshname2umlsid_sgl[m] = list(meshname2umlsid[m])[0]
        else:
            print(m, meshname2umlsid[m], sep='\t')

    print("one to one:", len(meshname2umlsid_sgl))
    stat_maps(meshname2umlsid_sgl)
    stat_assos(meshname2umlsid)
    write_assos(meshname2umlsid, "data/meshtermname2umlsid_hsdn.tsv")


def conv_meshname2umlsid_in_hsdn():
    hsdnsims = read_sims("data/ncomms5212-s5.txt", True)
    stat_sims(hsdnsims)

    meshid2umlsid = read_assos("data/meshtermname2umlsid_hsdn.tsv")
    stat_assos(meshid2umlsid)

    meshid2umlsid_one2one = {}
    for m in meshid2umlsid.keys():
        if len(meshid2umlsid[m]) == 1:
            meshid2umlsid_one2one[m] = list(meshid2umlsid[m])[0]
    print("number of meshids who have unique umlsid:", len(meshid2umlsid_one2one))

    hsdnumlssims = {}
    for m1 in hsdnsims.keys():
        if m1 in meshid2umlsid_one2one.keys():
            hsdnumlssims[meshid2umlsid_one2one[m1]] = {}
            for m2 in hsdnsims[m1].keys():
                if m2 in meshid2umlsid_one2one.keys():
                    hsdnumlssims[meshid2umlsid_one2one[m1]][meshid2umlsid_one2one[m2]] = hsdnsims[m1][m2]
    stat_sims(hsdnumlssims)
    write_sims(hsdnumlssims, "similarity_symptom_1445umlsid.tsv")

    hsdnumlssims = {}
    for m1 in hsdnsims.keys():
        if m1 in meshid2umlsid.keys():
            for um1 in meshid2umlsid[m1]:
                if um1 not in hsdnumlssims.keys():
                    hsdnumlssims[um1] = {}
            for m2 in hsdnsims[m1].keys():
                if m2 in meshid2umlsid.keys():
                    for um2 in meshid2umlsid[m2]:
                        for um1 in meshid2umlsid[m1]:
                            hsdnumlssims[um1][um2] = hsdnsims[m1][m2]
    stat_sims(hsdnumlssims)
    write_sims(hsdnumlssims, "similarity_symptom_1736umlsid.tsv")


def get_umls2keggpathways():
    mesh2pathways = read_assos('data/kegg_disease/mesh2pathway.tsv')
    stat_assos(mesh2pathways)

    uds = read_umls_disease_info(0)
    umlsterms = uds.getumlsdiseases()
    # allids_mrconso = allids_mrconsorrf("E:\\UMLS\\2016AA\\META\\MRCONSO.RRF")
    # print('mrconso id:', len(allids_mrconso))

    umls2mesh = {}
    for ut in umlsterms.keys():
        meshs = umlsterms[ut].getmeshids()
        if len(meshs) > 0:
            umls2mesh[ut] = set()
            umls2mesh[ut].update(meshs)
    # umlsdiseases_cal = set(umlsterms.keys())
    # for allid in allids_mrconso:
    #     umlsidset = allid['umls']
    #     for uid in umlsidset:
    #         if uid in umlsdiseases_cal and 'mesh' in allid.keys():
    #             if uid not in umls2mesh.keys():
    #                 umls2mesh[uid] = set()
    #             umls2mesh[uid].update(allid['mesh'])

    umls2pathway = {}
    for ut in umls2mesh.keys():
        pwtemp = set()
        for msh in umls2mesh[ut]:
            if msh in mesh2pathways.keys():
                pwtemp.update(mesh2pathways[msh])
        if len(pwtemp) > 0:
            umls2pathway[ut] = pwtemp
    stat_assos(umls2pathway)
    # write_assos(umls2pathway, 'data/kegg_disease/umls2keggpathway.tsv')


def get_doid2keggpathways():
    mesh2pathways = read_assos('data/kegg_disease/mesh2pathway.tsv')
    stat_assos(mesh2pathways)

    dionto = DiseaseOntology()
    dionto.readobofile('data/do/HumanDO.obo')
    doterms = dionto.getterms()

    doid2mesh = {}
    for d in doterms.keys():
        xrefs = doterms[d].getxrefs()
        meshset = set()
        for x in xrefs:
            if str(x).startswith("MeSH:") or str(x).startswith("MSH:"):
                meshset.add('mesh:' + str(x).split(':')[1])
        if len(meshset) > 0:
            doid2mesh[d] = meshset
    stat_assos(doid2mesh)

    doid2pathway = {}
    for d in doid2mesh.keys():
        pwset = set()
        for m in doid2mesh[d]:
            if m in mesh2pathways.keys():
                pwset.update(mesh2pathways[m])
        if len(pwset) > 0:
            doid2pathway[d] = pwset
    stat_assos(doid2pathway)
    write_assos(doid2pathway, 'data/kegg_disease/doid2keggpathway.tsv')


def convert_geneid_dgassos():
    dgassos = read_assos('data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab')
    stat_assos(dgassos)
    dgs = set()
    for d in dgassos.keys():
        dgs.update(dgassos[d])
    symbol2entrezid = mapping.geneida2geneidb('symbol', 'entrezgene', list(dgs))
    stat_assos(symbol2entrezid)
    dgassos_conv = mapping.convert_dict_values(dgassos, symbol2entrezid)
    stat_assos(dgassos_conv)
    # write_assos(dgassos_conv, 'data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd_entrezid.tab')
# -----------------------------------------------------------


# ---information---------------------------------------------
def disease_module_info():
    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome_rmslpe.tsv", False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    disease_genes = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv", 0)
    print("disease-gene associations: ", end="")
    stat_assos(disease_genes)

    disease_genes_ori = deepcopy(disease_genes)

    gvs = g.vs['name']
    dkeys = deepcopy(list(disease_genes.keys()))
    for d in dkeys:
        gs = deepcopy(disease_genes[d])
        for gene in gs:
            if gene not in gvs:
                disease_genes[d].remove(gene)
        if len(disease_genes[d]) < 5:
            del disease_genes[d]
    print("disease-gene associations: ", end="")
    stat_assos(disease_genes)

    d_density = similarity_module.density(disease_genes, g)
    print("density finished")
    d_lambdamodule = similarity_module.lambda_module(disease_genes, g)
    print("lambdamodule finished")
    d_avgdegree = similarity_module.avg_degree(disease_genes, g)
    print("avgdegree finished")
    d_avgsp = similarity_module.avg_shortestpath(disease_genes, g)
    print("avgsp finished")
    d_diameter = similarity_module.diameter(disease_genes, g)
    print("diameter finished")

    print("disease\tnumber of genes\tnumber of genes in interactome\tavg degree\t"
          "density\tavgsp\tdiameter\tlambdamodule")
    for d in disease_genes.keys():
        print(d, len(disease_genes_ori[d]), len(disease_genes[d]), d_avgdegree[d],
              d_density[d], d_avgsp[d], d_diameter[d], d_lambdamodule[d], sep="\t")


def disease_module_info_simplify():
    g = similarity_module.read_interactome("data/DataS1_interactome_rmslpe.tsv", False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    disease_genes = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv", 0.06)
    print("disease-gene associations: ", end="")
    stat_assos(disease_genes)

    gvs = set(g.vs['name'])
    disease_genes_filtered = {}
    for d in disease_genes.keys():
        gs = disease_genes[d].intersection(gvs)
        if len(gs) != 0:
            disease_genes_filtered[d] = gs
    print("disease-gene (in interactome) associations: ", end="")
    stat_assos(disease_genes_filtered)

    print("disease\tnumber of genes\tnumber of genes in interactome")
    for d in disease_genes.keys():
        print(d, len(disease_genes[d]), sep='\t', end="\t")
        if d in disease_genes_filtered.keys():
            print(len(disease_genes_filtered[d]))
        else:
            print(str(0))


def gene_neighbor_info():
    g = similarity_module.read_interactome("data/DataS1_interactome.tsv", False, False)
    print("number of vertices:", g.vcount())

    all_disease_genes = read_one_col("data/curated_gene_disease_associations.tsv", 2, True)
    all_disease_genes = set(all_disease_genes)
    print("number of disease genes:", len(all_disease_genes))
    all_disease_genes = all_disease_genes.intersection(set(g.vs['name']))
    print("number of disease genes in interactome:", len(all_disease_genes))

    for gene in all_disease_genes:
        nei = g.vs(g.neighbors(gene))['name']
        dn = all_disease_genes.intersection(nei)
        lennondn = len(nei)-len(dn)
        print(gene+"\t"+str(len(nei))+"\t"+str(len(dn))+"\t"+str(lennondn))

    all_non_disease_genes = set(g.vs['name']).difference(all_disease_genes)
    print("non disease genes info:")
    for gene in all_non_disease_genes:
        nei = g.vs(g.neighbors(gene))['name']
        dn = all_disease_genes.intersection(nei)
        lennondn = len(nei) - len(dn)
        print(gene + "\t" + str(len(nei)) + "\t" + str(len(dn)) + "\t" + str(lennondn))


def benchmarksetpairs_stats():
    # bmkpairs = read_assos("data/benchmarkset_funsim/ground_truth_70_disease_pairs_doid.tsv")
    dgassos = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    stat_assos(dgassos)
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
                                           False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())
    gnodes = g.vs['name']

    pairsrank = []
    with open('C:\\Users\\NP\\Desktop\\pairsrank.txt', mode='r') as rf:
        line = next(rf)
        heads = line.strip().split('\t')
        for fline in rf:
            words = fline.strip().split('\t')
            for i in range(2, len(words)):
                words[i] = float(words[i])
            pairsrank.append(words)

    method = 5
    print("----" + heads[method] + "--------")
    from operator import itemgetter
    pairsrank = sorted(pairsrank, key=itemgetter(method))
    # bmktuple = []
    # for d1 in bmkpairs.keys():
    #     for d2 in bmkpairs[d1]:
    #         bmktuple.append((d1, d2))
    for h in heads:
        print(h, end='\t')
    print()
    for pair in pairsrank:
        for p in pair:
            print(p, end='\t')
        a, b = pair[0], pair[1]
        print(len(dgassos[a]), len(dgassos[b]),
              len(dgassos[a].intersection(gnodes)), len(dgassos[b].intersection(gnodes)),
              len(dgassos[a].intersection(dgassos[b])), sep='\t')


def tissuespecificity_stats():
    expgenes = read_one_col('data/tissuespec_srep/srep35241-s4.txt', 1, True)
    print('number of experssion genes:', len(set(expgenes)))

    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())
    gvs = set(g.vs['name'])
    expnodes = gvs.intersection(expgenes)
    print(len(expnodes))

    disease2gene = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                      0.06, True, True)
    print("disease gene assos: ", end='')
    stat_assos(disease2gene)
# ------------------------------------------------------------


# ---hamaneh and rwr------------------------------------------
def disease_gene_network_hamaneh():
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
                                           False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])

    # disease2gene_entrez = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
    #                                                          0.06, True, True)
    disease2gene_rwrsidd = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    print("disease gene assos: ", end='')
    stat_assos(disease2gene_rwrsidd)

    dgassos_new = {}
    for d in disease2gene_rwrsidd.keys():
        dgleft = gvs.intersection(disease2gene_rwrsidd[d])
        if len(dgleft) >= 1:
            dgassos_new[d] = dgleft
    print("disease gene assos left: ", end='')
    stat_assos(dgassos_new)

    selfloop = set()
    with open("data/hamaneh/rwrsidd_hppin_withselfloop_hamaneh_dgnetwrok.tab", mode='w') as f:
        for e in g.es:
            ns = g.vs[e.source]['name']
            nt = g.vs[e.target]['name']
            if ns == nt:
                f.write(ns+'\t'+nt+'\t4.0\n')
                selfloop.add(ns)
            else:
                f.write(ns+'\t'+nt+"\t1.0\n")
                f.write(nt+'\t'+ns+'\t1.0\n')
        for g in gvs:
            if g not in selfloop:
                f.write(g+'\t'+g+'\t0.0\n')
        for d in dgassos_new.keys():
            f.write(d+'\t'+d+'\t0.0\n')
            for g in dgassos_new[d]:
                f.write(d+'\t'+g+'\t1.0\n')
                f.write(g + '\t' + d + '\t1.0\n')
    with open("data/hamaneh/rwrsidd_hppin_withselfloop_hamaneh_dgnetwork_nodes.tab", mode='w') as f:
        for g in gvs:
            f.write(g+'\n')
        for d in dgassos_new.keys():
            f.write(d+'\n')
    write_slist(list(dgassos_new.keys()), "data/hamaneh/rwrsidd_hppin_withselfloop_boundary.tab")


def get_disease_correlations_hamaneh():
    disease_included = set(read_one_col("data/hamaneh/rwrsidd_included_diseases", 1))
    print(type(disease_included), len(disease_included))
    disease_cors = read_sims("data/hamaneh/rwrsidd_hppin_withselfloop_hamaneh_disease_relations.txt",
                             True, '\t', 1, 2, 3)
    stat_sims(disease_cors)

    newdisease_cors = {}
    for d1 in disease_cors.keys():
        if d1 in disease_included:
            for d2 in disease_cors[d1].keys():
                if d2 in disease_included:
                    if ('DOID:'+d1) not in newdisease_cors.keys():
                        newdisease_cors['DOID:'+d1] = {}
                    newdisease_cors['DOID:'+d1]['DOID:'+d2] = disease_cors[d1][d2]
    stat_sims(newdisease_cors)
    # write_sims(newdisease_cors, "data/hamaneh/similarity_hamaneh_rwrsidd_hppinwsl.tsv")


def experiment():
    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv",
                                           False, False)
    # g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
    #                                        False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])
    # disease2gene_entrez = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
    #                                                          0.06, True, True)
    disease2gene_entrez = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd_entrezid.tab")
    stat_assos(disease2gene_entrez)
    dgassos_new = {}
    for d in disease2gene_entrez.keys():
        dgleft = gvs.intersection(disease2gene_entrez[d])
        if len(dgleft) >= 1:
            dgassos_new[d] = dgleft
    print("disease gene assos left: ", end='')
    stat_assos(dgassos_new)

    # # ---graphlet-------------------------------------------------------
    # gene2sig = experiments.read_gene_signature("data/test/graphlet_interactome_maxcc.tsv")
    # calgenes = set(gene2sig.keys())
    # dgassos_graphlet = {}
    # for d in disease2gene_entrez.keys():
    #     dgleft = calgenes.intersection(disease2gene_entrez[d])
    #     if len(dgleft) >= 1:
    #         dgassos_graphlet[d] = dgleft
    # print("disease gene assos left graphlet: ", end='')
    # stat_assos(dgassos_graphlet)
    # sim_gene2gene = experiments.sim_gene2gene_graphlet(gene2sig)
    # sim_d2d = experiments.sim_geneset2geneset(dgassos_graphlet, sim_gene2gene)
    # write_sims(sim_d2d, "outputs/similarity_experiments_graphlet_less_umls_dgcutoff006.tsv")
    # # ------------------------------------------------------------------

    # # ---bmc rwr------------------------------------------------------------
    sim_gene2geneset = read_simmatrix("data/test/rwr_geneset2genescore_rwrsidd_interactome.tsv")
    sim_d2d = experiments.sim_geneset2geneset_rwr(disease2gene_entrez, sim_gene2geneset)
    write_sims(sim_d2d, "outputs/similarity_rwr_rwrsidd_interactome.tsv")
    # # ----------------------------------------------------------------------

    # ---shortest path------------------------------------------------------
    # sps_norm = experiments.sim_gene2gene_shortestpath(g)

    # sim_d2d = experiments.sim_geneset2geneset(dgassos_new, sps_norm)
    # write_sims(sim_d2d,
    #            "outputs/similarity_experiments_shortestpath_transformed_less_disgenet_dgcutoff006_interactome.tsv")

    # sps_normdivide = experiments.sim_gene2gene_shortestpath(g, False)
    # sim_d2d = experiments.sim_geneset2geneset(disease2gene_entrez, sps_normdivide)
    # write_sims(sim_d2d, "outputs/similarity_experiments_shortestpath_divide_umls_dgcutoff006.tsv")
    # ----------------------------------------------------------------------


def rwr_bmc_ppi():
    ppi_hprd = {}
    with open("data/rwr_bmc_bioinfo/ppi/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", mode='r') as f:
        next(f)
        for line in f:
            words = line.strip().split('\t')
            p1 = words[0].strip()
            p2 = words[3].strip()
            if p1 != '-' and p2 != '-':
                if (p1 in ppi_hprd.keys() and p2 in ppi_hprd[p1]) or (p2 in ppi_hprd.keys() and p1 in ppi_hprd[p2]):
                    # print(words, "replicated!")
                    pass
                else:
                    if p1 not in ppi_hprd.keys():
                        ppi_hprd[p1] = set()
                    ppi_hprd[p1].add(p2)
    stat_network(ppi_hprd)
    # write_assos(ppi_hprd, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_hprd_withselfloop.tab")

    ppi_biogrid = {}
    with open("data/rwr_bmc_bioinfo/ppi/BIOGRID-ORGANISM-Homo_sapiens-3.4.140.tab2.txt", mode='r') as f:
        next(f)
        for line in f:
            words = line.strip().split('\t')
            if words[15].strip() == '9606' and words[16].strip() == '9606' and words[12].strip() == 'physical':
                p1 = words[7].strip()
                p2 = words[8].strip()
                if (p1 in ppi_biogrid.keys() and p2 in ppi_biogrid[p1]) or \
                        (p2 in ppi_biogrid.keys() and p1 in ppi_biogrid[p2]):
                    # print(words, "replicated!")
                    pass
                else:
                    if p1 not in ppi_biogrid.keys():
                        ppi_biogrid[p1] = set()
                    ppi_biogrid[p1].add(p2)
    stat_network(ppi_biogrid)
    # write_assos(ppi_biogrid, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_biogrid_withselfloop.tab")

    # ppi_homomint = {}
    # mg = mygene.MyGeneInfo()
    # uniprot2symbol = {}
    # with open("data/rwr_bmc_bioinfo/ppi/dump-homomint.txt", mode='r') as f:
    #     for line in f:
    #         words = line.strip().split('|')
    #         p1 = words[2].strip()
    #         p2 = words[5].strip()
    #         if p1 == words[1].strip():
    #             if p1 in uniprot2symbol.keys():
    #                 p1 = uniprot2symbol[p1]
    #             else:
    #                 p1 = ''
    #                 gmapping = mg.query(words[1].strip())
    #                 if gmapping['total'] == 1 and 'symbol' in gmapping['hits'][0]:
    #                     p1 = gmapping['hits'][0]['symbol']
    #                 uniprot2symbol[words[1].strip()] = p1
    #             print(line.strip(), p1)
    #         if p2 == words[4].strip():
    #             if p2 in uniprot2symbol.keys():
    #                 p2 = uniprot2symbol[p2]
    #             else:
    #                 p2 = ''
    #                 gmapping = mg.query(words[4].strip())
    #                 if gmapping['total'] == 1 and 'symbol' in gmapping['hits'][0]:
    #                     p2 = gmapping['hits'][0]['symbol']
    #                 uniprot2symbol[words[4].strip()] = p2
    #             print(line.strip(), p2)
    #         if p1 != '' and p2 != '':
    #             if (p1 in ppi_homomint.keys() and p2 in ppi_homomint[p1]) or \
    #                     (p2 in ppi_homomint.keys() and p1 in ppi_homomint[p2]):
    #                 # print(words, "replicated!")
    #                 pass
    #             else:
    #                 if p1 not in ppi_homomint.keys():
    #                     ppi_homomint[p1] = set()
    #                 ppi_homomint[p1].add(p2)
    # stat_network(ppi_homomint)
    # # write_assos(ppi_homomint, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_homomint1_withselfloop.tab")

    ppi_homomint = {}
    with open("data/rwr_bmc_bioinfo/ppi/dump-homomint.txt", mode='r') as f:
        for line in f:
            words = line.strip().split('|')
            p1 = words[1].strip()
            p2 = words[4].strip()
            if (p1 in ppi_homomint.keys() and p2 in ppi_homomint[p1]) or \
                    (p2 in ppi_homomint.keys() and p1 in ppi_homomint[p2]):
                # print(words, "replicated!")
                pass
            else:
                if p1 not in ppi_homomint.keys():
                    ppi_homomint[p1] = set()
                ppi_homomint[p1].add(p2)
    alluniprots = set()
    alluniprots |= set(ppi_homomint.keys())
    for p in ppi_homomint.keys():
        alluniprots |= ppi_homomint[p]
    uni2sym = mapping.geneida2geneidb('uniprot', 'symbol', list(alluniprots))
    ppi_homomint_t = {}
    for p in ppi_homomint.keys():
        if p in uni2sym.keys() and len(uni2sym[p]) == 1:
            ppi_homomint_t[list(uni2sym[p])[0]] = set()
            for q in ppi_homomint[p]:
                if q in uni2sym.keys() and len(uni2sym[q]) == 1:
                    ppi_homomint_t[list(uni2sym[p])[0]].add(list(uni2sym[q])[0])
    stat_network(ppi_homomint_t)
    # write_assos(ppi_homomint_t, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_homomint2_withselfloop.tab")

    ppi_intact = {}
    # interactiontypes = set(read_one_col("data/rwr_bmc_bioinfo/ppi/intact.txt", 11))
    # pprint(interactiontypes)
    taxidhuman = 'taxid:9606(human)|taxid:9606(Homo sapiens)'
    with open("data/rwr_bmc_bioinfo/ppi/intact.txt", mode='r', encoding='utf-8') as f:
        next(f)
        for line in f:
            words = line.strip().split('\t')
            if words[9].strip() == taxidhuman and words[10].strip() == taxidhuman:
                p1s = words[4].strip().split('|')
                p2s = words[5].strip().split('|')
                p1 = ''
                p2 = ''
                for p in p1s:
                    if p.startswith("uniprotkb:") and p.endswith("(gene name)"):
                        # here has 'EIF4G1 variant protein' exception
                        # here has exceptions such as 'MLL/AF6 fusion'
                        # here has '"WUGSC:' exception （没有处理此异常，需人工修改结果）
                        p1 = p.split(':')[1].split('(')[0].split(' ')[0]
                        break
                for p in p2s:
                    if p.startswith("uniprotkb:") and p.endswith("(gene name)"):
                        p2 = p.split(':')[1].split('(')[0].split(' ')[0]
                        break
                if p1 != '' and p2 != '':
                    if ((p1 in ppi_intact.keys() and p2 in ppi_intact[p1]) or
                            (p2 in ppi_intact.keys() and p1 in ppi_intact[p2])):
                        # print(words, "replicated!")
                        pass
                    else:
                        if p1 not in ppi_intact.keys():
                            ppi_intact[p1] = set()
                        ppi_intact[p1].add(p2)
    stat_network(ppi_intact)
    # write_assos(ppi_intact, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_intact_withselfloop.tab")


def rwr_bmc_ppi_combine():
    ppi_hprd = read_assos("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hprd_withselfloop.tab")
    print("hprd: ", end='')
    stat_network(ppi_hprd)
    ppi_biogrid = read_assos('data/rwr_bmc_bioinfo/ppi/rwr_ppi_biogrid_withselfloop.tab')
    print("biogrid: ", end='')
    stat_network(ppi_biogrid)
    ppi_intact = read_assos('data/rwr_bmc_bioinfo/ppi/rwr_ppi_intact_withselfloop.tab')
    print('intact: ', end='')
    stat_network(ppi_intact)
    ppi_homomint = read_assos('data/rwr_bmc_bioinfo/ppi/rwr_ppi_homomint2_withselfloop.tab')
    print('homomint: ', end='')
    stat_network(ppi_homomint)

    ppi_all = deepcopy(ppi_hprd)
    for p in ppi_biogrid.keys():
        for q in ppi_biogrid[p]:
            if not ((p in ppi_all.keys() and q in ppi_all[p]) or
                    (q in ppi_all.keys() and p in ppi_all[q])):
                if p not in ppi_all.keys():
                    ppi_all[p] = set()
                ppi_all[p].add(q)
    stat_network(ppi_all)
    for p in ppi_intact.keys():
        for q in ppi_intact[p]:
            if not ((p in ppi_all.keys() and q in ppi_all[p]) or
                    (q in ppi_all.keys() and p in ppi_all[q])):
                if p not in ppi_all.keys():
                    ppi_all[p] = set()
                ppi_all[p].add(q)
    stat_network(ppi_all)
    for p in ppi_homomint.keys():
        for q in ppi_homomint[p]:
            if not ((p in ppi_all.keys() and q in ppi_all[p]) or
                    (q in ppi_all.keys() and p in ppi_all[q])):
                if p not in ppi_all.keys():
                    ppi_all[p] = set()
                ppi_all[p].add(q)
    stat_network(ppi_all)
    # write_assos(ppi_all, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab")
    pks = list(ppi_all.keys())
    for p in pks:
        ppi_all[p].discard(p)
        if len(ppi_all[p]) == 0:
            del ppi_all[p]
    stat_network(ppi_all)
    # write_assos(ppi_all, "data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop.tab")


def rwr_bmc_dgassos():
    dg_sidd = {}
    with open("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.txt", mode='r') as f:
        next(f)
        for line in f:
            words = line.strip().split('\t')
            genes = words[3].strip().split('|')
            d = words[1].strip()
            dg_sidd[d] = set()
            for g in genes:
                gtemp = g.strip().split(" ")
                for gt in gtemp:
                    dg_sidd[d].add(gt.strip(";").strip())
    stat_assos(dg_sidd)
    # write_assos(dg_sidd, "data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")

    dg_disgenet = {}
    with open("data/rwr_bmc_bioinfo/dg/rwr_dgassos_disgenet.txt", mode='r') as f:
        next(f)
        for line in f:
            words = line.strip().split('\t')
            genes = words[3].strip().split('|')
            d = words[1].strip()
            dg_disgenet[d] = set()
            for g in genes:
                gtemp = g.strip().split(" ")
                for gt in gtemp:
                    dg_disgenet[d].add(gt.strip(";").strip())
    stat_assos(dg_disgenet)
    # write_assos(dg_disgenet, "data/rwr_bmc_bioinfo/dg/rwr_dgassos_disgenet.tab")


def rwr_dgfilter():
    dgassos = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd_entrezid.tab")
    stat_assos(dgassos)
    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv",
                                           False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])
    dgassos_new = {}
    for d in dgassos.keys():
        dgleft = gvs.intersection(dgassos[d])
        if len(dgleft) >= 1:
            dgassos_new[d] = dgleft
    print("disease gene assos left: ", end='')
    stat_assos(dgassos_new)
    write_assos(dgassos_new, 'data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd_entrezid_ininteractome.tab')
# -------------------------------------------------------------


# ---disease id mapping stats----------------------------------
def parse_morbidmap(morbidmap):
    mim = read_one_col(morbidmap, 1, True)
    mimnumber = set()
    for m in mim:
        words = str(m).strip().split(',')
        mtemp = words[len(words)-1].strip()
        mtemp = mtemp.split('(')[0].strip()
        if len(mtemp) == 6 and mtemp[0] in [str(i) for i in range(0, 10)]:
            mimnumber.add(mtemp)
    return mimnumber


def diseaseid_mapping_stats():
    uds = read_umls_disease_info(0)
    dionto = DiseaseOntology()
    dionto.readobofile('data/do/HumanDO.obo')
    hponto = HumanPhenotypeOntology()
    hponto.readobofile('data/hpo/hp.obo')

    allids_umls = []
    udsterms = uds.getumlsdiseases()
    for u in udsterms.keys():
        kids = dict()
        kids['umls'] = set()
        kids['umls'].add(udsterms[u].getid())

        meshids = udsterms[u].getmeshids()
        if len(meshids) > 0:
            kids['mesh'] = set()
            kids['mesh'].update(meshids)

        doids = udsterms[u].getdoids()
        if len(doids) > 0:
            kids['do'] = set()
            kids['do'].update(doids)

        omims = udsterms[u].getomimids()
        if len(omims) > 0:
            kids['omim'] = set()
            kids['omim'].update(omims)

        hpoids = udsterms[u].gethpoids()
        if len(hpoids) > 0:
            kids['hpo'] = set()
            kids['hpo'].update(hpoids)

        icd9cmids = udsterms[u].geticd9cmids()
        if len(icd9cmids) > 0:
            kids['icd9cm'] = set()
            kids['icd9cm'].update(icd9cmids)

        allids_umls.append(kids)
    print('umls id:', len(allids_umls))

    allids_do = []
    doterms = dionto.getterms()
    for d in doterms.keys():
        kids = dict()
        kids['do'] = set()
        kids['do'].add(doterms[d].getdoid())

        xrefs = doterms[d].getxrefs()
        for x in xrefs:
            if str(x).startswith("ICD9:") or str(x).startswith("ICD9CM:"):
                if 'icd9cm' not in kids.keys():
                    kids['icd9cm'] = set()
                kids['icd9cm'].add('icd9cm:' + str(x).split(':')[1])
            elif str(x).startswith("MeSH:") or str(x).startswith("MSH:"):
                if 'mesh' not in kids.keys():
                    kids['mesh'] = set()
                kids['mesh'].add('mesh:' + str(x).split(':')[1])
            elif str(x).startswith("UMLS_CUI:") or str(x).startswith("UMLS:"):
                if 'umls' not in kids.keys():
                    kids['umls'] = set()
                kids['umls'].add('umls:' + str(x).split(':')[1])
            elif str(x).startswith("OMIM:") or str(x).startswith("OMM:"):
                if 'omim' not in kids.keys():
                    kids['omim'] = set()
                kids['omim'].add('omim:' + str(x).split(':')[1])
            elif str(x).startswith("HP:"):
                if 'hpo' not in kids.keys():
                    kids['hpo'] = set()
                kids['hpo'].add('HP:' + str(x).split(':')[1])
        ncount = 0
        for kk in kids.keys():
            if len(kids[kk]) >= 1:
                ncount += 1
        if ncount > 1:
            allids_do.append(kids)
    print('doid:', len(allids_do))

    allids_hpo = []
    hpoterms = hponto.getterms()
    for h in hpoterms.keys():
        kids = dict()
        kids['hpo'] = set()
        kids['hpo'].add(hpoterms[h].gethpoid())

        xrefs = hpoterms[h].getxrefs()
        for x in xrefs:
            if str(x).startswith("MeSH:") or str(x).startswith("MESH:") or str(x).startswith("Mesh:"):
                if 'mesh' not in kids.keys():
                    kids['mesh'] = set()
                kids['mesh'].add('mesh:' + str(x).split(':')[1])
            if str(x).startswith("UMLS:") or str(x).startswith("UMLS_CUI:"):
                if 'umls' not in kids.keys():
                    kids['umls'] = set()
                kids['umls'].add('umls:' + str(x).split(':')[1])
            if str(x).startswith("DOID:"):
                if 'do' not in kids.keys():
                    kids['do'] = set()
                kids['do'].add('DOID:' + str(x).split(':')[1])
            if str(x).startswith("ICD-9:"):
                if 'icd9cm' not in kids.keys():
                    kids['icd9cm'] = set()
                kids['icd9cm'].add('icd9cm:' + str(x).split(':')[1])
        ncount = 0
        for kk in kids.keys():
            if len(kids[kk]) >= 1:
                ncount += 1
        if ncount > 1:
            allids_hpo.append(kids)
        # allids_hpo.append(kids)
    print('hpo id:', len(allids_hpo))

    allids_mrconso = allids_mrconsorrf("E:\\UMLS\\2016AA\\META\\MRCONSO.RRF")
    print('mrconso id:', len(allids_mrconso))

    allids = allids_do + allids_mrconso + allids_umls + allids_hpo
    allids = combine_allids_asumls(allids)
    print("allids original:", len(allids))
    # ---filter---
    meshids = set(read_one_col("data/mesh/meshtreehierarchy_C_F123.txt", 2, True, '\t', "GBK"))
    meshidscount = set()
    for m in meshids:
        meshidscount.add("mesh:" + m)
    print(len(meshidscount))
    omimids = parse_morbidmap('data/morbidmap.txt')
    omimidscount = set()
    for o in omimids:
        omimidscount.add("omim:" + o)
    omimids = set(read_one_col("data/diseasename_new.txt", 1))
    for o in omimids:
        omimidscount.add("omim:" + o)
    print(len(omimidscount))
    allidsnew = []
    for allid in allids:
        if 'umls' in allid.keys() and len(allid.keys()) == 2:
            if 'mesh' in allid.keys():
                for msh in allid['mesh']:
                    if msh in meshidscount:
                        allidsnew.append(allid)
                        break
            elif 'omim' in allid.keys():
                for mim in allid['omim']:
                    if mim in omimidscount:
                        allidsnew.append(allid)
                        break
            else:
                allidsnew.append(allid)
        else:
            allidsnew.append(allid)
    print("allids filtered:", len(allidsnew))
    # ------------
    umlsidcount = set()
    for allid in allidsnew:
        if 'umls' in allid.keys():
            umlsidcount.update(allid['umls'])
    print('umlsid:', len(umlsidcount))
    analyze_allids(allidsnew, False)

    idclasses = {'umls', 'do', 'omim', 'hpo', 'icd9cm', 'mesh'}

    # ----omim test---------------------------------------
    print("omim test:")
    omimids = set(read_one_col("data/diseasename_new.txt", 1))
    print(len(omimids))
    omimtestresult = idmapping_test(omimids, allidsnew)
    print(len(omimtestresult))
    print(omimids.difference(set(omimtestresult.keys())))
    omimmapcount = {}
    for idc in idclasses:
        omimmapcount[idc] = 0
    for mimid in omimtestresult.keys():
        for oidtype in omimtestresult[mimid].keys():
            omimmapcount[oidtype] += 1
    pprint(omimmapcount)
    omimidsp = set()
    for m in omimids:
        omimidsp.add('omim:' + m)
    allids_mim = []
    for allid in allidsnew:
        if 'omim' in allid.keys() and len(allid['omim'].intersection(omimidsp)) > 0:
            allids_mim.append(allid)
    print(len(allids_mim))
    notmap = set()
    for mim in omimtestresult.keys():
        if 'do' in omimtestresult[mim].keys():
            print(mim, omimtestresult[mim]['do'], omimtestresult[mim]['umls'], sep='\t')
        else:
            notmap.add(mim)
    for mim in notmap:
        print(mim, '', omimtestresult[mim]['umls'], sep='\t')
    mim2do = {}
    for mim in omimtestresult.keys():
        if mim not in notmap:
            mim2do[mim] = omimtestresult[mim]['do']
    stat_assos(mim2do)
    # umlsid_rel = set()
    # for mim in notmap:
    #     umlsid_rel.update(omimtestresult[mim]['umls'])
    # umlsid_ext = {}
    # for uid in umlsid_rel:
    #     umlsid_ext[uid] = set()
    # with open('E:\\UMLS\\2016AA\\META\\MRREL.RRF', mode='r') as rf:
    #     for line in rf:
    #         words = line.strip().split('|')
    #         umlsid1 = 'umls:' + words[0].strip()
    #         umlsid2 = 'umls:' + words[4].strip()
    #         if umlsid1 in umlsid_rel:
    #             umlsid_ext[umlsid1].add(umlsid2)
    #         if umlsid2 in umlsid_rel:
    #             umlsid_ext[umlsid2].add(umlsid1)
    # stat_assos(umlsid_ext)
    # umlsallid = {}
    # for allid in allidsnew:
    #     if 'umls' in allid.keys():
    #         for uid in allid['umls']:
    #             umlsallid[uid] = allid
    # # stat_maps(umlsallid)
    # for mim in notmap:
    #     dotemp = set()
    #     fstumlss = omimtestresult[mim]['umls']
    #     sndumlss = set()
    #     for uid in fstumlss:
    #         sndumlss.update(umlsid_ext[uid])
    #     for uid in sndumlss:
    #         if uid in umlsallid.keys():
    #             alltemp = umlsallid[uid]
    #             if 'do' in alltemp.keys():
    #                 dotemp.update(alltemp['do'])
    #     if len(dotemp) > 0:
    #         mim2do[mim] = set()
    #         mim2do[mim].update(dotemp)
    # stat_assos(mim2do)
    # for mim in mim2do:
    #     if mim not in notmap:
    #         print(mim, mim2do[mim])
    # for mim in mim2do:
    #     if mim in notmap:
    #         print(mim, mim2do[mim])
    # -----------------------------------------------------
    # --doid test------------------------------------------
    # print("doid test:")
    # doids = set(read_one_col("data/disease_137_name2doid.tsv", 2))
    # print(len(doids))
    # dotestresult = idmapping_test(doids, allids, itype='do', prefix=True)
    # print(len(dotestresult))
    # print(doids.difference(set(dotestresult.keys())))
    # domapcount = {}
    # for idc in idclasses:
    #     domapcount[idc] = 0
    # for did in dotestresult.keys():
    #     for oidtype in dotestresult[did].keys():
    #         domapcount[oidtype] += 1
    # pprint(domapcount)
    # allids_doid = []
    # for allid in allidsnew:
    #     if 'do' in allid.keys() and len(allid['do'].intersection(doids)) > 0:
    #         allids_doid.append(allid)
    # print(len(allids_doid))
    # notmap = set()
    # for did in dotestresult.keys():
    #     if 'omim' in dotestresult[did].keys():
    #         print(did, dotestresult[did]['omim'], dotestresult[did]['umls'], sep='\t')
    #     else:
    #         notmap.add(did)
    # for did in notmap:
    #     print(did, '', dotestresult[did]['umls'], sep='\t')
    # do2mim = {}
    # for did in dotestresult.keys():
    #     if did not in notmap:
    #         do2mim[did] = dotestresult[did]['omim']
    # stat_assos(do2mim)
    # -----------------------------------------------------
    # ---omim2mesh-----------------------------------------
    # omim2mesh = get_idmapping(allids)
    # stat_assos(omim2mesh)
    # # write_assos(omim2mesh, 'omim2mesh.tsv')
    # -----------------------------------------------------


def allids_mrconsorrf(filepath):
    mappings = []
    mtemp = {}
    idtypemap = {'MSH': 'mesh:', 'ICD9CM': 'icd9cm:',
                 'OMIM': 'omim:', 'HPO': ''}
    idtype2class = {'omim': 'omim', 'icd9cm': 'icd9cm',
                    'HP': 'hpo', 'mesh': 'mesh'}
    with open(filepath, mode='r', encoding='utf-8') as rf:
        for line in rf:
            words = line.strip().split('|')
            umlsid = 'umls:' + words[0].strip()
            # language = words[1].strip()
            otherid = words[13].strip()
            if otherid.startswith('MTHU'):
                otherid = ''
            othertype = words[11].strip()

            if othertype in idtypemap.keys() and otherid != '' and otherid != 'NOCODE':
                if umlsid not in mtemp.keys():
                    mtemp[umlsid] = set()
                mtemp[umlsid].add(idtypemap[othertype]+otherid)
    for umlsid in mtemp.keys():
        mids = dict()
        mids['umls'] = set()
        mids['umls'].add(umlsid)
        for oid in mtemp[umlsid]:
            idtype = str(oid).split(':')[0]
            if idtype in idtype2class.keys():
                idclass = idtype2class[idtype]
                if idclass not in mids.keys():
                    mids[idclass] = set()
                mids[idclass].add(oid)
        mappings.append(mids)
    return mappings


def idmapping_test(testids, allids, itype='omim', prefix=False):
    testidmap = {}
    for mids in allids:
        if itype in mids.keys():
            tids = mids[itype]
            for tid in tids:
                idtemp = tid
                if not prefix:
                    idtemp = str(idtemp).split(":")[1]
                if idtemp in testids:
                    if idtemp not in testidmap.keys():
                        testidmap[idtemp] = {}
                    for mk in mids.keys():
                        if mk not in testidmap[idtemp].keys():
                            testidmap[idtemp][mk] = set()
                        testidmap[idtemp][mk].update(mids[mk])
    return testidmap


def analyze_allids(allids, one2one):
    """
    analyze allids
    :param allids: a list of dicts
    :param one2one: True or False
    :return:
    """
    meshids = set(read_one_col("data/mesh/meshtreehierarchy_C_F123.txt", 2, True, '\t', "GBK"))
    meshidscount = set()
    for m in meshids:
        meshidscount.add("mesh:" + m)
    print('number of mesh ids we care about:', len(meshidscount))
    omimids = parse_morbidmap('data/morbidmap.txt')
    omimidscount = set()
    for o in omimids:
        omimidscount.add("omim:" + o)
    omimids = set(read_one_col("data/diseasename_new.txt", 1))
    for o in omimids:
        omimidscount.add("omim:" + o)
    print('number of omim ids we care about:', len(omimidscount))

    idtypes = ['umls', 'do', 'omim', 'hpo', 'icd9cm', 'mesh']
    stats_matrix = {}

    for i in range(0, len(idtypes)):
        keyid = idtypes[i]
        # ------------------------------
        stats_matrix[keyid] = {}
        stats_matrix[keyid][keyid] = '-'
        # ------------------------------
        for j in range(0, len(idtypes)):
            if i != j:
                valueid = idtypes[j]
                ijassos = {}
                for aid in allids:
                    if keyid in aid.keys() and valueid in aid.keys():
                        for kid in aid[keyid]:
                            if kid not in ijassos.keys():
                                ijassos[kid] = set()
                            ijassos[kid].update(aid[valueid])
                # ----特殊要求-----
                if keyid == 'mesh':
                    ks = list(ijassos.keys())
                    for k in ks:
                        if k not in meshidscount:
                            del ijassos[k]
                if valueid == 'mesh':
                    ks = list(ijassos.keys())
                    for k in ks:
                        ijassos[k] = ijassos[k].intersection(meshidscount)
                        if len(ijassos[k]) == 0:
                            del ijassos[k]
                if keyid == 'omim':
                    ks = list(ijassos.keys())
                    for k in ks:
                        if k not in omimidscount:
                            del ijassos[k]
                if valueid == 'omim':
                    ks = list(ijassos.keys())
                    for k in ks:
                        ijassos[k] = ijassos[k].intersection(omimidscount)
                        if len(ijassos[k]) == 0:
                            del ijassos[k]
                # ---------------
                if one2one:
                    invertij = invert_dict(ijassos)
                    iks = list(invertij.keys())
                    rmks = set()
                    for ik in iks:
                        if len(invertij[ik]) > 1:
                            rmks.update(invertij[ik])
                    ks = list(ijassos.keys())
                    for k in ks:
                        if len(ijassos[k]) > 1 or k in rmks:
                            del ijassos[k]
                # ---------------
                print(keyid, valueid, ' ', end='')
                stat_assos(ijassos)
                stats_matrix[keyid][valueid] = len(list(ijassos.keys()))
    # ----print stats 1----------------------
    for idt in idtypes:
        print("\t"+idt, end='')
    print()
    for id1 in idtypes:
        print(id1, end='')
        for id2 in idtypes:
            print("\t"+str(stats_matrix[id1][id2]), end='')
        print()
    countsum = {'umls': 52355, 'mesh': len(meshidscount),
                'omim': len(omimidscount), 'do': 6930, 'icd9cm': 22406, 'hpo': 11813}
    for idt in idtypes:
        print("\t"+idt, end='')
    print()
    for id1 in idtypes:
        print(id1, end='')
        for id2 in idtypes:
            if id1 == id2:
                print("\t" + str(stats_matrix[id1][id2]), end='')
            else:
                print("\t"+str(stats_matrix[id1][id2]/countsum[id1])[0:6], end='')
        print()
    # -------------------------------------
    # ----print stats 2 do2omim------------
    for allid in allids:
        if 'do' in allid.keys() or 'omim' in allid.keys():
            if 'umls' in allid.keys():
                print(allid['umls'], end='')
            else:
                print('', end='')
            print('\t', end='')
            if 'do' in allid.keys():
                print(allid['do'], end='')
            else:
                print('', end='')
            print('\t', end='')
            if 'omim' in allid.keys():
                print(allid['omim'])
            else:
                print('')
    # -------------------------------------


def get_idmapping(allids, idtype1='omim', idtype2='mesh'):
    mapassos = {}
    for allid in allids:
        if idtype1 in allid.keys() and idtype2 in allid.keys():
            type1ids = allid[idtype1]
            type2ids = allid[idtype2]
            for type1id in type1ids:
                if type1id not in mapassos.keys():
                    mapassos[type1id] = set()
                mapassos[type1id].update(type2ids)
    return mapassos


def combine_allids_asumls(allids):
    allidsnew = []
    alliddic = {}
    for allid in allids:
        if 'umls' in allid.keys():
            for umlsid in allid['umls']:
                if umlsid not in alliddic.keys():
                    alliddic[umlsid] = set()
                aks = set(allid.keys())
                aks.discard('umls')
                for ak in aks:
                    alliddic[umlsid].update(allid[ak])
        else:
            allidsnew.append(allid)

    idtype2class = {'omim': 'omim', 'icd9cm': 'icd9cm',
                    'HP': 'hpo', 'mesh': 'mesh', 'DOID': 'do', }
    for umlsid in alliddic.keys():
        mids = dict()
        mids['umls'] = set()
        mids['umls'].add(umlsid)
        for oid in alliddic[umlsid]:
            idtype = str(oid).split(':')[0]
            if idtype in idtype2class.keys():
                idclass = idtype2class[idtype]
                if idclass not in mids.keys():
                    mids[idclass] = set()
                mids[idclass].add(oid)
        allidsnew.append(mids)
    return allidsnew
# --------------------------------------------------------------


# ---cal disease sim for evaluation with disease gene/mirna prediction--------------------
def similarity_cal_predeva():
    # ---disease gene prediction--------------------------------------------
    dgassos = read_assos("data/test/geneOmimId_diseaseOmimId_in_ppi.txt",
                         False, "\t", 2, 1)
    stat_assos(dgassos)
    omim2hprd = read_mappings("D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\HPRD_ID_MAPPINGS.txt",
                              False, "\t", 6, 1)
    # omim2entrez = read_mappings("data/test/HPRD_ID_MAPPINGS.txt",
    #                             False, "\t", 6, 5)
    stat_maps(omim2hprd)
    dgassos = mapping.convert_dict_values(dgassos, omim2hprd)
    stat_assos(dgassos)

    g = similarity_module.read_interactome("D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\HPRD_ppi.txt",
                                           False, False)
    # g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    print(len(g.vs), len(g.es))

    dsim = similarity_module.similarity_cal_spmaxn(dgassos, g)

    # dsim = common_use.normalize_simdict(dsim)
    write_simmatrix(dsim,
                    "data/test/similarity_spmaxn_dgomim_hprd_matrix.tsv")
    # ---------------------------------------------------------------------------

    # ----disease mirna prediction-----------------------------------------------
    # d137order = read_one_col("data/test/disease_137_SIDD.txt", 1)
    # doname2doid = read_mappings("data/test/disease_137_name2doid.tsv", False)
    # dgassos = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    # # dgassos = read_assos("data/test/doid_result.txt", False, "\t", 2, 3)
    # dgassos_cal = {}
    # for dname in d137order:
    #     if doname2doid[dname] in dgassos.keys():
    #         dgassos_cal[dname] = dgassos[doname2doid[dname]]
    #
    # g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/normalize_ppi_symbol.txt", False, False)
    # print(len(g.vs), len(g.es))
    # # dsim = similarity_module.similarity_cal_spavgn(dgassos_cal, g)
    # sps_norm = experiments.sim_gene2gene_shortestpath(g)
    # dsim = experiments.sim_geneset2geneset(dgassos_cal, sps_norm)
    # dsim = common_use.normalize_simdict(dsim)
    # write_simmatrix(dsim, "data/test/similarity_spmax_sidd_humannet.tsv", True, d137order)
    #
    # g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab", False, False)
    # print(len(g.vs), len(g.es))
    # # dsim = similarity_module.similarity_cal_spavgn(dgassos_cal, g)
    # sps_norm = experiments.sim_gene2gene_shortestpath(g)
    # dsim = experiments.sim_geneset2geneset(dgassos_cal, sps_norm)
    # dsim = common_use.normalize_simdict(dsim)
    # write_simmatrix(dsim, "data/test/similarity_spmax_sidd_hppin.tsv", True, d137order)
    # -----------------------------------------------------------------------------------


def convert_similarity():
    from operator import itemgetter
    dgassos = read_assos("D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\"
                         "geneOmimId_diseaseOmimId_in_ppi.txt", False, '\t', 2, 1)
    stat_assos(dgassos)
    count = 0
    diseases = list(dgassos.keys())
    for i in range(0, len(diseases)-1):
        for j in range(i+1, len(diseases)):
            if len(dgassos[diseases[i]].intersection(dgassos[diseases[j]])) != 0:
                count += 1
    print('disease pairs:', count)
    origsim = read_simmatrix('D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\'
                             'similarity_spavgn_dgomim_hprd_matrix.tsv', True, False)
    stat_sims(origsim)
    simvalue = []
    for d1 in origsim.keys():
        for d2 in origsim[d1].keys():
            if d1 != d2:
                simvalue.append((d1, d2, origsim[d1][d2]))
    print('simvalue len:', len(simvalue))
    simvalue = sorted(simvalue, key=itemgetter(2), reverse=True)
    simvalue = simvalue[0:2*count]

    with open('D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\'
              'diseasepairs_spavgn_dgomim_hprd.tsv', mode='w') as wf:
        for simv in simvalue:
            wf.write(simv[0] + '\t' + simv[1] + '\n')

    dneighbors = {}
    for simv in simvalue:
        a, b = simv[0], simv[1]
        if a not in dneighbors.keys():
            dneighbors[a] = set()
        if b not in dneighbors.keys():
            dneighbors[b] = set()
        dneighbors[a].add(b)
        dneighbors[b].add(a)
    stat_assos(dneighbors)

    # ressim = {}
    # for i in range(0, len(diseases)):
    #     ressim[diseases[i]] = {}
    #     for j in range(i, len(diseases)):
    #         if i == j:
    #             ressim[diseases[i]][diseases[j]] = 1.0
    #         else:
    #             ressim[diseases[i]][diseases[j]] = 0.0
    #
    # for d1 in ressim.keys():
    #     for d2 in ressim[d1].keys():
    #         if d1 in dneighbors.keys() and d2 in dneighbors.keys():
    #             pass  # something like ecc
    # write_simmatrix(ressim, 'D:\\bioinformatics\\tools\\zrq\\PDGTR\\example\\'
    #                         'similarity_spavgn_ecc_dgomim_hprd_matrix.tsv')
# --------------------------------------------------------------


if __name__ == "__main__":
    # evaluation_groundtruth(evaluation_simfilepaths1, shortnames1, [gtpathlist1[1], ])
    # evaluation_70benchmarkset(evaluation_simfilepaths2, shortnames2, 0, 100,
    #                           'data/benchmarkset_funsim/ground_truth_70_disease_pairs_umlsid.tsv')
    # evaluation_validationpairs(evaluation_simfilepaths4, shortnames4, 100)
    diseaseid_mapping_stats()
    pass
