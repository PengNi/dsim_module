#! /usr/bin/env python3
from copy import deepcopy
from pprint import pprint
import similarity_module
import mapping
import experiments
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
from files import write_slist
from umls_disease import read_all_gene_disease_associations
from umls_disease import read_umls_disease_info
from similarity_coexpression import read_probeid_expfile
from similarity_coexpression import probeexp2geneexp
from similarity_coexpression import diseases_similarity_coexp
from gene_ontology import read_go_annotation_file
from gene_ontology import GeneOntology
from gene_ontology import invert_dict
from gene_ontology import add_implicit_annotations
from similarity_go import diseases_similarity_go
from evaluation import eva_readsims
from evaluation import eva_70benchmarkpairs
from evaluation import eva_ranking
from evaluation import eva_tprfpr
from evaluation import eva_tprfprs
from evaluation import eva_auc
from evaluation import eva_aucs
from evaluation import eva_groundtruth
from evaluation import eva_roc_getvalidationpairs
from evaluation import eva_get_validationpairsnlabels_comorbidity

namespaces = ("biological_process", "molecular_function", "cellular_component")


def evaluation_validationpairs():
    vpairsnlabels = eva_get_validationpairsnlabels_comorbidity("data/comorbidity/AllNet3_umlsid.net")

    pathlist = [  # 'outputs/similarity_icod_umls_dgcutoff006_triplet.tsv',
                'outputs/similarity_suntopo_disgenet_dgcutoff006_interactomemaxcc_triplet.tsv',
                'outputs/similarity_funsim_disgenet_dgcutoff006.tsv',
                'outputs/similarity_bog_disgenet_dgcutoff006_triplet.tsv',
                'outputs/similarity_hamaneh_interactomenumls_dgcuff006.tsv',
                'outputs/similarity_experiments_rwr_disgenet_dgcutoff006_interactome.tsv',
                'outputs/similarity_experiments_shortestpath_transformed_less_disgenet_dgcutoff006_interactome.tsv',
                'outputs/similarity_module5_disgenet_dgcutoff006_interactome.tsv'
                ]
    msims = eva_readsims(pathlist)
    print("------------scores and lable---------------------------------------------------")
    vpairsinfo = eva_roc_getvalidationpairs(msims, vpairsnlabels)
    d1 = list(vpairsinfo.keys())[0]
    d2 = list(vpairsinfo[d1].keys())[0]
    methodnames = list(vpairsinfo[d1][d2].keys())
    methodnames.remove('label')
    print("d1\td2\tlabel", end='')
    for m in methodnames:
        print("\t" + m, end='')
    print()
    for d1 in vpairsinfo.keys():
        for d2 in vpairsinfo[d1].keys():
            print(d1 + '\t' + d2 + '\t' + str(vpairsinfo[d1][d2]['label']), end='')
            for m in methodnames:
                print("\t" + str(vpairsinfo[d1][d2][m]), end='')
            print()
    print("------------tpr and fpr--------------------------------------------------------")
    vpairs_tprfpr = eva_tprfpr(eva_ranking(vpairsinfo))
    methodnames = list(vpairs_tprfpr.keys())
    for m in methodnames:
        print(str(m) + "_fpr\t" + m + "_tpr\t", end='')
    print()
    llen = len(vpairs_tprfpr[methodnames[0]])
    for i in range(0, llen):
        for name in methodnames:
            for j in range(0, 2):
                print(str(vpairs_tprfpr[name][i][j]) + "\t", end="")
        print()
    print("------------auc value----------------------------------------------------------")
    vpairs_auc = eva_auc(vpairs_tprfpr)
    for a in vpairs_auc.keys():
        print(a+'\t'+str(vpairs_auc[a]))


def evaluation_groundtruth():
    pathlist = ['outputs/similarity_icod_umls_dgcutoff006_triplet.tsv',
                'outputs/similarity_suntopo_umls_dgcutoff006_triplet.tsv',
                'outputs/similarity_funsim_umls_dgcutoff006.tsv',
                'outputs/similarity_bognew_umls_dgcutoff006_triplet.tsv',
                'outputs/similarity_hamaneh_interactomenumls_dgcuff006.tsv',
                'outputs/similarity_experiments_shortestpath_transformed_less_umls_dgcutoff006.tsv',
                'outputs/similarity_experiments_rwrzrq_umls_dgcutoff006.tsv',
                # 'outputs/similarity_module1_umls_dgcutoff006.tsv',
                'outputs/similarity_module5_umls_dgcutoff006.tsv'
                ]
    gtpathlist = ['outputs/similarity_go_bp_umls_dgcutoff006.tsv',
                  'outputs/similarity_go_cc_umls_dgcutoff006.tsv',
                  'outputs/similarity_go_mf_umls_dgcutoff006.tsv',
                  'outputs/similarity_coexp_umls_dgcutoff006.tsv',
                  'outputs/similarity_symptom_1445umlsid.tsv',
                  'outputs/similarity_symptom_1736umlsid.tsv']

    msims = eva_readsims(pathlist)

    for g in [0, 1, 2, 3, 4]:
        evagtres = eva_groundtruth(msims, gtpathlist[g])
        print(gtpathlist[g]+"----------------------------")
        mnames = list(evagtres.keys())
        sumtemp = {}
        for i in range(0, len(mnames)):
            if i == 0:
                print("rank"+"\t"+mnames[i], end='')
            else:
                print("\t"+mnames[i], end='')
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


def evaluation_70benchmarkset(times=1,
                              bmkpfile="data/benchmarkset_funsim/ground_truth_68_disease_pairs_umlsid.tsv"):
    pathlist = [  # 'outputs/similarity_icod_umls_dgcutoff006_triplet.tsv',
                'outputs/similarity_suntopo_disgenet_dgcutoff006_interactomemaxcc_triplet.tsv',
                'outputs/similarity_funsim_disgenet_dgcutoff006.tsv',
                'outputs/similarity_bog_disgenet_dgcutoff006_triplet.tsv',
                'outputs/similarity_hamaneh_interactomenumls_dgcuff006.tsv',
                'outputs/similarity_experiments_rwr_disgenet_dgcutoff006_interactome.tsv',
                'outputs/similarity_experiments_shortestpath_transformed_less_disgenet_dgcutoff006_interactome.tsv',
                'outputs/similarity_module5_disgenet_dgcutoff006_interactome.tsv',
                # 'outputs/similarity_funsim_rwrsidd.tsv',
                # 'outputs/similarity_hamaneh_rwrsidd_hppinwsl.tsv',
                # 'outputs/similarity_module5_rwrsidd_hppinwsl.tsv',
                # 'outputs/similarity_experiments_rwr_rwrsidd_hppinwsl.tsv',
                # 'outputs/similarity_experiments_shortestpath_transformed_less_rwrsidd_hppinwsl.tsv',
                # 'outputs/similarity_suntopo_rwrsidd_hppinwosl_triplet.tsv',
                # # 'outputs/similarity_icod_rwrsidd_hppinwosl_triplet.tsv',
                # 'outputs/similarity_bognew_rwrsidd_triplet.tsv'
                ]

    benchmarkpairs = read_assos(bmkpfile, header=False)
    stat_assos(benchmarkpairs)
    bmptuple = []
    for p in benchmarkpairs.keys():
        for q in benchmarkpairs[p]:
            bmptuple.append((p, q))
    msims = eva_readsims(pathlist)
    evaress = eva_70benchmarkpairs(msims, bmptuple, times)
    print("------------scores and lable---------------------------------------------------")
    t = 0
    for er in evaress:
        if t < 3:
            t += 1
            print("-----------------time------------------------------------------------------")
            d1 = list(er.keys())[0]
            d2 = list(er[d1].keys())[0]
            methodnames = list(er[d1][d2].keys())
            methodnames.remove('label')
            print("d1\td2\tlabel", end='')
            for m in methodnames:
                print("\t"+m, end='')
            print()
            for d1 in er.keys():
                for d2 in er[d1].keys():
                    print(d1+'\t'+d2+'\t' + str(er[d1][d2]['label']), end='')
                    for m in methodnames:
                        print("\t"+str(er[d1][d2][m]), end='')
                    print()
    print("-ranking disease pairs---------------------------------------------------------")
    t = 0
    for er in evaress:
        if t < 3:
            t += 1
            print("---------------time--------------------------------------------------------")
            ranktemp = eva_ranking(er)
            ms = list(ranktemp.keys())
            for m in ms:
                print(str(m) + "\t\t\t", end="")
            print()
            llen = len(ranktemp[ms[0]])
            for i in range(0, llen):
                for m in ms:
                    for j in range(0, 4):
                        print(str(ranktemp[m][i][j])+"\t", end="")
                print()
    print("-tpr--fpr----------------------------------------------------------------------")
    tpfprs = eva_tprfprs(evaress)
    t = 0
    for tpfpr in tpfprs:
        if t < 3:
            t += 1
            print("-----------------time------------------------------------------------------")
            methodnames = list(tpfpr.keys())
            for m in methodnames:
                print(str(m)+"_fpr\t"+m+"_tpr\t", end='')
            print()
            llen = len(tpfpr[methodnames[0]])
            for i in range(0, llen):
                for name in methodnames:
                    for j in range(0, 2):
                        print(str(tpfpr[name][i][j])+"\t", end="")
                print()
    print("---avg auc values------------------------------------------------------------------")
    aucs = eva_aucs(tpfprs)
    auclen = len(aucs)
    print("replicated times:", auclen)
    mns = list(aucs[0].keys())
    for mn in mns:
        avgauc = 0.0
        for auc in aucs:
            avgauc += auc[mn]
        print(str(mn)+'\t'+str(avgauc / auclen))


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
    disease2gene_entrez = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                             0.06, True, True)
    # disease2gene_entrez = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    print("disease gene assos: ", end='')
    stat_assos(disease2gene_entrez)

    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv", False, False)
    # g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
    #                                        False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    # sims = similarity_module.similarity_cal_module_1(disease2gene_entrez, g)
    # write_sims(sims, "similarity_module1_umls_dcutoff006.tsv")
    # sims = similarity_module.similarity_cal_module_2(disease2gene_entrez, g)
    # write_sims(sims, "similarity_module2_umls_dcutoff006.tsv")
    # sims = similarity_module.similarity_cal_module_3(disease2gene_entrez, g)
    # write_sims(sims, "similarity_module3_umls_dcutoff006.tsv")
    # sims = similarity_module.similarity_cal_module_4(disease2gene_entrez, g)
    # write_sims(sims, "similarity_module4_umls_dcutoff006.tsv")
    sims = similarity_module.similarity_cal_module_5(disease2gene_entrez, g)
    write_sims(sims, "outputs/similarity_module5_disgenet_dgcutoff006_interactome.tsv")


def do2umls_mapping():
    uds = read_umls_disease_info(0)
    do2umls = mapping.doid2umlsid(uds)
    stat_assos(do2umls)
    doid1 = read_one_col("data/ground_truth_70_disease_pairs_doid.tab", 1, True)
    doid2 = read_one_col("data/ground_truth_70_disease_pairs_doid.tab", 2, True)
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
    write_mappings(do2umls_map, "ground_truth_doid2umlsid_46diseases.tsv")
    with open("data/ground_truth_70_disease_pairs_doid.tab", mode='r') as rf:
        next(rf)
        with open("data/ground_truth_68_disease_pairs_umlsid.tsv", mode='w') as wf:
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


def diseaseidmapping_comorbidity():
    uds = read_umls_disease_info(0)
    icd92umls = mapping.icd9cmid2umlsid_alldigit(uds)
    stat_assos(icd92umls)
    pprint(icd92umls)
    icd9s = list(icd92umls.keys())
    for icd9 in icd9s:
        if len(icd9) != 3:
            del icd92umls[icd9]
    stat_assos(icd92umls)
    pprint(icd92umls)
    # write_assos(icd92umls, "data/disgenet/icd9_3digit2umlsid.tsv")


def conv_icd9id2umlsid_allnet():
    icd9id2umlsid = read_mappings("data/disgenet/icd9_3digit2umlsid.tsv")
    with open("data/comorbidity/AllNet3.net", mode='r') as rf:
        with open("data/comorbidity/AllNet3_umlsid.net", mode='w') as wf:
            for line in rf:
                words = line.strip().split("\t")
                if words[0] in icd9id2umlsid.keys() and words[1] in icd9id2umlsid.keys():
                    wf.write(icd9id2umlsid[words[0]] + '\t' + icd9id2umlsid[words[1]])
                    for i in words[2:]:
                        wf.write('\t'+i)
                    wf.write('\n')


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


def disease_module_info():
    g = similarity_module.read_interactome("data/DataS1_interactome_rmslpe.tsv", False, False)
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


def disease_gene_network_hamaneh():
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withselfloop.tab",
                                           False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])

    # disease2gene_entrez = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
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


def get_rwr_input():
    g = similarity_module.read_interactome("data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop.tab",
                                           False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])

    # disease2gene = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
    #                                                   0.06, True, True)
    disease2gene = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
    print("disease gene assos: ", end='')
    stat_assos(disease2gene)

    dgassos_new = {}
    for d in disease2gene.keys():
        dgleft = gvs.intersection(disease2gene[d])
        if len(dgleft) >= 1:
            dgassos_new[d] = dgleft
    print("disease gene assos left: ", end='')
    stat_assos(dgassos_new)

    write_assos(dgassos_new, "data/rwr_dgassos_sidd_inhppinwosl.tsv")
    # write_slist(list(gvs), "data/rwr_bmc_bioinfo/rwr_ppi_hppinwsl_nodes.tsv")
    # write_slist(list(dgassos_new.keys()), "data/rwr_bmc_bioinfo/rwr_dgassos_sidd_inhppinwsl_diseases.tsv")


def experiment():
    g = similarity_module.read_interactome("data/interactome_science/DataS1_interactome.tsv",
                                           False, False)
    print("number of vertices:", g.vcount())
    gvs = set(g.vs['name'])
    disease2gene_entrez = read_all_gene_disease_associations("data/disgenet/all_gene_disease_associations.tsv",
                                                             0.06, True, True)
    # disease2gene_entrez = read_assos("data/rwr_bmc_bioinfo/dg/rwr_dgassos_sidd.tab")
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
    sim_gene2geneset = read_simmatrix("data/test/rwr_geneset2genescore_disgenet_dgcutoff006_interactome.tsv")
    sim_d2d = experiments.sim_geneset2geneset_rwr(disease2gene_entrez, sim_gene2geneset)
    write_sims(sim_d2d, "outputs/similarity_experiments_rwr_disgenet_dgcutoff006_interactome.tsv")
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


if __name__ == "__main__":
    # evaluation_70benchmarkset(1000)
    evaluation_validationpairs()
