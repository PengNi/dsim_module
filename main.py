#! /usr/bin/env python3
import similarity_module
from files import read_one_col
from files import read_mappings
from files import read_assos
from files import stat_assos
from files import stat_maps
from files import write_assos
from files import write_mappings
from files import write_sims
from copy import deepcopy
from umls_disease import read_all_gene_disease_associations
from umls_disease import read_umls_disease_info
from mapping import termname2umlsid
from mapping import icd9cmid2umlsid_3digit
from mapping import doid2umlsid
from similarity_coexpression import read_probeid_expfile
from similarity_coexpression import probeexp2geneexp
from similarity_coexpression import diseases_similarity_coexp
from gene_ontology import read_go_annotation_file
from gene_ontology import GeneOntology
from gene_ontology import invert_dict
from gene_ontology import add_implicit_annotations
from similarity_go import diseases_similarity_go
from evaluation import eva_70benchmarkpairs


namespaces = ("biological_process", "molecular_function", "cellular_component")


def evaluation_70benchmarkset():
    pathlist = ['similarity_icod_umls_dcutoff006_triplet.tsv',
                'similarity_suntopo_umls_dcutoff006_triplet.tsv',
                'similarity_funsim_umls_dcutoff006.tsv',
                'similarity_bog_umls_dcutoff006_triplet.tsv',
                'simialrity_module1_umls_dcutoff006.tsv']

    benchmarkpairs = read_assos("data/ground_truth_68_disease_pairs_umlsid.tsv")
    stat_assos(benchmarkpairs)
    bmptuple = []
    for p in benchmarkpairs.keys():
        for q in benchmarkpairs[p]:
            bmptuple.append((p, q))
    evaress = eva_70benchmarkpairs(pathlist, bmptuple, 5)
    for er in evaress:
        print("-----------------time-----------------------------------------------------")
        print("d1\td2\ticod\tsuntopo\tfunsim\tbog\tmodule1\tlabel")
        for d1 in er.keys():
            for d2 in er[d1].keys():
                print(d1+"\t"+d2+'\t'+str(er[d1][d2][pathlist[0]])+'\t' +
                      str(er[d1][d2][pathlist[1]])+'\t' +
                      str(er[d1][d2][pathlist[2]])+'\t' +
                      str(er[d1][d2][pathlist[3]])+'\t' +
                      str(er[d1][d2][pathlist[4]])+'\t' +
                      str(er[d1][d2]['label']))


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
    disease2gene_entrez = read_all_gene_disease_associations("data/all_gene_disease_associations.tsv",
                                                             0.06, True, True)
    print("disease gene assos: ", end='')
    stat_assos(disease2gene_entrez)

    g = similarity_module.read_interactome("data/DataS1_interactome_rmslpe.tsv", False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    sims = similarity_module.similarity_cal_module_1(disease2gene_entrez, g)
    write_sims(sims, "simialrity_module1_umls_dcutoff006.tsv")


def do2umls_mapping():
    uds = read_umls_disease_info(0)
    do2umls = doid2umlsid(uds)
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
    icd92umls = icd9cmid2umlsid_3digit(uds)
    stat_assos(icd92umls)
    icd9s = list(icd92umls.keys())
    for icd9 in icd9s:
        if len(icd92umls[icd9]) > 1:
            del icd92umls[icd9]
    stat_assos(icd92umls)


def diseaseidmapping_hsdn():
    meshnames1 = read_one_col("ncomms5212-s5.txt", 1, True)
    meshnames2 = read_one_col("ncomms5212-s5.txt", 2, True)

    meshnames = set(meshnames1).union(meshnames2)
    print("meshnames needed to be map:", len(meshnames))
    umlsdiseases = read_umls_disease_info(0)
    meshname2meshid = read_mappings("data/MeshTreeHierarchy.csv", True, "\t", 3, 2)

    meshname2umlsid = termname2umlsid(meshnames, umlsdiseases, meshname2meshid)
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


if __name__ == "__main__":
    evaluation_70benchmarkset()
