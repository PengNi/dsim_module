#! /usr/bin/env python3
import similarity_module
from files import read_one_col
from files import read_mappings
from files import stat_assos
from copy import deepcopy
from umls_disease import read_all_gene_disease_associations
from umls_disease import read_umls_disease_info
from similarity_symptom import termnamemap2umlsid


def diseaseidmapping_hsdn():
    meshnames1 = read_one_col("ncomms5212-s5.txt", 1, True)
    meshnames2 = read_one_col("ncomms5212-s5.txt", 2, True)

    meshnames = set(meshnames1).union(meshnames2)
    print("meshnames needed to be map:", len(meshnames))
    umlsdiseases = read_umls_disease_info(0)
    meshname2meshid = read_mappings("MeshTreeHierarchy.csv", True, "\t", 3, 2)

    meshname2umlsid = termnamemap2umlsid(meshnames, umlsdiseases, meshname2meshid)
    print("mapped successfully:", len(meshname2umlsid))


def disease_module_info():
    g = similarity_module.read_interactome("DataS1_interactome_rmslpe.tsv", False, False)
    print("number of vertices:", g.vcount(), "number of edges:", g.ecount())

    disease_genes = read_all_gene_disease_associations("all_gene_disease_associations.tsv", 0)
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


def gene_neighbor_info():
    g = similarity_module.read_interactome("DataS1_interactome.tsv", False, False)
    print("number of vertices:", g.vcount())

    all_disease_genes = read_one_col("curated_gene_disease_associations.tsv", 2, True)
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
    diseaseidmapping_hsdn()
