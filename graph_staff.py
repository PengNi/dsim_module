#! /usr/bin/env python3
from igraph import Graph


def get_lcc(graphfile):
    with open(graphfile, 'r') as f:
        g = Graph.Read_Ncol(f, names=True, weights=False, directed=False)
    print('nodes:', len(g.vs))
    print('edges:', len(g.es))
    glcc = g.clusters(mode='WEAK').giant()
    print('lcc nodes:', len(glcc.vs))
    print('lcc edges:', len(glcc.es))
    return glcc


def write_lcc():
    glcc = get_lcc('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop.tab')
    # glccnodes = glcc.vs['name']
    # glccnodes2loc = {}
    # for i in range(0, len(glccnodes)):
    #     glccnodes2loc[glccnodes[i]] = i
    # with open('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop_maxcc_nodes_numbered.tab', 'w') as f:
    #     for gln in glccnodes:
    #         f.write(gln+'\t'+str(glccnodes2loc[gln])+'\n')
    #
    # with open('data/rwr_bmc_bioinfo/ppi/rwr_ppi_hppin_withoutselfloop_maxcc_edges_numbered.tab', 'w') as f:
    #     gles = glcc.es
    #     f.write(str(len(glcc.vs))+' '+str(len(glcc.es))+'\n')
    #     for e in gles:
    #         ns = glcc.vs[e.source]['name']
    #         nt = glcc.vs[e.target]['name']
    #         f.write(str(glccnodes2loc[ns])+' '+str(glccnodes2loc[nt])+'\n')


if __name__ == '__main__':
    write_lcc()
