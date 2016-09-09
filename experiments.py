#! /usr/bin/env python3
import math
import time

graphlet_weight = [1, 0.838444532557004, 0.838444532557004, 0.838444532557004, 0.743940642316373, 0.676889065114007,
                   0.743940642316373, 0.743940642316373, 0.676889065114007, 0.743940642316373, 0.676889065114007,
                   0.676889065114007, 0.676889065114007, 0.676889065114007, 0.743940642316373, 0.676889065114007,
                   0.582385174873377, 0.624879821261446, 0.676889065114007, 0.624879821261446, 0.582385174873377,
                   0.582385174873377, 0.676889065114007, 0.676889065114007, 0.676889065114007, 0.624879821261446,
                   0.546456463288587, 0.676889065114007, 0.582385174873377, 0.582385174873377, 0.546456463288587,
                   0.676889065114007, 0.582385174873377, 0.582385174873377, 0.582385174873377, 0.624879821261446,
                   0.582385174873377, 0.546456463288587, 0.546456463288587, 0.624879821261446, 0.546456463288587,
                   0.582385174873377, 0.546456463288587, 0.582385174873377, 0.624879821261446, 0.624879821261446,
                   0.582385174873377, 0.515333597671011, 0.546456463288587, 0.582385174873377, 0.582385174873377,
                   0.515333597671011, 0.582385174873377, 0.487881284632746, 0.624879821261446, 0.582385174873377,
                   0.676889065114007, 0.582385174873377, 0.582385174873377, 0.546456463288587, 0.515333597671011,
                   0.582385174873377, 0.582385174873377, 0.515333597671011, 0.546456463288587, 0.582385174873377,
                   0.546456463288587, 0.546456463288587, 0.515333597671011, 0.624879821261446, 0.582385174873377,
                   0.582385174873377, 0.676889065114007]


def read_gene_signature(filepath):
    gene2signature = {}
    with open(filepath, mode='r') as f:
        for line in f:
            words = line.strip().split('\t')
            gene2signature[words[0]] = [int(i) for i in words[1:]]
    return gene2signature


def sim_gene2gene_graphlet(gene2signature):
    weigetsum = sum(graphlet_weight)
    genenames = list(gene2signature.keys())
    result = {}
    for g in genenames:
        result[g] = {}
    for i in range(0, len(genenames)):
        for j in range(i, len(genenames)):
            u = gene2signature[genenames[i]]
            v = gene2signature[genenames[j]]
            sigsim = 0.0
            for s in range(0, 73):
                sigsim += graphlet_weight[s] * abs(math.log(u[s]+1)-math.log(v[s]+1)) / math.log(max(u[s], v[s])+2)
            sigsim = 1 - sigsim / weigetsum
            result[genenames[i]][genenames[j]] = sigsim
            result[genenames[j]][genenames[i]] = sigsim
        print("sim_gene2gene_graphlet(): ", i, "done..")
    return result


def sim_geneset2gene_max(g, gset, sim_gene2gene):
    result = 0.0
    for i in gset:
        if sim_gene2gene[g][i] > result:
            result = sim_gene2gene[g][i]
    return result


def sim_geneset2geneset_rwr(d2g, sim_geneset2gene):
    pass


def sim_geneset2geneset(d2g, sim_gene2gene):
    result = {}
    ds = list(d2g.keys())
    simgenes = set(sim_gene2gene.keys())
    d2g_sim = {}
    d2g_nsim = {}
    for d in ds:
        d2g_sim[d] = d2g[d].intersection(simgenes)
        d2g_nsim[d] = d2g[d].difference(simgenes)
    for i in range(0, len(ds)):
        result[ds[i]] = {}
        now = time.time()
        print("sim_geneset2geneset():", i, "dg len:", len(d2g[ds[i]]))
        for j in range(i, len(ds)):
            simsum = 0.0
            for m in d2g_sim[ds[i]]:
                simsum += sim_geneset2gene_max(m, d2g_sim[ds[j]], sim_gene2gene)
            for m in d2g_sim[ds[j]]:
                simsum += sim_geneset2gene_max(m, d2g_sim[ds[i]], sim_gene2gene)
            for m in d2g_nsim[ds[i]]:
                if m in d2g_nsim[ds[j]]:
                    simsum += 1
            for m in d2g_nsim[ds[j]]:
                if m in d2g_nsim[ds[i]]:
                    simsum += 1
            result[ds[i]][ds[j]] = simsum / (len(d2g[ds[i]]) + len(d2g[ds[j]]))
        print("---------------------------------------cost time:", str(time.time()-now))
    return result
