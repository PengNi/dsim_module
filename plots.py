#! /usr/bin/env python3
"""for plots"""
import files
import matplotlib
matplotlib.use('Agg')


def plot_roc(tprfpr, names2legends, names2aucs, savepath='roc.png', fformat='png'):
    """
    plot rocs based on fpr and tpr
    :param tprfpr: from method eva_tprfpr()
    :param names2legends: dict, names to legends
    :param names2aucs: dict, names to auc values
    :param savepath: path of saved file
    :param fformat: format of saved file
    :return: a plot
    """
    from matplotlib import pyplot as plt
    for tf in tprfpr.keys():
        line = tprfpr[tf]
        fprs = []
        tprs = []
        for dot in line:
            fprs.append(dot[0])
            tprs.append(dot[1])
        plt.plot(fprs, tprs, ls='-',
                 label=str(names2legends[tf]) + '(' + "%.3f" % names2aucs[tf] + ')')
    plt.plot([0, 1], [0, 1], ls='--', color='k')
    plt.title('ROC')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc='lower right')
    plt.savefig(savepath, format=fformat)
    plt.clf()


def plot_simshist(sims, title, bins=10, color='red', savepath='hist.png', fformat='png'):
    """
    plot histgram of sim values
    :param sims:
    :param title:
    :param bins:
    :param color:
    :param savepath:
    :param fformat:
    :return:
    """
    calpairs = {}
    for k in sims.keys():
        calpairs[k] = set()
        calpairs[k].update(set(sims[k].keys()))
        calpairs[k].discard(k)
    files.stat_assos(calpairs)

    simvalues = []
    for k1 in calpairs.keys():
        for k2 in calpairs[k1]:
            simvalues.append(sims[k1][k2])
    from matplotlib import pyplot as plt
    plt.hist(simvalues, bins=bins, color=color)
    plt.title(title)
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.yscale('log', nonposy='clip')
    plt.savefig(str(title + savepath), format=fformat)
    plt.clf()


def vennstats_twosets_ppi(ppi1, ppi2):
    """

    :param ppi1: ppi1, dict (key-value: string-set<string>)
    :param ppi2: ppi2, dict (key-value: string-set<string>)
    :return:
    """
    cppi = {}
    for p1 in ppi1.keys():
        if p1 not in cppi.keys():
            cppi[p1] = {}
        for p2 in ppi1[p1]:
            cppi[p1][p2] = 0
    for p1 in ppi2.keys():
        for p2 in ppi2[p1]:
            if not ((p1 in cppi.keys() and p2 in cppi[p1].keys()) or
                    (p2 in cppi.keys() and p1 in cppi[p2].keys())):
                if p1 not in cppi.keys():
                    cppi[p1] = {}
                cppi[p1][p2] = 0
    # ---scoring---
    for p1 in cppi.keys():
        for p2 in cppi[p1].keys():
            if (p1 in ppi1.keys() and p2 in ppi1[p1]) or (p2 in ppi1.keys() and p1 in ppi1[p2]):
                cppi[p1][p2] += 1
            if (p1 in ppi2.keys() and p2 in ppi2[p1]) or (p2 in ppi2.keys() and p1 in ppi2[p2]):
                cppi[p1][p2] += 2
    # ---classing---
    unique1, unique2, intersect12 = [], [], []
    unique1g, unique2g, intersect12g = set(), set(), set()
    for p1 in cppi.keys():
        for p2 in cppi[p1].keys():
            if cppi[p1][p2] == 1:
                unique1.append((p1, p2))
                unique1g.add(p1)
                unique1g.add(p2)
            elif cppi[p1][p2] == 2:
                unique2.append((p1, p2))
                unique2g.add(p1)
                unique2g.add(p2)
            else:
                intersect12.append((p1, p2))
                intersect12g.add(p1)
                intersect12g.add(p2)
    print('\tppi1\tppi2\tintersection')
    print('nodes\t' + str(len(unique1g)) + '\t' + str(len(unique2g)) + '\t' + str(len(intersect12g)))
    print('assos\t' + str(len(unique1)) + '\t' + str(len(unique2)) + '\t' + str(len(intersect12)))


def vennstats_twosets_dgassos(d2g1, d2g2):
    """

    :param d2g1: d2g assos 1
    :param d2g2: d2g assos 2
    :return: venn stats
    """
    alld2g = {}
    for d in d2g1.keys():
        alld2g[d] = {}
        for g in d2g1[d]:
            alld2g[d][g] = 0
    for d in d2g2.keys():
        for g in d2g2[d]:
            if d not in alld2g.keys():
                alld2g[d] = {}
            alld2g[d][g] = 0
    for d in alld2g.keys():
        for g in alld2g[d].keys():
            if d in d2g1.keys() and g in d2g1[d]:
                alld2g[d][g] += 1
            if d in d2g2.keys() and g in d2g2[d]:
                alld2g[d][g] += 2
    ud1, ud2, id12 = set(), set(), set()
    ug1, ug2, ig12 = set(), set(), set()
    ua1, ua2, ia12 = set(), set(), set()
    for d in alld2g.keys():
        for g in alld2g[d].keys():
            if alld2g[d][g] == 1:
                ud1.add(d)
                ug1.add(g)
                ua1.add((d, g))
            elif alld2g[d][g] == 2:
                ud2.add(d)
                ug2.add(g)
                ua2.add((d, g))
            else:
                id12.add(d)
                ig12.add(g)
                ia12.add((d, g))
    print('\td2g1\td2g2\tintersection')
    print('d\t' + str(len(ud1)) + '\t' + str(len(ud2)) + '\t' + str(len(id12)))
    print('g\t' + str(len(ug1)) + '\t' + str(len(ug2)) + '\t' + str(len(ig12)))
    print('a\t' + str(len(ua1)) + '\t' + str(len(ua2)) + '\t' + str(len(ia12)))
