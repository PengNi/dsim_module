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
    plt.plot([0, 1], [0, 1], ls='--')
    plt.title('roc')
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
