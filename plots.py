#! /usr/bin/env python3
"""for plots"""
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
