#! /usr/bin/env python3
from igraph import Graph


def read_interactome(interactomefile, weights, directed):
    """
    read a file into an igraph object using igraph read_Ncol method
    :param interactomefile: a file path of a "ncol" file contains interactome
    :param weights: same as "weights" in igraph object construction method:
    True, False, "auto", "if_present"
    :param directed: True or False
    :return: an igraph object
    """
    with open(interactomefile, mode='r') as f:
        g = Graph.Read_Ncol(f, names=True, weights=weights, directed=directed)
        return g


def density(dgassos, graph):
    """
    given a dict which contains a list of node groups, and an igraph object,
    return a dict which gives density of these node groups based on this igraph object.
    :param dgassos: a dict object which keys are module names and values are module nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are their density scores
    """
    ddensityscore = {}
    for name in dgassos:
        nodes = graph.vs.select(name_in=dgassos[name])
        ddensityscore[name] = graph.subgraph(nodes).density()
    return ddensityscore


def lambda_module(dgassos, graph):
    """
    get lambda value from the "lambda-module"
    :param dgassos: a dict object which keys are module names and values are module nodes sets
    :param graph: an igraph object
    :return: a dict which keys are module names and values are their lambda values
    """
    dlambda = {}
    for name in dgassos:
        nodes = graph.vs.select(name_in=dgassos[name])
        din = sum(graph.subgraph(nodes).degree())
        dout = sum(graph.degree(vertices=nodes))
        dout = dout-din
        dlambda[name] = din/dout
    return dlambda
