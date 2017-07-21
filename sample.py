#!/usr/bin/python
import sys,os
import networkx as nx
import random


G=nx.DiGraph()
f=open("breast_correlation_v2.txt")

for i in f:
        a=i.strip().split('\t')
        G.add_edge(a[0],a[1],corr=a[2],pval=a[3])


def count_successor(node, count_list):
        if len(G.successors(node)) == 0:
                return
        else:
                count_list += G.successors(node)
                for k in G.successors(node):
                        count_successor(k, count_list)

lists=[]
count_successor("A1BG",lists)
print lists



