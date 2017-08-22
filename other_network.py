import numpy as np
import pandas as pd
import networkx as nx
from sklearn.utils import shuffle
from network import *

def construct_graph(df, directed=True):
    """ Constructs networkx graph from network DataFrame.

    args:
        df (pandas DataFrame) : edge information. D
                                Dataframe must be in the format src/dest/attr1/...
        directed (bool) : deafuat True. If True, use DiGraph, otherwise use Graph.                

    returns:
        networkx graph

    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    # add edges to the graph
    for idx, row in df.iterrows():
        attr_dict = row.drop(['src', 'dest']).to_dict()
        G.add_edge(row.src, row.dest, attr_dict=attr_dict)
    ### END - for

    return G
### END - construct_graph

def do_centrality():
    df = pd.read_table('../HumanNet_all_uniq.txt', 
                    sep='\t', header=None, names=['src', 'dest'], index_col=None)
    G = construct_graph(df, directed=False)
    all_nodes = nx.nodes(G)
    in_degree, out_degree, closeness, between = centrality(G)
    central_df = pd.DataFrame(index=all_nodes)
    central_df['closeness'] = pd.Series(closeness)
    central_df['between'] = pd.Series(between)
    central_df.to_csv('../HumanNet_centrality.tsv', sep='\t')
### END - do_centrality

def main():
    df = pd.read_table('../HumanNet_all_uniq.txt', 
                sep='\t', header=None, names=['src', 'dest'], index_col=None)
    G = construct_graph(df, directed=False)
    all_nodes = nx.nodes(G)

    # centrality
    in_degree, out_degree, closeness, between = centrality(G)
    central_df = pd.DataFrame(index=all_nodes)
    central_df['closeness'] = pd.Series(closeness)
    central_df['between'] = pd.Series(between)
    central_df.to_csv('../HumanNet_centrality.tsv', sep='\t')
    
    # make essentiality dataframe
    node_df = pd.DataFrame(index=all_nodes)
    for c in ['c1', 'c3', 'c3_cs']:
        node_df.loc[:, c + 'ess'] = None
        for ess in ['Essen', 'Non']:
            with open('../Input_{0}_{1}.txt'.format(ess, c), 'r') as f:
                ls = f.read().split('\n')

                for gene in ls:
                    if gene in node_df.index:
                        if ess == 'Essen':
                            node_df.loc[gene, c + '_ess'] = 1
                        else:
                            node_df.loc[gene, c + '_ess'] = 0
                ### END - for gene
            ### END - with open
        ### END - for ess
    ### END - for c

    node_df.loc[:, 'normal_pred_ess'] = 0
    node_df.loc[:, 'tumor_pred_ess'] = 0
    for ess in ['N', 'T']:
        with open('../Gene_c3_top_case_3x_relu_1_{0}_0.5_3.tsv'.format(ess), 'r') as f:
            ls = f.read().split('\n')
            for gene in ls:
                if gene in node_df.index:
                    if ess == 'T':
                        node_df.loc[gene, 'tumor_pred_ess'] = 1
                    else:
                        node_df.loc[gene, 'normal_pred_ess'] = 1
                ### END - for gene
            ### END - with open
        ### END - for ess
    ### END - for c
                        
    node_df.apply(pd.to_numeric, errors='coerce', downcast='integer')
    node_df.to_csv('../HumanNet_essentiality.tsv', sep='\t')


if __name__ == "__main__":
    main()

