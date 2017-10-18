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
    df = pd.read_table('../HumanNet_all_uniq.tsv', 
                    sep='\t', header=None, names=['src', 'dest'], index_col=None)
    G = construct_graph(df, directed=False)
    print "constructed graph..."
    #all_nodes = nx.nodes(G)
    #in_degree, out_degree, closeness, between = centrality(G)
    #central_df = pd.DataFrame(index=all_nodes)
    #central_df['closeness'] = pd.Series(closeness)
    #central_df['between'] = pd.Series(between)
    # largest eigenvalue of the adjacency matrix
    central_df = pd.read_table('../HumanNet_centrality.tsv', sep='\t', header=0, index_col=0)
    max_eigenval = max(nx.adjacency_spectrum(G))
    print "max eigen value", max_eigenval
    central_df['katz'] = pd.Series(nx.katz_centrality(G))
    central_df.to_csv('../HumanNet_centrality_updated.tsv', sep='\t')
### END - do_centrality


def do_essentiality():
    df = pd.read_table('../PPI_vidal_all_uniq.tsv', 
            sep='\t', header=None, names=['src', 'dest'], index_col=None)
    G = construct_graph(df, directed=False)
    all_nodes = nx.nodes(G)

    # make essentiality dataframe
    node_df = pd.DataFrame(index=all_nodes)
    for c in ['c1', 'c3', 'c3_cs']:
        node_df.loc[:, c + '_ess'] = None
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
    node_df.to_csv('../ppi_essentiality.tsv', sep='\t')
### END - do_essentiality


def cal_distance(G, ess_df, title=""):
    path_df = pd.DataFrame()
    for col in ess_df.columns:
        print "running for " + col
        # make sets
        ess_set = set(ess_df[ess_df[col] == 1].index)
        non_set = set(ess_df[ess_df[col] == 0].index)
        pred_set = set(ess_df[ess_df['tumor_pred_ess'] == 1].index)
        
        # get paths
        ess_path = pairwise_distance(G, list(ess_set))
        non_path = pairwise_distance(G, list(non_set))
        ess_non_path = pairwise_distance(G, list(ess_set), another_list=list(non_set))
        ess_pred_path = pairwise_distance(G, list(ess_set), another_list=list(pred_set))
        
        # print path ratio
        print "{0}\tess\t{1}".format(col, path_exist_ratio(ess_path))
        print "{0}\tnon\t{1}".format(col, path_exist_ratio(non_path))
        print "{0}\tess_non\t{1}".format(col, path_exist_ratio(ess_non_path))
        print "{0}\tess_pred\t{1}".format(col, path_exist_ratio(ess_pred_path))

        # write path lengths to dataframe
        #path_list = [ess_path, non_path, ess_non_path, ess_pred_path]
        #name_list = ['ess', 'non', 'ess_non', 'ess_tumor_pred']
        #for path, name in zip(path_list, name_list):
        #    path = np.array(path)
        #    path = path[path != np.inf] # exclude infinite values
        #    path_df[col[:-3] + name] = pd.Series(path)
        ### END - for path, name
    ### END - for col

    #path_df.to_csv("../{0}_path_length_essentiality.tsv".format(title), sep='\t')
### END - cal_distance


def main():
    #df = pd.read_table('../HumanNet_all_uniq.txt', 
    #            sep='\t', header=None, names=['src', 'dest'], index_col=None)
    #G = construct_graph(df, directed=False)
    #all_nodes = nx.nodes(G)
    #ess_df = pd.read_table("../HumanNet_essentiality.tsv", sep='\t', header=0, index_col=0)
    #cal_distance(G, ess_df, title='HumanNet')
    do_centrality()
### END - main

if __name__ == "__main__":
    main()

