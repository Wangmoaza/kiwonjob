import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def print_stat(network_df):
	""" This method prints statistics of the network based on pandas DataFrame.
	Appropriate data frame should have fields src, dest.src

	args:
		network_df (pd.DataFrame) : data frame containig network information

	"""
    src_nodes = set(np.unique(network_df.src.values))
    dest_nodes = set(np.unique(network_df.dest.values))
    print "number of all nodes: " + str(len(src_nodes | dest_nodes))
    print "number of src nodes: " + str(np.unique(network_df.src.values).shape[0])
    print "number of dest nodes: " + str(np.unique(network_df.dest.values).shape[0])
    print "number of src_nodes - dest_nodes: " + str(len(src_nodes - dest_nodes))
    print "number of dest_nodes - src_nodes: " + str(len(dest_nodes - src_nodes))
    print "number of all edges: " + str(network_df.shape[0])
### END - print_stat


def mark_clusters():
	""" Make columns that mark the clusters (c1, c2, c3, f) and their intersections (c1f1, c1f2, ...)
	and return the new network DataFrame

	returns:
		network DataFrame containing cluster information

	"""
    network = pd.read_table("../breast_cancer_corr_essential_marked", sep='\t', header=0, index_col=0)

    # mark case (c) clusters
    for c in range(1, 4):
        df_c = pd.read_table('../clusters/mini_c{0}.txt'.format(c), sep='\t', header=None, names=['gene'])
        df_c = df_c['gene'].str.split(':', expand=True)
        df_c.columns = ['gene', 'cell']
        print "df_c" + str(c) 
        
        
        diff_set = set(df_c.gene.values) - set(network.src.values)
        
        if len(diff_set) != 0:
            for gene in diff_set:
                row = pd.Series(np.zeros((network.shape[1], ), dtype=object), index=network.columns)
                row['src'] = gene
                row['dest'] = None
                network = network.append(row, ignore_index=True)
            ### END - for
        ### END - if
        
        print "difference " + str(len(diff_set))
        
        c_col = np.zeros(network.shape[0], dtype=np.int32)
        for idx, row in network.iterrows():
            if row.src in df_c.gene.values:
                c_col[idx] = 1
        ### END - for

        network['c{0}'.format(c)] = c_col
    ### END - for
    
    # mark feature (f) clusters
    network['f'] = np.zeros(network.shape[0], dtype=np.int32)  # make new column

    for f in range(1, 6):
        df_f = pd.read_table('../clusters/mini_f{0}.txt'.format(f), sep='\t', header=None, names=['gene', 'cluster'])
        diff_set = set(df_f.gene.values) - set(network.src.values)
        
        if len(diff_set) != 0:
            for gene in diff_set:
                row = pd.Series(np.zeros((network.shape[1], ), dtype=object), index=network.columns)
                row['src'] = gene
                row['dest'] = None
                network = network.append(row, ignore_index=True)
        	### END - for
        ### END - if
        print "difference " + str(len(diff_set))
        
        for idx, row in network.iterrows():
            if row.src in df_f.gene.values:
                network.at[idx, 'f'] = f
        ### END - for
    ### END - for
    
    # mark intersection of case and feature clusters
    for c in range(1, 4):
        for f in range(1, 6):
            idx = network[(network['c{0}'.format(c)] == 1) | (network['f'] == f)].index
            col = np.zeros(network.shape[0], dtype=np.int32)
            col[idx] = 1
            network['c{0}f{1}'.format(c, f)] = col
        ### END - for
    ### END - for
    return network
### END - mark_clusters


def construct_node_df(G):
	""" Constructs node information DataFrame from edge information DataFrame.

	args:
		G (nx.DiGraph): graph corresponding to the data frame for removing the isolated nodes

	returns:
		pd.DataFrame containing nnetwork node information

	"""
	df = pd.read_table('../breast_cancer_cluster_marked_cf', sep='\t', header=0, index_col=0)
	df.drop_duplicates(subset='src', inplace=True)
	df.drop(['dest', 'r', 'p'], axis=1, inplace=True)
	df.set_index('src', inplace=True)
	# remove node that are not in the main connected component
	to_be_removed = isolated_nodes(G)
	df.drop(list(to_be_removed), axis=0, inplace=True)
	return df
### END - construct_node_df


def construct_graph(set_cluster=True):
	""" Constructs networkx graph from network DataFrame. 
	If set_cluster is True, cluster information is added to the nodes.

	args:
		set_cluster (bool) : whether to put cluster information to nodes

	returns:
		networkx constructed from edge information data frame

	"""

	# network_df should also have clustered information
	df = pd.read_table('../breast_cancer_cluster_marked_cf', sep='\t', header=0, index_col=0)
	ess_dic, c1_dic, c2_dic, c3_dic, f_dic = {}, {}, {}, {}, {}
	G = nx.DiGraph()
	# add edges to the graph
	for idx, row in df.iterrows():
		G.add_edge(row.src, row.dest, r=row.r, p=row.p)
		if set_cluster:
			ess_dic[row.src] = row.src_ess
			c1_dic[row.src] = row.c1
			c2_dic[row.src] = row.c2
			c3_dic[row.src] = row.c3
			f_dic[row.src] = row.f
	### END - for
	if set_cluster:
		nx.set_node_attributes(G, 'ess', ess_dic)
		nx.set_node_attributes(G, 'c1', c1_dic)
		nx.set_node_attributes(G, 'c2', c2_dic)
		nx.set_node_attributes(G, 'c3', c3_dic)
		nx.set_node_attributes(G, 'f', f_dic)

	try:
		G.remove_node(np.nan)
	except NetworkXError:
		print "node requested to remove does not exists"

	# remove node that are not in the main conneted_component
	to_be_removed = isolated_nodes(G)
	for n in to_be_removed:
		G.remove_node(n)

	return G
### END - construct_graph


def isolated_nodes(G):
	""" Returns the nodes that are not in the largest weakly connected component.

	args:
		G (nx.DiGraph) : graph

	returns:
		set of nodes that are not in the largest weakly connected component

	"""
	to_be_removed = set()
	idx = 0
	for comp in nx.weakly_connected_components(G):
		if idx != 0:
			to_be_removed.union(comp)
		idx += 1
	return to_be_removed
### END - isolated_nodes


def centrality(G):
	""" Calculates the in-degree, out-degree, closeness, betweenness centrality
	for the given graph.

	args:
		G (nx.DiGraph) : input graph

	returns:
		tuple of in-degree, out-degree, closeness, betweenness centrality dictionaries
		that the keys are nodes.

	"""
	print "calculating in_degree centrality..."
	in_degree = nx.in_degree_centrality(G)
	print "calculating out_degree centrality..."
	out_degree = nx.out_degree_centrality(G)
	print "calculating closeness centrality..."
	closeness = nx.closeness_centrality(G)
	print "calculating betweenness centrality..."
	betweenness = nx.betweenness_centrality(G)

	return in_degree, out_degree, closeness, betweenness
### END - centrality

def main():
	G = construct_graph()
	in_degree, out_degree, closeness, betweenness = centrality(G)
	node_df = construct_node_df(G)
	node_df['in_central'] = pd.Series(in_degree)
	node_df['out_central'] = pd.Series(out_degree)
	node_df['close'] = pd.Series(closeness)
	node_df['between'] = pd.Series(betweenness)
	node_df.to_csv('breast_cancer_node_df', sep='\t')
### END - main

if __name__ == "__main__":
	main()




