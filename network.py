import numpy as np
import pandas as pd
import networkx as nx
from sklearn.utils import shuffle


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
	#print df.head()

	# append nodes that exist only in dest (leaf nodes)
	src_nodes = set(np.unique(df.src.values))
	dest_nodes = set(np.unique(df.dest.values))
	leaves = dest_nodes - src_nodes
	leaves.discard(np.nan)
	leaves_len = len(leaves)

	leaves_dic = {'src' : list(leaves),
	               'c1' : np.zeros(leaves_len, dtype=np.int32),
	               'c2' : np.zeros(leaves_len, dtype=np.int32),
	               'c3' : np.zeros(leaves_len, dtype=np.int32),
	               'f'  : np.zeros(leaves_len, dtype=np.int32)}

	for c in range(1, 4):
		for f in range(1, 6):
			leaves_dic['c{0}f{1}'.format(c, f)] = np.zeros(leaves_len, dtype=np.int32)

	leaves_df = pd.DataFrame(leaves_dic)

	df.drop_duplicates(subset='src', inplace=True)
	df.drop(['dest', 'r', 'p'], axis=1, inplace=True)
	df = df.append(leaves_df, ignore_index=True) # append leaf nodes
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
	for the given graph. If the graph is undirected, return empty dictionary for indegree
	and outdegree centrality.

	args:
		G (nx.DiGraph) : input graph

	returns:
		tuple of in-degree, out-degree, closeness, betweenness centrality dictionaries
		that the keys are nodes.

	"""
	in_degree, out_degree = {}, {}

	if G.is_directed():
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


def pairwise_distance(G, node_list, another_list=None):
	""" Calculate pairwise shortest path length for given nodes (node_list) in graph G. 
	If a path exists in either way, use that path length. Every edge is considered as length 1.
	If no path exists, use np.inf as length. If another_list is given, calculate distances
	between the nodes in node_list and the nodes in another_list.

	args:
		G (NetworkX graph) : directed or undirected graph
		node_list (list-like obj) : list of nodes to calculate pairwise shortest length
		another_list (list-like obj) : default None.

	returns:
		list of shortest path lengts in unrolled fasion

	"""

	if another_list is None:
		ls = list(node_list)
	else:
		ls = list(set(node_list) & set(another_list))
		node_list = list(node_list)
		another_list = list(another_list)

	# calculate within-list pairwise distance
	list_len = len(ls)
	path_len_list = []
	for i in range(list_len):
		for k in range(i+1, list_len):
			n1 = ls[i]
			n2 = ls[k]
			try:
				length = nx.shortest_path_length(G, source=n1, target=n2, weight=None)
			except nx.NetworkXNoPath:
				try:
					length = nx.shortest_path_length(G, source=n2, target=n1, weight=None)
				except nx.NetworkXNoPath:
					length = np.inf
			path_len_list.append(length)
		### END - for k
	### END - for i

	# calculate two seperate list pairwise distance
	# for instance, if list 1 comprises of A, B, C and list 2 a, b, C
	# calculate A-a, B-a, C-a, A-b, B-b, C-b distance
	if another_list is not None:
		for n1 in node_list:
			for n2 in another_list:
				# if both are in two list's intersection, distance is already calculated
				if (n1 in ls) and (n2 in ls):
					continue
				try:
					length = nx.shortest_path_length(G, source=n1, target=n2, weight=None)
				except nx.NetworkXNoPath:
					try:
						length = nx.shortest_path_length(G, source=n2, target=n1, weight=None)
					except nx.NetworkXNoPath:
						length = np.inf
				path_len_list.append(length)
			### END - for n2
		### END - for n1
	### END - if
	return path_len_list
### pairwise_distance

def cal_distance(G, node_df, n_cases=4, n_features=6, random=False):
	""" Calcualtes pairwise distance fo each cluster.
	args:
		G (NetworkX Graph) : graph containing nodes in node_df
		node_df (Pandas DataFrame) : node_df should contain src, cf clusters (ex. c1f1, c1f2) columns
		n_cases (int) : default=4. number of cases for clusters
		n_features (int) : default=6, number of features for clusters
		random (bool) : if True, select random nodes which the count is equal to 
		                 the corresponding cluster
	
	returns:
		a dataframe that contains pairwise distance for each cluster
	"""
	df = node_df
	all_nodes = df.src.values
	len_df = pd.DataFrame()
	for c in range(1, n_cases):
		for f in range(1, n_features):
			cluster = df[df['c{0}f{1}'.format(c, f)] == 1]['src'].values
			if not random:
				len_df['c{0}f{1}'.format(c, f)] = pd.Series(pairwise_distance(G, list(cluster)))
			else:
				random_nodes = list(shuffle(all_nodes, n_samples=cluster.shape[0]))
				len_df['rand_c{0}f{1}'.format(c, f)] = pd.Series(pairwise_distance(G, random_nodes))
		### END - for f
	### END - for c
	return len_df
### END - compare_distance


def path_exist_ratio(len_list):
	len_arr = np.array(len_list)
	return len_arr[len_arr != np.inf].shape[0] / float(len(len_list))
### END - path_exist_ratio


def construct_ratio_df(G, node_df, n_cases=4, n_features=6, n_repeat=1000, verbose=False):
	df = node_df
	all_nodes = df.src.values
	ratio_df = pd.DataFrame()
	if verbose:
		print "construct_ratio_df..."
	for c in range(1, n_cases):
		for f in range(1, n_features):
			if verbose:
				print "running for c{0}f{1}...".format(c, f)
			ratio_list = []
			for n in range(n_repeat):
				cluster = df[df['c{0}f{1}'.format(c, f)] == 1]['src'].values
				random_nodes = list(shuffle(all_nodes, n_samples=cluster.shape[0], random_state=n))
				len_list = pairwise_distance(G, random_nodes)
				ratio_list.append(path_exist_ratio(len_list))
			### END - for n
			ratio_df['rand_c{0}f{1}'.format(c, f)] = pd.Series(ratio_list)
		### END - for f
	### END - for c
	return ratio_df
### END - construct_ratio_df


def set_essentiality(G):
	""" Sets essentiality attribute to nodes. The attribute value is either 0, 1, None.

	args:
		G (NetworkX Graph)

	"""
	df = pd.read_table('../all_node_essentiality.tsv', sep='\t', header=0, index_col=0)
	for ess in df.columns:
		nx.set_node_attributes(G, ess, df[ess].T.to_dict())
### END - set_essentiality


def ess_path_ratio():
	""" Calculates path existence ratio for essential and non-essential gene group.
	Both directed and undirected paths are calculated.
	"""
	ess_df = pd.read_table('../all_node_essentiality.tsv', sep='\t', header=0, index_col=0)
	G_direct = construct_graph(set_cluster=False)
	G_undirect = G_direct.to_undirected()
	#set_essentiality(G)
	for direct in ['directed', 'undirected']:
		for col in ess_df.columns:
			if direct == 'undirected':
				G = G_undirect
			else:
				G = G_direct
			print "running {0} {1}...".format(direct, col)
			ess_path = pairwise_distance(G, ess_df[ess_df[col] == 1].index)
			non_path = pairwise_distance(G, ess_df[ess_df[col] == 0].index)
			print "{0}\tess\t{1}\t{2}".format(col, direct, path_exist_ratio(ess_path))
			print "{0}\tnon\t{1}\t{2}".format(col, direct, path_exist_ratio(non_path))
		### END - for col
	### END - for direct
### END - ess_path_ratio


def main():
	#G = construct_graph()
	#print len(G)
	#in_degree, out_degree, closeness, betweenness = centrality(G)
	#print len(closeness.keys())
	#node_df = construct_node_df(G)
	#node_df['in_central'] = pd.Series(in_degree)
	#node_df['out_central'] = pd.Series(out_degree)
	#node_df['close'] = pd.Series(closeness)
	#node_df['between'] = pd.Series(betweenness)
	#node_df.to_csv('breast_cancer_all_node_df', sep='\t')
	ess_df = pd.read_table('../all_node_essentiality.tsv', sep='\t', header=0, index_col=0)
	G_direct = construct_graph(set_cluster=False)
	#G_undirect = G_direct.to_undirected()
	path_df = pd.read_table('../directed_path_length_essentiality.tsv', sep='\t', header=0, index_col=0)
	#set_essentiality(G)
	
	# make sets
	ess_set = set(ess_df[ess_df['c3_cs_ess'] == 1].index)
	non_set = set(ess_df[ess_df['c3_cs_ess'] == 0].index)
	pred_set = set(ess_df[ess_df['tumor_pred_ess'] == 1].index)

	print "for ess_non"
	print "only in ess " + str(len(ess_set - non_set))
	print "only in non " + str(len(non_set - ess_set))
	print "intersection " + str(len(ess_set & non_set))

	print "for ess_pred"
	print "only in ess " + str(len(ess_set - pred_set))
	print "only in pred " + str(len(pred_set - ess_set))
	print "intersection " + str(len(ess_set & pred_set))
	
	#print "running for {0}...".format(col)
	ess_non_path = pairwise_distance(G_direct, list(ess_set), another_list=list(non_set))
	ess_non_path = np.array(ess_non_path)
	print "ess_non_path existence ratio " + str(path_exist_ratio(ess_non_path))
	ess_non_path = ess_non_path[ess_non_path != np.inf]  # exclude infinite values
	path_df['c3_cs_ess_non'] = pd.Series(ess_non_path)
	
	ess_pred_path = pairwise_distance(G_direct, list(ess_set), another_list=list(pred_set))
	ess_pred_path = np.array(ess_pred_path)
	print "ess_pred_path existence ratio " + str(path_exist_ratio(ess_pred_path))
	ess_pred_path = ess_pred_path[ess_pred_path != np.inf]  # exclude infinite values
	path_df['c3_cs_ess_tumor_pred'] = pd.Series(ess_pred_path)
	
	path_df.to_csv('../directed_path_length_essentiality.tsv', sep='\t')
	#if 'pred' not in col:
	#	non_path = pairwise_distance(G_direct, ess_df[ess_df[col] == 0].index)
	#	non_path = np.array(non_path)
		#non_path = non_path[non_path != np.inf]  # exclude infinite values
	#	path_df[col[:-3] + 'non'] = pd.Series(non_path)
	### END - for col
	#path_df.to_csv('../undirected_path_length_essentiality.tsv', sep='\t')
	#ratio_df = construct_ratio_df(G, df, verbose=True)
	#ratio_df.to_csv('../rand_path_exist_ratio', sep='\t')
	#len_df = compare_distance()
### END - main

#def mark_essentiality(gene_df, ess_file, non_file):

if __name__ == "__main__":
	main()




