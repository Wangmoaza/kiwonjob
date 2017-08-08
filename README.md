# kiwonjob
Helping Kiwon

## file description

#### files from Kiwon

* *files in clusters folder* : gene list of each clusters using hierarchical clustering of the heat map. c stands for case, f stands for features. mini_ files are artificial clusters made by Kiwon.
* *breast_correlation_v2.txt* : edge information of network format src/dest/corr/p-value
* *deep_input_c3_pert_v2LogN* : cut off 3
* *deep_input_c3_pert_v2LogN_top* : contains only features with high variance
* *deep_input_c3_pert_v2LogN_top_cs_3x* : plus case selection
* *Gene_c3_top_case_3x_relu_1_T_0.5_3.tsv* : normal predicted essential genes
* *Gene_c3_top_case_3x_relu_1_N_0.5_3.tsv* : tumor predicted essential genes

#### files made by me

* *gene_essential_cnt* : gene essentiality is based on deep_input_c3_pert_v2LogN file.
*  *breast_cancer_node_df* : node data frame of network. cluster marked, centrality calculated.
*  *breast_cancer_node_df_central_ess_marked* : src / centralities... / c1_ess/ c3_ess / c3_cs_ess
*  *breast_cancer_all_nodes_df* : breast_cancer_nodes_df does not include leaf nodes in the network. This file does.
* *rand_path_exist_ttest_result* : one-sample t-test result of randomly chosen nodes (n_samples=1000) corresponding to cluster nodes
* *all_node_centrality.tsv* : in-degree, out-degree, closeness, betweeness centrality of all nodes
* *all_node_essentiality.tsv* : c1, c3, c3_cs, tumor predicted, normal predicted essentiality of all nodes in the network
* *directed_path_length_essentiality.tsv* : directed pairwise distances within group (c1, c3, c3_cs, normal predicted, tumor predicted). Excluded infinite values.
* *directed_path_length_essentiality_for_r.tsv* : converted dataframe in the layout group, distance. For easy use in R
* *undirected_path_length_essentiality.tsv* : directed pairwise distances within group (c1, c3, c3_cs, normal predicted, tumor predicted). Included infinite values.
* *undirected_path_length_essentiality_for_r.tsv* : equivelent of directed version.

## lab note by date

#### 2017-07-21

Visualized cluster c2f4 with cytoscape

#### 2017-07-24

* labeled clusters to breast_cancer_cluster_labeled_cf file.
* manually visualizing each clusters using cytoscape

#### 2017-07-25

* centrality calculation uisng networkx functions
* data saved in file breast_cancer_node_df

#### 2017-07-26

* drew box plot for centralities in two groups (essential, non-essential)

#### 2017-07-27

* centrality group divided by c1, c3, c3_cs cutoff (by Input data files)

#### 2017-07-31

* pairwise distance calculation something wrong. need to check random selection along with cluster nodes
* node_df missing many dest nodes. need to check this again.

#### 2017-08-02

* fixed node_df to include all nodes (leaf nodes too). saved in file breast_cancer_all_nodes_df
* fixed pairwise distance calculation

#### 2017-08-03

* one-sample t-test on path existance ratio of randomly chosen nodes corresponding to each cluster. All significantly different. Saved in the file rand_path_exist_ttest_result.txt.

#### 2017-08-07

* calculating path existence ratio grouped by essentiality (c1, c3, c3_ess, normal_predicted, tumor_predicted)

#### 2017-08-08

* calculated undirected path length within each group.
* calculated directed path length within each group.
* plotted (boxplot, frequency polygon) the results. 

#### TODO

* clueGO gene ontology plugin for cytoscape: map clustered or essential genes to gene ontology, or KEGG pathway
    - find a module for mapping
    - 1. [KGML pathway](http://biopython.org/DIST/docs/api/Bio.KEGG.KGML.KGML_pathway-module.html): biopython
    - 2. [Pathview](http://pathview.r-forge.r-project.org/): bioconductor for R
    - 3. [clueGO](http://apps.cytoscape.org/apps/cluego): cytoscape plugin

* pairwise distance calculation for clustered genes compared with randomly selected genes (same count)
