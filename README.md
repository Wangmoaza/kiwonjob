# kiwonjob
Helping Kiwon

## file description

#### files from Kiwon

* *files in clusters folder* : gene list of each clusters using hierarchical clustering of the heat map. c stands for case, f stands for features. mini_ files are artificial clusters made by Kiwon.
* *breast_correlation_v2.txt* : edge information of network format src/dest/corr/p-value
* *deep_input_c3_pert_v2LogN* : cut off 3
* *deep_input_c3_pert_v2LogN_top* : contains only features with high variance
* *deep_input_c3_pert_v2LogN_top_cs_3x* : plus case selection

#### files made by me

* *gene_essential_cnt* : gene essentiality is based on deep_input_c3_pert_v2LogN file.
*  *breast_cancer_node_df* : node data frame of network. cluster marked, centrality calculated.
*  

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


