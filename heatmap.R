data <- read.table("deep_input_c3_pert_v2LogN_top_cs_3x", header=T, sep='\t', row.names = "Gene.Cell" )

# remove T.F feature
wolabel <- data
wolabel$T.F <- NULL
df <- scale(wolabel)

# uncolored branches
Heatmap(wolabel, name = "uncaled expression",
        column_title = "Features",
        column_title_gp = gpar(fontsize=10),
        row_title = "Gene:Cell",
        row_title_gp = gpar(fontsize=10),
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 2),
        cluster_columns = FALSE,
        clustering_distance_rows = "pearson",
        row_dend_width = unit("20", "mm")
        )

# colored branches

# T/F annotation
annot_df <- data.frame(T.F = data$T.F, cell = )