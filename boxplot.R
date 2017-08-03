df <- read.table("breast_cancer_node_df_central_ess_marked", sep='\t', header=T)
#ess_genes <-read.table("gene_essential_cnt", sep='\t', header=T, row.names='gene')

library(ggplot2)

df$c1_ess <- factor(df$c1_ess, levels=c(0, 1))


# box plot
# x axis: group divided by essentiality
# y axis: centrality
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$close))
p <- p + geom_boxplot() + coord_cartesian(ylim=c(0,0.005))
p <- p + scale_x_discrete(name= "Group")
p

# violin plot
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$close))
p <- p + geom_violin() + coord_cartesian(ylim=c(0,0.01))
p <- p + scale_x_discrete(name= "Group")
p

# log scale box plot
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$close), log="y")
p <- p + geom_boxplot() + coord_cartesian(ylim=c(0,0.001))
p <- p + scale_x_discrete(name= "Group")
p

# two sample t-test
t.test(df[df$c1_cs_ess == 1, ]$close, df[df$c3_cs_ess == 0, ]$close)


df[df$c1_ess == 1,]
