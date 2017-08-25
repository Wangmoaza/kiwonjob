setwd("/home/hwangbo/kiwonjob/network_analysis")
library(ggplot2)


# esesentiality dataframe
ess_df <- read.table("HumanNet_essentiality.tsv", sep='\t', header=T)
# centrality dataframe
cen_df <- read.table("HumanNet_centrality.tsv", sep='\t', header=T)

df <- merge(ess_df, cen_df, by='X')

df$c1_ess <- factor(df$c1_ess, levels=c(0, 1))
df$c3_ess <- factor(df$c3_ess, levels=c(0, 1))
df$c3_cs_ess <- factor(df$c3_cs_ess, levels=c(0, 1))


# box plot
# x axis: group divided by essentiality
# y axis: centrality
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$closeness))
p <- p + geom_boxplot()
p <- p + scale_x_discrete(name= "c1 essentiality")
p

# violin plot
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$close))
p <- p + geom_violin() + coord_cartesian(ylim=c(0,0.01))
p <- p + scale_x_discrete(name= "Group")
p

# log scale box plot
p <- ggplot(data=df, aes(x=df$c1_ess, y=df$between), log="y")
p <- p + geom_boxplot() + coord_cartesian(ylim=c(0,0.001))
p <- p + scale_x_discrete(name= "Group")
p

# two sample t-test
t.test(df[df$c1_ess == 1, ]$closeness, df[df$c1_ess == 0, ]$closeness)
t.test(df[df$c3_ess == 1, ]$closeness, df[df$c3_ess == 0, ]$closeness)
t.test(df[df$c3_cs_ess == 1, ]$closeness, df[df$c3_cs_ess == 0, ]$closeness)
t.test(df[df$c1_ess == 1, ]$between, df[df$c1_ess == 0, ]$between)
t.test(df[df$c3_ess == 1, ]$between, df[df$c3_ess == 0, ]$between)
t.test(df[df$c3_cs_ess == 1, ]$between, df[df$c3_cs_ess == 0, ]$between)


# Wilcoxon rank-sum test
wilcox.test(df[df$c1_ess == 1, ]$closeness, df[df$c1_ess == 0, ]$closeness)
wilcox.test(df[df$c3_ess == 1, ]$closeness, df[df$c3_ess == 0, ]$closeness)
wilcox.test(df[df$c3_cs_ess == 1, ]$closeness, df[df$c3_cs_ess == 0, ]$closeness)
wilcox.test(df[df$c1_ess == 1, ]$between, df[df$c1_ess == 0, ]$between)
wilcox.test(df[df$c3_ess == 1, ]$between, df[df$c3_ess == 0, ]$between)
wilcox.test(df[df$c3_cs_ess == 1, ]$between, df[df$c3_cs_ess == 0, ]$between)

