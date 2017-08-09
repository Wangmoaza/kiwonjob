library(ggplot2)
setwd('/home/hwangbo/kiwonjob/src')

df <- read.table("../undirected_path_length_essentiality_for_r.tsv", sep='\t', header=T)
df_c <- read.table("../undirected_path_length_essentiality.tsv", sep='\t', header=T)

df$X <- NULL
df_c$X <- NULL

attach(df)

df_c1 <- subset(df, group == 'c1_ess' | group == 'c1_non')
df_c3 <- subset(df, group == 'c3_ess' | group == 'c3_non')
df_c3_cs <- subset(df, group == 'c3_cs_ess' | group == 'c3_cs_non')
df_ess_non <- subset(df, group == 'c3_cs_ess' | group == 'c3_cs_ess_non' |  group == 'c3_cs_ess_tumor_pred')

# draw frequency plot
p <- ggplot(df_ess_non, aes(distance, colour=group)) + 
  geom_freqpoly(binwidth=1) +
  labs(x = '(undirected) pairwise distance', title='Pairwise Distances within Group')
p

# boxplot
q <- ggplot(df, aes(group, distance)) +
  geom_boxplot()
q

# t.test(df_c$c3_cs_ess_non, df_c$c3_cs_ess_tumor_pred, alternative='greater')
# t.test(df_c$c1_ess, df_c$c1_non, alternative='less')

