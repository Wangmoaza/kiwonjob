library(ggplot2)
setwd('/home/hwangbo/kiwonjob/')

df <- read.table("../directed_path_length_essentiality_for_r.tsv", sep='\t', header=T)
df_c <- read.table("../directed_path_length_essentiality.tsv", sep='\t', header=T)

df$X <- NULL
df_c$X <- NULL

attach(df)

df_c1 <- subset(df, group == 'c1_ess' | group == 'c1_non')
df_c3 <- subset(df, group == 'c3_ess' | group == 'c3_non')
df_c3_cs <- subset(df, group == 'c3_cs_ess' | group == 'c3_cs_non')

# draw frequency plot
p <- ggplot(df_c3_cs, aes(distance, colour=group)) + 
  geom_freqpoly(binwidth=1) +
  labs(x = '(directed) pairwise distance', title='Pairwise Distances within Group')
p

# boxplot
q <- ggplot(df, aes(group, distance)) +
  geom_boxplot(varwidth = TRUE)
q

t.test(df_c$c1_ess, df_c$c1_non)
