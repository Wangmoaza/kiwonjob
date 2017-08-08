library(ggplot2)
df <- read.table("rand_path_exist_ratio", sep='\t', header=T)

# cluster ratio values
c1f1 <- 0.0552475247525
c1f2 <- 0.066091954023
c1f3 <- 0.0772277227723
c1f4 <- 0.0841514726508
c1f5 <- 0.0681188118812
c2f1 <- 0.0473684210526
c2f2 <- 0.619047619048
c2f3 <- 0.0337381916329
c2f4 <- 0.363636363636
c2f5 <- 0.0774193548387
c3f1 <- 0.050099009901
c3f2 <- 0.0626959247649
c3f3 <- 0.0669306930693
c3f4 <- 0.0828185769847
c3f5 <- 0.0643564356436

# histogram
p <- qplot(df$rand_c1f5, geom="histogram")
p <- p + geom_vline(xintercept = c1f5, colour="red")
p

  
t.test(df$rand_c1f1, mu=c1f1)
t.test(df$rand_c1f2, mu=c1f2)
t.test(df$rand_c1f3, mu=c1f3)
t.test(df$rand_c1f4, mu=c1f4)
t.test(df$rand_c1f5, mu=c1f5)
t.test(df$rand_c2f1, mu=c2f1)
t.test(df$rand_c2f2, mu=c2f2)
t.test(df$rand_c2f2, mu=c2f2)
t.test(df$rand_c2f3, mu=c2f3)
t.test(df$rand_c2f4, mu=c2f4)
t.test(df$rand_c2f5, mu=c2f5)
t.test(df$rand_c3f1, mu=c3f1)
t.test(df$rand_c3f2, mu=c3f2)
t.test(df$rand_c3f3, mu=c3f3)
t.test(df$rand_c3f4, mu=c3f4)
t.test(df$rand_c3f5, mu=c3f5)


