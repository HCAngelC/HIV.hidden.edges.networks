#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
#Date: 2025/08/18
##############################################################
# Input: dataframe containing clones and distances to histone marks.
# Object: The R squared of models fitting different degrees (max. 5).
##############################################################

# functions

curve.fitting.polynomial.regression <- function(df) {
  
  df.bs <- sample_n(df, 12, replace = T)
  
  #fit polynomial regression models up to degree 5
  x = df.bs[,3]
  y = df.bs[,2]

  fit1 <- lm(y~x, data = df.bs)
  fit2 <- lm(y~poly(x, 2, raw = T), data = df.bs)
  fit3 <- lm(y~poly(x, 3, raw = T), data = df.bs)
  fit4 <- lm(y~poly(x, 4, raw = T), data = df.bs)
  fit5 <- lm(y~poly(x, 5, raw = T), data = df.bs)
  
  Rsq1 <- summary(fit1)$adj.r.squared
  Rsq2 <- summary(fit2)$adj.r.squared
  Rsq3 <- summary(fit3)$adj.r.squared
  Rsq4 <- summary(fit4)$adj.r.squared
  Rsq5 <- summary(fit5)$adj.r.squared
  
  df.out <- data.frame(Rsq = c(Rsq1, Rsq2, Rsq3, Rsq4, Rsq5), degree = c("1", "2", "3", "4", "5"))
  
  return(df.out)
}

decompose.df.list <- function(df.list) {
  df.merge <- data.frame()
  
  for (df in df.list) {
    df.merge <- rbind(df.merge, df)
  }
  
  return(df.merge)
}

plot.polynomial.regression <- function(df) {

   x = df[,3]
  y = df[,2]
  
  fit1 <- lm(y~x, data = df)
  fit2 <- lm(y~poly(x, 2, raw = T), data = df)
  fit3 <- lm(y~poly(x, 3, raw = T), data = df)
  fit4 <- lm(y~poly(x, 4, raw = T), data = df)
  #fit5 <- lm(y~poly(x, 5, raw = T), data = df)
  
  ggplot(df, aes(x = x, y = y))+geom_point()+theme_bw()+geom_line(aes(x = x, y = predict(fit1, list(x = x))), data = df, color = "red")+geom_line(aes(x = x, y = predict(fit2, list(x = x))), data = df, color = "orange")+geom_line(aes(x = x, y = predict(fit3, list(x = x))), data = df, color = "blue")+geom_line(aes(x = x, y = predict(fit4, list(x = x))), data = df, color = "purple")+xlab("Histone")+ylab("Number of hidden edges")+ylim(0, 21)+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
}

# Execution

df.sinpro.network.IS.protein.edge.count.histone <- merge(df.sinpro.network.IS.protein.edge.count, IS.komp, by = "sample") %>% dplyr::rename(group = group.x)
df.sinpro.network.IS.protein.edge.count.histone$group.y <- NULL

df.edge.sub.H3K27ac <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K27ac, group)
df.Rsq.H3K27ac <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K27ac)}))) %>% dplyr::mutate(histone = "H3K27ac")

df.edge.sub.H3K4me3 <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K4me3, group)
df.Rsq.H3K4me3 <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K4me3)}))) %>% dplyr::mutate(histone = "H3K4me3")

df.edge.sub.H3K4me1 <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K4me1, group)
df.Rsq.H3K4me1 <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K4me1)}))) %>% dplyr::mutate(histone = "H3K4me1")

df.edge.sub.H3K36me3 <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K36me3, group)
df.Rsq.H3K36me3 <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K36me3)}))) %>% dplyr::mutate(histone = "H3K36me3")

df.edge.sub.H3K9me3 <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K9me3, group)
df.Rsq.H3K9me3 <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K9me3)}))) %>% dplyr::mutate(histone = "H3K9me3")

df.edge.sub.H3K27me3 <- df.sinpro.network.IS.protein.edge.count.histone %>% dplyr::select(sample, edge.sub, d.H3K27me3, group)
df.Rsq.H3K27me3 <- decompose.df.list(as.data.frame(replicate(500, {curve.fitting.polynomial.regression(df.edge.sub.H3K27me3)}))) %>% dplyr::mutate(histone = "H3K27me3")

df.edge.sub.pool <- bind_rows(df.Rsq.H3K27ac, df.Rsq.H3K4me3, df.Rsq.H3K4me1, df.Rsq.H3K36me3, df.Rsq.H3K9me3, df.Rsq.H3K27me3)

Richtung.histone <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K9me3", "H3K27me3")
df.edge.sub.pool$histone <- factor(df.edge.sub.pool$histone, levels = Richtung.histone)

Verglichung.degree <- list(c("1", "2"), c("2", "3"), c("3", "4"), c("4", "5"), c("1", "3"), c("2", "4"), c("3", "5"), c("1", "4"), c("2", "5"), c("1", "5"))

ggplot(df.edge.sub.pool, aes(x = degree, y = Rsq, fill = histone))+geom_boxplot(color = "black")+facet_wrap(.~histone, nrow = 2)+theme_bw()+xlab("polynomial regression degree")+ylab("R squared")+scale_fill_manual(values = c(rep("orangered", 4), rep("blue4", 2)))+stat_compare_means(comparisons = Verglichung.degree, label = "p.signif")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 0, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
