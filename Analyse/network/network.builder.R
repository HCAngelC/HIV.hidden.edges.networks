#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2025/08/23
##############################################################
# Input: df containing B-HIVE HIV integration sites at the coordinates of all features.
# Object: Construct a network with less than 5 attributes.
          # - edge: correlation coefficients (Pearson test).
##############################################################

# function

require("rxtatix")
require("igraph")

network.builder <- function(df) {
  df$id <- c(1:nrow(df))
  node <- df %>% dplyr::select(id, latent)
  
  df.tpm <- df[,c(26, 22:24, 27)]
  df.tpm$id <- as.character(df.tpm$id)
  df.tpm.t <- t(df.tpm)
  colnames(df.tpm.t) <- df.tpm.t[1,]
  df.tpm.t <- as.data.frame(df.tpm.t[-1,])
  df.tpm.t.nu <- as.data.frame(sapply(df.tpm.t, as.numeric))
  df.tpm.mx <- cor_test(df.tpm.t.nu)
  
  df.tpm.mx.out.p_sel <- df.tpm.mx %>% dplyr::filter(p < 0.05 | p == 0.05) %>% dplyr::select(var1, var2, cor)

  df.tpm.mx.out.p_sel <- df.tpm.mx.out.p_sel %>% dplyr::filter(var1 != var2) #remove self contact
          
  df.tpm.mx.out.p_sel$var1 <- as.integer(df.tpm.mx.out.p_sel$var1)
  df.tpm.mx.out.p_sel$var2 <- as.integer(df.tpm.mx.out.p_sel$var2)
  df.tpm.mx.out.p_sel$cor <- as.integer(df.tpm.mx.out.p_sel$cor)

# igraph
df.igraph <- graph_from_data_frame(d = df.tpm.mx.out.p_sel, vertices = node, directed = F)

#assign colors for types
Farbe <- c("green3", "grey75")
Farbe.net <- Farbe[as.numeric(as.factor(V(df.igraph)$latent))]

plot(df.igraph, vertex.size = 7, vertex.label = NA,  edge.arrow.size = 0, vertex.color = Farbe.net)
}
