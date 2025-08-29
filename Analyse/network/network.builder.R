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

#' @title network.builder
#'
#' @param df A dataframe containing B-HIVE HIV integration sites at the coordinates of all features.
#' Assign series ID number to the df.
#' Generate the edge list containing two columns (ID and dominant).
#' 
#' Select required features used to construct the network.
#' Convert series ID number as character.
#' Swap the columns and the rows in the dataframe.
#' Use the firs first row (series ID number) as the column name in the matrix.
#' Convert the matrix to the datafame, excluding the first row (series ID number).
#' Convert the entire dataframe to be numeric.
#' Compute Pearson correlation coefficient between 2 adjacent vertices.
#' 
#' Select the significant correlation coefficients (the p value equal or smaller than 0.05) with the value greater than 0.
#' 
#' Remove self-connected edges.
#' 
#' Convert var1 as integer.
#' Convert var2 as integer.
#' Convert correlation coefficients as integer.
#' 
#' #igraph
#' Create graph (d, A data frame containing a symbolic edge list in the first two columns; vertices, A data frame with vertex metadata, or NULL.)
#' 
#' #assign colors for types.
#' Assign color codes.
#' Assign color codes to the attribute "dominant".
#' 
#' Plot the graph
network.builder <- function(df) {
  df$id <- c(1:nrow(df))
  node <- df %>% dplyr::select(id, dominant)
  
  df.tpm <- df[,c(26, 22:24, 27)]
  df.tpm$id <- as.character(df.tpm$id)
  df.tpm.t <- t(df.tpm)
  colnames(df.tpm.t) <- df.tpm.t[1,]
  df.tpm.t <- as.data.frame(df.tpm.t[-1,])
  df.tpm.t.nu <- as.data.frame(sapply(df.tpm.t, as.numeric))
  df.tpm.mx <- cor_test(df.tpm.t.nu)
  
  df.tpm.mx.out.p_sel <- df.tpm.mx %>% dplyr::filter(p < 0.05 | p == 0.05) %>% dplyr::select(var1, var2, cor) %>% dplyr::filter(cor > 0)

  df.tpm.mx.out.p_sel <- df.tpm.mx.out.p_sel %>% dplyr::filter(var1 != var2) #remove edges originated from self nodes.
          
  df.tpm.mx.out.p_sel$var1 <- as.integer(df.tpm.mx.out.p_sel$var1)
  df.tpm.mx.out.p_sel$var2 <- as.integer(df.tpm.mx.out.p_sel$var2)
  df.tpm.mx.out.p_sel$cor <- as.integer(df.tpm.mx.out.p_sel$cor)

# igraph
df.igraph <- graph_from_data_frame(d = df.tpm.mx.out.p_sel, vertices = node, directed = F)

#assign colors for types
Farbe <- c("green3", "grey75")
Farbe.net <- Farbe[as.numeric(as.factor(V(df.igraph)$dominant))]

plot(df.igraph, vertex.size = 7, vertex.label = NA,  edge.arrow.size = 0, vertex.color = Farbe.net)
}
