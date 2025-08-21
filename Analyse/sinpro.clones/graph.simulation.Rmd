#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
#Date: 2025/08/20
##############################################################
# Input: igraph matrix.
# Object: 
##############################################################

library("igraph")
library("gplots")
library("graphsim")
library("scales")

# function
graph.simulation <- function(igraph) {

svg("/path_to_file/activating.graph.svg", height = 6, width = 6)
state <- rep(1, length(E(igraph)))
plot_directed(igraph, state=state, cex.node=2, cex.arrow=4, arrow_clip = 0)
dev.off()

#plot relationship matrix
#Here we plot the number of edges in the shortest paths between each pair of nodes in the graph (as an integrer value). Relative to the “diameter” (length of the longest shortest path between any 2 nodes), we can show which genes are more similar or different based on the graph structure.

svg("/home/labadmin/Documents/Arbeitplatz/Projekten/IS_Genome_Property/Etage_II_I/Illustracja/graph.simulation/distance.mx.svg", height = 6, width = 6)
short.path <- shortest.paths(igraph)
short.path[is.infinite(short.path)] <- 0

diam <- diameter(igraph)
relative_dist <- (1 + diam - short.path)/diam

heatmap.2(relative_dist, scale = "none", trace = "none", col = colorpanel(100, "grey75", "blue"), colsep = 1:length(V(igraph)), rowsep = 1:length(V(igraph)))
dev.off()

#plot sigma matrix
svg("/path_to_file/graph.simulation/sigma.mx.svg", height = 6, width = 6)
sigma_mat <- make_sigma_mat_graph(igraph,
               cor = 0.8, comm = TRUE)
heatmap.2(sigma_mat, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "blue"),
          colsep = 1:length(V(igraph)),
          rowsep = 1:length(V(igraph)))
dev.off()

#simulated data
expr <- generate_expression(100, igraph, cor = 0.8, mean = 0,
          comm = FALSE, dist = FALSE, absolute = FALSE, state = state)
svg("/path_to_file/expr.simulated.expression.svg", height = 6, width = 6)

#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:length(V(igraph)),
          rowsep = 1:length(V(igraph)))
dev.off()

svg("/path_to_file/expr.simulated.correlation.svg", height = 10, width = 10)

#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(IS.igraph)),
          rowsep = 1:length(V(IS.igraph)), labCol = "")
dev.off()
}
