#' @title Retrieve.hidden.edges
#'
#' @param dataframe containing epigenomic and provirus transcription-related attributes.
#' select columns including barcode. H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, PHA.sum, VOR.sum, DMSO.sum, expr, tpm (RNA-seq)
#' 
#' aggregation of dataframe based on barcode sequences. Compute mean values associated with identical barcode sequences.
#' Assign series ID numbers across rows.
#' 
#' epigenomic attributes  
#' select attributes related to the epigenomic network, including id, H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3.
#' Run the function "correlation.matrix" to create the correlation matix (vertex = individual barcode (single virus)).
#' Run the function "count.edges" to compute a total number of edges per vertex, forming a dataframe. Change column names.
#' 
#' provirus transcriptional attributes measured at the protein level
#' select attributes related to the transcriptional network, including id, PHA.sum, VOR.sum, DMSO.sum, expr.
#' Run the function "correlation.matrix" to create the correlation matix (vertex = individual barcode (single virus)).
#' Run the function "count.edges" to compute a total number of edges per vertex, forming a dataframe. Change column names.
#' 
#' Merge two dataframes possessing edge counts by the column "var1".
#' Convert NA to 0.
#' Compute the mean of correlation and the difference in edge counts from idnetical vertices.
#' 
#' @return A dataframe.
Retrieve.hidden.edges <- function(df) {
  df.tmp <- df %>% dplyr::select(brcd, H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, PHA.sum, VOR.sum, DMSO.sum, expr, tpm)
  
  #' keep all barcodes unique
  df.tmp.agg <- aggregate(. ~ brcd, df.tmp, mean) 
  df.tmp.agg$id <- c(1:nrow(df.tmp.agg))
  
  #' epigenomic attributes
  df.epi <- df.tmp.agg %>% dplyr::select(id, H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  mx.epi <- correlation.matrix(df.epi)
  df.edge.epi <- count.edges(mx.epi) %>% dplyr::rename(cor.epi = cor, count.epi = count)
  
  #' provirus transcriptional attributes measured at the protein level
  df.tx <- df.tmp.agg %>% dplyr::select(id, PHA.sum, VOR.sum, DMSO.sum, expr)
  mx.tx <- correlation.matrix(df.tx)
  df.edge.tx <- count.edges(mx.tx) %>% dplyr::rename(cor.tx = cor, count.tx = count)
  
  df.merge <- dplyr::full_join(df.edge.epi, df.edge.tx, by = "var1")
  df.merge[is.na(df.merge)] <- 0
  df.merge <- df.merge %>% dplyr::mutate(cor.mean = (cor.epi + cor.tx)/2, delta.edge = (count.epi - count.tx))
  
  return(df.merge)
}

#' @title correlation.matrix
#'
#' @param dataframe containing epigenomic and provirus transcription-related attributes.
#' #assign series ID numbers across each row.
#' #assign the edge list (node).
#' 
#' Convert numeric ID as character.
#' Swap columns and rows of a dataframe.
#' Replace column names by series ID.
#' Convert the matrix to a dataframe and remove the first row.
#' Convert all values as numeric.
#' Compute correlation matrix (Pearson).
#' 
#' Select the correlation coefficients with the p-value equal or smaller than 0.05 and the correlation coefficients greater than 0.
#' 
#' Remove self-connection of edges.
#' 
#' Convert var1 as integer.
#' Convert var2 as integer.
#' Convert correlation coefficients as integer.
#' 
#' @return A dataframe.
correlation.matrix <- function(df) {
  #df$id <- c(1:nrow(df))
  #node <- df %>% dplyr::select(id, dominant)
  
  df$id <- as.character(df$id)
  df.t <- t(df)
  colnames(df.t) <- df.t[1,]
  df.t <- as.data.frame(df.t[-1,])
  df.t.nu <- as.data.frame(sapply(df.t, as.numeric))
  df.mx <- cor_test(df.t.nu)
  
  df.mx.out.p_sel <- df.mx %>% dplyr::filter(p < 0.05 | p == 0.05) %>% dplyr::select(var1, var2, cor) %>% dplyr::filter(cor > 0)
  
  df.mx.out.p_sel <- df.mx.out.p_sel %>% dplyr::filter(var1 != var2)
  
  df.mx.out.p_sel$var1 <- as.integer(df.mx.out.p_sel$var1)
  df.mx.out.p_sel$var2 <- as.integer(df.mx.out.p_sel$var2)
  df.mx.out.p_sel$cor <- as.integer(df.mx.out.p_sel$cor)
  
  return(df.mx.out.p_sel)
}

#' @title count.edges
#' 
#' @param dataframe containing correlation coefficients between 2 adjacent vertices.
#' Retrieve the first column (var1) and corresponding correlation coefficients; assign the count number 1.
#' Compute the mean of correlation coefficients by grouping var1; sum the count number by grouping var1.
#' 
#' @return A dataframe.
count.edges <- function(cor.mx) {
  cor.mx.var1 <- cor.mx %>% dplyr::select(var1, cor) %>% dplyr::mutate(count = 1)
  cor.mx.var1.agg <- merge(aggregate(cor ~ var1, cor.mx.var1, mean), aggregate(count ~ var1, cor.mx.var1, sum), by = c("var1"))
  
  return(cor.mx.var1.agg)
}
