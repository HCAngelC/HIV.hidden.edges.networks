#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
#Date: 2025/08/11
##############################################################
# Input: (1) data frame containing Jurkat methylation loci (top 25% and bottom 25%); (2) data frame containing all features associated with B-HIVE integration sites IS).
# Object: Pairing of top 25% and bottom 25% methylation loci to HIV IS.
##############################################################

#R functions
merge.by.position <- function(df.bhive, df.methyl) {
  df.merge <- merge(df.bhive, df.methyl, by = "gene_name")
  return(df.merge)
}

pair.to.methylation <- function(df.bhive, df.methyl) {
  
  df.bhive.wo.X <- df.bhive %>% dplyr::filter(chr.num != "X")
  df.bhive.wo.X$chr.num <- as.numeric(df.bhive.wo.X$chr.num)
  df.methyl.wo.X <- df.methyl %>% dplyr::filter(chr.num != "X")
  df.methyl.wo.X$chr.num <- as.numeric(df.methyl.wo.X$chr.num)
  
  df.bhive.X <- df.bhive %>% dplyr::filter(chr.num == "X")
  df.methyl.X <- df.methyl %>% dplyr::filter(chr.num == "X")
  
  df.pool <- data.frame()
  k <- 1
  
  for(i in 1:nrow(df.bhive.wo.X)) {
    for(j in 1:nrow(df.methyl.wo.X)) {
      if(k > 22) {
        break
      } else if (df.bhive.wo.X[i,20] == k & df.methyl.wo.X[j,9] == k) {
        df.merge <- merge.by.position(df.bhive.wo.X, df.methyl.wo.X)
      }
      df.pool <- rbind(df.pool, df.merge)
      k = k + 1
    }
    
  }
  df.merge.X <- merge.by.position(df.bhive.X, df.methyl.X)
    
  df.final <- rbind(df.pool, df.merge.X)
    
  Richtung.chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
  
  df.final$chr.x <- factor(df.final$chr.x, levels = Richtung.chr)
   
  return(df.final)
}
