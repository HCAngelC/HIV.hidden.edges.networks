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
