#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2025/08/08
##############################################################
# Input: df containing B-HIVE HIV integration sites at the coordinates of all features.
# Object: The area under the curve (AUC) of logistic regression models demonstrates the classification power on the prediction of latent (GFP-negative_ and non-latent (GFP-positive) IS (a.k.a. prediction of molecular microenvironment of latent reservoirs).
##############################################################

require(pROC)

#Functions
## Scenario "a", features: (1) the distance of HIV IS to histone marks.
Random.selection.bhive.latent.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 5)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id")
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[2] <- "x1"
  colnames(df_boostrap_label)[3] <- "x2"
  colnames(df_boostrap_label)[4] <- "x3"
  colnames(df_boostrap_label)[5] <- "x4"
  colnames(df_boostrap_label)[6] <- "x5"
  
  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

## Scenario "b", features: (1) the distance of HIV IS to histone marks, (2) RNA barcode counts corresponding to reactivation drugs.
Random.selection.bhive.latent.drug.1h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent, PHA.sum, VOR.sum, DMSO.sum) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 1)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id") %>% unique()
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[5] <- "x1"

  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ PHA.sum + VOR.sum + DMSO.sum + x1, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

## Scenario "c", features: (1) the distance of HIV IS to histone marks, (2) Endogenous gene expression.
Random.selection.bhive.latent.rnaseq.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent, tpm) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 5)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id")
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[3] <- "x1"
  colnames(df_boostrap_label)[4] <- "x2"
  colnames(df_boostrap_label)[5] <- "x3"
  colnames(df_boostrap_label)[6] <- "x4"
  colnames(df_boostrap_label)[7] <- "x5"
  
  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ tpm + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

## Scenario "d", features: (1) the distance of HIV IS to histone marks, (2) RNA barcode counts corresponding to reactivation drugs, (3) provirus expression (barcode counts).
Random.selection.bhive.latent.drug.expr.3h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent, PHA.sum, VOR.sum, DMSO.sum, expr) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 3)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id") %>% unique()
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[6] <- "x1"
  colnames(df_boostrap_label)[7] <- "x2"
  colnames(df_boostrap_label)[8] <- "x3"

  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ PHA.sum + VOR.sum + DMSO.sum + expr + x1 + x2 + x3, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

## Scenario "e", features: (1) the distance of HIV IS to histone marks, (2) RNA barcode counts corresponding to reactivation drugs, (3) Endogenous gene expression.
Random.selection.bhive.latent.drug.rnaseq.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent, PHA.sum, VOR.sum, DMSO.sum, tpm) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 5)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id") %>% unique()
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[6] <- "x1"
  colnames(df_boostrap_label)[7] <- "x2"
  colnames(df_boostrap_label)[8] <- "x3"
  colnames(df_boostrap_label)[9] <- "x4"
  colnames(df_boostrap_label)[10] <- "x5"
  
  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ PHA.sum + VOR.sum + DMSO.sum + tpm + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

## Scenario "f", features: (1) the distance of HIV IS to histone marks, (2) RNA barcode counts corresponding to reactivation drugs, (3) provirus expression (barcode counts), (4) Endogenous gene expression.
Random.selection.bhive.latent.drug.expr.rnaseq.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(latent, PHA.sum, VOR.sum, DMSO.sum, expr, tpm) %>% dplyr::mutate(id = nrow(df))
  df_input <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  
  df_boostrap <- dplyr::select(df_input, sample(seq_len(ncol(df_input)), size = 5)) %>% dplyr::mutate(id = nrow(df))
  df_boostrap_label <- merge(df_provirus, df_boostrap, by = "id") %>% unique()
  df_boostrap_label$id <- NULL
  
  colnames(df_boostrap_label)[7] <- "x1"
  colnames(df_boostrap_label)[8] <- "x2"
  colnames(df_boostrap_label)[9] <- "x3"
  colnames(df_boostrap_label)[10] <- "x4"
  colnames(df_boostrap_label)[11] <- "x5"
  
  # 2. Create training and test samples
   Sample <-sample(c(TRUE, FALSE), nrow(df_boostrap_label), replace=TRUE, prob=c(0.8,0.2))
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Logistic Regression Model
   model <- glm(latent ~ PHA.sum + VOR.sum + DMSO.sum + expr + tpm + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$latent, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}
