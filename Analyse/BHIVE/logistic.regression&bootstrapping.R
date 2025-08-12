#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2025/08/08
##############################################################
# Input: df containing B-HIVE HIV integration sites at the coordinates of all features.
# Object: The area under the curve (AUC) of logistic regression models demonstrates the classification power on the prediction of latent PHA-specific and VOR-specific MMEs.
##############################################################

require(pROC)

#Functions
1. FUN::logistic regression::dominant versus histone mark
Random.selection.bhive.dominant.histone.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = nrow(df))
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
  
   df_boostrap_label$dominant.level <- relevel(df_boostrap_label$dominant, ref = "DMSO")
   
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Multinominal Logistic Regression Model
   model <- glm(dominant.level ~ x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$dominant.level, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
  }


2. FUN::logistic regression::dominant versus histone mark + rnaseq
Random.selection.bhive.dominant.histone.rnaseq.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(dominant, tpm) %>% dplyr::mutate(id = nrow(df))
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
  
   df_boostrap_label$dominant.level <- relevel(df_boostrap_label$dominant, ref = "DMSO")
   
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Multinominal Logistic Regression Model
   model <- glm(dominant.level ~ tpm + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$dominant.level, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}
  
3. FUN::logistic regression::dominant versus histone mark + rnaseq
Random.selection.bhive.dominant.histone.rnaseq.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(dominant, tpm) %>% dplyr::mutate(id = nrow(df))
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
  
   df_boostrap_label$dominant.level <- relevel(df_boostrap_label$dominant, ref = "DMSO")
   
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Multinominal Logistic Regression Model
   model <- glm(dominant.level ~ tpm + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$dominant.level, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}

4. FUN::logistic regression::histone marks + methylation
Random.selection.bhive.dominant.methyl.5h <- function(df) {
  
  df_provirus <- df %>% dplyr::select(dominant, methyl.percent) %>% dplyr::mutate(id = nrow(df))
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
  
   df_boostrap_label$dominant.level <- relevel(df_boostrap_label$dominant, ref = "DMSO")
   
   train <- df_boostrap_label[Sample, ]
   test <- df_boostrap_label[!Sample, ] 

  # 3. Fit the Multinominal Logistic Regression Model
   model <- glm(dominant.level ~ methyl.percent + x1 + x2 + x3 + x4 + x5, na.action=na.exclude, family = "binomial", data = train)
  
   # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # create roc curve
   roc <- multiclass.roc(test$dominant.level, prediction)
   
   auc <- auc(roc)
   
   #print(auc)
   return(auc)
}
