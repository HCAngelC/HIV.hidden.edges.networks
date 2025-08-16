#!/usr/bin/env Rscript
rm(list=ls())
##############################################################
# Author: Heng-Chang Chen
# Date: 2025/08/18
##############################################################
# Input: df containing B-HIVE HIV integration sites at the coordinates of all features.
# Object: Calculation of sensitivity, specificity, F1 score and AUC for each model.
##############################################################

require(pROC)
require(caret)

logistic.model.confusion <- function(df) {
  # 2. Create training and test samples
  Sample <-sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.8,0.2))
  train <- df[Sample, ]
  test <- df[!Sample, ] 
 
  # 3. Fit the Multinominal Logistic Regression Model
   model <- glm(dominant ~ H3K27ac + H3K36me3 + H3K4me1 + H3K4me3 + H3K79me3 + H3K9me3 + tpm, na.action=na.exclude, family = "binomial", data = train)

  # 4 Prediction
   prediction <- predict(model, test, type = "response")
   
   # 5 convert dominant from "PHA" and "VOR" to 1's and 0's
    test$dominant<- ifelse(test$dominant =="PHA", 1, 0)
    
    test$predict.glm <- ifelse(prediction > 0.5, "1", "0")
    
   # 6 create confusion matrix
    confusion.mx <- confusionMatrix(as.factor(test$dominant), as.factor(test$predict.glm))
    
    sensitivity <- sensitivity(as.factor(test$dominant), as.factor(test$predict.glm))
    #The “true positive rate” – the percentage of individuals the model correctly predicted would default.
    
    specificity <- specificity(as.factor(test$dominant), as.factor(test$predict.glm))
    #The “true negative rate” – the percentage of individuals the model correctly predicted would not default.
    
    df.confusion.mx <- confusion.mx$table
    
    True.Pos <- df.confusion.mx[1,1]
    Fal.Pos <- df.confusion.mx[1,2]
    Fal.Neg <- df.confusion.mx[2,1]
    
    Precision <- True.Pos/(True.Pos + Fal.Pos)
    Recall <- True.Pos/(True.Pos + Fal.Neg)
    
    F1 <- 2*((Precision*Recall)/(Precision + Recall))
    
    # create roc curve
    roc <- multiclass.roc(test$dominant, prediction)
   
    auc <- auc(roc)
   
   
    df.output <- data.frame(label = c("Sensitivity", "Specificity", "F1", "AUC"), value = c(sensitivity, specificity, F1, auc))
    
    return(df.output)
}

model.a.confusion.mx <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))

  df.model.a.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  df.model.a.variable <- as.data.frame(scale(df.model.a.variable)) # Scale df
  df.model.a.variable <- df.model.a.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.a <- cbind(df.label, df.model.a.variable, by = "id")
  df.model.a$id <- NULL
  df.model.a <- df.model.a %>% dplyr::filter(dominant != "DMSO")
  df.model.a.output <- logistic.model.confusion(df.model.a) %>% dplyr::mutate(model = "a")
  
  return(df.model.a.output)
}

model.b.confusion.mx <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))

  # model b
  df.model.b.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, expr)
  df.model.b.variable <- as.data.frame(scale(df.model.b.variable)) # Scale df
  df.model.b.variable <- df.model.b.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.b <- cbind(df.label, df.model.b.variable, by = "id")
  df.model.b$id <- NULL
  df.model.b <- df.model.b %>% dplyr::filter(dominant != "DMSO")
  df.model.b.output <- logistic.model.confusion(df.model.b) %>% dplyr::mutate(model = "b")
}

model.c.confusion.mx <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))
 
  # model c
  df.model.c.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, tpm)
  df.model.c.variable <- as.data.frame(scale(df.model.c.variable)) # Scale df
  df.model.c.variable <- df.model.c.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.c <- cbind(df.label, df.model.c.variable, by = "id")
  df.model.c$id <- NULL
  df.model.c <- df.model.c %>% dplyr::filter(dominant != "DMSO")
  df.model.c.output <- logistic.model.confusion(df.model.c) %>% dplyr::mutate(model = "c")
}

model.d.confusion.mx <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))
  
  # model d
  #df.model.d.variable <- df.methyl %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, methyl.percent)
  #df.model.d.variable <- as.data.frame(scale(df.model.d.variable)) # Scale df
  #df.model.d.variable <- df.model.d.variable %>% dplyr::mutate(id = 1:nrow(df.methyl))
  #df.model.d <- cbind(df.label.methyl, df.model.d.variable, by = "id")
  #df.model.d$id <- NULL
  #df.model.d <- df.model.d %>% dplyr::filter(dominant != "DMSO")
  #df.model.d.output <- logistic.model.confusion(df.model.d) %>% dplyr::mutate(model = "d")
}
