require(neuralnet)

NN.model.accuracy <- function(df) {
  # 2. Create training and test samples
  Sample <-sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.8,0.2))
  train <- df[Sample, ]
  test <- df[!Sample, ] 
 
  # 3. Fit the Multinominal Logistic Regression Model
   model <- neuralnet::neuralnet(dominant ~ H3K27ac + H3K36me3 + H3K4me1 + H3K4me3 + H3K79me3 + H3K9me3, data = train)

  # 4 Prediction
   prediction <- predict(model, test)
   
# 5. calculate accuracy
## To check the accuracy, we have to first convert actual categorical values  into numerical ones and compare them with predicted values. As a result, we will receive a list of boolean values. 

##We can use the `sum` function to find the number of `TRUE` values and divide it  by the total number of samples to get the accuracy. 
  
  check = as.numeric(test$dominant) == max.col(prediction)
  accuracy = (sum(check)/nrow(test))*100
  
  return(accuracy)
}

model.a.NN <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))

  df.model.a.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3)
  df.model.a.variable <- as.data.frame(scale(df.model.a.variable)) # Scale df
  df.model.a.variable <- df.model.a.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.a <- cbind(df.label, df.model.a.variable, by = "id")
  df.model.a$id <- NULL
  df.model.a.output <- logistic.model.confusion(df.model.a) %>% dplyr::mutate(model = "a")
  
  return(df.model.a.output)
}

model.b.NN <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))

  # model b
  df.model.b.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, expr)
  df.model.b.variable <- as.data.frame(scale(df.model.b.variable)) # Scale df
  df.model.b.variable <- df.model.b.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.b <- cbind(df.label, df.model.b.variable, by = "id")
  df.model.b$id <- NULL
  df.model.b.output <- logistic.model.confusion(df.model.b) %>% dplyr::mutate(model = "b")
}

model.c.NN <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))
 
  # model c
  df.model.c.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, tpm)
  df.model.c.variable <- as.data.frame(scale(df.model.c.variable)) # Scale df
  df.model.c.variable <- df.model.c.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.c <- cbind(df.label, df.model.c.variable, by = "id")
  df.model.c$id <- NULL
  df.model.c.output <- logistic.model.confusion(df.model.c) %>% dplyr::mutate(model = "c")
}

model.d.NN <- function(df) {
  
  df.label <- df %>% dplyr::select(dominant) %>% dplyr::mutate(id = 1:nrow(df))
  
  # model d
  df.model.d.variable <- df %>% dplyr::select(H3K27ac, H3K36me3, H3K4me1, H3K4me3, H3K79me3, H3K9me3, methyl.percent)
  df.model.d.variable <- as.data.frame(scale(df.model.d.variable)) # Scale df
  df.model.d.variable <- df.model.d.variable %>% dplyr::mutate(id = 1:nrow(df))
  df.model.d <- cbind(df.label, df.model.d.variable, by = "id")
  df.model.d$id <- NULL
  df.model.d.output <- logistic.model.confusion(df.model.d) %>% dplyr::mutate(model = "d")
}
