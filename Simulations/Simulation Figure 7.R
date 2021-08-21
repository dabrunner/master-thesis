# Supplemental Material  -  Master Thesis  -  Daniel Brunner  -  Aug 23, 2021
#############################################################################
# This file contains the code for figure 7

setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions v5.R")

library(dplyr)
library(MASS)
library(caret)
library(tensorflow)
library(keras)

methods <- c("One-Hot", "Lasso",
             "GLMM", "Low Rank", "Entity Embedding", 
             "Group Lasso", "Sparse-Group Lasso", "Fused Lasso")

# Number of levels for each simulation (apply for each categorical column)
values <- c(seq(3,10,1), seq(12,50,2))

# Simulation
t <- Sys.time()
sim <- matrix(NA, ncol=length(methods)+1, nrow=length(values))
for (j in values){
  
  #### Generate Data Set ####
  n= 300   # Observations
  cat = 3  # Categorical Variables
  cont = 3 # Continuous Variables
 
  levels <- c(rep(j,cat)) # Levels determined in row 20

  rho = 0.3 # Correlation
  SNR = 3   # Signal-To-Noise Ratio
  
  betas <- betas_sim(cat, cont, levels, 0, 0)

  m = 250
  res <- matrix(NA, nrow = m, ncol = length(methods))
  colnames(res) <- methods
  
  start.time <- Sys.time()
  for (i in 1:m){
    data <- simulation3(n = n, levels = levels, cat = cat,cont = cont,
                        rho = rho,SNR = SNR, split = 2/3, betas = betas)
    train.X <- data$train.x
    train.Y <- data$train.y
    test.X <- data$test.x
    test.Y <- data$test.y

    train_data <- data.frame("Y" = train.Y, train.X)
    test_data <- data.frame("Y" = test.Y, test.X)
    
    # linear Regression
    lin_one_hot <- get_mse_linear_regression(train_data, test_data)
    lin_one_hot
    # Lasso
    mse_lasso <- get_mse_lasso(train_data, test_data)
    
    
    ##### Encoding Methods ####
    # Regularized target encoding           
    GLMM <- GLMM_encode(train_data, test_data, 5)
    GLMM_train <- GLMM$train
    GLMM_test <- GLMM$test
    colnames(GLMM_train) <- colnames(GLMM_test)
    
    lin_GLMM <- get_mse_linear_regression(GLMM_train, GLMM_test)
    
    # low rank label
    categorical <- colnames(train_data)[grepl("C", colnames(train_data))]
    most_lvls <- which.max(sapply(train_data, nlevels))
    label_train <- data.frame(sapply(train_data[,-most_lvls], as.numeric),train_data[,most_lvls])
    label_test <- data.frame(sapply(test_data[,-most_lvls], as.numeric),test_data[,most_lvls])
    most_lvls <- colnames(train_data)[most_lvls]
    colnames(label_train)[ncol(label_train)] <- most_lvls
    colnames(label_test)[ncol(label_test)] <- most_lvls
    num_components <- num_comp(label_train, most_lvls, "Y", "linear_regression", 5)
    low_rank_train <- low_rank_enc(label_train, most_lvls, num_components)
    low_rank_test <- low_rank_enc(label_test, most_lvls, num_components)
    
    lin_low_rank_label <- get_mse_linear_regression(low_rank_train, low_rank_test)
    
    
    # entity embedding
    #reduction <- Embed_Size(train_data, folds = 2)
    embed_data <- Embedding_encode(train_data, test_data, dimensionality_reduction = "very strong", batchsize = 1, epochs = 5)
    embed_train <- embed_data$train
    embed_test <- embed_data$test
    
    lin_embed <- get_mse_linear_regression(embed_train, embed_test)
    
    ####### Regression Methods ######
    
    #### Group Lasso
    mse_grplasso <- get_mse_group_lasso(train_data = train_data, test_data = test_data)
    
    #### Sparse-Group Lasso
    mse_sgl <- get_mse_sparse_group_lasso_seagull(train_data = train_data, test_data = test_data)
    
    # Fused Lasso
    mse_fused <- get_mse_fused_lasso(train_data = train_data, test_data = test_data, ordering=F)
    
    
    #####################
    # Results
    results <- c(lin_one_hot, mse_lasso, 
                 lin_GLMM, lin_low_rank_label, lin_embed,
                 mse_grplasso, mse_sgl, mse_fused)
    
    
    res[i,] <- results
  }
  time.diff  <- Sys.time() - start.time
  time.diff
  
  
  sim[which(values==j),] <- c(j, colMeans(res))
}
Sys.time()-t


sim <- as.data.frame(sim)
colnames(sim) = c("testvariable", methods)


#saveRDS(sim, "./simulation figure 7.rds")
#sim <- readRDS("./simulations/simulation figure 7.rds")

# Plot the MSE for different setups
library(ggplot2)
library(tidyr)

# Create data frame for plot
df <- cbind(sim[,1], sim[,2:9]/sim[,2]) # Relative MSE
colnames(df)[1] <- "testvariable"
df <- df %>%
  gather(key = variable, value="value", -testvariable)
df <- df[-which(df[,3]>1.5),]

# Order legend
df$variable <- factor(df$variable, levels = methods)

# Create Plot
ggplot(df, aes(x=testvariable, y=value))+
  geom_line(aes(color=variable), size=0.4)+
  labs(x="Number of feature levels", y="Relative MSE", colour = "Method") +
  theme_minimal(base_family = "Palatino")

ggsave("plot_MA.jpg", units = c("cm"), height = 1.1*7.44, width = 1.1*17.8, dpi = 150)





