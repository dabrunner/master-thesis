# Supplemental Material  -  Master Thesis  -  Daniel Brunner  -  Aug 23, 2021
#############################################################################
# This file contains the code for the simulations

# Readme: the simulation has to be started two times, because keras will crash the first time 

rm(list = ls(all = TRUE)) 
# Load the packages
library(tidyverse)
library(nnls)
library(Matrix)
library(matrixStats)
library(keras)
library(tensorflow)
library(MASS)
library(SimDesign)
library(caret)
library(gglasso)
library(genlasso)
library(glmnet)
library(seagull)
library(mlrCPO)


# Working directory
setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions.R")


methods <- c("One-Hot", "Lasso",
             "GLMM", "Low Rank", "Entity Embedding", 
             "Group Lasso", "Sparse-Group Lasso", "Fused Lasso",
             "Super Learner")

#### Generate Data Set ####
n= 300      # Observations
cat = 4     # Number of categorical variables
cont = 4    # Number of continuous variables
rho = 0.3   # correlation
SNR = 3     # signal-to-noise ratio

# Levels of each categorical variable
# levels <- c(80)
# levels <- c(40,15,10,5)
# levels <- c(30, 15, 10, 5)
levels <- c(8, 5, 4, 3)
# levels <- rep(3,cat)

# (generate betas with desired variance across or within groups)
betas <- beta_sel(cat, cont, levels,0,0) #(very inefficient and slow function)
beta_var(betas, levels)  # description of betas


runs = 5  # Number of Runs of the Montecarlo-Simulation
res <- matrix(NA, nrow = runs, ncol = length(methods))
colnames(res) <- methods
weights <- matrix(NA, nrow = runs, ncol = 6)

# Simulation
set.seed(0)
start.time <- Sys.time()
for (i in 1:runs){
print(paste("##########  Run", i, "of", runs, " ##########"))
  
# Simulated Data
  data <- simulation(n = n, levels = levels, cat = cat,cont = cont,
                      rho = rho,SNR = SNR, split = 2/3, betas = betas)
  train.X <- data$train.x
  train.Y <- data$train.y
  test.X <- data$test.x
  test.Y <- data$test.y

  train_data <- data.frame("Y" = train.Y, train.X)
  test_data <- data.frame("Y" = test.Y, test.X)

  
# REAL DATA
  # data_split <- split_data(data, 31, 2/3)
  # 
  # train.X <- data_split$train.x
  # train.Y <- data_split$train.y
  # test.X <- data_split$test.x
  # test.Y <- data_split$test.y
  # 
  # train_data <- data.frame("Y" = train.Y, train.X)
  # test_data <- data.frame("Y" = test.Y, test.X)

  
  
main_effects <- rep(TRUE, ncol(train_data))

# GLMM
glmm_ols = create_method("glmm_ols", x_select = main_effects, name = "GLMM OLS")
# low_rank_encoding
low_rank_ols = create_method("low_rank_ols", x_select = main_effects, name = "Low Rank OLS")
# entity embedding
embedding_ols = create_method("embedding_ols", x_select = main_effects,  name = "Embedding OLS")
# Group Lasso
group_lasso = create_method("group_lasso", x_select = main_effects,  name = "Group Lasso")
# Sparse-Group Lasso
sparse_group_lasso = create_method("sparse_group_lasso", x_select = main_effects,  name = "Sparse-Group Lasso")
# Fused Lasso
fused_lasso = create_method("fused_lasso", x_select = main_effects, name = "Fused Lasso")


# Run the ensemble
m = ensemble(list(glmm_ols,low_rank_ols,embedding_ols,group_lasso,sparse_group_lasso,fused_lasso),
             train_data, test_data,nfolds=5,quiet=F)


# linear Regression
mse_one_hot <- get_mse_linear_regression(train_data, test_data)

# Lasso
mse_lasso <- get_mse_lasso(train_data, test_data)


# Here are the weights that each estimator receives
weights[i, ] <- m$nnls_weights
colnames(weights) <- c("glmm", "low rank", "embedding", "group lasso", "sparse group lasso", "fused lasso")

mse_methods <- colMeans((m$fit_full$predictions - test_data$Y)^2)
mse_superlearner <- get_mse(m$ensemble, test_data$Y)
mse_full <- c(mse_one_hot, mse_lasso, mse_methods, mse_superlearner)

res[i,] <- mse_full
}
time <- Sys.time() - start.time
time

colnames(res) <- methods


as.data.frame(colMeans(weights))
as.data.frame(colMeans(res))

# saveRDS(res, "./simulations/xxx.rds")

