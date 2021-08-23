rm(list = ls(all = TRUE)) 

# Load packages
library(grf)
library(tidyverse)
library(hdm)
library(glmnet)
library(nnls)
library(Matrix)
library(matrixStats)
library(gglasso)
library(genlasso)
library(caret)
library(psych)

# Change to paths you are using

# This source is just imported for the functions to simulate data
source("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit/Abgabe/master-thesis/Code/Sim_functions.R")

source("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit/Abgabe/master-thesis/Ensemble Learner/utils_ensemble_brunner.R")
source("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit/Abgabe/master-thesis/Ensemble Learner/ensemble_brunner.R")
source("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit/Abgabe/master-thesis/Ensemble Learner/ml_wrapper_brunner.R")
source("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit/Abgabe/master-thesis/Ensemble Learner/plasso.R")


# The group lasso and the fused lasso are integrated into the ensemble learner.
# Simulated as well as real data is prepared to test the ensemble learner.
# FYI: Some functions from utils_ensemble have been extended.



set.seed(0)
#####################
# SIMULATED DATA ####
#####################
# Set up the data
n = 750
cat = 4
cont = 4
rho = 0.3
SNR = 3

levels <- c(8, 6, 5, 4)

#betas <- betas_sim(cat, cont, levels, 0.5, 0.5) # poor function to generate betas
betas <- c(1.1,-0.11,0.54,0.4,-0.13,1.09,1,1.41,2.14,2.03,2.27,1.54,
           2.52,0,0.41,0.68,1.31,1.68,1.89,1.23,0.9,1.3,0.8,-0.1)
beta_var(betas, levels)

# Same simulation as from my thesis 
data <- simulation(n = n, levels = levels, cat = cat,cont = cont,
                   rho = rho,SNR = SNR, betas = betas,  test = F)
data.x <- data$train.x
data.y <- data$train.y
# Column names "C0x" indicate categorical variables and "X0x" continuous variables

# Categorical data are already factors from the simulation. For that example we temporarily transform them to numeric
data.x <- sapply(data.x, as.numeric)

# Create design matrix
factor = c("C01", "C02", "C03", "C04")  # indicate factors
int = colnames(data.x)
log = c("X02")
poly = c("X03","X04")
# This function was extended to also (one-hot)-encode categorical variables
data.x = design_matrix(data.x, factor=factor,int=int,int_o=2,poly=poly,poly_o=4,log=log)  # in utils_ensemble.R
# Since we expand the matrix in one single design matrix, we later can extract
# the group assignment of the columns of a model matrix (group <- attr(., "assign"))

dim(data.x)
# This function was slightly adapted to also adjust the group membership to the cleaned data
data.x = data_screen(data.x, bin_cut = 0.05, corr_cut = 0.95)
dim(data.x)

# group membership after eliminating redundant variables
group <- attr(data.x, "assign")

# Split the data
index <- sample(nrow(data.x), 2/3*nrow(data.x))
X <- data.x[index,]
Y <- data.y[index,]
X_NEW <- data.x[-index,]
Y_VAL <- data.y[-index,]

# Set the main effects 
main_effects <- rep(F, ncol(X))
main_effects[c(1:23, 27)] <- T
colnames(X)[1:40]
# Check the main effects again (may change because other variables are eliminated while sampling the data; 
# would end with an error in "group lasso low", because there will be gaps in the vector indicating group membership) 

# Information about the group-membership of encoded categorical variables is given to the 
# respective Algorithms through the variable "group" in the "args"-list. Groups not only 
# refer to groups of dummy variables, but also to polynomials and interactions.

# Please jump down to the implementation of the ensemble learner. The code to be skipped
# is the preparation for the pension dataset with one categorical variable and can be used 
# in a separate run

################'
# REAL DATA ####
################'

set.seed(5555)
# Get data that are stored in the hdm package
data(pension)
# We reduce the sample size, so we will have more feature levels compared to
# the sample size, to highlight the strength of the group lasso & fused lasso
# And that the function runs a little faster
pension <- pension[sample(nrow(pension), 1000),]

# For illustration, we combine the dummy columns (nohs hs smcol col), that indicate
# educational degree to a categorical variable by reverse dummy encoding.
pen_names <- colnames(pension)
# This function is nothing special, but does the job in that case
reverse_dummy<- function(data, levels, intercept){ 
  lev <- levels
  data_new <- c()
  for(i in 1:length(lev)){
    temp <- as.data.frame(data)[,1:lev[i]]
    rev <- matrix(0, nrow=nrow(as.data.frame(data)))
    for(j in 1:ncol(as.data.frame(temp))){
      rev <- rev + (j)*as.data.frame(temp)[,j]
    }
    data <- as.data.frame(data)[,-c(1:ncol(as.data.frame(temp)))]
    data_new <- cbind(data_new, rev)
  }
  colnames(data_new) <-  paste("X", formatC(1:length(levels), width=2, flag=0), sep = "")
  data_rev <- cbind(data_new, data)
  for(i in 1:length(lev)){
    data_rev <- as.data.frame(data_rev)
  }
  return(data_rev)
}
pension <- reverse_dummy(pension, c(rep(1,23), 4, rep(1,5), 7, 5))
colnames(pension) <- c(pen_names[1:23], "degree", pen_names[28:32], "i", "a")

# Treatment (not used)
W = pension$e401

# Outcome
Y = pension$net_tfa

# Create covariates
X = model.matrix(~ 0 + age + inc + fsize + educ + marr + twoearn + db + pira + hown + male + degree, data = pension)
factor = c("degree")  # Indicate the factors (categorical variables)
poly = c("age","inc","fsize","educ")
log = c("age","fsize","educ")
# This function was extended to also (one-hot)-encode categorical variables
X = design_matrix(X,factor=factor,int=colnames(X),int_o=2,poly=poly,poly_o=4,log=log)  # in utils_ensemble.R
# Since we expand the matrix in one single design matrix, we later can extract
# the group assignment of the columns of a model matrix (group <- attr(., "assign"))

dim(X)
# This function was slightly adapted to also adjust the group membership to the cleaned data
X = data_screen(X, bin_cut = 0.05, corr_cut = 0.975)
dim(X)
# group membership after eliminating redundant variables
group <- attr(X, "assign")

# Split the data
index <- sample(nrow(X), 2/3*nrow(X))
X_NEW <- X[-index,]
X <- X[index,]
Y_VAL <- Y[-index]
Y <- Y[index]

# Some methods like RF only use the main effects, while Lasso might also be supplied with interactions polynomials...
main_effects = rep(FALSE,ncol(X))
main_effects[c(1:12,16,20)] = TRUE
# Check the main effects again (may change because other variables are eliminated while sampling the data; ends with an error in "group lasso low", because there will be gaps in the vector indicating group membership)
colnames(X)[1:25]

# Information about the group-membership of encoded categorical variables is given to the respective
# Algorithms through the variable "group" in the "args"-list. Groups not only refer to groups of dummy
# variables, but also to polynomials and interactions!


#######################'
# ENSEMBLE LEARNER ####
#######################'

# We define the different methods of the ensemble for the outcome
# OLS with only main effects (low-dimensional) or also interactions (high-dimensional)
ols_low = create_method("ols",x_select = main_effects, name="OLS low")
ols_high = create_method("ols", name="OLS high")
# Ridge with only main effects or also interactions
ridge_low = create_method("ridge", x_select = main_effects,name="Ridge low")
ridge_high = create_method("ridge",name="Ridge high")
# Post-Lasso with only main effects or also interactions
plasso_low = create_method("plasso",x_select = main_effects,name="Post-Lasso low")
plasso_high = create_method("plasso",name="Post-Lasso high")
# Random Forest with and without honesty
forest =  create_method("forest_grf",x_select = main_effects,name="Forest",args=list(tune.parameters = "all",honesty=FALSE))
forest_hon =  create_method("forest_grf",x_select = main_effects,name="Forest honest",args=list(tune.parameters = "all"))
# Group Lasso with only main effects or also interactions
group_lasso_low <- create_method("group_lasso", x_select = main_effects, name = "Group Lasso low", args=list(group=group[main_effects])) 
group_lasso_high <- create_method("group_lasso", name = "Group Lasso high", args = list(group = group, eps=exp(250))) # higher eps makes the function faster, surprisingly with no prediction loss (according to tests). If the eps is under some threshold, the function almost freezes. 
# Fused Lasso with only main effects or also interactions
fused_lasso_low <- create_method("fused_lasso", x_select = main_effects, name = "Fused Lasso low")
fused_lasso_high <- create_method("fused_lasso", name = "Fused Lasso high") # May be slow when there are too many columns

# Run the ensemble
m = ensemble(list(ols_low,ols_high,ridge_low,ridge_high,plasso_low,plasso_high,forest,forest_hon,
                  group_lasso_low, group_lasso_high, fused_lasso_low, fused_lasso_high),
             X,Y,nfolds=5,quiet=F,xnew=X_NEW)


# RESULTS ####
# These are the predictions of the ensemble
hist(m$ensemble)

# Here are the predictions of all components
summary(m$fit_full$predictions)

# Here are the weights that each estimator receives
m$nnls_weights

# MSE of single Predictions
round(colMeans((m$fit_full$predictions - Y_VAL)^2),3)
# MSE of Super Learner
round(mean((m$ensemble - Y_VAL)^2),3)

