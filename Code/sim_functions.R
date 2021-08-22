# Supplemental Material  -  Master Thesis  -  Daniel Brunner  -  Aug 23, 2021
#############################################################################
#This file contains all the functions used for the master thesis


# Performance Measure ####
# mean-squared error 
get_mse <- function (pred, obs) {
  return(round(mean((pred - obs)^2)))
}

#####
# Data Cleaning ####
# This function bins numerical values e.g. years 
bin <- function(column, start, end, bin, bin_min= FALSE, bin_max = FALSE){
  # Bin all values between start and end into bins of size "bin"
  for (i in 1:((end-start)/bin)){
    small <- seq(start,end,bin)[i]
    large <- seq(start,end,bin)[i+1]
    column[which(column>=small & column<large)]<- mean(small,large)
  }
  column[which(column==large)]<- mean(small,large)
  # Bins all values > end
  if(bin_max == TRUE){
    column[which(column>end)]<- end
  }
  # Bins all values < start
  if(bin_min == TRUE){
    column[which(column<start)]<- start
  }
  return(column)
}
# Puts factor with less categories than a threshold into a "others" category
others_factor <- function(data, categorical, threshold){
  for(j in categorical){
    lvl <- c()
    for(l in 1:nlevels(data[,j])){
      lvl <- c(lvl,sum(data[,j] == levels(data[,j])[l]))
    }
    names(lvl) <- levels(data[,j])
    if(any(lvl < threshold)){
      l <- list()
      for (i in names(which(lvl <= threshold))){
        l1 <- list(i = 'OTHERS')
        names(l1) <- i
        l <- c(l, l1)
      }
      data[,j] <- recode_factor(data[,j], !!!l)
      data[,j] <- droplevels(data[,j])
    }
  }
  return(data)
}
# Creates a list of all factors with the number of each level
summarise_cat <- function(data, categorical){
  z = 1
cat_list <- list()
for (i in categorical){
  lvl <- c()
  for(l in 1:nlevels(data[,i])){
    lvl <- c(lvl,sum(data[,i] == levels(data[,i])[l]))
  }
  names(lvl) <- levels(data[,i])
  cat_list[[z]] <- lvl
  z <- z+1
}
  return(cat_list)
}
# Computes the normalized Shannon Entropy and creates a plot
entropy_cat <- function(data, categorical){
  require(SwissR)
  count <- summarise_cat(data, categorical)
    ent <- lapply(count, function(x){normalizedShannonEntropy(rep(1:length(x),x))})
    ent <- unlist(ent)
    names(ent) <- names(data[,categorical])
    graphics::boxplot(x = ent, horizontal = T, ylim = c(0,1), xlim = c(0.75,1.25), frame=F, axes=T)
}

#####
# Simulate Betas ####
# Simulates betas with predefined variance among and within groups
betas_sim <- function(cat, cont, levels, var_among_groups = 1, var_within_groups = 0){
  p <- cat + cont
  if(var_among_groups > 0.5 & var_within_groups > 0.5){mu=0}else{mu=1}
  group_means <- round(rnorm(p+1, mu, sqrt(var_among_groups)),1)
  group_means <- sample(group_means, size =p+1, replace = FALSE)
  
  groups <- c(1,(levels-1),rep(1,length(group_means[-1])-length(levels)))
  betas <- c()
  for (i in 1:length(groups)){
    xx <- rnorm(groups[i], group_means[i], if(groups[i]==1){0}else{sqrt(var_within_groups)})
    if(groups[i]!=1){
      yy <- xx-mean(xx)
      zz <- yy + group_means[i]
    }else{zz=xx}
    betas <- c(betas, zz)
  }
  return(round(betas,2))
}
# Loops through beta_sim to find exact betas with betas_sim (very slow and poor function)
beta_sel <- function(cat, cont, levels, var_among_groups, var_within_groups){
  if(var_among_groups>2 & var_within_groups>2)stop("procedure would go too long")
  if(var_among_groups + var_within_groups > 4)stop("procedure would go too long")
  x <- 0
  y <- 0
  if(var_among_groups==0 & var_within_groups==0){
    betas <- betas_sim(cat,cont,levels,var_among_groups,var_within_groups)
  }
  if(var_among_groups!=0 & var_within_groups==0){
    while (x != var_among_groups){
      betas <- betas_sim(cat, cont, levels, var_among_groups,var_within_groups)
      x <- beta_var(betas, levels)$`variance among groups`
    }
  }
  if(var_among_groups==0 & var_within_groups!=0){
      x <- 1000
    while(min(x) <= 0.925 * var_within_groups | max(x) >= 1.075 * var_within_groups){
        betas <- betas_sim(cat, cont, levels, var_among_groups, var_within_groups)
        x <- as.numeric(na.omit(beta_var(betas, levels)$`variance within groups`))
    }
  }
  if(var_among_groups!=0 & var_within_groups!=0){
    y <- 0
    # while (x != var_within_groups | y != var_among_groups){
    #   betas <- betas_sim(cat, cont, levels, var_among_groups, var_within_groups)
    #   x <- beta_var(betas, levels)$`average variance within groups`
    #   y <- beta_var(betas, levels)$`variance among groups`
    # }
    x = 0
    while (min(y) <= 0.75 * var_among_groups | max(y) >= 1.25 * var_among_groups){
      while(min(x) <= 0.75 * var_within_groups | max(x) >= 1.25 * var_within_groups){
      betas <- betas_sim(cat, cont, levels, var_among_groups, var_within_groups)
      x <- as.numeric(na.omit(beta_var(betas, levels)$`variance within groups`))
      y <- beta_var(betas, levels)$`variance among groups`
      }
    }
    
  }
  return(betas)
}
# Some functions to create different vectors of betas (not used in the simulation)
betas_zeros <- function(cat, cont, levels, zeros, group_distance=0){
  p <- cat + cont
  group_means <- matrix(c(sample(c(-1, 1), size = ceiling((1-zeros) * p-1), replace = TRUE), rep(0, floor(zeros * p))), ncol = 1)
  group_means <- sample(group_means, size = p-1, replace = FALSE)
  group_means <- c(1,1, group_means)
  
  groups <- c(1,(levels-1),rep(1,length(group_means[-1])-length(levels)))
  betas <- c()
  for (i in 1:length(groups)){
    xx <- rnorm(groups[i], group_means[i], sqrt(group_distance))
    betas <- c(betas, xx)
  }
  return(betas)
}
betas_quantile_within <- function(cat, cont, levels){
  lev <- c(1,levels-1, rep(1,cont))
  betas <- list()
  for(l in 1:length(lev)){
    betas[[l]] <- if(lev[l]==1){0.5}else{seq(0,1,1/(lev[l]-1))}
  }
  betas <- unlist(betas)
  return(betas)
}
betas_quantile_among <- function(cat, cont, levels){
  p <- cat + cont + 1
  beta_groups <- sample(seq(0,1,1/(p-1)))
  
  lev <- c(1,levels-1, rep(1,cont))
  
  betas <- list()
  for (l in 1:p){
    betas[[l]] <- rep(beta_groups[l], lev[l])
  }
  betas <- unlist(betas)
  return(betas)
}
# Computes the variance of betas among and within groups
beta_var <- function(betas, levels){
  lev <- c(1,levels-1, rep(1,length(betas)-1-sum(levels-1)))
  beta_groups <- list()
  for(g in 1:length(lev)){
    beta_groups[[g]] <- betas[1:lev[g]]
    betas <- betas[-c(1:lev[g])]
  }
  betas <- NULL
  
  group_means <- unlist(lapply(beta_groups, mean))
  var_among_groups <- var(group_means)
  var_within_groups <- unlist(lapply(beta_groups, var))
  ave_var_within_groups <- mean(var_within_groups, na.rm=TRUE)
  near_zero_groups <- sum(abs(group_means) < 0.2 * median(abs(group_means)))
  
  return(list("group_means"=round(group_means,2), 
              "variance within groups"=round(var_within_groups,2), 
              "average variance within groups"= round(ave_var_within_groups,2), 
              "variance among groups"=round(var_among_groups,2),
              "groups near zero" = near_zero_groups))
}

# Generate Data ####
# Nice function to reverse the function "model.matrix()" (not used)
reverse_OHE <- function(data, levels, intercept){
  if(all(data[,1] == rep(1, nrow(data)))){data <- data[,-1]}
  lev <- levels - 1
  data_new <- c()
  for(i in 1:length(lev)){
    temp <- as.data.frame(data)[,1:lev[i]]
    rev <- matrix(0, nrow=nrow(as.data.frame(data)))
    for(j in 1:ncol(as.data.frame(temp))){
      rev <- rev + (j+1)*as.data.frame(temp)[,j]
    }
    rev[which(rev==0)] <- 1
    data <- as.data.frame(data)[,-c(1:ncol(as.data.frame(temp)))]
    data_new <- cbind(data_new, rev)
  }
  colnames(data_new) <-  paste("C", formatC(1:length(levels), width=2, flag=0), sep = "")
data_rev <- cbind(data_new, data)
for(i in 1:length(lev)){
  data_rev <- as.data.frame(data_rev)
  data_rev[,i] <- as.factor(data_rev[,i])
}
return(data_rev)
}


# Function for the simulation study
simulation <- function(n=500, levels=c(10, rep(3,9)), cat=10, cont=10, 
                        rho=0.3, SNR=3, 
                        test = TRUE, split = 2/3, 
                        betas = rep(1, 38)){
  data <- generate_categorical_data(n, levels, cat, cont, rho)
  Y <- make_Y(data, betas, SNR)
  if(test == TRUE){data <- split_data(data, Y, split)} else {data <- list("train.y" = Y,"train.x" = data)}
  return(data)
}
# Generating the explanatory variables
generate_categorical_data <- function(n, levels, cat, cont, rho){
  require(MASS)
  require(caret)
  require(SimDesign)
  # Check if there are levels for each categorical variable
  if (cat != length(levels)) stop("Length of Level-Vector must match the Number of Categorical Predictors")
  
  p = cat + cont
  
  sigma <- matrix(rho, ncol=p, nrow=p)
  diag(sigma) <- 1

  # categorical variables are drawn from a mvnd like continuous variables
  data <- SimDesign::rmvnorm(n = n, mean = rep(0, p), sigma = sigma)%>% 
    magrittr::set_colnames(c(if(cat!=0)(paste("C", formatC(1:cat, width=2, flag=0), sep = "")),
                             if(cont!=0)(paste("X", formatC(1:cont, width=2, flag=0), sep = ""))))
  data <- data %>%
    as.data.frame()
  
  # and then the values of the variables are uniformely binned according to the levels
  if(cat != 0){
    for(i in 1:cat){
      data[,i] <- ntile(data[,i], levels[i]) %>% 
        as.factor()
      # remove the order inside the categorical variable
      data[,i] <- factor(data[,i], levels = sample(levels(data[,i])))
    }
  }
  return(data)
}
# Generate the output 
make_Y <- function(data, betas, SNR){
  data.x <- model.matrix(~., data)
  
  # error with signal-to-noise ratio
  sigma.epsilon <- sqrt(var(data.x %*% betas)/SNR)
  err <- rnorm(nrow(data), 0, sigma.epsilon)
  
  data.y <- data.x %*% betas + err
  data.y <- data.frame("Y" = data.y)
  
  return(data.y)
}
# Split data into a training and a test set
split_data <- function(data, Y, split){
  if(is.numeric(Y)==T){
    data <- data[sample(nrow(data)), ]
    sep_Y <- Y
    Y <- as.data.frame(data[,sep_Y])
    data <- data[,-sep_Y]
  }
  strat_categ <- which.max(sapply(data, nlevels))
  
  # Stratified split over the category with the most levels
  index <- createDataPartition(data[,strat_categ], p = 2/3, list = FALSE)
  
  train.x <- data[index,]
  train.y <- Y[index]
  test.x <- data[-index,]
  test.y <- Y[-index]

  return(list("train.x"=train.x,"train.y"=train.y, "test.x"=test.x, "test.y"=test.y))
}

#####
# Regression methods ####
# The functions directly compute the MSE of the prediction
# The output column  has to be named Y

# Used Methods ####
# Linear Regression implemented with glm
get_mse_linear_regression <- function(train_data, test_data){
  options(warn = -1)
  train_data <- train_data %>% as.data.frame()
  test_data <- test_data %>% as.data.frame()
  
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) == "Y")
  
  test.x <- test_data[,-sep_Y]
  test.y <- test_data[,sep_Y]
  
  # had to add that bc of cv
  test.x <- model.matrix(~., as.data.frame(test.x))[,-1]
  train_data <- model.matrix(~.,train_data)[,-1]
  
  # fit linear model
  model <- glm(Y~., as.data.frame(train_data), family = 'gaussian')
  
  # new data
  test.yhat <- predict(model, newdata = as.data.frame(test.x))
  
  mse <- get_mse(test.y, as.numeric(test.yhat))
  options(warn=1)
  return(mse)
}
# Lasso implemented with glm
get_mse_lasso <- function(train_data, test_data){
  library(glmnet)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x<- as(model.matrix(~., train_data[,-sep_Y]), 'dgCMatrix')[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- as(model.matrix(~., test_data[,-sep_Y]), 'dgCMatrix')[,-1]
  test.y <- test_data[,sep_Y]
  
  # cross-validated model
  lasso.cv <- cv.glmnet(x = train.x, y = train.y, nfolds = 5, alpha = 1, type.measure = "mse", )
  
  # new data
  test.yhat <- predict(lasso.cv, newx = test.x, s = lasso.cv$lambda.1se)
  
  mse <- get_mse(test.y, test.yhat)
  return(mse)
}
# Group Lasso implemented with gglasso
get_mse_group_lasso <- function(train_data, test_data, already_encoded = FALSE, group = NULL){
  library(gglasso)
  if  (already_encoded == TRUE & is.null(group))(stop("If the categorical variables are already encoded, the variable GROUP must be provided."))
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  if (already_encoded == FALSE){
    # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
    train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
    train.y <- train_data[,sep_Y]
    
    # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
    test.x <- model.matrix(~., test_data[,-sep_Y], )[,-1]
    test.y <- test_data[,sep_Y]
    
    group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  } else {
    train.x <- as.matrix(train_data[, -sep_Y])
    train.y <- train_data[, sep_Y]
    test.x <- as.matrix(test_data[, -sep_Y])
    test.y <- test_data[, sep_Y]
  }
  
  # cross-validated model
  gglasso.cv <- cv.gglasso(x = train.x, y = train.y, group = group, eps = exp(-5))
  
  # new data
  test.yhat <- predict(gglasso.cv, newx = data.matrix(test.x), s = gglasso.cv$lambda.min)
  
  mse <- get_mse(test.y, test.yhat)
  return(mse)
}
# Sparse Group Lasso implemented with seagull
get_mse_sparse_group_lasso_seagull <- function(train_data, test_data, real_data=TRUE){
  library(seagull)
  if (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y])[,-1]
  test.y <- test_data[,sep_Y]
  
  #groups
  group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  
  #### Dumb workaround --> with real data, R crashes when 0<alpha<1, while finding the
  ####                     ideal lambda_max. So, we approximate it manually with alpha = 0
  ####                     and alpha = 1.
  ####                     fun fact: Works also better in the simulation
  
  if (real_data == TRUE){
    max_lambda <- max(seagull(y = train.y, Z = train.x, groups = group, alpha = 1, rel_acc = 0.1)[["lambda"]])
  }
  # cross-validated model
  require(caret)
  flds <- createFolds(train.y, k = 3, list = F, returnTrain = FALSE)
  grid <- seq(0.05,0.95,0.3)
  mse.cv <- matrix(NA, ncol = length(grid), nrow = length(unique(flds)))
  colnames(mse.cv) <- grid
  
  # Cross-Validation of alpha
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    for (i in 1:length(grid)) {
      sgl <- seagull(y = try.cv, Z = trx.cv, groups = group, alpha = grid[i], max_lambda = max_lambda, rel_acc = 0.1, xi = 0.001)
      y.cv <- as.matrix(tex.cv) %*% t(as.matrix(sgl$random_effects[,]))
      mse <- min(colMeans((tey.cv - y.cv)^2))
      mse.cv[j, i] <- mse
    }
  }
  alpha <- grid[which.min(colMeans(mse.cv))] # alpha = 0 -> Group Lasso | alpha = 1 -> Lasso
  
  
  # CV for lambda
  lambda.cv <- c()
  flds <- createFolds(train.y, k = 5, list = F, returnTrain = FALSE) 
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    sgl <- seagull(y = try.cv, Z = trx.cv, groups = group, alpha = alpha, max_lambda = max_lambda, rel_acc= 0.05, xi = 0.001) #rel_acc = 0.05
    # compute mse
    mse <- colMeans((as.matrix(tex.cv) %*% as.matrix(t(sgl$random_effects)) - tey.cv)^2)
    
    lambda.cv <-c(lambda.cv, which.min(mse))
  }
  lambda <- round(median(lambda.cv),0)
  
  sgl <- seagull(y = train.y, Z = train.x, groups = group, alpha = alpha, max_lambda = max_lambda, rel_acc = 0.001, xi = 0.001) #rel_acc = 0.001
  
  # new data
  test.yhat <- as.matrix(test.x) %*% as.matrix(sgl$random_effects[lambda,])
  mse <- get_mse(test.y, test.yhat)
  
  print(paste("Alpha of Seagull-Lasso:", alpha))
  return(mse)
} 
# Fused lasso implemented with genlasso
get_mse_fused_lasso <- function(train_data, test_data, ordering = FALSE){
  library(genlasso)
  options(warn=-1)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y])[,-1]
  test.y <- test_data[,sep_Y]
  
  if (ordering == TRUE){
    corr <- cor(train.x)
    dis <- dist(corr)
    clust <- hclust(dis, method = "ward.D")
    index <- clust$order
    
    train.x <- train.x[, index]
    test.x <- test.x[, index]
  }
  # cross-validation of lambda
  require(caret)
  flds <- createFolds(train.y, k = 5, list = F, returnTrain = FALSE)
  lambda.cv <- c()
  
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    model <- fusedlasso1d(y=try.cv, X = trx.cv, approx = T, minlam = exp(-4), maxsteps = 200) #exp(-4), maxsteps=200
    pred <- predict(model, Xnew = tex.cv)
    mse <- as.numeric(colMeans((tey.cv - pred$fit)^2))
    lambda.cv <- c(lambda.cv, model$lambda[which.min(mse)])
  } 
  lambda <- median(lambda.cv)

  # Model
  fused <- fusedlasso1d(y=train.y, X = train.x, minlam = lambda/2)
  
  # Prediction
  test.yhat <- predict(fused, Xnew = test.x, lambda = lambda)
  
  # Error routine, because sometimes, the predictions are 10^30, when lambda
  # is under a certain threshold 
  if(mean(abs(test.yhat$fit)) > 10 * mean(abs(pred$fit))){
    i <- 1.1
    while (mean(abs(test.yhat$fit)) > 10*mean(abs(pred$fit))){
      test.yhat <- predict(fused, Xnew = test.x, lambda = i*lambda)
      i <- i + 0.1
    }}
  mse <- get_mse(test.y, test.yhat$fit)
  options(warn=1)
  return(mse)
}
# Not used Methods ####
# Ridge implemented with glm
get_mse_ridge <- function(train_data, test_data){
  library(glmnet)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- as(model.matrix(~., train_data[,-sep_Y]), 'dgCMatrix')[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- as(model.matrix(~., test_data[,-sep_Y]), 'dgCMatrix')[,-1]
  test.y <- test_data[,sep_Y]
  
  # cross-validated model
  ridge.cv <- cv.glmnet(x = train.x, y = train.y, nfolds = 10, alpha = 0, type.measure = "mse")
  
  # new data
  test.yhat <- predict(ridge.cv, newx = test.x, s = ridge.cv$lambda.1se)
  
  mse <- get_mse(test.y, test.yhat)
  
  return(mse)
}
# Elastic Net implemented with glm
get_mse_elastic_net <- function(train_data, test_data){
  library(glmnet)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- as(model.matrix(~., train_data[,-sep_Y]), 'dgCMatrix')[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- as(model.matrix(~., test_data[,-sep_Y]), 'dgCMatrix')[,-1]
  test.y <- test_data[,sep_Y]
  
  # cross-validated model
  require(caret)
  flds <- createFolds(train.y, k = 4, list = F, returnTrain = FALSE)
  grid <- seq(0.1,0.9,0.1)
  mse.cv <- matrix(NA, ncol = length(grid), nrow = length(unique(flds)))
  colnames(mse.cv) <- grid
  
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    for (i in grid) {
      elnet.cv <- cv.glmnet(trx.cv, try.cv, type.measure="mse", alpha=i,family="gaussian", nfolds = 5)
      test.yhat <- predict(elnet.cv, newx = tex.cv, s = elnet.cv$lambda.1se)
      mse <- get_mse(tey.cv, test.yhat)
      
      mse.cv[j, 10*i] <- mse
    }
  }
  
  alpha <- grid[which.min(colMeans(mse.cv))] # alpha = 1 -> Lasso | alpha = 0 -> Ridge
  elnet.cv <- cv.glmnet(train.x, train.y, type.measure="mse", alpha=alpha,family="gaussian", nfolds = 10)
  
  # new data
  test.yhat <- predict(elnet.cv, newx = test.x, s = elnet.cv$lambda.1se)
  
  mse <- get_mse(test.y, test.yhat)
  
  print(paste("Alpha of EN-Lasso:", alpha))
  return(mse)
}
# Modified Group Lasso (weights are set to one)
get_mse_modified_group_lasso <- function(train_data, test_data){
  library(gglasso)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y], )[,-1]
  test.y <- test_data[,sep_Y]
  
  
  group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  
  # cross-validated model
  gglasso.cv <- cv.gglasso(x = train.x, y = train.y, group = group, pf = rep(1, length(unique(group))), eps = exp(-5))
  cbind(coef(gglasso.cv, s = c(gglasso.cv$lambda.1se, gglasso.cv$lambda.min)), betas)
  
  lambda <- gglasso.cv$lambda.1se
  
  # new data
  test.yhat <- predict(gglasso.cv, newx = data.matrix(test.x), s = lambda)
  
  mse <- get_mse(test.y, test.yhat)
  return(mse)
}
# Sparse Group Lasso implemented wit SGL (very slow!)
get_mse_sparse_group_lasso <- function(train_data, test_data, alpha = 0.95){
  # Very slow. Preferring the seagull implementation of sparse-group lasso
  library(SGL)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  # colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y])[,-1]
  test.y <- test_data[,sep_Y]
  
  data.train <- list(x = train.x, y = train.y)
  data.test <- list(x = test.x, y = test.y)
  
  group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  
  # cross-validated model
  sgl.cv <- cvSGL(data = data.train, index = group, type = 'linear', nfold = 2, alpha = alpha, verbose = F)
  lam <- which.min(sgl.cv[["lldiff"]])
  
  sgl <- SGL(data = data.train, index = group, type = 'linear', alpha = alpha)
  
  # new data
  test.yhat <- predictSGL(sgl, newX = data.matrix(test.x), lam = lam)
  
  mse <- get_mse(test.y, test.yhat)
  return(mse)
}
# LightGBM
get_mse_lightgbm <- function(train_data, test_data){
  library(lightgbm)
  if (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # Create LGB Dataset
  varnames <- 
    train_data %>% colnames()
  
  test.x <- as.matrix(sapply(test_data, as.numeric))[,-sep_Y]
  train.x <- as.matrix(sapply(train_data, as.numeric))[,-sep_Y]
  
  train.y <- train_data[,sep_Y]
  test.y <- test_data[,sep_Y]
  
  # Define Categorical Variables
  categoricals.vec <- varnames[which(sapply(train_data,class) == "factor")]
  
  
  # Build dataset for analysis
  lgb.train = lgb.Dataset(data=train.x, label=train.y, free_raw_data = F, categorical_feature = categoricals.vec)
  
  # Grid_search for hyperparameters
  grid_search <- expand.grid(objective = "regression",
                             metric = "mse",
                             Depth = c(1,2,3,4), #2
                             lambda_l1 = c(0,10), #0
                             lambda_l2 = c(0,10), #0
                             Ffraction = c(0.5), #1
                             min_gts = c(0), #0
                             min_data_in_leaf = c(0), #0,5
                             bagging_fraction = c(0.3),
                             bagging_freq = c(0,1),
                             min_data_per_group = c(100), #10,100
                             max_cat_threshold = c(32),
                             cat_smooth = c(0),
                             cat_l2 = c(10),
                             learning_rate = c(0.03),#0.03
                             verbosity = -1,
                             boosting = c("gbdt") # or dart
  )
  
  model <- list()
  perf <- numeric(nrow(grid_search))
  iter <- numeric(nrow(grid_search))
  
  start  <- Sys.time()
  for (i in 1:nrow(grid_search)) {
    model[[i]] <- lgb.cv(list(objective = "regression",
                              metric = "mse",
                              lambda_l1 = grid_search[i, "lambda_l1"],
                              lambda_l2 = grid_search[i, "lambda_l2"],
                              max_depth = grid_search[i, "Depth"],
                              feature_fraction = grid_search[i, "Ffraction"],
                              min_data_in_leaf = grid_search[i, "min_data_in_leaf"],
                              bagging_fraction = grid_search[i, "bagging_fraction"],
                              bagging_freq = if(grid_search[i, "boosting"] !="goss"){grid_search[i, "bagging_freq"]} else {0},
                              min_data_per_group = grid_search[i, "min_data_per_group"],
                              max_cat_threshold = grid_search[i, "max_cat_threshold"],
                              cat_smooth = grid_search[i, "cat_smooth"],
                              cat_l2 = grid_search[i, "cat_l2"],
                              min_gain_to_split = grid_search[i, "min_gts"]
    ),
    data = lgb.train, nfold = 5, learning_rate = grid_search[i, "learning_rate"],
    num_leaves = 2^(grid_search[i, "Depth"]), nrounds = 2500, early_stopping_rounds = 250, 
    eval = "mse", categorical_feature = categoricals.vec, boosting = grid_search[i, "boosting"], 
    verbose = -1
    )
    perf[i] <- model[[i]]$best_score
    iter[i] <- model[[i]]$best_iter
    cat(round(i/nrow(grid_search)*100,2), "%", " ", sep = "")
  }
  diff <- Sys.time()-start
  
  best_tuning <- grid_search[which.min(perf),]
  best_iter <-  iter[which.min(perf)]
  
  lgb.grid <- as.list(c(sapply(best_tuning[,1:2], as.character), best_tuning[-c(1,2,15)]))
  
  # Train final model
  lgb.model = lgb.train(params = lgb.grid, data = lgb.train, learning_rate = best_tuning[1,'learning_rate'],
                        num_leaves = 2^best_tuning[1, "Depth"], nrounds = best_iter, verbose=-1,  eval = "mse", categorical_feature = categoricals.vec)
  

  pred <- predict(lgb.model, test.x)
  
  mse <- get_mse(test.y, pred)
  
  return(mse)
}
# XGBoost
get_mse_XGBoost <- function(train_data, test_data){
  library(xgboost)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  # make Sparse Matrices
  colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  matrix_train <- as(model.matrix(~.,train_data[,-sep_Y])[,-1], "dgCMatrix")
  
  colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  matrix_test <- as(model.matrix(~.,test_data[,-sep_Y])[,-1], "dgCMatrix")
  
  train.Y <- train_data[,sep_Y]
  test.Y <- test_data[,sep_Y]
  
  # Set up Data
  trainxg <- xgboost::xgb.DMatrix(matrix_train, label = train.Y)
  
  
  grid_search <- expand.grid(booster = c("gbtree"),
                             objective = "reg:squarederror",
                             eval_metric = "rmse",
                             eta = c(0.01), #learning rate
                             gamma = c(0),
                             lambda = c(1),
                             alpha = c(0),
                             min_child_weight = 1, 
                             max_depth = c(2, 4, 6), 
                             subsample = c(0.5, 0.75), 
                             colsample_bytree = c(0.6, 0.8, 1))
  
  model <- list()
  perf <- numeric(nrow(grid_search))
  iter <- numeric(nrow(grid_search))
  
  for (i in 1:nrow(grid_search)) {
    model[[i]] <- xgb.cv(data = trainxg,
                         params = list(
                           booster = "gbtree",
                           objective = "reg:squarederror",
                           eval_metric = "rmse",
                           eta = 0.01, 
                           gamma = 0,
                           lambda = 1,
                           alpha = 0,
                           min_child_weight = 1, 
                           max_depth = grid_search[i, "max_depth"],
                           subsample = grid_search[i, "subsample"],
                           colsample_bytree = grid_search[i, "colsample_bytree"]
                         ),
                         nfold = 5, nrounds = 2000, early_stopping_rounds = 200, verbose = 0
    )
    perf[i] <- min(model[[i]]$evaluation_log$test_rmse_mean)
    iter[i] <- which(model[[i]]$evaluation_log$test_rmse_mean == min(model[[i]]$evaluation_log$test_rmse_mean))
    cat(round(i/nrow(grid_search)*100,2), "%", " ", sep = "")
  }

  best_tuning <- grid_search[which.min(perf), ]
  best_iter <- iter[which.min(perf)]
  
  params <- as.list(c(sapply(best_tuning[,1:3], as.character), best_tuning[-c(1:3)]))
  
  # Train the Model
  model <- xgboost(data = trainxg, params = params, nrounds = best_iter, verbose = 0)
  pred <- predict(model, matrix_test)
  mse <- get_mse(pred, test.Y)
  
  return(mse)
}
# Tuned Forest
get_mse_forest <- function(train_data, test_data){
  train_data <- as.data.frame(train_data)
  test_data <- as.data.frame(test_data)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  train.X <-sapply(train_data[,-sep_Y], as.numeric)
  train.Y <- train_data[,sep_Y]
  
  test.X <- sapply(test_data[,-sep_Y], as.numeric)
  test.Y <- test_data[,sep_Y]
  
  forest <-  grf::regression_forest(X = train.X,Y = train.Y, tune.parameters = c('mtry', 'sample.fraction', 'min.node.size', 'alpha', 'imbalance.penalty', 'honesty.fraction', 'honesty.prune.leaves'))
  
  grid <- c(seq(100, 1000, 100), seq(1500, 3000, 500))
  
  model <- list()
  perf <- numeric(length(grid))
  
  for (i in 1:length(grid)){
    model[[i]] <- grf::regression_forest(train.X,train.Y, 
                                         sample.fraction = forest$tuning.output$params[['sample.fraction']],
                                         mtry = forest$tuning.output$params[['mtry']],
                                         min.node.size = forest$tuning.output$params[['min.node.size']],
                                         honesty.fraction = forest$tuning.output$params[['honesty.fraction']],
                                         honesty.prune.leaves = forest$tuning.output$params[['honesty.prune.leaves']],
                                         alpha = forest$tuning.output$params[['alpha']],
                                         imbalance.penalty = forest$tuning.output$params[['imbalance.penalty']], 
                                         compute.oob.predictions = T,
                                         num.trees = grid[i])
    perf[i] <- get_mse(model[[i]]$predictions, train.Y)
  }
  
  model <- model[[which.min(perf)]]
  
  prediction <- predict(model, test.X)
  mse <- get_mse(t(prediction), test.Y)
  return(mse)
}
# Local Linear Forest (Favorite Prediction Algorithm :-)
get_mse_llforest <- function(train_data, test_data){
  library(grf)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  train.x <- sapply(train_data[,-sep_Y], as.numeric)
  train.y <- train_data[,sep_Y]
  
  test.x <- sapply(test_data[,-sep_Y], as.numeric)
  test.y <- test_data[,sep_Y]
  
  # Model
  LLforest <- ll_regression_forest(train.x, train.y)
  
  # Fined the optimal lambda for the ridge penalty in the prediction
  tuned.lambda <- tune_ll_regression_forest(LLforest)
  
  # Prediction
  test.yhat <- predict(LLforest, newdata = test.x, ll.lambda = tuned.lambda$lambda.min)$predictions
  
  mse <- get_mse(test.y, test.yhat)
  return(mse)
  
}
# Prediction algorithm from the paper "Sparse Modeling of categorical explanatory data - Gertheiss"
get_mse_gertheiss <- function(train_data, test_data){
  require(gvcm.cat)

  sep_Y <- which(grepl("Y", colnames(train_data)))
  categorical <- which(grepl("C", colnames(train_data)))
  continuous <- which(grepl("X", colnames(train_data)))
  
  
  p <- paste0(paste0(" + p(", colnames(train_data)[categorical], collapse = ", 'L1')"), ", 'L1')")
  c <- paste0(" + ", paste0( colnames(train_data)[continuous]), collapse = "")
  f <- as.formula(paste("Y ~ 0", p, c))
  
  
  m1 <- gvcm.cat(f, train_data, tuning = list(Lambda = TRUE))
  
  prediction <- predict.gvcm.cat(m1, test_data[,-sep_Y])
  
  mse <- get_mse(test_data[,1], prediction$fit)
  
  return(mse)
}

#####
# Encoding methods ####
# Low Rank Encoding
low_rank_enc <- function(data, categorical, num_comp){
  Y <- data[, which(colnames(data) == "Y")]
  train.X <- data[, -which(colnames(data) == "Y")]
  X <- train.X[,-which(colnames(train.X) %in% categorical)]
  G <- train.X[, which(colnames(train.X) %in% categorical)]
  
  if (is.null(dim(X)) == TRUE){
    X <- as.matrix(X)
    colnames(X) <- "X1"
  }
  
  if (length(categorical) != 1){
    X_aug <- X
    for (i in categorical){
      X_aug_help <- X_aug
      
      # means encoding
      p <- dim(X)[2]
      CM <- as.matrix(aggregate(X, list(G[,i]), mean)[, 2:(p + 1)])
      
      # low rank encoding
      decomp <- svd(CM)
      num_components <- min(dim(X)[2], ncol(decomp$u), num_comp)
      CM <- as.matrix(decomp$u[, 1:num_components])
      CM <- as.matrix(CM)
 
      # Augment original matrix
      X_aug <- cbind(X_aug, CM[as.integer(droplevels(G[,i])), ])
      colnames(X_aug) <- c(X, paste(paste(i, "ENC", sep=""), 1:dim(CM)[2], sep =""))
    }
    } else {
    X_aug <- X
    
   # plot(X[,1],X[,2], xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
    # means encoding
    p <- dim(X)[2]
    CM <- as.matrix(aggregate(X, list(G), mean)[, 2:(p + 1)])
 #   plot(CM[,1],CM[,2], xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
    
    # low rank encoding
    decomp <- svd(CM)
    num_components <- min(dim(X)[2], ncol(decomp$u), num_comp)
    CM <- as.matrix(decomp$u[, 1:num_components])
    CM <- as.matrix(CM)
    
    # Augment original matrix
    X_aug <- cbind(X_aug, CM[as.integer(droplevels(G)), ])
    colnames(X_aug) <- c(colnames(X), paste(paste(categorical, "ENC", sep=""), 1:dim(CM)[2], sep =""))
    }
  
  encoded_data <- cbind(X_aug, Y)
  return(encoded_data)
}
# Cross-Validation of Low Rank Encoding (Adapted from: Johannemann et al., 2019)
num_comp <- function(data, categorical, response, model, folds=5, quiet = TRUE){
  options(warn=-1)
  r <- NULL
  attempt <- 1
  while( is.null(r) && attempt <= 10 ) {
    attempt <- attempt + 1
    try(
      r <- suppressWarnings(num_comp_low_rank(data, categorical, response, model, folds, quiet))
      ) 
  } 
  if(is.null(r)){r <- 1}
  options(warn=1)
  return(r)
}
num_comp_low_rank <- function(data, categorical, response, model, folds=5, quiet = TRUE){
  cv_vals <- c(1:15)
  cv_mses <- c()
  randomized_df2 <- data[sample(nrow(data)), ]
  rownames(randomized_df2) <- NULL
  require(caret)
  fold_cat2  <- createFolds(randomized_df2[,1], k = folds, list = T, returnTrain = F)
  for (ii in 1:length(cv_vals)) {
    if(quiet==FALSE) {print(paste("Running CV...", ii, "of", length(cv_vals)))}
    cv_mse <- c()
    for (jj in 1:folds) {
      testIndexes2 <- fold_cat2[[jj]]
      testData2 <- randomized_df2[testIndexes2, ]
      trainData2 <- randomized_df2[-testIndexes2, ]
      
      trainData2enc <- low_rank_enc(data = trainData2,categorical = categorical, num_comp = cv_vals[ii])
      testData2enc <-  low_rank_enc(data = testData2,categorical = categorical, num_comp = cv_vals[ii])
      
      cv_mse <- c(cv_mse, switch(model, 
                                 'linear_regression' = get_mse_linear_regression(trainData2enc, testData2enc),
                                 'lasso' = get_mse_lasso(trainData2enc, testData2enc),
                                 'ridge' = get_mse_ridge(trainData2enc, testData2enc),
                                 'elastic_net' = get_mse_elastic_net(trainData2enc, testData2enc),
                                 'group_lasso' = get_mse_group_lasso(trainData2enc, testData2enc),
                                 'modified_group_lasso' = get_mse_modified_group_lasso(trainData2enc, testData2enc),
                                 'sparse_group_lasso' = get_mse_sparse_group_lasso_seagull(trainData2enc, testData2enc, real_data = F),
                                 'fused_lasso' = get_mse_fused_lasso(trainData2enc, testData2enc),
                                 'forest' = get_mse_forest(trainData2enc, testData2enc),
                                 'forest_untuned' = get_mse_forest_notune(trainData2enc, testData2enc),
                                 'llforest' = get_mse_llforest(trainData2enc, testData2enc), 
                                 'light_gbm' = get_mse_lightgbm(trainData2enc, testData2enc),
                                 'xgboost' = get_mse_XGBoost(trainData2enc, testData2enc)
      )
      )
      
    }
    cv_mses <- c(cv_mses, mean(cv_mse))
    
    # If there is no improvement over the last 3 iterations, break loop
    if (length(cv_mses) > 3){
      if((cv_mses[length(cv_mses)]>=cv_mses[length(cv_mses)-1]) & (cv_mses[length(cv_mses)-1]>=cv_mses[length(cv_mses)-2]) & (cv_mses[length(cv_mses)-2]>=cv_mses[length(cv_mses)-3])){
        if(quiet==FALSE){print("CV stopped. No improvement over last iterations")}
        break
      }
    }
    
  }
  mx <- which(cv_mses == min(cv_mses))[1]
  num_components <- cv_vals[mx]
  
  return(num_components)
}
# Regularized Target encoding with glmm  (Source: Pargent et al., 2021)
GLMM_encode <- function(train_data, test_data, n.folds=5){#Pargent et al. (2021)
  flagNewLvls = function(fac, lvls) {
    char = as.character(fac)
    char[!(char %in% lvls | is.na(char))] = "..new..level.."
    factor(char)
  }
  fitLmer = function(feature, target) {
    # performance tips from https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
    nlopt <- function(par, fn, lower, upper, control) {
      .nloptr <<- res <- nloptr::nloptr(par, fn, lb = lower, ub = upper, 
                                        opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                                    maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
      list(par = res$solution,
           fval = res$objective,
           conv = if (res$status > 0) 0 else res$status,
           message = res$message)
    }
    
    args = list(formula = y ~ 1 + (1 | lvl),
                data = data.frame(lvl = feature, y = target),
                na.action = na.omit,
                control = lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
    
    mod = do.call(lme4::lmer, args)
    coefs = coef(mod)$lvl
    lvls = rownames(coefs)
    coefs = coefs[,1]
    names(coefs) = lvls
    intercept = unname(lme4::fixef(mod))
    # replace missing coefs with intercept value
    coefs[is.na(coefs)] = intercept
    # save intercept value for new levels during prediction
    coefs = c(coefs, ..new..level.. = intercept)
    coefs
  }
  
  # Split data to do the encoding
  sep_Y <- which(colnames(train_data)=='Y')
  cat <- which(grepl("C", colnames(train_data)))
  cont <- which(grepl("X", colnames(train_data)))
  
  train.cont <- train_data[,cont]
  test.cont <- test_data[,cont]
  
  train.cat <- train_data[,cat]
  test.cat <- test_data[,cat]
  
  if (is.null(ncol(train.cat))){
    train.cat <- train.cat %>% as.data.frame() %>% magrittr::set_colnames("C1")
    test.cat <- test.cat %>% as.data.frame() %>% magrittr::set_colnames("C1")
  }
  
  train.Y <- train_data[,sep_Y]
  test.Y <- test_data[,sep_Y]
  
  
  # for prediction, use complete encoding model
  control = lapply(train.cat, function(col) {
    fitLmer(col, train.Y)
  })
  
  # if n.folds == 1 use only the complete model in training
  if (n.folds == 1) {
    train.cat = lapply(colnames(train.cat), function(cname) {
      as.numeric(control[[cname]][as.character(train.cat[[cname]])])
    })
  }
  
  # else use n.folds encoding models in crossvalidation scheme in training
  if(n.folds > 1){
    library(mlrCPO)
    rinst = makeResampleInstance(makeResampleDesc("CV", iters = n.folds), 
                                 size = nrow(train.cat))
    # list with n.folds models for each data column
    mod_list = lapply(colnames(train.cat), function(cname) {
      lapply(rinst$train.inds, function(inds) {
        fitLmer(train.cat[[cname]][inds], train.Y[inds])
      })
    })
    names(mod_list) = names(train.cat)
    # list with replaced values in n.folds for each data column
    num_vals_list = lapply(colnames(train.cat), function(cname) {
      lapply(seq_len(n.folds), function(fold) {
        # add new level for all values not present during training
        fac = flagNewLvls(train.cat[[cname]][rinst$test.inds[[fold]]],
                          names(mod_list[[cname]][[fold]]))
        mod = mod_list[[cname]][[fold]]
        
        as.numeric(mod[as.character(fac)])
      })
    })
    names(num_vals_list) = names(mod_list)
    
    # recombine replaced values from n.folds for each data column
    train.cat[] = lapply(seq_along(train.cat), function(cnum) {
      unlist(num_vals_list[[cnum]])[order(unlist(rinst$test.inds))]
    })
  }
  
  # For the test data, use the complete model
  test.cat[] = lapply(colnames(test.cat), function(cname) {
    as.numeric(control[[cname]][as.character(test.cat[[cname]])])
  })
  
  # merge data back together
  train.cat <- train.cat %>% as.data.frame() %>% magrittr::set_colnames(colnames(train_data)[cat])
  train.Y <- train.Y %>% as.data.frame() %>% magrittr::set_colnames("Y")
  train_enc <- data.frame(train.Y, train.cat, train.cont)
  
  
  test.cat <- test.cat %>% as.data.frame() %>% magrittr::set_colnames(colnames(test_data)[cat])
  test.Y <- test.Y %>% as.data.frame() %>% magrittr::set_colnames("Y")
  test_enc <- data.frame(test.Y, test.cat, test.cont)
  
  # with columns with a lot of levels, na's may happen
  na_col <- apply(is.na(test_enc),2, any)
  if(any(na_col)==TRUE){
    for(na in which(na_col)){
      test_enc[which(is.na(test_enc[,na])), na] <- mean(test_enc[,na],na.rm=T)
    }
  }
  
  return(list('train' = train_enc, 'test' = test_enc))
}
# Entity Embedding
Embedding_encode <- function(train_data, test_data, dimensionality_reduction = "standard", epochs=15, batchsize = 1, quiet = TRUE){
  # Bachsize 5, epochs = 40
  require(dplyr)
  require(keras)
  require(tensorflow)
  suppressPackageStartupMessages(require(purrr))
  
  if(quiet==TRUE){verbose = 0} else {verbose =1}
  
  x_train <- train_data %>% dplyr::select(-Y)
  x_test <- test_data %>% dplyr::select(-Y)
  y_train <- train_data[, "Y"]
  y_test <- test_data[, "Y"] 
  
  categorical <- colnames(train_data %>% select_if(is.factor))
  
  # We have to set the factor stings to numeric strings to use keras 
  # In order that the test data has the same new levels, although there are missing (dropped) or new levels in the test data, we create a mapping
  for(cc in categorical){
    oldLevels <- levels(x_train[,cc])
    newLevels <- 1:length(levels(x_train[,cc]))
    levels(x_train[,cc]) <- newLevels
    levels(x_test[,cc]) <- tidyr::replace_na(newLevels[match(levels(x_test[,cc]), oldLevels)],0)
  }
  
  # scaling (when scaling, you should also put scale(train.Y) in the model fitting step)
  train_means <- colMeans(x_train[sapply(x_train, is.double)]) %>% unname()
  train_sds <- apply(x_train[sapply(x_train, is.double)], 2, sd)  %>% unname()
  train_sds[train_sds == 0] <- 0.000001
  
  x_train[sapply(x_train, is.double)] <- sweep(
    x_train[sapply(x_train, is.double)], 2, train_means) %>%
    sweep(2, train_sds, "/")
  x_test[sapply(x_test, is.double)] <- sweep(
    x_test[sapply(x_test, is.double)], 2, train_means) %>%
    sweep(2, train_sds, "/")
  
  
  # number of levels per factor, required to specify input dimensionality for
  # the embedding layers
  # n_levels_in <- purrr::map(x_train %>% select_if(is.factor), compose(length, base::levels)) %>%
  #   unlist()
  n_levels_in <- sapply(x_train %>% select_if(is.factor), nlevels)
  
  # output dimensionality for the embedding layers, need +1 because Python is 0-based
  #n_levels_out <- n_levels_in %>% sqrt() %>% trunc() %>% `+`(1)
  n_levels_out <- switch(dimensionality_reduction, 
                         "weak" = n_levels_in%/%2 + 1,
                         "standard" = n_levels_in^(1/2) %>% trunc() %>% `+` (1),
                         "strong" = n_levels_in^(1/3) %>% trunc() %>% `+` (1),
                         "very strong" = n_levels_in^(1/4) %>% trunc() %>% `+` (1)
  )
  
  # each embedding layer gets its own input layer
  cat_inputs <- purrr::map(n_levels_in, function(l) layer_input(shape = 1)) %>%
    unname()
  
  # construct the embedding layers, connecting each to its input
  embedding_layers <- vector(mode = "list", length = length(cat_inputs))
  for (i in 1:length(cat_inputs)) {
    embedding_layer <-  cat_inputs[[i]] %>% 
      layer_embedding(input_dim = n_levels_in[[i]] + 1, 
                      output_dim = n_levels_out[[i]], name = paste("Embedding", categorical[i], sep = "")) %>%
      layer_flatten()
    embedding_layers[[i]] <- embedding_layer
  }
  
  # create a single input and a dense layer for the numeric data
  quant_input <- layer_input(shape = length(which(sapply(x_train, is.double))))
  
  quant_dense <- quant_input %>% layer_dense(ncol(train_data)-length(categorical)-1)
  
  intermediate_layers <- list(embedding_layers, list(quant_dense)) %>% flatten()
  inputs <- list(cat_inputs, list(quant_input)) #%>% flatten()
  
  
  output <- layer_concatenate(intermediate_layers) %>%
    layer_dense(units = 512, activation = "relu") %>% #256
    # layer_dense(units = 64, activation = "relu", kernel_initializer = "normal") %>% #256
    layer_dense(units = 1)
  
  model <- keras_model(inputs, output)
  model %>% compile(loss = "mean_squared_error", optimizer = optimizer_adam(lr=0.0002, decay = exp(-5))) #optimizer_adam(lr=0.001)
  
  
  # Each layer has to be feed into the model as list element
  xcont <- as.matrix(x_train %>% select_if(is.double))
  xcat <- apply(as.matrix(sapply(x_train[,categorical], as.integer)), 2, list)
  train_input <- list(xcat, xcont)
  
  hist <- model %>% fit(x = train_input, y = scale(y_train), epochs = epochs, 
                        batch_size = batchsize, verbose = verbose)
  
  embeddings <- list()
  for (zz in 1:length(cat_inputs)){
    embedding <- get_layer(model, paste("Embedding", categorical[zz], sep=""))$get_weights() %>%
      as.data.frame()
    colnames(embedding) <- paste(categorical[zz], "ENC", 1:ncol(embedding), sep = "")
    embedding$name <- c("0", levels(x_train[,categorical[zz]]))
    embeddings[[zz]] <- embedding
  }
  
  train_enc <- data.frame("ID" =  1:length(y_train), "Y"=y_train, x_train)
  test_enc <- data.frame("ID" = 1:length(y_test), "Y"=y_test, x_test)
  for (j in 1:length(categorical)){
    train_enc <- merge(train_enc, embeddings[[j]], by.x=categorical[j], by.y="name",sort = F)[,-1]
    test_enc <- merge(test_enc, embeddings[[j]], by.x=categorical[j], by.y="name", sort=F)[,-1]
  }
  train_enc <- train_enc[order(train_enc$ID),-1]
  test_enc <- test_enc[order(test_enc$ID),-1]
  get_mse_linear_regression(train_enc, test_enc)
  return(list("train"=train_enc, "test"=test_enc))
}
# Cross-Validation of dimensionality reduction (not used)
Embed_Dim_Reduction <- function(data, model="linear_regression", folds=5){
  cv_vals <- c("weak", "standard", "strong", "very strong")
  cv_mses <- c()
  randomized_df2 <- data[sample(nrow(data)), ]
  rownames(randomized_df2) <- NULL
  require(caret)
  fold_cat2  <- createFolds(randomized_df2[,1], k = folds, list = T, returnTrain = F)
  for (ii in 1:length(cv_vals)) {
    print("Running CV...")
    cv_mse <- c()
    for (jj in 1:folds) {
      testIndexes2 <- fold_cat2[[jj]]
      testData2 <- randomized_df2[testIndexes2, ]
      trainData2 <- randomized_df2[-testIndexes2, ]
      
      Data2enc <- Embedding_encode(trainData2, testData2, dimensionality_reduction = cv_vals[ii])
      trainData2enc <- Data2enc$train
      testData2enc <- Data2enc$test
      
      cv_mse <- c(cv_mse, switch(model, 
                                 'linear_regression' = get_mse_linear_regression(trainData2enc, testData2enc),
                                 'lasso' = get_mse_lasso(trainData2enc, testData2enc),
                                 'ridge' = get_mse_ridge(trainData2enc, testData2enc),
                                 'fused_lasso' = get_mse_fused_lasso(trainData2enc, testData2enc),
                                 'group_lasso' = get_mse_group_lasso(trainData2enc, testData2enc),
                                 'light_gbm' = get_mse_lightgbm(trainData2enc, testData2enc),
                                 'sparse_group_lasso' = get_mse_sparse_group_lasso_seagull(trainData2enc, testData2enc, real_data = F),
                                 'elastic_net' = get_mse_elastic_net(trainData2enc, testData2enc),
                                 'forest' = get_mse_forest(trainData2enc, testData2enc),
                                 'forest_untuned' = get_mse_forest_notune(trainData2enc, testData2enc),
                                 'llforest' = get_mse_llforest(trainData2enc, testData2enc), 
                                 'stand_group_lasso' = get_mse_stand_group_lasso(trainData2enc, testData2enc),
                                 'modified_group_lasso' = get_mse_modified_group_lasso(trainData2enc, testData2enc),
                                 'xgboost' = get_mse_XGBoost(trainData2enc, testData2enc)
      )
      )
      
    }
    cv_mses <- c(cv_mses, mean(cv_mse))
    
  }
  mx <- which(cv_mses == min(cv_mses))[1]
  
  
  
  reduction <- cv_vals[mx]
  
  return(reduction)
}

#####
# Ensemble by M.Knaus with my methods ####
# UTILS_ENSEMBLE ####
prep_cf_mat = function(n,cf,cl=NULL) {
  
  if (cf == 1) cf_mat = matrix(rep(1,n),ncol=1)
  
  if (!is.null(cl)) {
    fold = ntile(runif(length(unique(cl))),cf)
    fold = factor(fold[match(cl,unique(cl))])
    cf_mat = model.matrix(~0+fold)
  }
  else {
    fold = factor(ntile(runif(n),cf))
    cf_mat = model.matrix(~0+fold)
  }
  colnames(cf_mat) = sprintf("CF %d",1:cf)
  
  return(cf_mat==1)
}

ensemble = function(ml,
                    train_data, test_data, 
                    nfolds=5, 
                    weights=FALSE, quiet=TRUE) {
  # Matrix to store the cross-validated predictions
  fit_cv = matrix(NA,nrow(train_data),length(ml))
  colnames(fit_cv) = sprintf("Method%s",seq(1:length(ml)))
  for (i in 1:length(ml)) {
    if (!is.null(ml[[i]]$name)) colnames(fit_cv)[i] = ml[[i]]$name
  }
  
  # Get CV folds
  cvf = prep_cf_mat(nrow(train_data),nfolds) # see utils_ensemble.R
  
  # Loop over different folds
  if (length(ml) > 1) {
    for (i in 1:nfolds) {
      # Define training and test sample for this fold
      fold = cvf[,i]
      tr = train_data[!fold,]
      te = train_data[fold,]
      # Get predictions of all methods for test sample
      fit_cv[fold,] = ensemble_core(ml,tr, te,quiet=quiet)$predictions
    }
    # Cross-validated MSE
    fit_cv[is.na(fit_cv)] = mean(train_data$Y) # e.g. glmnet produced sometimes NaN for logistic Ridge
    mse_cv = colMeans((c(train_data$Y) - fit_cv)^2)
    
    # Calculate weights each method receives
    nnls_weights = nnls::nnls(fit_cv,train_data$Y)$x
    nnls_weights = nnls_weights / sum(nnls_weights)
    
    # Run all methods on the full sample
    fit_full = ensemble_core(ml,train_data, test_data,weights=weights,quiet=quiet)
    best = fit_full$predictions[,which.min(mse_cv)]
    ensemble = fit_full$predictions %*% nnls_weights
    
    # Calculate Smoothing matrix if weights=TRUE
    w = NULL
    if (isTRUE(weights)) {
      w = matrix(0,nrow(test_data),nrow(train_data))
      for (i in 1:length(ml)) {
        w = w + nnls_weights[i] * fit_full$weights[[i]]
      }
      w = Matrix::Matrix(w,sparse=T)
    }
    
    colnames(fit_full$predictions) = names(mse_cv) = names(nnls_weights) = colnames(fit_cv)
  }
  else { # in case only one method defined, no cross-validation needed
    fit_full = ensemble_core(ml,train_data,test_data,weights=weights,quiet=quiet)
    ensemble = best = fit_full$predictions
    w = nnls_weights = mse_cv = fit_cv = NULL
    if (isTRUE(weights)) w = fit_full$weights[[1]]
  }
  
  output = list("ensemble" = ensemble,"best" = best,"fit_full" = fit_full,"weights" = w,
                "nnls_weights" = nnls_weights, "mse_cv" = mse_cv, "fit_cv" = fit_cv)
  class(output) = "ensemble"
  return(output)
}
ensemble_core = function(ml,
                         tr, te,
                         weights=FALSE,quiet=TRUE) {
  # Initialize objects to be filled
  fit_mat = matrix(NA,nrow(te),length(ml))
  weights_list = vector("list",length(ml))
  
  # Loop over specified methods
  for (i in 1:length(ml)) {
    wrapper = paste0(ml[[i]]$method,"_fit")
    if (isFALSE(quiet)) print(wrapper)
    
    # Check whether subset of variables specified and/or additional arguments are defined and run method
    if (is.null(ml[[i]]$x_select) & length(ml[[i]]$args) == 0)          prediction = do.call(wrapper,list(train_data=tr,test_data=te, weights=weights))
    else if (is.null(ml[[i]]$x_select) & !(length(ml[[i]]$args) == 0))  prediction = do.call(wrapper,list(train_data=tr,test_data=te,args=ml[[i]]$args, weights=weights))
    else if (!is.null(ml[[i]]$x_select) & length(ml[[i]]$args) == 0)    prediction = do.call(wrapper,list(train_data=tr[,ml[[i]]$x_select],test_data=te[,ml[[i]]$x_select]))
    else                                                                prediction = do.call(wrapper,list(train_data=tr[,ml[[i]]$x_select],test_data=te[,ml[[i]]$x_select],args=ml[[i]]$args, weights=weights))
    # Get predictions
    # if (is.null(ml[[i]]$x_select))  temp = do.call(paste0("predict.",wrapper),list(fit,x=x_tr,y=y_tr,xnew=x_te,weights=weights))
    # else                            temp = do.call(paste0("predict.",wrapper),list(fit,x=x_tr[,ml[[i]]$x_select],y=y_tr,xnew=x_te[,ml[[i]]$x_select],weights=weights))
    # 
    fit_mat[,i] = prediction
    #  weights_list[[i]] = temp$weights
    weights_list <- NULL
  }
  
  return(list("predictions" = fit_mat, "weights" = weights_list))
}


create_method = function(method,
                         x_select=NULL,
                         args=list(),
                         name=NULL) {
  
  if (!(is.character(method) & length(method) == 1)) stop("Provide single string to define method.")
  if (!(any(method == c("glmm_ols", "low_rank_ols", "embedding_ols", "group_lasso", "sparse_group_lasso", "fused_lasso")))
  ) stop("Provide one of these options c(\"glmm_ols\",\"low_rank_ols\",\"embedding_ols\",\"group_lasso\",\"sparse_group_lasso\",\"fused_lasso\").")
  if (!(is.null(args) | is.list(args))) stop("Provide either NULL or list for args.")
  if (!(is.null(x_select) | is.logical(x_select))) stop("Provide either NULL or logical for x_select.")
  if (!((is.character(name) & length(name) == 1) | is.null(name))) stop("Provide single string to name method.")
  
  return(list(method=method,args=args,x_select=x_select,name=name,weights=weights))
}
# ML-WRAPPER ####
glmm_ols_fit <- function(train_data, test_data, n.folds=5){
  flagNewLvls = function(fac, lvls) {
    char = as.character(fac)
    char[!(char %in% lvls | is.na(char))] = "..new..level.."
    factor(char)
  }
  fitLmer = function(feature, target) {
    # performance tips from https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
    nlopt <- function(par, fn, lower, upper, control) {
      .nloptr <<- res <- nloptr::nloptr(par, fn, lb = lower, ub = upper, 
                                        opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                                    maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
      list(par = res$solution,
           fval = res$objective,
           conv = if (res$status > 0) 0 else res$status,
           message = res$message)
    }
    
    args = list(formula = y ~ 1 + (1 | lvl),
                data = data.frame(lvl = feature, y = target),
                na.action = na.omit,
                control = lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
    
    mod = do.call(lme4::lmer, args)
    coefs = coef(mod)$lvl
    lvls = rownames(coefs)
    coefs = coefs[,1]
    names(coefs) = lvls
    intercept = unname(lme4::fixef(mod))
    # replace missing coefs with intercept value
    coefs[is.na(coefs)] = intercept
    # save intercept value for new levels during prediction
    coefs = c(coefs, ..new..level.. = intercept)
    coefs
  }
  
  # Split data to do the encoding
  sep_Y <- which(colnames(train_data)=='Y')
  cat <- which(grepl("C", colnames(train_data)))
  cont <- which(grepl("X", colnames(train_data)))
  
  train.cont <- train_data[,cont]
  test.cont <- test_data[,cont]
  
  train.cat <- train_data[,cat]
  test.cat <- test_data[,cat]
  
  if (is.null(ncol(train.cat))){
    train.cat <- train.cat %>% as.data.frame() %>% magrittr::set_colnames("C1")
    test.cat <- test.cat %>% as.data.frame() %>% magrittr::set_colnames("C1")
  }
  
  train.Y <- train_data[,sep_Y]
  test.Y <- test_data[,sep_Y]
  
  
  # for prediction, use complete encoding model
  control = lapply(train.cat, function(col) {
    fitLmer(col, train.Y)
  })
  
  # if n.folds == 1 use only the complete model in training
  if (n.folds == 1) {
    train.cat = lapply(colnames(train.cat), function(cname) {
      as.numeric(control[[cname]][as.character(train.cat[[cname]])])
    })
  }
  
  # else use n.folds encoding models in crossvalidation scheme in training
  if(n.folds > 1){
    library(mlrCPO)
    rinst = makeResampleInstance(makeResampleDesc("CV", iters = n.folds), 
                                 size = nrow(train.cat))
    # list with n.folds models for each data column
    mod_list = lapply(colnames(train.cat), function(cname) {
      lapply(rinst$train.inds, function(inds) {
        fitLmer(train.cat[[cname]][inds], train.Y[inds])
      })
    })
    names(mod_list) = names(train.cat)
    # list with replaced values in n.folds for each data column
    num_vals_list = lapply(colnames(train.cat), function(cname) {
      lapply(seq_len(n.folds), function(fold) {
        # add new level for all values not present during training
        fac = flagNewLvls(train.cat[[cname]][rinst$test.inds[[fold]]],
                          names(mod_list[[cname]][[fold]]))
        mod = mod_list[[cname]][[fold]]
        
        as.numeric(mod[as.character(fac)])
      })
    })
    names(num_vals_list) = names(mod_list)
    
    # recombine replaced values from n.folds for each data column
    train.cat[] = lapply(seq_along(train.cat), function(cnum) {
      unlist(num_vals_list[[cnum]])[order(unlist(rinst$test.inds))]
    })
  }
  
  # For the test data, use the complete model
  test.cat[] = lapply(colnames(test.cat), function(cname) {
    as.numeric(control[[cname]][as.character(test.cat[[cname]])])
  })
  
  # merge data back together
  train.cat <- train.cat %>% as.data.frame() %>% magrittr::set_colnames(colnames(train_data)[cat])
  train.Y <- train.Y %>% as.data.frame() %>% magrittr::set_colnames("Y")
  train_enc <- data.frame(train.Y, train.cat, train.cont)
  
  test.cat <- test.cat %>% as.data.frame() %>% magrittr::set_colnames(colnames(test_data)[cat])
  test.Y <- test.Y %>% as.data.frame() %>% magrittr::set_colnames("Y")
  test_enc <- data.frame(test.Y, test.cat, test.cont)
  
  # with columns with a lot of levels, na's may happen
  na_col <- apply(is.na(test_enc),2, any)
  if(any(na_col)==TRUE){
    for(na in which(na_col)){
      test_enc[which(is.na(test_enc[,na])), na] <- mean(test_enc[,na],na.rm=T)
    }
  }
  
  fit <- glm(Y~., as.data.frame(train_enc), family = 'gaussian')
  
  # new data
  test.yhat <- predict(fit, newdata = as.data.frame(test_enc))
  return(test.yhat)
}
group_lasso_fit <- function(train_data, test_data, already_encoded = FALSE, group = NULL){
  library(gglasso)
  if  (already_encoded == TRUE & is.null(group))(stop("If the categorical variables are already encoded, the variable GROUP must be provided."))
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  if (already_encoded == FALSE){
    colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
    train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
    train.y <- train_data[,sep_Y]
    
    colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
    test.x <- model.matrix(~., test_data[,-sep_Y], )[,-1]
    test.y <- test_data[,sep_Y]
    
    group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  } else {
    train.x <- as.matrix(train_data[, -sep_Y])
    train.y <- train_data[, sep_Y]
    test.x <- as.matrix(test_data[, -sep_Y])
    test.y <- test_data[, sep_Y]
  }
  
  # cross-validated model
  gglasso.cv <- cv.gglasso(x = train.x, y = train.y, group = group, eps = exp(-4), nfolds = 10)
  #cbind(coef(gglasso.cv, s = c(gglasso.cv$lambda.1se, gglasso.cv$lambda.min)), trueBeta)
  
  # new data
  test.yhat <- predict(gglasso.cv, newx = data.matrix(test.x), s = gglasso.cv$lambda.min)
  
  return(test.yhat)
}
embedding_ols_fit <- function(train_data, test_data, dimensionality_reduction = "very strong", epochs=5, batchsize = 1){
  require(dplyr)
  require(keras)
  require(tensorflow)
  suppressPackageStartupMessages(require(purrr))
  
  x_train <- train_data %>% dplyr::select(-Y)
  x_test <- test_data %>% dplyr::select(-Y)
  y_train <- train_data[, "Y"]
  y_test <- test_data[, "Y"] 
  
  categorical <- colnames(train_data %>% select_if(is.factor))
  
  # We have to set the factor stings to numeric strings to use keras 
  # In order that the test data has the same new levels, although there are missing (dropped) or new levels in the test data, we create a mapping
  for(cc in categorical){
    oldLevels <- levels(x_train[,cc])
    newLevels <- 1:length(levels(x_train[,cc]))
    levels(x_train[,cc]) <- newLevels
    levels(x_test[,cc]) <- tidyr::replace_na(newLevels[match(levels(x_test[,cc]), oldLevels)],0)
  }
  
  # scaling (when scaling, you should also put scale(train.Y) in the model fitting step)
  train_means <- colMeans(x_train[sapply(x_train, is.double)]) %>% unname()
  train_sds <- apply(x_train[sapply(x_train, is.double)], 2, sd)  %>% unname()
  train_sds[train_sds == 0] <- 0.000001
  
  x_train[sapply(x_train, is.double)] <- sweep(
    x_train[sapply(x_train, is.double)], 2, train_means) %>%
    sweep(2, train_sds, "/")
  x_test[sapply(x_test, is.double)] <- sweep(
    x_test[sapply(x_test, is.double)], 2, train_means) %>%
    sweep(2, train_sds, "/")
  
  
  # number of levels per factor, required to specify input dimensionality for
  # the embedding layers
  # n_levels_in <- purrr::map(x_train %>% select_if(is.factor), compose(length, base::levels)) %>%
  #   unlist()
  n_levels_in <- sapply(x_train %>% select_if(is.factor), nlevels)
  
  # output dimensionality for the embedding layers, need +1 because Python is 0-based
  #n_levels_out <- n_levels_in %>% sqrt() %>% trunc() %>% `+`(1)
  n_levels_out <- switch(dimensionality_reduction, 
                         "weak" = n_levels_in%/%2 + 1,
                         "standard" = n_levels_in^(1/2) %>% trunc() %>% `+` (1),
                         "strong" = n_levels_in^(1/3) %>% trunc() %>% `+` (1),
                         "very strong" = n_levels_in^(1/4) %>% trunc() %>% `+` (1)
  )
  
  # each embedding layer gets its own input layer
  cat_inputs <- purrr::map(n_levels_in, function(l) layer_input(shape = 1)) %>%
    unname()
  
  # construct the embedding layers, connecting each to its input
  embedding_layers <- vector(mode = "list", length = length(cat_inputs))
  for (i in 1:length(cat_inputs)) {
    embedding_layer <-  cat_inputs[[i]] %>% 
      layer_embedding(input_dim = n_levels_in[[i]] + 1, 
                      output_dim = n_levels_out[[i]], name = paste("Embedding", categorical[i], sep = "")) %>%
      layer_flatten()
    embedding_layers[[i]] <- embedding_layer
  }
  
  # create a single input and a dense layer for the numeric data
  quant_input <- layer_input(shape = length(which(sapply(x_train, is.double))))
  
  quant_dense <- quant_input %>% layer_dense(ncol(train_data)-length(categorical)-1)
  
  intermediate_layers <- list(embedding_layers, list(quant_dense)) %>% flatten()
  inputs <- list(cat_inputs, list(quant_input)) #%>% flatten()
  
  
  output <- layer_concatenate(intermediate_layers) %>%
    layer_dense(units = 512, activation = "relu") %>% #256
    # layer_dense(units = 64, activation = "relu", kernel_initializer = "normal") %>% #256
    layer_dense(units = 1)
  
  model <- keras_model(inputs, output)
  model %>% compile(loss = "mean_squared_error", optimizer = optimizer_adam(lr=0.0002, decay = exp(-5))) #optimizer_adam(lr=0.001)
  
  
  # Each layer has to be feed into the model as list element
  xcont <- as.matrix(x_train %>% select_if(is.double))
  xcat <- apply(as.matrix(sapply(x_train[,categorical], as.integer)), 2, list)
  train_input <- list(xcat, xcont)
  
  hist <- model %>% fit(x = train_input, y = (y_train), epochs = epochs, 
                        batch_size = batchsize, verbose = 0)
  
  embeddings <- list()
  for (zz in 1:length(cat_inputs)){
    embedding <- get_layer(model, paste("Embedding", categorical[zz], sep=""))$get_weights() %>%
      as.data.frame()
    colnames(embedding) <- paste(categorical[zz], "ENC", 1:ncol(embedding), sep = "")
    embedding$name <- c("0", levels(x_train[,categorical[zz]]))
    embeddings[[zz]] <- embedding
  }
  
  train_enc <- data.frame("ID" =  1:length(y_train), "Y"=y_train, x_train)
  test_enc <- data.frame("ID" = 1:length(y_test), "Y"=y_test, x_test)
  for (j in 1:length(categorical)){
    train_enc <- merge(train_enc, embeddings[[j]], by.x=categorical[j], by.y="name",sort = F)[,-1]
    test_enc <- merge(test_enc, embeddings[[j]], by.x=categorical[j], by.y="name", sort=F)[,-1]
  }
  train_enc <- train_enc[order(train_enc$ID),-1]
  test_enc <- test_enc[order(test_enc$ID),-1]
  
  fit <- glm(Y~., train_enc, family = 'gaussian')
  
  # new data
  test.yhat <- predict(fit, newdata = test_enc)
  
  return(test.yhat)
}
low_rank_ols_fit <- function(train_data, test_data){
  most_lvls <- which.max(sapply(train_data, nlevels))
  label_train <- data.frame(sapply(train_data[,-most_lvls], as.numeric),train_data[,most_lvls])
  label_test <- data.frame(sapply(test_data[,-most_lvls], as.numeric),test_data[,most_lvls])
  most_lvls <- colnames(train_data)[most_lvls]
  colnames(label_train)[ncol(label_train)] <- most_lvls
  colnames(label_test)[ncol(label_test)] <- most_lvls
  num_components <- num_comp(label_train, most_lvls, "Y", "linear_regression", 5)
  low_rank_train <- low_rank_enc(label_train, most_lvls, num_components)
  low_rank_test <- low_rank_enc(label_test, most_lvls, num_components)
  
  fit <- glm(Y~., as.data.frame(low_rank_train), family = 'gaussian')
  
  # new data
  test.yhat <- predict(fit, newdata = as.data.frame(low_rank_test))
  return(test.yhat)
}
sparse_group_lasso_fit <- function(train_data, test_data, real_data=TRUE){
  library(seagull)
  if (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y])[,-1]
  test.y <- test_data[,sep_Y]
  
  #groups
  group <- attr(model.matrix(~., train_data[,-sep_Y]), "assign")[-1]
  
  #### Dumb workaround --> with real data, R crashes when 0<alpha<1, while finding the
  ####                     ideal lambda_max. So, we approximate it manually with alpha = 0
  ####                     and alpha = 1.
  ####                     fun fact: Works also better in the simulation
  
  if (real_data == TRUE){
    max_lambda <- max(seagull(y = train.y, Z = train.x, groups = group, alpha = 0, rel_acc = 0.1, xi = 0.001)[["lambda"]])
    #Lambda1 <- max(seagull(y = train.y, Z = train.x, groups = group, alpha = 1, rel_acc = 0.1)[["lambda"]])
    #max_lambda <- max(Lambda0, Lambda1)
  }
  
  
  # cross-validated model
  require(caret)
  flds <- createFolds(train.y, k = 5, list = F, returnTrain = FALSE)
  grid <- seq(0.05,0.95,0.3)
  mse.cv <- matrix(NA, ncol = length(grid), nrow = length(unique(flds)))
  colnames(mse.cv) <- grid
  
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    for (i in 1:length(grid)) {
      sgl <- seagull(y = try.cv, Z = trx.cv, groups = group, alpha = grid[i], max_lambda = max_lambda, rel_acc = 0.01, xi = 0.001)
      y.cv <- as.matrix(tex.cv) %*% t(as.matrix(sgl$random_effects))
      mse <- min(colMeans((tey.cv - y.cv)^2))
      mse.cv[j, i] <- mse
    }
  }
  alpha <- grid[which.min(colMeans(mse.cv))] # alpha = 0 -> Group Lasso | alpha = 1 -> Lasso
  
  
  # CV for lambda
  lambda.cv <- c()
  flds <- createFolds(train.y, k = 10, list = F, returnTrain = FALSE)
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    sgl <- seagull(y = try.cv, Z = trx.cv, groups = group, alpha = alpha, max_lambda = max_lambda, rel_acc = 0.005, xi = 0.001)
    # compute mse
    mse <- colMeans((as.matrix(tex.cv) %*% as.matrix(t(sgl$random_effects[,])) - tey.cv)^2)
    
    lambda.cv <-c(lambda.cv, which.min(mse))
  }
  lambda <- round(median(lambda.cv),0)
  
  sgl <- seagull(y = train.y, Z = train.x, groups = group, alpha = alpha, max_lambda = max_lambda, rel_acc = 0.00001, xi = 0.001)
  
  # new data
  test.yhat <- as.matrix(test.x) %*% as.matrix(sgl$random_effects[lambda,])
  
  print(paste("Alpha of Seagull-Lasso:", alpha))
  return(test.yhat)
} 
fused_lasso_fit <- function(train_data, test_data, order = T){
  library(genlasso)
  if  (!("Y" %in% colnames(train_data)))(stop("There is no column Y in the data"))
  if  (!("Y" %in% colnames(test_data)))(stop("There is no column Y in the data"))  
  sep_Y <- which(colnames(test_data) %in% "Y")
  
  
  colnames(train_data)[grepl("C", colnames(train_data))] <- paste(colnames(train_data)[grepl("C", colnames(train_data))], "ENC", sep = "")
  train.x <- model.matrix(~., train_data[,-sep_Y])[,-1]
  train.y <- train_data[,sep_Y]
  
  colnames(test_data)[grepl("C", colnames(test_data))] <- paste(colnames(test_data)[grepl("C", colnames(test_data))], "ENC", sep = "")
  test.x <- model.matrix(~., test_data[,-sep_Y])[,-1]
  test.y <- test_data[,sep_Y]
  
  if (order == TRUE){
    corr <- cor(train.x)
    dis <- dist(corr)
    clust <- hclust(dis, method = "complete") #ward.D
    index <- clust$order
    
    train.x <- train.x[, index]
    test.x <- test.x[, index]
  }
  
  ###### Lambda via CV
  # cross-validated model
  require(caret)
  flds <- createFolds(train.y, k = 10, list = F, returnTrain = FALSE)
  lambda.cv <- c()
  
  
  for (j in 1:length(unique(flds))){
    trx.cv <- train.x[which(flds != j),]
    try.cv <- train.y[which(flds != j)]
    tex.cv <- train.x[which(flds == j),]
    tey.cv <- train.y[which(flds == j)]
    
    model <- fusedlasso1d(y=try.cv, X = trx.cv, approx = T, maxsteps = 2000)
    pred <- predict(model, Xnew = tex.cv)
    mse <- as.numeric(colMeans((tey.cv - pred$fit)^2))
    
    lambda.cv <- c(lambda.cv, model$lambda[which.min(mse)])
  } 
  
  lambda <- median(lambda.cv)
  
  #######
  fused <- fusedlasso1d(y=train.y, X = train.x, minlam = lambda/2)
  
   test.yhat <- predict(fused, Xnew = test.x, lambda = lambda)
  
  if(mean(abs(test.yhat$fit)) > 10 * mean(abs(pred$fit))){
    i <- 1.1
    while (mean(abs(test.yhat$fit)) > 10*mean(abs(pred$fit))){
      test.yhat <- predict(fused, Xnew = test.x, lambda = i*lambda)
      i <- i + 0.1
    }}
  
  test.yhat <- predict(fused, Xnew = test.x, lambda = lambda)$fit
  return(test.yhat)
}

#####