#' Calculates arithmetic mean.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns list containing mean and number of observations
#'
#' @keywords internal
#'
mean_fit = function(x,y) {
  
  mean = mean(y)
  output = list("mean"=mean,"n"=nrow(x))
  output
}

#' Predicts arithmetic mean and provides prediction weights if required.
#' @param mean_fit Output of \code{\link{mean_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.mean_fit = function(mean_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) fit = rep(mean_fit$mean,nrow(x))
  else fit = rep(mean_fit$mean,nrow(xnew))

  if (isTRUE(weights)) w = matrix(1 / length(y),nrow(xnew),nrow(x))
  else w = NULL

  list("prediction"=fit,"weights"=w)
}


#' Calculates OLS fit.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns OLS coefficients
#'
#' @keywords internal
#'
ols_fit = function(x,y) {
  
  x = cbind(rep(1,nrow(x)),x)
  ols_coef = lm.fit(x,y)$coefficients
  ols_coef
}


#' Prediction based on OLS and provides prediction weights if required.
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.ols_fit = function(ols_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  x = cbind(rep(1,nrow(x)),x)
  xnew = cbind(rep(1,nrow(xnew)),xnew)

  # Remove variables that were dropped due to collinearity
  x = x[,!is.na(ols_fit)]
  xnew = xnew[,!is.na(ols_fit)]

  # Calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x),tol=2.225074e-308) %*% t(x)
  fit = hat_mat %*% y

  if (weights==FALSE) hat_mat = NULL

  list("prediction"=fit,"weights"=hat_mat)
}


#' This function estimates cross-validated ridge regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
#' @keywords internal
#'
ridge_fit = function(x,y,args=list()) {
  
  ridge = do.call(cv.glmnet,c(list(x=x,y=y,alpha=0),args))
  ridge
}

#' Prediction based on Ridge and provides prediction weights if required.
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.ridge_fit = function(ridge_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  fit = predict(ridge_fit,newx=xnew,type="response")

  if (weights==FALSE) hat_mat = NULL
  else {
    # Get covariate matrices
    n = nrow(x)
    p = ncol(x)
    x = scale(x)
    x = cbind(rep(1,nrow(x)),x)
    xnew = scale(xnew)
    xnew = cbind(rep(1,nrow(xnew)),xnew)

    # Calculate hat matrix, see also (https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat)
    hat_mat = xnew %*% solve(crossprod(x) + ridge_fit$lambda.min  * n / sd(y) * diag(x = c(0, rep(1,p)))) %*% t(x)
    fit = hat_mat %*% y
  }

  list("prediction"=fit,"weights"=hat_mat)
}


#' This function estimates cross-validated Post-Lasso based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "plasso"
#'
#' @keywords internal
#'
plasso_fit = function(x,y,args=list()) {
  
  plasso = do.call(plasso,c(list(x=x,y=y),args))
  plasso
}


#' Prediction based on Post-Lasso and provides prediction weights if required.
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.plasso_fit = function(plasso_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  x = add_intercept(x)
  xnew = add_intercept(xnew)

  # Fitted values for post lasso
  nm_act = names(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl])[which(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl] != 0)]

  xact = x[,nm_act,drop=F]
  xactnew = xnew[,nm_act,drop=F]

  # Remove potentially collinear variables
  coef = lm.fit(xact,y)$coefficients
  xact = xact[,!is.na(coef)]
  xactnew = xactnew[,!is.na(coef)]

  hat_mat = xactnew %*% solve(crossprod(xact),tol=2.225074e-308) %*% t(xact)
  fit_plasso = hat_mat %*% y
  if (weights==FALSE) hat_mat = NULL

  list("prediction"=fit_plasso,"weights"=hat_mat)
}


#' Calculates Random Forest fit using the \code{\link{grf}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{regression_forest}}
#' @import grf
#'
#' @return An object with S3 class "regression_forest"
#'
#' @keywords internal
#'
forest_grf_fit = function(x,y,args=list()) {
  rf = do.call(regression_forest,c(list(X=x,Y=y),args))
  rf
}


#' Prediction based on Random Forest and provides prediction weights if required.
#' @param forest_grf_fit Output of \code{\link{regression_forest}} or \code{\link{forest_grf_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.forest_grf_fit = function(forest_grf_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x

  fit = predict(forest_grf_fit,newdata=xnew)$prediction

  if (weights==TRUE) w = get_forest_weights(forest_grf_fit,newdata=xnew)
  else w = NULL

  list("prediction"=fit,"weights"=w)
}




#' This function estimates cross-validated lasso regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
#' @keywords internal
#'
lasso_fit = function(x,y, args=list()) {

  lasso = do.call(cv.glmnet,c(list(x=x,y=y),args))
  lasso
}


#' Prediction based on Lasso Forest and provides prediction weights if required.
#' @param lasso_fit Output of \code{\link{glmnet}} or \code{\link{lasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not available for Lasso, only for Post-Lasso}
#'
#' @keywords internal
#'
predict.lasso_fit = function(lasso_fit,x,y,xnew=NULL,weights=FALSE) {
  f = function() stop("No weighted representation of Lasso available.",call.=FALSE)
  if (isTRUE(weights)) f()
  if (is.null(xnew)) xnew = x
  
  
  fit = predict(lasso_fit,newx=xnew,type="response",s="lambda.min")

  list("prediction"=fit,"weights"="No weighted representation of Lasso available.")
}


#' This function estimates cross-validated group lasso regression based on the \code{\link{gglasso}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param one_hot_encoded boolean indicating if x is one-hot encoded
#' @param group vector indicating group memberships (same number for each covariate in a group)
#' @param args List of arguments passed to  \code{\link{gglasso}}
#' @import gglasso
#'
#' @return An object with S3 class "gglasso"
#'
#' @keywords internal
#'
group_lasso_fit <- function(x, y, args=list()){
 
  # cross-validated model  (the parameter "groups" is feed into the function through args)
  gglasso.cv <- do.call(cv.gglasso, c(list(x = x, y = y), args))
  gglasso.cv
}


#' Prediction based on group lasso
#' 
#' @param group_lasso_fit Output of \code{\link{gglasso}} or \code{\link{group_lasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not available for group lasso}
#'
#' @keywords internal
#'
predict.group_lasso_fit = function(group_lasso_fit,x,y,xnew=NULL,weights=FALSE) {
  f = function() stop("No weighted representation of Group Lasso available.",call.=FALSE)
  if (isTRUE(weights)) f()
  if (is.null(xnew)) xnew = x
  
  # predictions of the group lasso
  fit <- predict(group_lasso_fit, newx = xnew, s = group_lasso_fit$lambda.min)
  
  list("prediction"=fit,"weights"="No weighted representation of Lasso available.")
}


#' This function estimates cross-validated fused lasso regression based on the \code{\link{genlasso}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param one_hot_encoded boolean indicating if the predictors should be ordered (hierarchical clustering).
#' @param ordering vector indicating group memberships (same number for each covariate in a group)
#' @param args List of arguments passed to  \code{\link{genlasso}}
#' @import genlasso
#'
#' @return An object with S3 class "genlasso"
#'
#' @keywords internal
#'
fused_lasso_fit <- function(x, y, args=list()){
  
  fused <- do.call(fusedlasso1d, c(list(y=y, X = x), args))
  fused
}


#' Prediction based on fused lasso
#' 
#' @param fused_lasso_fit Output of \code{\link{genlasso}} or \code{\link{fused_lasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#' @import genlasso caret
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not available for fused lasso}
#'
#' @keywords internal
#'
predict.fused_lasso_fit = function(fused_lasso_fit,x,y,xnew=NULL,weights=FALSE) {
  f = function() stop("No weighted representation of Fused Lasso available.",call.=FALSE)
  if (isTRUE(weights)) f()
  if (is.null(xnew)) xnew = x
  
  # Finding an appropriate penalization parameter by cross-validation
  flds <- createFolds(y, k = 5, list = F, returnTrain = FALSE)
  lambda.cv <- c()

  for (j in 1:length(unique(flds))){
    trx.cv <- x[which(flds != j),]
    try.cv <- y[which(flds != j)]
    tex.cv <- x[which(flds == j),]
    tey.cv <- y[which(flds == j)]

    model <- fusedlasso1d(y=try.cv, X = trx.cv)
    pred <- predict(model, Xnew = tex.cv)$fit

    mse <- as.numeric(colMeans((tey.cv - pred)^2))
    lambda.cv <- c(lambda.cv, model$lambda[which.min(mse)])
  }
  lambda <- mean(lambda.cv)

  # Fit the model with the cross-validated lambda
   fit <- predict(fused_lasso_fit, Xnew = xnew, lambda = lambda)$fit

  
   # In some rare cases, the predictions are +/-10^30, when lambda is under a certain threshold 
   # To avoid that weird bug, the predictions are compared to them from the cross-validation.
   # If they are ten times larger, the lambda is increased 10 Percent, until the predictions
   # are "normal". 
  if(mean(abs(fit)) > 10 * mean(abs(pred))){
    i <- 1.1
    while (mean(abs(fit)) > 10*mean(abs(pred))){
      fit <- predict(fused_lasso_fit, Xnew = xnew, lambda = i*lambda)$fit
      i <- i + 0.1
    }}
  
  list("prediction"=fit,"weights"="No weighted representation of Fused Lasso available.")
}

