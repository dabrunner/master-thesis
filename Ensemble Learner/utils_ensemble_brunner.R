#' #' Creates the methods to be used in \code{\link{ensemble}}
#' #'
#' #' @param method Choose method from currently c("mean","ols",ridge","plasso",forest_grf","lasso")
#' #' @param x_select Optional logical vector of length equal to the number of columns of the covariate matrix
#' #' indicating which variables should be used by this method. E.g. tree-based methods usually should not be provided
#' #' with the interactions that Lasso is using
#' #' @param args Optional list containing the additional arguments that should be passed to the underlying method
#' #' @param name Optional sting naming the method, otherwise named as Method1, Method2, ...
#' #'
#' #' @return Method that can be passed to \code{\link{ensemble}}
#' 
#' create_method = function(method,
#'                          x_select=NULL,
#'                          args=list(),
#'                          name=NULL) {
#'   
#'   if (!(is.character(method) & length(method) == 1)) stop("Provide single string to define method.")
#'   if (!(any(method == c("mean","ols","ridge","plasso","forest_grf","lasso", "group_lasso", "fused_lasso")))
#'   ) stop("Provide one of these options c(\"mean\",\"ols\",\"ridge\",\"plasso\",\"forest_grf\",\"lasso\",\"group_lasso\",\"fused_lasso\").")
#'   if (!(is.null(args) | is.list(args))) stop("Provide either NULL or list for args.")
#'   if (!(is.null(x_select) | is.logical(x_select))) stop("Provide either NULL or logical for x_select.")
#'   if (!((is.character(name) & length(name) == 1) | is.null(name))) stop("Provide single string to name method.")
#'   
#'   list(method=method,args=args,x_select=x_select,name=name,weights=weights)
#' }



#' This function extends the covariate matrix with interactions, polynomials and logarithms
#'
#' @param data Covariate matrix to be extended
#' @param int Vector of strings with variable names to be interacted
#' @param int_o Order of interactions (default 2)
#' @param poly Vector of strings with variable names for which polynomials should be added
#' @param poly_o Order of polynomials to be added (default 2)
#' @param log Vector of strings with variable names for which the logarithm should be added
#'
#' @import stats matrixStats
#'
#' @return Extended covariate matrix
#' 
#' @keywords internal
#' 
design_matrix = function(data,
                         factor = NULL,
                         int=NULL,
                         int_o=2,
                         poly=NULL,
                         poly_o=2,
                         log=NULL) {
  


  # Part for factors (new)
  if(!is.null(factor)){
    int[int %in% factor] = paste0("factor(", factor, ")", sep="")
    factor = paste0("+factor(",paste0(factor,collapse=")+factor("),")")
  }
  
  
  # Part for interactions
  if (!is.null(int)) {
    int = paste0("+(",paste0(int,collapse="+"),")^",toString(int_o))
  }
  
  # Part for polynomials
  if (!is.null(poly)) {
    poly = paste0("+poly(",paste0(poly,collapse=paste0(",",toString(poly_o),",raw=TRUE)+poly(")),",",toString(poly_o),",raw=TRUE)",
                   "-(",paste0(paste0(poly,collapse="+")),")")
  }
  
  # Part for logs
  if (!is.null(log)) {
    # Check whether some variables can't be logged because not positive
    if(length(log)>1){
      ind_neg = colMins(data[,log])<=0} else {ind_neg = min(data[,log])<=0}
    if (sum(ind_neg)>0) {
      cat("\n Following variables not modified to be logged because of non-positive values:",paste(log[ind_neg]),"\n" )
      log = log[!ind_neg]
    }
    if (identical(log, character(0)) == FALSE) log = paste0("+log(",paste0(log,collapse=")+log("),")") else log = NULL
  }
  
  # Combine the three parts
  fmla = as.formula(paste("~0",int,poly,log))
  
  # Generate matrix
  data = model.matrix(fmla,data=as.data.frame(data))
  # Clean variable names to make sense
  colnames(data) = gsub("poly\\(","",colnames(data))
  colnames(data) = gsub(paste0(", ",toString(poly_o),", raw = TRUE)"),"",colnames(data))
  colnames(data) = gsub("log\\(","ln_",colnames(data))
  # Maybe the factors should be named differently (Right now it is the same as for polynomials)
  colnames(data) = gsub("factor\\(","",colnames(data))
  colnames(data) = gsub("\\)","",colnames(data))
  
  return(data)
}


#' This function takes a matrix of data and removes (i) variables without variation
#' (ii) dummy variables where one group is nearly empty (optional in one of both treatment groups), and
#' (iii) redundant (highly correlated variables)
#'
#' @param data Covariate matrix to be screened
#' @param bin_cut Cut-off for binary variables to be considered empty (default 0.01)
#' @param corr_cut Cut-off for too high correlation (default 0.99)
#' @param treat Optional binary treatment indicator to check binary variables in each treatment arm
#' @param print If TRUE, summary printed
#'
#' @importFrom  psych describe 
#' @import matrixStats
#'
#' @return Screened covariate matrix
#' 
#' @keywords internal
#' 
data_screen = function(data,
                       bin_cut=0.01,
                       corr_cut=0.99,
                       treat=NULL,
                       print=TRUE) {
  
  
  # Save the group assignment to eliminate redundant columns as well. 
  group <- attr(data, "assign")
  
  ## Kick out variables with no variation
  # Identify the names
  nm_del = colnames(data)[colSds(data) == 0]
  # Describe identified variables
  if (print==TRUE) {
    cat("\n\n Variables with no variation:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  # Remove identified variables (from data and groups)
  if( identical(nm_del, character(0)) == FALSE) group <- group[-which(colnames(data) %in% nm_del)]
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]
  
 
  
  ## Remove dummy variables lower than threshold in one of the two treatment groups
  # Identify dummies
  bin = apply(data,2,function(x) { all(x %in% 0:1) })

  # Calculate means of all variables and check whether they are potentially close to 0 or 1
  if (is.null(treat)) {
    mean = colMeans(data)
    bel_cut = (mean<bin_cut | mean > (1-bin_cut))
  } else {
    mean1 = colMeans(data[treat==1,])
    mean0 = colMeans(data[treat==0,])
    bel_cut = (mean1<bin_cut | mean1 > (1-bin_cut) | mean0<bin_cut | mean0 > (1-bin_cut))
  }
  
  # Identify names that are binary and close to 0 and 1
  nm_del = colnames(data)[bin & bel_cut]
  if (print==TRUE) {
    cat("\n\n Dummy variables close to 0 or 1:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  
  # Remove identified variables (from data & groups)
  if( identical(nm_del, character(0)) == FALSE) group <- group[-which(colnames(data) %in% nm_del)]
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]
  
  ## Remove all redundant (nearly perfectly correlated) variables
  # Calculate correlation matrix and consider only upper diagonal
  cor = (abs(cor(data))>corr_cut)
  cor[lower.tri(cor, diag=TRUE)] = FALSE
  
  # Identify names of redundant variables
  nm_del = colnames(cor)[colSums(cor)>0]
  
  if (print==TRUE) {
    cat("\n\n Variables (nearly) perfectly correlated:",nm_del,"\n\n")
    if (identical(nm_del, character(0)) == FALSE) print(describe(data[,nm_del],fast=TRUE))
  }
  
  # Remove identified variables (from data & groups)
  if( identical(nm_del, character(0)) == FALSE) group <- group[-which(colnames(data) %in% nm_del)]
  if (identical(nm_del, character(0)) == FALSE) data = data[,colSums(cor)==0]
  
  
  # The group membership must not have any steps larger than 1. This code eliminates gaps, if present
  for(i in length(group):2){
    if(group[i]-1>group[i-1]) group[i:length(group)] = group[i:length(group)]-(group[i]-group[i-1]-1)
  }
  # The new group membership is assigned to the model matrix
 attr(data, "assign") <- group
  
  return(data)
}


#' Matrix for cross-fitting indicators
#'
#' Creates matrix of binary fold indicators (n x # cross-folds)
#'
#' @param n Number of of observations
#' @param cf Number of cross-fitting folds
#' @param cl Optional vector of cluster variable if cross-fitting should account for clusters.
#'
#' @importFrom dplyr ntile
#' @import stats
#'
#' @return n times # cross-folds matrix of binary fold indicators
#'
#' @export

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









