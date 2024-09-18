#scale and center the variables
#df: dataset to scale,
#y: target variable (not to scale)
#x: regressors (to be scaled)


#' @title scaling_df_for_pye
#'
#' @description function to scale and center variables, not applied to the
#' intercept, dummys and non-regressors variables.
#'
#' @param n_folds number of fold of the cross validation
#' @param df the input dataset
#' @param x regressors to be scaled. The function is not applied to the
#' intercept, dummys and non-regressors variables.
#' @param y the target variable, to not scale (we are in a classification
#' environment)
#'
#' @return the scaled database
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' X <- simMicroarrayData_cov02_dim50_covariates$X
#' y <- simMicroarrayData_cov02_dim50_covariates$y
#' C <- simMicroarrayData_cov02_dim50_covariates$C
#'
#' df <- scaling_df_for_pye (df=df, X=colnames(df[, names(df) %in% c(X,C)]), y="y")
#' print(df)
#'
#' @export
scaling_df_for_pye <- function(df, X, y){

  #(not applyied to the intercept, dummys and non-regressors variables)
  to_not_scale = c(y, names(df[,!(colnames(df) %in% c(X))]), names(which(sapply(X, function(z) length(unique(df[,z])))==2)))

  original_var_mean = colMeans(df[!(colnames(df) %in% to_not_scale)])
  original_var_stdev = sapply(df[!(colnames(df) %in% to_not_scale)], stats::sd)
  df_scaled = df
  df_scaled[,!(colnames(df) %in% to_not_scale)] <- as.data.frame(scale(df[,!(colnames(df) %in% to_not_scale)], center = TRUE, scale = TRUE))
  #scaling--DONE

  return(list(df_scaled=df_scaled, original_mu=original_var_mean, original_stdev=original_var_stdev))
}
