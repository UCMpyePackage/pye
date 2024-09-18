#' @title calibrate_lambda_max
#'
#' @description Function to calibrate the starting value of lambda to make
#' cross-validation more effective. In particular, this function aim to solve
#' the problem that every model has its own "best" range of lambdas at which
#' the parameters are different from 0 and performs better. The function stops
#' where the number of elements of var_to_check is greater then 0.
#'
#' @param function_to_run Function that needs to be cross-validated. It needs to
#' have a have just one parameter to change (theoretically called lambda),
#' changing which the number of selected variables might change.
#' @param var_to_check Name of the variable (a vector) to check if, reducing
#' lambda, the number of its elements increases.
#' @param lambda_start Starting value of lambda. If left big might take longer
#' (default is 5)
#' @param lambda_min Minimum value of lambda to try (default is 0.1)
#' @param n_min_var Minimum number of variables at which the process stops
#' @param factor Number to multiply or divide the considered lambda in the
#' step. The higher, the slower to find the optimal value. Default is 2.
#'
#' @return The first encountered value of the parameter (lambda inside the
#' function) at which the number of elements of var_to_check is greater then 0.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' #create the wrapper
#' penalty <- "L1"
#' wrapper <- function(lambda){
#'    return(pye_KS_estimation(df=df, penalty=penalty, trace=2, lambda=lambda,
#'    c_zero_fixed=FALSE, max_iter=5))
#' }
#' #find the lambda max
#' var_to_check <- paste0("betas_hat_", penalty)
#' lambda_max <- calibrate_lambda_max (function_to_run=wrapper,
#'                var_to_check=var_to_check, lambda_start = 5, lambda_min=0.1)
#' print(lambda_max)
#'
#' @export
calibrate_lambda_max <- function(function_to_run, var_to_check="betas_hat",
                                 lambda_start = 5, lambda_inf=1e-20, n_min_var=0,
                                 factor=2){
  cat("\n Starting calibrating lambda MAX \n")
  fun <- match.fun(function_to_run)
  function_to_run <- function(lambda) fun(lambda)
  lambda0 <- lambda_start
  lambda_max <- lambda_start
  #first evaluation
  eval <- getElement(function_to_run(lambda=lambda0), var_to_check)
  result0 <- length(eval[eval!=0])
  convergence <- FALSE
  #till convergence
  while ((convergence == FALSE) & (lambda0 >= lambda_inf)) {
    #if lambda is too high
    if (result0 <= n_min_var){
      lambda1 <- lambda0/factor
      eval <- getElement(function_to_run(lambda=lambda1), var_to_check)
      result1 <- length(eval[eval!=0])
      lambda_max <- lambda1*5 #for sure the max lambda is at least this
      if (result1 > n_min_var){
        convergence <- TRUE
      } else {
        lambda0 <- lambda1
        result0 <- result1
      }
      #if lambda is too low
    } else {
      lambda1 <- lambda0*factor
      eval <- getElement(function_to_run(lambda=lambda1), var_to_check)
      result1 <- length(eval[eval!=0])
      if (result1 <=n_min_var){
        convergence <- TRUE
        lambda_max <- lambda1
      } else {
        if (lambda1 >= 10){
          convergence <- TRUE
          lambda_max <- lambda1
          cat("I stop since lambda_max is >= 10 \n")
        } else {
          lambda0=lambda1
          result0=result1
        }
      }
    }
  }
  if(lambda0 < lambda_inf) {
    convergence <- TRUE
    lambda_max = length(lambda0)/100000
	warning("calibrate_lambda_max did not converge to a valid lambda max parameter.
         We set lambda_max=", lambda_max, "as a general value.")
  }
  return(lambda_max)
}

#-------------------------------------------------------------------------------

#' @title calibrate_lambda_min
#'
#' @description Function to calibrate the last used value of lambda, to make
#' cross-validation more effective. In particular, this function aim to solve
#' the problem that every model has its own "best" range of lambdas at which
#' the parameters are different from 0 and performs better. The function stops
#' where the number of elements of var_to_check is greater then 0.
#'
#' @param function_to_run Function that needs to be cross-validated. It needs to
#' have a have just one parameter to change (theoretically called lambda),
#' changing which the number of selected variables might change.
#' @param var_to_check Name of the variable (a vector) to check if, reducing
#' lambda, the number of its elements increases.
#' @param lambda_start Starting value of lambda. If left too small it might take
#' longer (default is 0.0001)
#' @param max_var Maximum number of variable to accept the result of the
#' function (default is 100.000)
#' @param factor Number to multiply or divide the considered lambda in the
#' step. The higher, the slower to find the optimal value. Default is 2.
#'
#' @return The first encountered value of the parameter (lambda inside the
#' function) at which the number of elements of var_to_check is lower then
#' max_var.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' #create the wrapper
#' penalty <- "L1"
#' wrapper <- function(lambda){
#'    return(pye_KS_estimation(df=df, penalty=penalty, trace=2, lambda=lambda,
#'    c_zero_fixed=FALSE, max_iter=5))
#' }
#' #find the lambda min
#' var_to_check <- paste0("betas_hat_", penalty)
#' #accept maximum 100 betas
#' max_var <- 100
#' lambda_min <- calibrate_lambda_min (function_to_run=wrapper,
#'   var_to_check=var_to_check, lambda_start = 0.0001, max_var = max_var)
#' print(lambda_min)
#'
#' @export

calibrate_lambda_min <- function(function_to_run, var_to_check="betas_hat",
                                 lambda_start = 0.0001, lambda_inf=1e-20, max_var = 1000,
                                 factor=2){
  cat("\n Starting calibrating lambda MIN \n")
  fun <- match.fun(function_to_run)
  function_to_run <- function(lambda) fun(lambda)
  lambda0 <- lambda_start
  lambda_min <- lambda_start
  #first evaluation with lambda=0
  eval0 <- getElement(function_to_run(lambda=0), var_to_check)
  result <- min(length(eval0[eval0!=0])-1, max_var) # min among the maximum number of betas -1 (to select) and max_var
  if (result==0){stop("The number of parameter of the model without penalty is zero! Something is wrong!")}
  eval1 <- getElement(function_to_run(lambda=lambda0), var_to_check)
  result0 <- length(eval1[eval1!=0])
  convergence <- FALSE
  #till convergence
  while ((convergence == FALSE) & (lambda0 >= lambda_inf)) {
    #if lambda is too high
    if ((result0 >= result)){
      lambda1 <- lambda0*factor
      eval2 <- getElement(function_to_run(lambda=lambda1), var_to_check)
      result1 <- length(eval2[eval2!=0])
      if (result1 < result){
        convergence <- TRUE
        lambda_min <- lambda1/5
      } else {
        lambda0 <- lambda1
        result0 <- result1
      }
      #if lambda is too low
    } else {
      lambda1 <- lambda0/factor
      eval2 <- getElement(function_to_run(lambda=lambda1), var_to_check)
      result1 <- length(eval2[eval2!=0])
      if (result1 >= result){
        convergence <- TRUE
        lambda_min <- lambda1
      } else {
        lambda0 <- lambda1
        result0 <- result1
      }
    }
  }
  if(lambda0 < lambda_inf) {
    convergence <- TRUE
    lambda_min = length(lambda0)/10000000
	  warning("calibrate_lambda_min did not converge to a valid lambda min parameter.
         We set lambda_min=", lambda_min, "as a general value.")
  }
  return(lambda_min)
}
