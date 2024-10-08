% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/longpye_kernel_smoothing_with_covariates_same_betas.R
\name{longpye_KS_compute_cv_same_betas}
\alias{longpye_KS_compute_cv_same_betas}
\title{longpye_KS_compute_cv_same_betas}
\usage{
longpye_KS_compute_cv_same_betas(
  n_folds,
  df,
  X = names(df[, !(names(df) \%in\% c(y))]),
  y = "y",
  id = "id",
  t = "t",
  lambda,
  trace = 1,
  alpha = 0.5,
  a1 = 3.7,
  a2 = 3,
  penalty = "L1",
  regressors_betas = NULL,
  seed = 1,
  used_cores = 1,
  trend = "monotone",
  delta = 1e-05,
  max_alpha = 10000,
  kernel = "gaussian",
  beta_start_input = NULL,
  max_iter = 20,
  min_alpha = 1e-09,
  convergence_error = 1e-07,
  stepsizeShrink = 0.5,
  beta_start_default = NULL,
  scaling = FALSE,
  c_zero_fixed = FALSE,
  measure_to_select_lambda = "ccr"
)
}
\arguments{
\item{n_folds}{number of fold of the cross validation}

\item{df}{the input dataset}

\item{X}{regressors to consider in the estimation. It can be of type
dataframe, containing also the same name of the regressors included in df,
of just a vector of character. Default is all not present in y and C}

\item{y}{the target variable. It can be only binomial 0,1. It can be of type
dataframe, containing also the same name of the same target variable included
in df, or just a character. Default is "y".}

\item{lambda}{the penalization parameter of the regressors X}

\item{trace}{2:visualize all the steps, 1:visualize just the result,
0:visualize nothing. Default is 1}

\item{alpha}{parameter for the Elastic-Net penalization term. Default is 0.5}

\item{a1}{parameter for the SCAD and MCP penalization term. Default is 3.7}

\item{a2}{parameter for the MCP penalization term. Default is 3.0}

\item{penalty}{the considered penalty. To be chosen among L12, L1, EN, SCAD
and MCP. Default is "L1"}

\item{regressors_betas}{a vector containing the real betas (if known). Default
is NULL, i.e. we do not know the real regressors}

\item{seed}{fix the seed. Default is 1}

\item{used_cores}{number of core used for the parallelization of the
process. if equal to 1, then no parallelization is adopted. Default is 1.}

\item{trend}{if "monotone", mmAPG is used, if "nonmonotone", mnmAPG is used.
Default is "monotone"}

\item{delta}{parameter for the convergence condition of the optimization
algorithm. Default is 1e-5}

\item{max_alpha}{maximum value of the step-parameter alpha. Default is 1000}

\item{kernel}{the kernel type to use for the estimation of the density
function (tested only for "gaussian").  Default is "gaussian"}

\item{beta_start_input}{vector of a specific starting point for betas.
Default is NULL, i.e. no input vector}

\item{max_iter}{maximum number of iterations in the algorithms mmAPG and
mnmAPG. Default is 10000}

\item{min_alpha}{minimum value of the step-parameter alpha. Default is 1e-12}

\item{convergence_error}{error to accept for considering the algorithm
converged. Default is 1e-5}

\item{stepsizeShrink}{parameter to adjust the step-size in the backtracking
line-search, in the optimization of pye. Taking values between 0 and 1,
the closer to 1, the more accurate the estimation will be, the longer it
will take and viceversa. Default is 0.8}

\item{beta_start_default}{set the default starting point of betas. If
"zeros", it starts with a vector of all zeros, if "corr" it starts with the
value of the correlation of every regressor with the target variable y.
Default is "zeros"}

\item{scaling}{if TRUE, the dataset is scaled. FALSE otherwise, Default is
FALSE.}

\item{c_zero_fixed}{if TRUE the estimation process considers c, the cut-off
point, as fixed and equal to zero, to reduce the complexity of the estimation.
If FALSE c can vary and be different from zero, estimated by pye. Default
is FALSE}

\item{measure_to_select_lambda}{the measure used to select lambda in the
cross-validation process. Default is "ccr", i.e. correct classification rate}

\item{C}{covariate variables. It can be of type dataframe, containing the
same covariates included in df, or just a vector of character. Default is NULL}
}
\value{
a list containing the optimal value of lambda to estimate betas and
c, the value of the main accuracy measure for all the folds.
}
\description{
function to perform the cross-validation to select the best
value of lambda for the estimation of betas and c (and
possibly gammas) using pye (and possibly covYI).
}
\examples{
library(pye)
cols <- 2000
cols_cov <- 20
seed=1
simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
X <- simMicroarrayData_cov02_dim50_covariates$X
y <- simMicroarrayData_cov02_dim50_covariates$y
C <- simMicroarrayData_cov02_dim50_covariates$C
regressors_betas<-simMicroarrayData_cov02_dim50_covariates$nregressors
regressors_gammas<-simMicroarrayData_cov02_dim50_covariates$ncovariates
penalty <- "SCAD"
c <- 0
prox_penalty = get(paste0("proximal_operator_", penalty))
trend = "monotone" #or "nonmonotone"
pye_starting_point = "zeros"
alpha = 0.5
c_zero_fixed <- FALSE
used_cores <- 1
used_penalty_pye <- c("L1", "MCP") #c("L12", "L1", "EN", "SCAD", "MCP")
max_iter <- 10
n_folds <- 3

#pye Gaussian (and others) Kernel Smooth
for (p in used_penalty_pye){

  name <- paste0("param_estimate_PYE_KS_covYI_", p)

  #wrapper of the function
  wrapper_lambda <- function(lambda){
    return(pye_KS_estimation(df=df, X=X, penalty=p, trend=trend, trace=1,
      beta_start_default=pye_starting_point, beta_start_input=NULL,
      lambda=lambda, alpha=alpha, a1=3.7, a2=3,
      regressors_betas=regressors_betas, c_zero_fixed=c_zero_fixed,
      max_iter=max_iter))
  }

  #lambda to max and min
  lambda_max <- calibrate_lambda_max (function_to_run=wrapper_lambda,
    var_to_check=paste0("betas_hat_",p), lambda_start = 5)
  lambda_min <- calibrate_lambda_min (function_to_run=wrapper_lambda,
    var_to_check=paste0("betas_hat_",p), lambda_start = 0.0001, max_var=100)

 cat("\n Pye penalty:", p, "; lambda_max=", lambda_max, "; lambda_min=",
  lambda_min, "\n")

  #create a suited lambda
  lambda <- create_lambda(n=3, lmax=lambda_max, lmin=lambda_min)
  lambda <- as.numeric(formatC(lambda, format = "e", digits = 9))

  #start cross-validation
  assign(name, pye_KS_compute_cv(penalty=p, df=df,
    X=names(df[,!(names(df) \%in\% c(y,C))]), y=y, C=C, trace=1,
    beta_start_default=pye_starting_point, trend = trend,
    beta_start_input=NULL, n_folds=n_folds, lambda=lambda,
    alpha=alpha, a1=3.7, a2=3,
    regressors_betas=regressors_betas, regressors_gammas=regressors_gammas,
    used_cores=used_cores, kernel="gaussian", c_zero_fixed=c_zero_fixed,
    penalty_g=p, max_iter=max_iter, max_iter_g=max_iter))

  #take the best lambda per measures (in case of same measure for diff lambdas,
  #we take the one associated to less betas)
  assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_yi"), get(name)$lambda_hat_pye_KS_yi)
  assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_auc"), get(name)$lambda_hat_pye_KS_auc)
  assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_ccr"), get(name)$lambda_hat_pye_KS_ccr)
  assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_gm"), get(name)$lambda_hat_pye_KS_gm)
  assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_pye"), get(name)$lambda_hat_pye_KS_pye)
}


}
