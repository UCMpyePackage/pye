% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmAPG2_and_mnmAPG2_algorithms.R
\name{mmAPG2}
\alias{mmAPG2}
\title{mmAPG}
\usage{
mmAPG2(
  x0,
  g0,
  delta_fx_betas,
  delta_fx_gammas,
  proxx_betas,
  proxx_gammas,
  Fx,
  lambda = NULL,
  tau = NULL,
  penalty = NULL,
  fold = NULL,
  stepsizeShrink = 0.8,
  max_alpha = 10000,
  min_alpha = 1e-10,
  delta = 1e-05,
  trace = 2,
  seed = 1,
  max_iter = 10000,
  convergence_error = 1e-07,
  zeros_stay_zeros_from_iteration = 20,
  max.print = 10
)
}
\arguments{
\item{x0}{starting point of the regressor vector (to encourage a fast
convergence to a sparse result we suggest to use the zero vector as starting
point)}

\item{g0}{starting point of the covariate vector (to encourage a fast
convergence to a sparse result we suggest to use the zero vector as starting
point)}

\item{delta_fx_betas}{gradient of betas}

\item{delta_fx_gammas}{gradient of gammas}

\item{proxx_betas}{proximal operator related to gx fo the betas}

\item{proxx_gammas}{proximal operator related to gx fo the gammas}

\item{Fx}{fx + gx}

\item{lambda}{the penalization parameter of the regressors X, related to gx,
just for the visualization in the trace report. Default is NULL}

\item{tau}{the penalization parameter of the covariates C, related to gx,
just for the visualization in the trace report. Default is NULL}

\item{penalty}{penalty parameter, related to gx, just for the visualization
in the trace report. Default is NULL}

\item{fold}{fold in which you are doing the cross-validation, just for the
visualization in the trace report. Default is NULL}

\item{stepsizeShrink}{parameter to adjust the step-size in the backtracking
line-search, in the optimization of pye. Taking values between 0 and 1,
the closer to 1, the more accurate the estimation will be, the longer it
will take and viceversa. Default is 0.8}

\item{max_alpha}{maximum value of the step-parameter alpha. Default is 1000}

\item{min_alpha}{minimum value of the step-parameter alpha. Default is 1e-12}

\item{delta}{parameter for the convergence condition of the optimization
algorithm. Default is 1e-5}

\item{trace}{2:visualize all the steps, 1:visualize just the result,
0:visualize nothing. Default is 2}

\item{seed}{fix the seed. Default is 1}

\item{max_iter}{maximum number of iterations in the algorithms mmAPG and
mnmAPG. Default is 500}

\item{convergence_error}{error to accept for considering the algorithm
converged. Default is 1e-5}

\item{zeros_stay_zeros_from_iteration}{the number of iteration from which
each parameter that reached zero in the estimation cannot change anymore.
This is done to preserve  the sparsity of the solution. Default is 5}

\item{max.print}{number of elements to show if printing the results. Default
is 10}
}
\value{
a list with x1 the estimated vector of parameters, tot_iters the
total number of iterations and backtrack_iters the number of backtracking
iterations.
}
\description{
Inspired by the paper Li and Lin, 2015, "Accelerated Proximal
Gradient Methods for Nonconvex Programming", this is the monotonone
optimization algorithm used in PYE. It is a modified version of the original
of Li and Lin: Changes have been proposed in the initialization of the
algorithm and in the selection part. With respect of this second point:
we did not give the possibility to come back to deleted variables, i.e. the
selection is not reversible
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
ID <- rownames(df)
const <- rep(1, length(ID)) #the constant for the coefficients gammas
df1 <- cbind(ID, df[,(names(df) \%in\% c(y,X))], const, df[,(names(df) \%in\% C)])
C1 <- c("const", C) #adding "const" in C
penalty <- "L12"
alpha <- 0.5
lambda <- 0.1
tau <- 0.05
max_iter <- 10
kernel <- "gaussian"
if (penalty == "SCAD") {a=3.7} else {a=3.0}
prox_penalty <- get(paste0("proximal_operator_", penalty))

#wrappers
delta_fx_betas <- function(betas, gammas){
  result <- -cPYE_KS(df=df1[,!(names(df1) \%in\% c("ID"))], X=X, y=y, C=C1,
    betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
    alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
  result <- result[names(result) \%in\% X]
  return(result)
}

delta_fx_gammas <- function(betas, gammas){
  result <- -cPYE_KS(df=df1[,!(names(df1) \%in\% c("ID"))], X=X, y=y, C=C1,
    betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
    alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
  result <- result[names(result) \%in\% C1]
  return(result)
}

proxx_betas <- function(x, eta){
  result <- prox_penalty(betas=x, lambda=eta*lambda, alpha=alpha, a=a)
  return(result)
}

proxx_gammas <- function(x, eta){
  result <- prox_penalty(betas=x, lambda=eta*tau, alpha=alpha, a=a)
  return(result)
}

Fx <- function(x,z){
  result <- -getElement(cPYE_KS(df=df1[,!(names(df1) \%in\% c("ID"))], X=X,
    y=y, C=C1, betas=x, gammas=z, lambda=lambda, tau=tau, kernel=kernel,
    alpha=alpha, a1=3.7, a2=3.0, penalty=penalty), paste0("cPYE_", penalty))
  return(result)
}
#starting points:
x0 <- c(rep(0, length(X)))
names(x0) <- c(X)
g0 <- c(0, rep(0, length(C))) #the first zero is the constant term
names(g0) <- C1

estim <- mmAPG2(x0=x0, g0=g0, delta_fx_betas=delta_fx_betas,
  delta_fx_gammas=delta_fx_gammas, proxx_betas=proxx_betas,
  proxx_gammas=proxx_gammas, Fx=Fx, lambda=lambda, tau=tau, penalty=penalty,
  max_iter=max_iter, trace=2, seed=1)
print(estim)
}
