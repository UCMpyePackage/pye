#' @title mmAPG
#'
#' @description Inspired by the paper Li and Lin, 2015, "Accelerated Proximal
#' Gradient Methods for Nonconvex Programming", this is the monotonone
#' optimization algorithm used in PYE. It is a modified version of the original
#' of Li and Lin: Changes have been proposed in the initialization of the
#' algorithm and in the selection part. With respect of this second point:
#' we did not give the possibility to come back to deleted variables, i.e. the
#' selection is not reversible
#'
#' @param x0 starting point of the regressor vector (to encourage a fast
#' convergence to a sparse result we suggest to use the zero vector as starting
#' point)
#' @param g0 starting point of the covariate vector (to encourage a fast
#' convergence to a sparse result we suggest to use the zero vector as starting
#' point)
#' @param delta_fx_betas gradient of betas
#' @param delta_fx_gammas gradient of gammas
#' @param proxx_betas proximal operator related to gx fo the betas
#' @param proxx_gammas proximal operator related to gx fo the gammas
#' @param Fx fx + gx
#' @param lambda the penalization parameter of the regressors X, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param tau the penalization parameter of the covariates C, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param penalty penalty parameter, related to gx, just for the visualization
#' in the trace report. Default is NULL
#' @param fold fold in which you are doing the cross-validation, just for the
#' visualization in the trace report. Default is NULL
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 2
#' @param seed fix the seed. Default is 1
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param zeros_stay_zeros_from_iteration the number of iteration from which
#' each parameter that reached zero in the estimation cannot change anymore.
#' This is done to preserve  the sparsity of the solution. Default is 5
#' @param max.print number of elements to show if printing the results. Default
#' is 10
#'
#' @return a list with x1 the estimated vector of parameters, tot_iters the
#' total number of iterations and backtrack_iters the number of backtracking
#' iterations.
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
#' ID <- rownames(df)
#' const <- rep(1, length(ID)) #the constant for the coefficients gammas
#' df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], const, df[,(names(df) %in% C)])
#' C1 <- c("const", C) #adding "const" in C
#' penalty <- "L12"
#' alpha <- 0.5
#' lambda <- 0.1
#' tau <- 0.05
#' max_iter <- 100
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx_betas <- function(betas, gammas){
#'   result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X, y=y, C=C1,
#'     betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
#'   result <- result[names(result) %in% X]
#'   return(result)
#' }
#'
#' delta_fx_gammas <- function(betas, gammas){
#'   result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X, y=y, C=C1,
#'     betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
#'   result <- result[names(result) %in% C1]
#'   return(result)
#' }
#'
#' proxx_betas <- function(x, eta){
#'   result <- prox_penalty(betas=x, lambda=eta*lambda, alpha=alpha, a=a)
#'   return(result)
#' }
#'
#' proxx_gammas <- function(x, eta){
#'   result <- prox_penalty(betas=x, lambda=eta*tau, alpha=alpha, a=a)
#'   return(result)
#' }
#'
#' Fx <- function(x,z){
#'   result <- -getElement(cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X,
#'     y=y, C=C1, betas=x, gammas=z, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty), paste0("cPYE_", penalty))
#'   return(result)
#' }
#' #starting points:
#' x0 <- c(rep(0, length(X)))
#' names(x0) <- c(X)
#' g0 <- c(0, rep(0, length(C))) #the first zero is the constant term
#' names(g0) <- C1
#'
#' estim <- mmAPG2(x0=x0, g0=g0, delta_fx_betas=delta_fx_betas,
#'   delta_fx_gammas=delta_fx_gammas, proxx_betas=proxx_betas,
#'   proxx_gammas=proxx_gammas, Fx=Fx, lambda=lambda, tau=tau, penalty=penalty,
#'   max_iter=max_iter, trace=2, seed=1)
#' print(estim)
#' @export

#2 steps Monotone APG with line search
#for simplicity we test only with the same penalty betas and c but it could be possible to test
#the result also with different penalties for betas and c.
mmAPG2 <- function(x0, g0, delta_fx_betas, delta_fx_gammas, proxx_betas, proxx_gammas, Fx,
                   lambda=NULL, tau=NULL, penalty=NULL, fold=NULL, stepsizeShrink=0.8,
                   max_alpha=10000, min_alpha=1e-10, delta=1e-5, trace=2, seed=1, max_iter=10000,
                   convergence_error=1e-7, zeros_stay_zeros_from_iteration=20, max.print=10){

  start_time <- Sys.time()
  #dont print all long results
  options(max.print = max.print)

  t1 <- 1
  t0 <- 0

  z0 <- x1 <- x0
  h0 <- g1 <- g0

  #counter
  i=1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)){

    y1 <- x1 + (t0/t1)*(z0 - x1) + ((t0-1)/t1)*(x1 - x0)
    yy1 <- g1 + (t0/t1)*(h0 - g1) + ((t0-1)/t1)*(g1 - g0) #for gammas (the same in all the operations)
    #keep not considering the lastly deleted betas
    if (i>zeros_stay_zeros_from_iteration){
      y1[x1==0] <- 0 #we leave c the possibility vary
      yy1[-1][g1[-1]==0] <- 0 #we leave to the const the possibility vary
    }

    if (i==1){ #the first iteration z0 does not exists
      s1 <- z0
      ss1 <- h0
    } else {
      s1 <- z0 - y0
      ss1 <- h0 - yy0
    }
    dfx_z0 <- delta_fx_betas(z0, h0)
    dfx_h0 <- delta_fx_gammas(z0, h0)
    if (i==1){
      r1 <- dfx_z0
      rr1 <- dfx_h0
    }else{
      dfx_y0 <- delta_fx_betas(y0, yy0)
      dfx_yy0 <- delta_fx_gammas(y0, yy0)
      r1 <- dfx_z0 - dfx_y0
      rr1 <- dfx_h0 - dfx_yy0
    }
    if (sum(s1)==0){
      alpha_y=1
    } else {
      #print((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1))))
      alpha_y = min(max_alpha, abs((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }
    if (sum(ss1)==0){
      alpha_yy=1
    } else {
      #print((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1))))
      alpha_yy = min(max_alpha, abs((t(ss1)%*% t(t(ss1)))/(t(ss1)%*%t(t(rr1)))), na.rm=T) #abs since sometimes might be negative
    }
    # alternatively: alpha_y = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_y)){stop("alpha_y is NaN -> Go and discover why!")}
    if (is.nan(alpha_yy)){stop("alpha_yy is NaN -> Go and discover why!")}


    if (i==1){
      s1 <- - x0
      ss1 <- - g0
    } else {
      s1 <- v0 - x0
      ss1 <- vv0 - g0
    }
    dfx_x0 <- delta_fx_betas(x0, g0)
    dfx_g0 <- delta_fx_gammas(x0, g0)
    if (i==1){
      r1 <- - dfx_x0
      rr1 <- - dfx_g0
    }else{
      dfx_v0 <- delta_fx_betas(v0, vv0)
      dfx_vv0 <- delta_fx_gammas(v0, vv0)
      r1 <- dfx_v0 - dfx_x0
      rr1 <- dfx_vv0 - dfx_g0
    }
    if (sum(s1)==0){
      alpha_x = 1
    } else {
      #print((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1))))
      alpha_x = min(max_alpha, abs((t(s1)%*%t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }
    if (sum(ss1)==0){
      alpha_g = 1
    } else {
      #print((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1))))
      alpha_g = min(max_alpha, abs((t(ss1)%*%t(t(ss1)))/(t(ss1)%*%t(t(rr1)))), na.rm=T) #abs since sometimes might be negative
    }
    # alternatively: alpha_x = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_x)){stop("alpha_x is NaN -> Go and discover why!")}
    if (is.nan(alpha_g)){stop("alpha_g is NaN -> Go and discover why!")}

    #cat("alpha_y:", alpha_y)
    #cat("; alpha_x:", alpha_x,"\n")

    dfx_y1 <- delta_fx_betas(y1, yy1)
    dfx_yy1 <- delta_fx_gammas(y1, yy1)
    Fx_yyy1 <- Fx(y1, yy1)
    #line search 1
    while (TRUE) {

      z1 <- proxx_betas(y1 - alpha_y*dfx_y1, alpha_y)
      h1 <- proxx_gammas(yy1 - alpha_yy*dfx_yy1, alpha_yy)
      if (i>zeros_stay_zeros_from_iteration){
        z1[y1==0] <- 0
        h1[-1][yy1[-1]==0] <- 0
      }
      #print(z1[z1!=0])

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_zh1 <- Fx(z1, h1)
      barrier <- Fx_yyy1 - delta*norm(c(z1-y1, h1-yy1), type="2")^2
      condition <- (round(Fx_zh1, 10) <= round(barrier, 10))

      #cat("Fx_zh1", Fx_zh1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_y <- alpha_y * stepsizeShrink
      alpha_yy <- alpha_yy * stepsizeShrink
    }

    dfx_x1 <- delta_fx_betas(x1, g1)
    dfx_g1 <- delta_fx_gammas(x1, g1)
    Fx_xg1 <- Fx(x1, g1)
    #line search 2
    while (TRUE) {

      v1 <- proxx_betas(x1 - alpha_x*dfx_x1, alpha_x)
      vv1 <- proxx_gammas(g1 - alpha_g*dfx_g1, alpha_g)
      if (i>zeros_stay_zeros_from_iteration){
        v1[x1==0] <- 0
        vv1[-1][g1[-1]==0] <- 0
      }

      #print(v1[v1!=0])

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_vvv1 <- Fx(v1, vv1)
      barrier <- Fx_xg1 - delta*norm(c(v1-x1, vv1-g1), type="2")^2
      condition <- (round(Fx_vvv1, 10) <= round(barrier, 10))

      #cat("Fx_vvv1", Fx_vvv1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_x <- alpha_x * stepsizeShrink
      alpha_g <- alpha_g * stepsizeShrink
    }

    #sum the loops of the line search
    totalBacktracks <- totalBacktracks + backtrackCount

    #results and one step forward
    t0 <- t1
    t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
    z0 <- z1
    h0 <- h1
    v0 <- v1
    vv0 <- vv1
    y0 <- y1
    yy0 <- yy1
    x0 <- x1
    g0 <- g1
    Fx_zh1 <- Fx(z1, h1)
    Fx_vvv1 <- Fx(v1, vv1)

    #cat("Fx_zh1",Fx_zh1, "\n")
    #cat("Fx_vvv1",Fx_vvv1, "\n")
    if(Fx_zh1 <= Fx_vvv1){
      x1 <- z1
      g1 <- h1
      Fx_xg1 <- Fx_zh1
    } else {
      x1 <- v1
      g1 <- vv1
      Fx_xg1 <- Fx_vvv1
    }

    #test
    if (i>zeros_stay_zeros_from_iteration){
      if(sum(x1)==0){
        x1 <- x1[1]
      } else {
        x1 <- x1[x1!=0]
      }
      if(sum(g1)==0){
        g1 <- g1[1]
      } else {
        g1 <- g1[g1!=0]
      }
      if (length(x0) != length(x1)){exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
      x0 <- x0[names(x0) %in% names(x1)]
      v0 <- v0[names(v0) %in% names(x1)]
      v1 <- v1[names(v1) %in% names(x1)]
      z0 <- z0[names(z0) %in% names(x1)]
      z1 <- z1[names(z1) %in% names(x1)]
      y0 <- y0[names(y0) %in% names(x1)]
      g0 <- g0[names(g0) %in% names(g1)]
      vv0 <- vv0[names(vv0) %in% names(g1)]
      vv1 <- vv1[names(vv1) %in% names(g1)]
      h0 <- h0[names(h0) %in% names(g1)]
      h1 <- h1[names(h1) %in% names(g1)]
      yy0 <- yy0[names(yy0) %in% names(g1)]
    }
    #end-test

    #print the partial results
    if (trace == 2){
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming 2 (APG2), the monotone version. \n")
      if (length(fold)!=0){
        cat("Fold:", fold, "; \n")
      }
      visualize_x1 = x1[which(x1!=0)]
      visualize_g1 = g1[which(g1!=0)]
      cat("iter:", i, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), ifelse(length(tau)!=0,paste0("; tau:", tau),""), "; alpha_y:", alpha_y, "; alpha_yy:", alpha_yy, "; alpha_x:", alpha_x, "; alpha_g:", alpha_g, "; F(x1):", Fx_xg1, "; \n")
      print(visualize_x1)
      print(visualize_g1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i==1){
      Fx_xg1_best <- Fx_xg1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_xg1_best - Fx_xg1 < convergence_error){
        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4){break}
      } else {
        Fx_xg1_best <- Fx_xg1
        exit_if_5_times_the_same <- 0
      }
    }

    #counters
    i=i+1

    if (i==max_iter){
	  if (trace %in% c(1,2)){
        cat("The maximum number of iteration has been reached. \n")
	  }
    }

    #if all alphas are very small, the new increment will be very poor, so we consider the algorithm converged
    if (sum(alpha_x, alpha_g, alpha_y, alpha_y) < min_alpha){break}
  }

  #print the partial results
  if (trace %in% c(2)){
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming 2 (APG2), the monotone version. \n")
    if (length(fold)!=0){
      cat("Fold:", fold, "; \n")
    }
    visualize_x1 = x1[which(x1!=0)]
    visualize_g1 = g1[which(g1!=0)]
    cat("iter:", i, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), ifelse(length(tau)!=0,paste0("; tau:", tau),""), "; alpha_y:", alpha_y, "; alpha_yy:", alpha_yy, "; alpha_x:", alpha_x, "; alpha_g:", alpha_g, "; F(xg1):", Fx_xg1, "; \n")
    print(visualize_x1)
    print(visualize_g1)
    cat("\n")
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(2)){
    print(estimation_time)
  }

  return(list(x1=x1, g1=g1, tot_iters = i, backtrack_iters = totalBacktracks, estimation_time=estimation_time))
}


#-------------------------------------------- NON-MONOTONE VERSION -------------------------------------------------



#' @title mnmAPG2
#'
#' @description Inspired by the paper Li and Lin, 2015, "Accelerated Proximal
#' Gradient Methods for Nonconvex Programming", this is the nonmonotonone
#' optimization algorithm used in PYE. It is a modified version of the original
#' of Li and Lin: changes have been proposed in the initialization of the
#' algorithm and in the selection part. With respect of this second point:
#' we did not give the possibility to come back to deleted variables, i.e. the
#' selection is not reversible
#'
#' @param x0 starting point (to encourage a fast convergence to a sparse result
#' we suggest to use the zero vector as starting point)
#' @param g0 starting point of the covariate vector (to encourage a fast
#' convergence to a sparse result we suggest to use the zero vector as starting
#' point)
#' @param delta_fx_betas gradient of betas
#' @param delta_fx_gammas gradient of gammas
#' @param proxx_betas proximal operator related to gx fo the betas
#' @param proxx_gammas proximal operator related to gx fo the gammas
#' @param Fx fx + gx
#' @param lambda the penalization parameter of the regressors X, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param tau the penalization parameter of the covariates C, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param penalty penalty parameter, related to gx, just for the visualization
#' in the trace report. Default is NULL
#' @param fold fold in which you are doing the cross-validation, just for the
#' visualization in the trace report. Default is NULL
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param eta parameter that control the degree of nonmonotonicity.
#' Li and Lin, 2015 set it to 0.8, as we do, i.e. default is 0.8
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 2
#' @param seed fix the seed. Default is 1
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param zeros_stay_zeros_from_iteration the number of iteration from which
#' each parameter that reached zero in the estimation cannot change anymore.
#' This is done to preserve  the sparsity of the solution. Default is 5
#' @param max.print number of elements to show if printing the results. Default
#' is 10
#'
#' @return a list with x1 the estimated vector of parameters, tot_iters the
#' total number of iterations and backtrack_iters the number of backtracking
#' iterations.
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
#' ID <- rownames(df)
#' const <- rep(1, length(ID)) #the constant for the coefficients gammas
#' df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], const, df[,(names(df) %in% C)])
#' C1 <- c("const", C) #adding "const" in C
#' penalty <- "L12"
#' alpha <- 0.5
#' lambda <- 0.1
#' tau <- 0.05
#' max_iter <- 100
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx_betas <- function(betas, gammas){
#'   result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X, y=y, C=C1,
#'     betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
#'   result <- result[names(result) %in% X]
#'   return(result)
#' }
#'
#' delta_fx_gammas <- function(betas, gammas){
#'   result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X, y=y, C=C1,
#'     betas=betas, gammas=gammas, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty)$gr_yi
#'   result <- result[names(result) %in% C1]
#'   return(result)
#' }
#'
#' proxx_betas <- function(x, eta){
#'   result <- prox_penalty(betas=x, lambda=eta*lambda, alpha=alpha, a=a)
#'   return(result)
#' }
#'
#' proxx_gammas <- function(x, eta){
#'   result <- prox_penalty(betas=x, lambda=eta*tau, alpha=alpha, a=a)
#'   return(result)
#' }
#'
#' Fx <- function(x,z){
#'   result <- -getElement(cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X,
#'     y=y, C=C1, betas=x, gammas=z, lambda=lambda, tau=tau, kernel=kernel,
#'     alpha=alpha, a1=3.7, a2=3.0, penalty=penalty), paste0("cPYE_", penalty))
#'   return(result)
#' }
#' #starting points:
#' x0 <- c(rep(0, length(X)))
#' names(x0) <- c(X)
#' g0 <- c(0, rep(0, length(C))) #the first zero is the constant term
#' names(g0) <- C1
#'
#' estim <- mnmAPG2(x0=x0, g0=g0, delta_fx_betas=delta_fx_betas,
#'   delta_fx_gammas=delta_fx_gammas, proxx_betas=proxx_betas,
#'   proxx_gammas=proxx_gammas, Fx=Fx, lambda=lambda, tau=tau, penalty=penalty,
#'   max_iter=max_iter, trace=2, seed=1)
#' print(estim)
#' @export

#Nonmonotone APG with line search
mnmAPG2 <- function(x0, g0, delta_fx_betas, delta_fx_gammas, proxx_betas, proxx_gammas, Fx,
                    lambda=NULL, tau=NULL, penalty=NULL, fold=NULL, stepsizeShrink=0.8,
                    eta=0.8, max_alpha=10000, min_alpha=1e-10, delta=1e-5, trace=2, seed=1,
                    max_iter=10000, convergence_error=1e-7, zeros_stay_zeros_from_iteration=20,
                    max.print=10){

  start_time <- Sys.time()
  #dont print all long results
  options(max.print = max.print)

  #initializations
  t1 <- 1
  t0 <- 0
  z0 <- x1 <- x0
  h0 <- g1 <- g0
  c1 <- Fx(x1, g1)
  q1 <- 1
  alpha_x <-0

  #counter
  i=1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)){

    y1 <- x1 + (t0/t1)*(z0 - x1) + ((t0-1)/t1)*(x1 - x0)
    yy1 <- g1 + (t0/t1)*(h0 - g1) + ((t0-1)/t1)*(g1 - g0)
    #keep not considering the lastly deleted betas
    if (i>zeros_stay_zeros_from_iteration){
      y1[x1==0] <- 0 #we leave c the possibility vary
      yy1[-1][g1[-1]==0] <- 0 #we leave to the const the possibility vary
    }

    if (i==1){
      s1 <- y1
      ss1 <- yy1
    } else {
      s1 <- y1 - y0
      ss1 <- yy1 - yy0
    }
    dfx_y1 <- delta_fx_betas(y1, yy1)
    dfx_yy1 <- delta_fx_gammas(y1, yy1)

    if (i==1){
      r1 <- dfx_y1
      rr1 <- dfx_yy1
    }else{
      dfx_y0 <- delta_fx_betas(y0, yy0)
      dfx_yy0 <- delta_fx_gammas(y0, yy0)
      r1 <- dfx_y1 - dfx_y0
      rr1 <- dfx_yy1 - dfx_yy0
    }
    if (sum(s1)==0){
      alpha_y = 1
    } else {
      alpha_y = min(max_alpha, abs((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }

    if (sum(ss1)==0){
      alpha_yy=1
    } else {
      alpha_yy = min(max_alpha, abs((t(ss1)%*% t(t(ss1)))/(t(ss1)%*%t(t(rr1)))), na.rm=T) #abs since sometimes might be negative
    }

    # alternatively: alpha_y = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_y)){stop("alpha_y is NaN -> Go and discover why!")}
    if (is.nan(alpha_yy)){stop("alpha_yy is NaN -> Go and discover why!")}

    Fx_yyy1 <- Fx(y1, yy1)
    #line search 1
    while (TRUE) {

      z1 <- proxx_betas(y1 - alpha_y*dfx_y1, alpha_y)
      h1 <- proxx_gammas(yy1 - alpha_yy*dfx_yy1, alpha_yy)
      if (i>zeros_stay_zeros_from_iteration){
        z1[y1==0] <- 0
        h1[-1][yy1[-1]==0] <- 0
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_zh1 <- Fx(z1, h1)
      barrier <- Fx_yyy1 - delta*norm(c(z1-y1, h1-yy1), type="2")^2
      condition <- (round(Fx_zh1, 10) <= round(barrier, 10))

      #cat("Fx_zh1", Fx_zh1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_y <- alpha_y * stepsizeShrink
      alpha_yy <- alpha_yy * stepsizeShrink
    }

    #second condition
    barrier <- c1 - delta*norm(c(z1-y1, h1-yy1), type="2")^2
    condition <- (round(Fx_zh1, 10) <= round(barrier, 10))

    #cat("Fx_zh1", Fx_zh1, "\n")
    #cat("barrier", barrier, "\n")
    #cat("condition", condition, "\n")

    if (condition){

      x0 <- x1
      g0 <- g1
      x1 <- z1
      g1 <- h1
      Fx_xg1 <- Fx_zh1

    } else {

      #set the stepsize alpha_x
      if (i==1){
        s1 <- x1
        ss1 <- g1
      } else {
        s1 <- x1 - y0
        ss1 <- g1 - yy0
      }
      dfx_x1 <- delta_fx_betas(x1, g1)
      dfx_g1 <- delta_fx_gammas(x1, g1)
      if (i==1){
        r1 <- dfx_x1
        rr1 <- dfx_g1
      }else{
        dfx_y0 <- delta_fx_betas(y0, yy0)
        dfx_yy0 <- delta_fx_gammas(y0, yy0)
        r1 <- dfx_x1 - dfx_y0
        rr1 <- dfx_g1 - dfx_yy0
      }
      if (sum(s1)==0){
        alpha_x = 1
      } else {
        alpha_x = min(max_alpha, abs((t(s1)%*%t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
      }

      if (sum(s1)==0){
        alpha_g = 1
      } else {
        alpha_g = min(max_alpha, abs((t(ss1)%*%t(t(ss1)))/(t(ss1)%*%t(t(rr1)))), na.rm=T) #abs since sometimes might be negative
      }
      # alternatively: alpha_x = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
      if (is.nan(alpha_x)){stop("alpha_x is NaN -> Go and discover why!")}
      if (is.nan(alpha_g)){stop("alpha_x is NaN -> Go and discover why!")}

      Fx_xg1 <- Fx(x1, g1)
      #line search 2
      while (TRUE) {

        v1 <- proxx_betas(x1 - alpha_x*dfx_x1, alpha_x)
        vv1 <- proxx_gammas(g1 - alpha_g*dfx_g1, alpha_g)
        if (i>zeros_stay_zeros_from_iteration){
          v1[x1==0] <- 0
          vv1[-1][g1[-1]==0] <- 0
        }

        backtrackCount <- backtrackCount + 1

        #condition
        Fx_vvv1 <- Fx(v1, vv1)
        barrier <- c1 - delta*norm(c(v1-x1, vv1-g1), type="2")^2
        condition <- (round(Fx_vvv1, 10) <= round(barrier, 10))

        #cat("Fx_vvv1", Fx_vvv1, "\n")
        #cat("barrier", barrier, "\n")
        #cat("condition", condition, "\n")

        if (condition == TRUE) {break}
        #if not, update alpha
        alpha_x <- alpha_x * stepsizeShrink
        alpha_g <- alpha_g * stepsizeShrink
      }

      x0 <- x1
      g0 <- g1
      Fx_zh1 <- Fx(z1, h1)
      Fx_vvv1 <- Fx(v1, vv1)

      #cat("Fx_zh1",Fx_zh1, "\n")
      #cat("Fx_vvv1",Fx_vvv1, "\n")

      if(Fx_zh1 <= Fx_vvv1){
        x1 <- z1
        g1 <- h1
        Fx_xg1 <- Fx_zh1
      } else {
        x1 <- v1
        g1 <- vv1
        Fx_xg1 <- Fx_vvv1
      }
    }

    #sum the loops of the line search
    totalBacktracks <- totalBacktracks + backtrackCount

    #results and one step forward
    t0 <- t1
    t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
    z0 <- z1
    h0 <- h1
    y0 <- y1
    yy0 <- yy1
    q0 <- q1
    q1 <- eta*q0 + 1
    c1 <- (eta*q0*c1 + Fx_xg1)/q1

    #test
    if (i>zeros_stay_zeros_from_iteration){
      if(sum(x1)==0){
        x1 <- x1[1]
      } else {
        x1 <- x1[x1!=0]
      }
      if(sum(g1)==0){
        g1 <- g1[1]
      } else {
        g1 <- g1[g1!=0]
      }
      if (length(x0) != length(x1)){exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
      x0 <- x0[names(x0) %in% names(x1)]
      z0 <- z0[names(z0) %in% names(x1)]
      z1 <- z1[names(z1) %in% names(x1)]
      y0 <- y0[names(y0) %in% names(x1)]
      g0 <- g0[names(g0) %in% names(g1)]
      h0 <- h0[names(h0) %in% names(g1)]
      h1 <- h1[names(h1) %in% names(g1)]
      yy0 <- yy0[names(yy0) %in% names(g1)]
    }
    #end-test

    #print the partial results
    if (trace == 2){
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming 2 (APG2), the nonmonotone version. \n")
      if (length(fold)!=0){
        cat("Fold:", fold, "; \n")
      }
      visualize_x1 = x1[which(x1!=0)]
      visualize_g1 = g1[which(g1!=0)]
      cat("iter:", i, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), ifelse(length(tau)!=0,paste0("; tau:", tau),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_xg1, "; \n")
      print(visualize_x1)
      print(visualize_g1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i==1){
      Fx_xg1_best <- Fx_xg1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_xg1_best - Fx_xg1 < convergence_error){
        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4){break}
      } else {
        Fx_xg1_best <- Fx_xg1
        exit_if_5_times_the_same <- 0
      }
    }

    #counters
    i=i+1

    if (i==max_iter){
	  if (trace %in% c(1,2)){
        cat("The maximum number of iteration has been reached. \n")
	  }
    }

    #if both the two alphas are very small, the new increment will be very poor, so we consider the algorithm converged
    if ((sum(alpha_x, alpha_y) < min_alpha) || ((0 < alpha_x) & (alpha_x < min_alpha))){break}
  }

  #print the partial results
  if (trace %in% c(1,2)){
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming 2 (APG2), the nonmonotone version. \n")
    if (length(fold)!=0){
      cat("Fold:", fold, "; \n")
    }
    visualize_x1 = x1[which(x1!=0)]
    visualize_g1 = g1[which(g1!=0)]
    cat("iters:", i, "; backtraking iters:", totalBacktracks, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), ifelse(length(tau)!=0,paste0("; tau:", tau),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(xg1):", Fx_xg1, "; \n")
    print(visualize_x1)
    print(visualize_g1)
    cat("\n")
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
  }

  return(list(x1=x1, g1=g1, tot_iters = i, backtrack_iters = totalBacktracks, estimation_time=estimation_time))
}

