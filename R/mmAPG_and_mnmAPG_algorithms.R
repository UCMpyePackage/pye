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
#' @param x0 starting point (to encourage a fast convergence to a sparse result
#' we suggest to use the zero vector as starting point)
#' @param c_pos position of the constant term. Equal it to NULL if you don't use
#' any constant.
#' @param delta_fx gradient of fx
#' @param proxx proximal operator related to gx
#' @param Fx fx + gx
#' @param lambda the penalization parameter of the regressors X, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param penalty penalty parameter, related to gx, just for the visualization
#' in the trace report. Default is NULL
#' @param fold fold in which you are doing the cross-validation, just for the
#' visualization in the trace report. Default is NULL
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 2
#' @param seed fix the seed. Default is 1
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
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
#' ID <- rownames(df)
#' df1 <- cbind(ID, df[,(names(df) %in% c(y,X))])
#' penalty <- "L12"
#' c_zero_fixed <- TRUE
#' lambda <- 0.1
#' max_iter <- 100
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx <- function(x){
#'  if (c_zero_fixed==TRUE){
#'    x[names(x) == "c"]<-0
#'   }
#'   result <- -pye_KS(df=df1[,names(df1)!="ID"],, X=X, y=y,
#'     betas=x[!(names(x) == "c")], lambda=lambda, c=x[(names(x) == "c")],
#'     kernel=kernel, alpha=0.5, a1=3.7, a2=3, penalty=penalty)$gr_yi
#'   if (c_zero_fixed==TRUE){
#'     result[names(result) == "c"] <- 0
#'   }
#'   return(result)
#' }
#'
#' proxx <- function(x, eta){
#'  if (c_zero_fixed==TRUE){
#'     x[names(x) == "c"]<-0
#'   }
#'   result <- c(prox_penalty(betas=x[!(names(x) == "c")], lambda=eta*lambda,
#'     alpha=0.5, a=a), x[(names(x) == "c")])
#'   return(result)
#' }
#'
#' Fx <- function(x){
#'   if (c_zero_fixed==TRUE){
#'     x[(names(x) == "c")]<-0
#'   }
#'   result <- -getElement(pye_KS(df=df1[,names(df1)!="ID"], X=X, y=y,
#'     betas=x[!(names(x) == "c")], lambda=lambda, c=x[(names(x) == "c")],
#'     kernel=kernel, alpha=0.5, a1=3.7, a2=3, penalty=penalty),
#'     paste0("pye_", penalty))
#'   return(result)
#' }
#'
#' #starting point:
#' x0 <- c(0, rep(0, length(X))) #the first zero is the constant term
#' names(x0) <- c("c",X)
#'
#' estim <- mmAPG(x0=x0, delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda,
#' penalty=penalty, max_iter=max_iter, trace=2)
#' print(estim)
#' @export

#Monotone APG with line search
mmAPG <- function(x0, c_pos=length(x0), delta_fx, proxx, Fx, lambda=NULL, penalty=NULL, fold=NULL, stepsizeShrink=0.8,
                                 max_alpha=10000, min_alpha=1e-10, delta=1e-5, trace=2, seed=1, max_iter=10000, #max_iter=500, the old max_iter
                                 convergence_error=1e-7, zeros_stay_zeros_from_iteration=20, max.print=10){

  start_time <- Sys.time()
  #dont print all long results
  options(max.print = max.print)

  if(length(names(x0))==0){stop("x0 needs to have a name vector corresponding to the variable names.")}

  preparing_the_solution <- rep(0, length(x0))
  names(preparing_the_solution) <- names(x0)

  t1 <- 1
  t0 <- 0
  z0 <- x1 <- x0

  #counter
  i=1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)){

    y1 <- x1 + (t0/t1)*(z0 - x1) + ((t0-1)/t1)*(x1 - x0)
    #keep not considering the lastly deleted betas
    zeros <- FALSE
    if (i>zeros_stay_zeros_from_iteration){
      c_proxy <- c_pos
      if (length(c_pos)==0){
        c_proxy <- 100000000000000000000000000
      }
      if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros){
      if (length(c_pos)==0){
        y1[x1==0] <- 0
      } else {
        y1[-c_pos][x1[-c_pos]==0] <- 0 #we leave c the possibility vary
      }
    }

    if (i==1){ #the first iteration y0 does not exists
      s1 <- z0
    } else {
      s1 <- z0 - y0
    }
    dfx_z0 <- delta_fx(z0)
    if (i==1){
      r1 <- dfx_z0
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_z0))<0.00000001){
        dfx_z0 <- dfx_z0*(10^(min(round(abs(log10(abs(dfx_z0))+1)))))
        }
    }else{
      dfx_y0 <- delta_fx(y0)
      r1 <- dfx_z0 - dfx_y0
    }
    if (sum(s1)==0){
      alpha_y=1
    } else {
      alpha_y <- min(max_alpha, abs((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }
    # alternatively: alpha_y = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_y)){stop("alpha_y is NaN -> Go and discover why!")}

    if (i==1){
      s1 <- - x0
    } else {
      s1 <- v0 - x0
    }
    dfx_x0 <- delta_fx(x0)
    if (i==1){
      r1 <- - dfx_x0
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_x0))<0.00000001){
        dfx_x0 <- dfx_x0*(10^(min(round(abs(log10(abs(dfx_x0))+1)))))
      }
    }else{
      dfx_v0 <- delta_fx(v0)
      r1 <- dfx_v0 - dfx_x0
    }
    if (sum(s1)==0){
      alpha_x = 1
    } else {
      #print((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1))))
      alpha_x = min(max_alpha, abs((t(s1)%*%t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }
    # alternatively: alpha_x = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_x)){stop("alpha_x is NaN -> Go and discover why!")}

    #cat("alpha_y:", alpha_y)
    #cat("; alpha_x:", alpha_x,"\n")

    dfx_y1 <- delta_fx(y1)
    if(i==1){
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_y1))<0.00000001){
        dfx_y1 <- dfx_y1*(10^(min(round(abs(log10(abs(dfx_y1))+1)))))
      }
    }
    Fx_y1 <- Fx(y1)
    #line search 1
    while (TRUE) {

      z1 <- proxx(y1 - alpha_y*dfx_y1, alpha_y)
      zeros <- FALSE
      if (i>zeros_stay_zeros_from_iteration){
        c_proxy <- c_pos
        if (length(c_pos)==0){
          c_proxy <- 100000000000000000000000000
        }
        if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
          zeros <- TRUE
        }
      }
      if (zeros){
        if (length(c_pos)==0){
          z1[y1==0] <- 0
        } else {
          z1[-c_pos][y1[-c_pos]==0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_z1 <- Fx(z1)
      barrier <- Fx_y1 - delta*norm(z1 - y1, type="2")^2
      condition <- (round(Fx_z1, 10) <= round(barrier, 10))

      #cat("Fx_z1", Fx_z1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_y <- alpha_y*stepsizeShrink
    }

    #If all betas are 0, c goes to 0
    #I put it here, out of the above loop, because otherwise it interfere with the optimization
    if (length(c_pos)!=0){
      if (sum(z1[-c_pos])==0) {z1[c_pos] <- 0}
    }


    dfx_x1 <- delta_fx(x1)
    if(i==1){
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_x1))<0.00000001){
        dfx_x1 <- dfx_x1*(10^(min(round(abs(log10(abs(dfx_x1))+1)))))
      }
    }
    Fx_x1 <- Fx(x1)
    #line search 2
    while (TRUE) {

      v1 <- proxx(x1 - alpha_x*dfx_x1, alpha_x)
      zeros <- FALSE
      if (i>zeros_stay_zeros_from_iteration){
        c_proxy <- c_pos
        if (length(c_pos)==0){
          c_proxy <- 100000000000000000000000000
        }
        if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
          zeros <- TRUE
        }
      }
      if (zeros){
        if (length(c_pos)==0){
          v1[x1==0] <- 0
        } else {
          v1[-c_pos][x1[-c_pos]==0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_v1 <- Fx(v1)
      #Fx_x1 <- Fx(x1)
      barrier <- Fx_x1 - delta*norm(v1 - x1, type="2")^2
      condition <- (round(Fx_v1, 10) <= round(barrier, 10))

      #cat("Fx_v1", Fx_v1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_x <- alpha_x * stepsizeShrink
    }

    #If all betas are 0, c goes to 0
    #I put it here, out of the above loop, because otherwise it interfere with the optimization
    if (length(c_pos)!=0){
      if (sum(v1[-c_pos])==0) {v1[c_pos]=0}
    }

    #sum the loops of the line search
    totalBacktracks <- totalBacktracks + backtrackCount

    #results and one step forward
    t0 <- t1
    t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
    z0 <- z1
    v0 <- v1
    y0 <- y1
    x0 <- x1
    Fx_z1 <- Fx(z1)
    Fx_v1 <- Fx(v1)

    #cat("Fx_z1",Fx_z1, "\n")
    #cat("Fx_v1",Fx_v1, "\n")
    if(Fx_z1 <= Fx_v1){
      x1 <- z1
      Fx_x1 <- Fx_z1
    } else {
      x1 <- v1
      Fx_x1 <- Fx_v1
    }

    #test
    zeros <- FALSE
    if (i>zeros_stay_zeros_from_iteration){
      c_proxy <- c_pos
      if (length(c_pos)==0){
        c_proxy <- 100000000000000000000000000
      }
      if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros){
      if (length(c_pos)==0){
        if(sum(x1)==0){
          x1 <- x1[1]
        } else {
          x1 <- x1[x1!=0]
        }
        if (length(x0) != length(x1)){exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
        x0 <- x0[names(x0) %in% names(x1)]
        v0 <- v0[names(v0) %in% names(x1)]
        v1 <- v1[names(v1) %in% names(x1)]
        z0 <- z0[names(z0) %in% names(x1)]
        z1 <- z1[names(z1) %in% names(x1)]
        y0 <- y0[names(y0) %in% names(x1)]
      }
    }
    #end-test

    #cat("\n z1",z1, " Fx_z1:", Fx_z1, "\n")
    #cat("\n v1",v1, " Fx_v1:", Fx_v1, "\n")
    #cat("\n x1",x1, "\n")

    #print the partial results
    if (trace == 2){
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the monotone version. \n")
      if (length(fold)!=0){
        cat("Fold:", fold, "; \n")
      }
      if (length(c_pos)==0){
        visualize_x1 <- x1[which(x1!=0)]
      } else {
        visualize_x1 <- c(x1[c_pos], x1[-c_pos][which(x1[-c_pos]!=0)]) #we leave c the possibility vary
      }
      cat("iter:", i, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
      print(visualize_x1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i==1){
      Fx_x1_best <- Fx_x1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_x1_best - Fx_x1 < convergence_error){
        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4){
          reason_for_exit <- "min_change"
          if (trace %in% c(1,2)){
            cat("mmAPG converged because the convergence error is below the threshold \n")
          }
          break
          }
      } else {
        Fx_x1_best <- Fx_x1
        exit_if_5_times_the_same <- 0
      }
    }

    #if both the two alphas are very small, the new increment will be very poor, so we consider the algorithm converged
    if (sum(alpha_x, alpha_y) < min_alpha){
      reason_for_exit <- "min_alpha"
      if (trace %in% c(1,2)){
        cat("mmAPG converged because the alpha is below the threshold \n")
      }
      break
    }

    #counters
    i=i+1

    if (i==max_iter){
      reason_for_exit <- "max_iter"
  	  if (trace %in% c(1,2)){
          cat("mnmAPG converged because the maximum number of iteration has been reached. \n")
  	  }
    }
  }

  preparing_the_solution[names(preparing_the_solution) %in% names(x1)] <- x1
  x1 <- preparing_the_solution

  #print the partial results
  if (trace %in% c(2)){
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the monotone version. \n")
    if (length(fold)!=0){
      cat("Fold:", fold, "; \n")
    }
    if (length(c_pos)==0){
      visualize_x1 = x1[which(x1!=0)]
    } else {
      visualize_x1 = c(x1[c_pos], x1[-c_pos][which(x1[-c_pos]!=0)] ) #we leave c the possibility vary
    }
    cat("iters:", i, "; backtraking iters:", totalBacktracks, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
    print(visualize_x1)
    cat("\n")
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(2)){
    print(estimation_time)
  }

  return(list(x1=x1, tot_iters=i, backtrack_iters=totalBacktracks, estimation_time=estimation_time,
              reason_for_exit=reason_for_exit))
}


#-------------------------------------------- NON-MONOTONE VERSION -------------------------------------------------



#' @title mnmAPG
#'
#' @description Inspired by the paper Li and Lin, 2015, "Accelerated Proximal
#' Gradient Methods for Nonconvex Programming", this is the monotonone
#' optimization algorithm used in PYE. It is a modified version of the original
#' of Li and Lin: Changes have been proposed in the initialization of the
#' algorithm and in the selection part. With respect of this second point:
#' we did not give the possibility to come back to deleted variables, i.e. the
#' selection is not reversible
#'
#' @param x0 starting point (to encourage a fast convergence to a sparse result
#' we suggest to use the zero vector as starting point)
#' @param c_pos position of the constant term. Equal it to NULL if you don't use any constant.
#' @param delta_fx gradient of fx
#' @param proxx proximal operator related to gx
#' @param Fx fx + gx
#' @param lambda the penalization parameter of the regressors X, related to gx,
#' just for the visualization in the trace report. Default is NULL
#' @param penalty penalty parameter, related to gx, just for the visualization
#' in the trace report. Default is NULL
#' @param fold fold in which you are doing the cross-validation, just for the
#' visualization in the trace report. Default is NULL
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 2
#' @param seed fix the seed. Default is 1
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param eta parameter that control the degree of nonmonotonicity
#' (Li and Lin, 2015 set it to 0.8, as we do)
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
#' ID <- rownames(df)
#' df1 <- cbind(ID, df[,(names(df) %in% c(y,X))])
#' penalty <- "L12"
#' c_zero_fixed <- TRUE
#' lambda <- 0.1
#' max_iter <- 100
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx <- function(x){
#'  if (c_zero_fixed==TRUE){
#'    x[names(x) == "c"]<-0
#'   }
#'   result <- -pye_KS(df=df1[,names(df1)!="ID"],, X=X, y=y,
#'   betas=x[!(names(x) == "c")], lambda=lambda, c=x[(names(x) == "c")],
#'   kernel=kernel, alpha=0.5, a1=3.7, a2=3, penalty=penalty)$gr_yi
#'   if (c_zero_fixed==TRUE){
#'     result[length(result)] <- 0
#'   }
#'   return(result)
#' }
#'
#' proxx <- function(x, eta){
#'  if (c_zero_fixed==TRUE){
#'     x[names(x) == "c"] <- 0
#'   }
#'   result <- c(prox_penalty(betas=x[!(names(x) == "c")], lambda=eta*lambda,
#'   alpha=0.5, a=a), x[(names(x) == "c")])
#'   return(result)
#' }
#'
#' Fx <- function(x){
#'   if (c_zero_fixed==TRUE){
#'     x[(names(x) == "c")]<-0
#'   }
#'   result <- -getElement(pye_KS(df=df1[,names(df1)!="ID"],
#'     betas=x[!(names(x) == "c")], lambda=lambda, c=x[(names(x) == "c")],
#'     kernel="gaussian", alpha=0.5, a1=3.7, a2=3, penalty=penalty),
#'     paste0("pye_", penalty))
#'   return(result)
#' }
#'
#' #starting point:
#' x0 <- c(0, rep(0, length(X))) #the first zero is the constant term
#' names(x0) <- c("c",X)
#'
#' estim <- mnmAPG(x0=x0, delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda,
#'   penalty=penalty, max_iter=max_iter, trace=2)
#' print(estim)
#' @export

#Nonmonotone APG with line search
mnmAPG <- function(x0, c_pos=length(x0), delta_fx, proxx, Fx, lambda=NULL, penalty=NULL, fold=NULL, stepsizeShrink=0.8,
                                    eta=0.8, max_alpha=10000, min_alpha=1e-10, delta=1e-5, trace=2, seed=1, max_iter=10000,
                                    convergence_error=1e-7, zeros_stay_zeros_from_iteration=20, max.print=10){

  start_time <- Sys.time()
  #dont print all long results
  options(max.print = max.print)

  if(length(names(x0))==0){stop("x0 needs to have a name vector corresponding to the variable names.")}

  preparing_the_solution <- rep(0, length(x0))
  names(preparing_the_solution) <- names(x0)

  #initializations
  t1 <- 1
  t0 <- 0
  z0 <- x1 <- x0
  c1 <- Fx(x1)
  q1 <- 1
  alpha_x <- 0

  #counter
  i=1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)){

    y1 <- x1 + (t0/t1)*(z0 - x1) + ((t0-1)/t1)*(x1 - x0)
    #keep not considering the lastly deleted betas
    zeros <- FALSE
    if (i>zeros_stay_zeros_from_iteration){
      c_proxy <- c_pos
      if (length(c_pos)==0){
        c_proxy <- 100000000000000000000000000
      }
      if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros){
      if (length(c_pos)==0){
        y1[x1==0] <- 0
      } else {
        y1[-c_pos][x1[-c_pos]==0] <- 0 #we leave c the possibility vary
      }
    }

    if (i==1){
      s1 <- y1
    } else {
      s1 <- y1 - y0
    }
    dfx_y1 <- delta_fx(y1)
    if (i==1){
      r1 <- dfx_y1
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_y1))<0.00000001){
        dfx_y1 <- dfx_y1*(10^(min(round(abs(log10(abs(dfx_y1))+1)))))
      }
    }else{
      dfx_y0 <- delta_fx(y0)
      r1 <- dfx_y1 - dfx_y0
    }
    if (sum(s1)==0){
      alpha_y = 1
    } else {
      alpha_y = min(max_alpha, abs((t(s1)%*% t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
    }
    # alternatively: alpha_y = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
    if (is.nan(alpha_y)){stop("alpha_y is NaN -> Go and discover why!")}

    Fx_y1 <- Fx(y1)
    #line search 1
    while (TRUE) {

      z1 <- proxx(y1 - alpha_y*dfx_y1, alpha_y)
      if (i>3){
        if (length(c_pos)==0){
          z1[y1==0] <- 0
        } else {
          z1[-c_pos][y1[-c_pos]==0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_z1 <- Fx(z1)
      barrier <- Fx_y1 - delta*norm(z1 - y1, type="2")^2
      condition <- (round(Fx_z1, 10) <= round(barrier, 10))

      #cat("Fx_z1", Fx_z1, "\n")
      #cat("barrier", barrier, "\n")
      #cat("condition", condition, "\n")

      if (condition == TRUE) {break}
      #if not, update alpha
      alpha_y <- alpha_y * stepsizeShrink
    }

    #If all betas are 0, c goes to 0
    #I put it here, out of the above loop, because otherwise it interfere with the optimization
    if (length(c_pos)!=0){
      if (sum(z1[-c_pos])==0) {z1[c_pos]=0}
    }

    #second condition
    barrier <- c1 - delta*norm(z1 - y1, type="2")^2
    condition <- (round(Fx_z1, 10) <= round(barrier, 10))

    #cat("Fx_z1", Fx_z1, "\n")
    #cat("barrier", barrier, "\n")
    #cat("condition", condition, "\n")

    if (condition){

      x0 <- x1
      x1 <- z1
      y0 <- y1
      Fx_x1 <- Fx_z1

      if(i==1){
        best_x1 <- z1
        best_Fx_x1 <- Fx_z1
      } else if(round(Fx_z1, 10) <= round(best_Fx_x1, 10)) {
        best_x1 <- z1
        best_Fx_x1 <- Fx_z1
        cat("\n I am updating the best x1: \n")
        print(best_x1)
        cat("\n")
      }

    } else {

      #set the stepsize alpha_x
      if (i==1){
        s1 <- x1
      } else {
        s1 <- x1 - y0
      }
      dfx_x1 <- delta_fx(x1)
      if (i==1){
        r1 <- dfx_x1
        #to try not to stuck to the initial zero point when the derivative is too low
        if (sum(abs(dfx_x1))<0.00000001){
          dfx_x1 <- dfx_x1*(10^(min(round(abs(log10(abs(dfx_x1))+1)))))
        }
      }else{
        dfx_y0 <- delta_fx(y0)
        r1 <- dfx_x1 - dfx_y0
      }
      if (sum(s1)==0){
        alpha_x = 1
      } else {
        alpha_x = min(max_alpha, abs((t(s1)%*%t(t(s1)))/(t(s1)%*%t(t(r1)))), na.rm=T) #abs since sometimes might be negative
      }
      # alternatively: alpha_x = (t(s1)%*%t(t(r1)))/(t(r1)%*%t(t(r1)))
      if (is.nan(alpha_x)){stop("alpha_x is NaN -> Go and discover why!")}

      #Fx_x1 <- Fx(x1)
      #line search 2
      while (TRUE) {

        v1 <- proxx(x1 - alpha_x*dfx_x1, alpha_x)
        zeros <- FALSE
        if (i>zeros_stay_zeros_from_iteration){
          c_proxy <- c_pos
          if (length(c_pos)==0){
            c_proxy <- 100000000000000000000000000
          }
          if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
            zeros <- TRUE
          }
        }
        if (zeros){
          if (length(c_pos)==0){
            v1[x1==0] <- 0
          } else {
            v1[-c_pos][x1[-c_pos]==0] <- 0 #we leave c the possibility vary
          }
        }

        backtrackCount <- backtrackCount + 1

        #condition
        Fx_v1 <- Fx(v1)
        barrier <- c1 - delta*norm(v1 - x1, type="2")^2
        condition <- (round(Fx_v1, 10) <= round(barrier, 10))

        #cat("Fx_v1", Fx_v1, "\n")
        #cat("barrier", barrier, "\n")
        #cat("condition", condition, "\n")

        if (condition == TRUE) {break}
        #if not, update alpha
        alpha_x <- alpha_x * stepsizeShrink
      }

      #If all betas are 0, c goes to 0
      #I put it here, out of the above loop, because otherwise it interfere with the optimization
      if (length(c_pos)!=0){
        if (sum(v1[-c_pos])==0) {v1[c_pos]=0}
      }

      #sum the loops of the line search
      totalBacktracks <- totalBacktracks + backtrackCount

      #results and one step forward
      t0 <- t1
      t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
      z0 <- z1
      y0 <- y1
      q0 <- q1
      q1 <- eta*q0 + 1
      c1 <- (eta*q0*c1 + Fx_x1)/q1
      x0 <- x1
      Fx_z1 <- Fx(z1)
      Fx_v1 <- Fx(v1)

      #cat("Fx_z1",Fx_z1, "\n")
      #cat("Fx_v1",Fx_v1, "\n")

      if(Fx_z1 <= Fx_v1){
        x1 <- z1
        Fx_x1 <- Fx_z1
      } else {
        x1 <- v1
        Fx_x1 <- Fx_v1
      }

      if(i==1){
        best_x1 <- x1
        best_Fx_x1 <- Fx_x1
      } else if(round(Fx_x1, 10) <= round(best_Fx_x1, 10)) {
        best_x1 <- x1
        best_Fx_x1 <- Fx_x1
        cat("\n I am updating the best x1: \n")
        print(best_x1)
        cat("\n")
      }
    }

    #test
    zeros <- FALSE
    if (i>zeros_stay_zeros_from_iteration){
      c_proxy <- c_pos
      if (length(c_pos)==0){
        c_proxy <- 100000000000000000000000000
      }
      if(!((sum(x1[-c_proxy])==0) & (i < 4))){ #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros){
      if (length(c_pos)==0){
        if(sum(x1)==0){
          x1 <- x1[1]
        } else {
          x1 <- x1[x1!=0]
        }
        if (length(x0) != length(x1)){exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
        x0 <- x0[names(x0) %in% names(x1)]
        v1 <- v1[names(v1) %in% names(x1)]
        z0 <- z0[names(z0) %in% names(x1)]
        z1 <- z1[names(z1) %in% names(x1)]
        y0 <- y0[names(y0) %in% names(x1)]
      }
    }
    #end-test

    #print the partial results
    if (trace == 2){
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the nonmonotone version. \n")
      if (length(fold)!=0){
        cat("Fold:", fold, "; \n")
      }
      if (length(c_pos)==0){
        visualize_x1 = x1[which(x1!=0)]
      } else {
        visualize_x1 = c(x1[c_pos], x1[which(x1[-c_pos]!=0)]) #we leave c the possibility vary
      }
      cat("iter:", i, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
      print(visualize_x1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i==1){
      Fx_x1_best <- Fx_x1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_x1_best - Fx_x1 < convergence_error){

        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4){
            cat("mnmAPG converged because the convergence error is below the threshold \n")
            break
          }
      } else {
        Fx_x1_best <- Fx_x1
        exit_if_5_times_the_same <- 0
      }
    }

    #if both the two alphas are very small, the new increment will be very poor, so we consider the algorithm converged
    if ((sum(alpha_x, alpha_y) < min_alpha) || ((0 < alpha_x) & (alpha_x < min_alpha))){
      if (trace %in% c(1,2)){
        cat("mnmAPG converged because alpha is below the threshold \n")
      }
      break
    }

    #counters
    i=i+1

    if (i==max_iter){
  	  if (trace %in% c(1,2)){
        cat("mnmAPG converged because the maximum number of iteration has been reached. \n")
  	  }
    }
  }

  #take the best solution
  x1 <- best_x1
  Fx_x1 <- best_Fx_x1
  preparing_the_solution[names(preparing_the_solution) %in% names(x1)] <- x1
  x1 <- preparing_the_solution

  #print the partial results
  if (trace %in% c(1,2)){
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the nonmonotone version. \n")
    if (length(fold)!=0){
      cat("Fold:", fold, "; \n")
    }
    if (length(c_pos)==0){
      visualize_x1 = x1[which(x1!=0)]
    } else {
      visualize_x1 = c(x1[c_pos], x1[which(x1[-c_pos]!=0)]) #we leave c the possibility vary
    }
    cat("iters:", i, "; backtraking iters:", totalBacktracks, ifelse(length(penalty)!=0, paste0("; penalty:", penalty),""), ifelse(length(lambda)!=0,paste0("; lambda:", lambda),""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
    print(visualize_x1)
    cat("\n")
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
  }

  return(list(x1=x1, tot_iters=i, backtrack_iters=totalBacktracks, estimation_time=estimation_time))
}







#The next part of the code is internal and aims to test the differences in term of time and precision of using
#diffent values of stepsizeShrink and min_alpha, to better understand the cost-benefit of the choice if this params
test_parameters <- function(df, X=names(df[,!(names(df) == y)]) , y="y", stepsizeShrink= c(0.5,0.6,0.7,0.8),
                            min_alpha=c(1e-7, 1e-10, 1e-12), lambda=0.5, kernel="gaussian", alpha=0.5,
                            a1=3.7, a2=3, penalty="L1"){

  if (class(X) != "character"){stop("X can only be of class character")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X))])

  betas1_initial_zeros <- rep(0, ncol(df1)-1)
  names <- c(colnames(df1[,3:length(df1)]),"c")

  betas_start <- betas1_initial_zeros
  names(betas_start) <- names

  delta_fx <- function(x){
    result <- -pye_KS(df=df1, X=X, y=y, betas=x[-length(x)], lambda=lambda, c=x[length(x)], kernel=kernel,alpha=alpha,
                      a1=a1, a2=a2, penalty=penalty)$gr_yi
    return(result)
  }

  if (penalty == "SCAD") {a=a1} else {a=a2}
  prox_penalty = get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  proxx <- function(x, eta){
    result <- c(prox_penalty(betas=x[-length(x)], lambda=eta*lambda, alpha=alpha, a=a), x[length(x)])
    return(result)
  }

  Fx <- function(x){
    result <- -getElement(pye_KS(df=df1, X=X, y=y, betas=x[-length(x)], lambda=lambda, c=x[length(x)], kernel=kernel,
                                 alpha=alpha, a1=a1, a2=a2, penalty=penalty),paste0("pye_", penalty))
    return(result)
  }

  # lunch the test
  for (iii in stepsizeShrink){
    for (iiii in min_alpha){
      #monotone
      assign(paste0("mmAPG_test_", iii, "_", iiii), mmAPG(x0=betas_start, c_pos=length(betas_start), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda,
                                                   penalty="L1", stepsizeShrink=iii, min_alpha=iiii, trace=1, seed=1, max_iter=100))
      assign(paste0("pye_KS_value_mm_", iii, "_", iiii), pye_KS(df=df1, X=X, y=y, betas=get(paste0("mmAPG_test_", iii, "_", iiii))$x1[-length(get(paste0("mmAPG_test_", iii, "_", iiii))$x1)], lambda=lambda,
                                                                c=get(paste0("mmAPG_test_", iii, "_", iiii))$x1[length(get(paste0("mmAPG_test_", iii, "_", iiii))$x1)], kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty))

      #nonmonotone
      assign(paste0("mnmAPG_test_", iii, "_", iiii), mnmAPG(x0=betas_start, c_pos=length(betas_start), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda,
                                                          penalty="L1", stepsizeShrink=iii, min_alpha=iiii, trace=1, seed=1, max_iter=100))
      assign(paste0("pye_KS_value_mnm_", iii, "_", iiii), pye_KS(df=df1, X=X, y=y, betas=get(paste0("mnmAPG_test_", iii, "_", iiii))$x1[-length(get(paste0("mnmAPG_test_", iii, "_", iiii))$x1)], lambda=lambda,
                                                             c=get(paste0("mnmAPG_test_", iii, "_", iiii))$x1[length(get(paste0("mnmAPG_test_", iii, "_", iiii))$x1)], kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty))

    }
  }

  result <- lapply(stepsizeShrink, function(x) lapply(min_alpha, function(y) list(get(paste0("mmAPG_test_", x, "_", y)), get(paste0("pye_KS_value_mm_", iii, "_", iiii)),get(paste0("mnmAPG_test_", x, "_", y)), get(paste0("pye_KS_value_mnm_", iii, "_", iiii)))))
  return(result)
}

