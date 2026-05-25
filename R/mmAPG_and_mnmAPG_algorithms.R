#' @title Monotone Accelerated Proximal Gradient (APG) method
#'
#' @description Implements the Monotone Accelerated Proximal Gradient (APG)
#'   method, inspired by Li and Lin (2015), for the optimization problem
#'   within the Penalized Youden index Estimator (pye) framework.
#'
#'   This variant is tailored for pye with key modifications:
#'   1. Initialization: Adjusted for sparse starting points.
#'   2. Selection Rule: Implements a non-reversible variable selection.
#'      Once a parameter is set to zero (after a warm-up phase), it is
#'      permanently excluded from the active set to enforce sparsity.
#'
#' @param x0 A named numeric vector representing the starting point for the
#'   optimization. It is highly recommended to use the zero vector
#'   to encourage a sparse solution.
#' @param c_pos The index position of the constant term (intercept) in the
#'   \eqn{x0} vector. Use \eqn{NULL} if no constant is included.
#'   Default is \eqn{NULL}.
#' @param delta_fx A function that computes the gradient of the smooth
#'   (loss) component, \eqn{f(x)}, of the objective function.
#' @param proxx A function that computes the proximal operator related
#'   to the non-smooth (penalty) component, \eqn{g(x)}.
#' @param Fx A function that computes the full objective function,
#'   \eqn{F(x) = f(x) + g(x)}.
#' @param lambda The penalization parameter (\eqn{\lambda}) related to \eqn{g(x)}.
#'   Used primarily for tracing and reporting. Default is \eqn{NULL}.
#' @param penalty The type of penalty (e.g., "L1", "SCAD"). Used for tracing
#'   and reporting. Default is \eqn{NULL}.
#' @param fold An optional numeric fold number, typically used when the
#'   function is called within a cross-validation loop. Default is \eqn{NULL}.
#' @param stepsizeShrink The shrinking factor for the step-size \eqn{\alpha}
#'   in the backtracking line search. Must be in $(0, 1)$. A value closer
#'   to 1 increases accuracy but slows convergence. Default is \eqn{0.8}.
#' @param max_alpha The maximum value allowed for the step-size \eqn{\alpha}.
#'   Default is \eqn{10000}.
#' @param min_alpha The minimum value allowed for the step-size \eqn{\alpha}.
#'   If \eqn{\alpha_x + \alpha_y} falls below this threshold, the algorithm
#'   is considered converged. Default is \eqn{1e-10}.
#' @param delta The convergence criterion parameter used in the line-search
#'   condition: \eqn{F(z_k) \le F(y_k) - \delta \|z_k - y_k\|_2^2}.
#'   Default is \eqn{1e-5}.
#' @param trace An integer to control output verbosity:
#'   \eqn{2} = print details for every step;
#'   \eqn{1} = print only final result and convergence message;
#'   \eqn{0} = no output. Default is \eqn{2}.
#' @param seed Numeric seed for reproducibility. Default is \eqn{1}.
#' @param max_iter The maximum number of iterations. Default is \eqn{10000}.
#' @param convergence_error The absolute tolerance used for the 'minimum
#'   change' convergence criterion: convergence is assumed if the objective
#'   function changes by less than this value for 5 consecutive iterations.
#'   Default is \eqn{1e-7}.
#' @param zeros_stay_zeros_from_iteration The iteration number after which
#'   any coefficient that becomes zero is permanently fixed to zero.
#'   This enforces sparsity. Default is \eqn{20}.
#' @param max.print The number of non-zero elements to display in the trace
#'   output. Default is \eqn{10}.
#'
#' @return A \eqn{list} containing the following elements:
#' \item{x1}{The final estimated vector of parameters, including the constant.}
#' \item{tot_iters}{The total number of main algorithm iterations performed.}
#' \item{backtrack_iters}{The total number of backtracking steps executed.}
#' \item{estimation_time}{The total time taken for the estimation (in minutes).}
#' \item{reason_for_exit}{The convergence reason ("max_iter", "min_change",
#'   or "min_alpha").}
#'
#' @references
#' Li, C., & Lin, Z. (2015). Accelerated Proximal Gradient Methods for
#' Nonconvex Programming. *Advances in Neural Information Processing
#' Systems*, *28*.
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' ID <- rownames(df)
#' df1 <- cbind(ID, df[, c(y, X)], drop = FALSE)
#' penalty <- "L12"
#' c_zero_fixed <- TRUE
#' lambda <- 0.1
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx <- function(x) {
#'   if (c_zero_fixed == TRUE) {
#'     x[names(x) == "c"] <- 0
#'   }
#' 
#'   result <- -pye_KS(df = df1[, names(df1) != "ID", drop = FALSE],
#'     X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")],
#'     lambda = lambda, c = x[(names(x) == "c")], kernel = kernel,
#'     alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty)$gr_yi
#'   if (c_zero_fixed == TRUE) {
#'     result[names(result) == "c"] <- 0
#'   }
#'   return(result)
#' }
#'
#' proxx <- function(x, eta) {
#'   if (c_zero_fixed == TRUE) {
#'     x[names(x) == "c"] <- 0
#'   }
#'   result <- c(prox_penalty(betas = x[!(names(x) == "c")], lambda = eta * lambda,
#'     alpha = 0.5, a = a), x[(names(x) == "c")])
#'   return(result)
#' }
#'
#' Fx <- function(x) {
#'   if (c_zero_fixed == TRUE) {
#'     x[(names(x) == "c")] <- 0
#'   }
#'   result <- -getElement(pye_KS(df = df1[, names(df1) != "ID", drop = FALSE],
#'     X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")],
#'     lambda = lambda, c = x[(names(x) == "c")], kernel = kernel,
#'     alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty), paste0("pye_", penalty))
#'   return(result)
#' }
#'
#' #starting point:
#' x0 <- c(rep(0, length(X)), 0) #the first zero is the constant term
#' names(x0) <- c(X, "c")
#' position_c <- which(names(x0) == "c")
#'
#' estim <- mmAPG(x0 = x0, c_pos = position_c, delta_fx = delta_fx, proxx = proxx, Fx = Fx, 
#'      lambda = lambda, penalty = penalty, max_iter = 8, trace = 2)
#' print(estim)
#'
#' @export

#Monotone APG with line search
mmAPG <- function(x0, c_pos = NULL, delta_fx, proxx,
                  Fx, lambda = NULL, penalty = NULL,
                  fold = NULL, stepsizeShrink = 0.8,
                  max_alpha = 10000, min_alpha = 1e-10,
                  delta = 1e-5, trace = 2, seed = 1,
                  max_iter = 10000, convergence_error = 1e-7,
                  zeros_stay_zeros_from_iteration = 20, max.print = 10) {

  start_time <- Sys.time()
  
	if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
  # Set max.print temporarily
  old_options <- options(max.print = max.print)
  on.exit(options(old_options))

  if (length(names(x0)) == 0) {stop("x0 needs to have a name vector corresponding to the variable names.")}

  preparing_the_solution <- rep(0, length(x0))
  names(preparing_the_solution) <- names(x0)

  t1 <- 1
  t0 <- 0
  z0 <- x1 <- x0

  #counter
  i <- 1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)) {

    y1 <- x1 + (t0 / t1) * (z0 - x1) + ((t0 - 1) / t1) * (x1 - x0)
    #keep not considering the lastly deleted betas
    zeros <- FALSE
    if (i > zeros_stay_zeros_from_iteration) {
      c_proxy <- c_pos
      if (length(c_pos) == 0) {
        c_proxy <- 100000000000000000000000000
      }
      if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros) {
      if (length(c_pos) == 0) {
        y1[x1 == 0] <- 0
      } else {
        y1[-c_pos][x1[-c_pos] == 0] <- 0 #we leave c the possibility vary
      }
    }

    if (i == 1) { #in the first iteration y0 does not exists
      s1 <- z0
    } else {
      s1 <- z0 - y0
    }
    dfx_z0 <- delta_fx(z0)
    if (i == 1) {
      r1 <- dfx_z0
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_z0)) < 0.00000001) {
        dfx_z0 <- dfx_z0 * (10^(min(round(abs(log10(abs(dfx_z0)) + 1)))))
        }
    } else {
      dfx_y0 <- delta_fx(y0)
      r1 <- dfx_z0 - dfx_y0
    }
    if (sum(s1) == 0) {
      alpha_y <- 1
    } else {
      alpha_y <- min(max_alpha, abs((t(s1) %*% t(t(s1))) / (t(s1) %*% t(t(r1)))), na.rm = TRUE) #abs since sometimes might be negative
    }
    # alternatively: alpha_y = (t(s1) %*% t(t(r1))) / (t(r1) %*% t(t(r1)))
    if (is.nan(alpha_y)) {stop("alpha_y is NaN -> Go and discover why!")}

    if (i == 1) {
      s1 <- -x0
    } else {
      s1 <- v0 - x0
    }
    dfx_x0 <- delta_fx(x0)
    if (i == 1) {
      r1 <- -dfx_x0
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_x0)) < 0.00000001) {
        dfx_x0 <- dfx_x0 * (10^(min(round(abs(log10(abs(dfx_x0)) + 1)))))
      }
    } else {
      dfx_v0 <- delta_fx(v0)
      r1 <- dfx_v0 - dfx_x0
    }
    if (sum(s1) == 0) {
      alpha_x <- 1
    } else {
      #print((t(s1) %*% t(t(s1))) / (t(s1) %*% t(t(r1))))
      alpha_x <- min(max_alpha, abs((t(s1) %*% t(t(s1))) / (t(s1) %*% t(t(r1)))), na.rm = TRUE) #abs since sometimes might be negative
    }
    # alternatively: alpha_x = (t(s1) %*% t(t(r1))) / (t(r1) %*% t(t(r1)))
    if (is.nan(alpha_x)) {stop("alpha_x is NaN -> Go and discover why!")}

    #cat("alpha_y:", alpha_y)
    #cat("; alpha_x:", alpha_x, "\n")

    dfx_y1 <- delta_fx(y1)
    if (i == 1) {
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_y1)) < 0.00000001) {
        dfx_y1 <- dfx_y1 * (10^(min(round(abs(log10(abs(dfx_y1)) + 1)))))
      }
    }
    Fx_y1 <- Fx(y1)
    #line search 1
    while (TRUE) {

      z1 <- proxx(y1 - alpha_y * dfx_y1, alpha_y)
      zeros <- FALSE
      if (i > zeros_stay_zeros_from_iteration) {
        c_proxy <- c_pos
        if (length(c_pos) == 0) {
          c_proxy <- 100000000000000000000000000
        }
        if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
          zeros <- TRUE
        }
      }
      if (zeros) {
        if (length(c_pos) == 0) {
          z1[y1 == 0] <- 0
        } else {
          z1[-c_pos][y1[-c_pos] == 0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_z1 <- Fx(z1)
      barrier <- Fx_y1 - delta * norm(z1 - y1, type = "2")^2
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
    if (length(c_pos) != 0) {
      if (sum(z1[-c_pos]) == 0) {z1[c_pos] <- 0}
    }


    dfx_x1 <- delta_fx(x1)
    if (i == 1) {
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_x1)) < 0.00000001) {
        dfx_x1 <- dfx_x1 * (10^(min(round(abs(log10(abs(dfx_x1)) + 1)))))
      }
    }
    Fx_x1 <- Fx(x1)
    #line search 2
    while (TRUE) {

      v1 <- proxx(x1 - alpha_x * dfx_x1, alpha_x)
      zeros <- FALSE
      if (i > zeros_stay_zeros_from_iteration) {
        c_proxy <- c_pos
        if (length(c_pos) == 0) {
          c_proxy <- 100000000000000000000000000
        }
        if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
          zeros <- TRUE
        }
      }
      if (zeros) {
        if (length(c_pos) == 0) {
          v1[x1 == 0] <- 0
        } else {
          v1[-c_pos][x1[-c_pos] == 0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_v1 <- Fx(v1)
      #Fx_x1 <- Fx(x1)
      barrier <- Fx_x1 - delta * norm(v1 - x1, type = "2")^2
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
    if (length(c_pos) != 0) {
      if (sum(v1[-c_pos]) == 0) {v1[c_pos] <- 0}
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
    if (Fx_z1 <= Fx_v1) {
      x1 <- z1
      Fx_x1 <- Fx_z1
    } else {
      x1 <- v1
      Fx_x1 <- Fx_v1
    }

    #test
    zeros <- FALSE
    if (i > zeros_stay_zeros_from_iteration) {
      c_proxy <- c_pos
      if (length(c_pos) == 0) {
        c_proxy <- 100000000000000000000000000
      }
      if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros) {
      if (length(c_pos) == 0) {
        if (sum(x1) == 0) {
          x1 <- x1[1]
        } else {
          x1 <- x1[x1 != 0]
        }
        if (length(x0) != length(x1)) {exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
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
    if (trace == 2) {
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the monotone version. \n")
      if (!is.null(fold)) {
        cat("Fold:", fold, "; \n")
      }
      if (length(c_pos) == 0) {
        visualize_x1 <- x1[which(x1 != 0)]
      } else {
        visualize_x1 <- c(x1[c_pos], x1[-c_pos][which(x1[-c_pos] != 0)]) #we leave c the possibility vary
      }
      cat("iter:", i, ifelse(length(penalty) != 0, paste0("; penalty:", penalty), ""), ifelse(length(lambda) != 0, paste0("; lambda:", lambda), ""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
      print(visualize_x1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i == 1) {
      Fx_x1_best <- Fx_x1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_x1_best - Fx_x1 < convergence_error) {
        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4) {
          reason_for_exit <- "min_change"
          if (trace %in% c(1, 2)) {
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
    if (sum(alpha_x, alpha_y) < min_alpha) {
      reason_for_exit <- "min_alpha"
      if (trace %in% c(1, 2)) {
        cat("mmAPG converged because alpha is below the threshold min_alpha. \n")
      }
      break
    }

    #counters
    i <- i + 1

    if (i == max_iter) {
      reason_for_exit <- "max_iter"
      if (trace %in% c(1, 2)) {
          cat("mmAPG converged because the maximum number of iteration has been reached. \n")
      }
    }
  }

  preparing_the_solution[names(preparing_the_solution) %in% names(x1)] <- x1
  x1 <- preparing_the_solution

  # print the partial results
  if (trace == 2) {
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the monotone version. \n")
    if (!is.null(fold)) {
      cat("Fold:", fold, "; \n")
    }
    if (length(c_pos) == 0) {
      visualize_x1 <- x1[which(x1 != 0)]
    } else {
      visualize_x1 <- c(x1[c_pos], x1[-c_pos][which(x1[-c_pos] != 0)] ) #we leave c the possibility vary
    }
    cat("iters:", i, "; backtraking iters:", totalBacktracks, ifelse(length(penalty) != 0, paste0("; penalty:", penalty), ""), ifelse(length(lambda) != 0, paste0("; lambda:", lambda), ""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
    print(visualize_x1)
    cat("\n")
  }

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  return(list(x1 = x1,
              tot_iters = i,
              backtrack_iters = totalBacktracks,
              estimation_time = estimation_time,
              reason_for_exit = reason_for_exit))
}


#-------------------------------------------- NON-MONOTONE VERSION -------------------------------------------------



#' @title Non-Monotone Accelerated Proximal Gradient (APG) method
#'
#' @description Implements the Non-Monotone Accelerated Proximal Gradient
#'   (APG) method, inspired by Li and Lin (2015), for the optimization
#'   problem within the Penalized Youden index Estimator (pye) framework.
#'
#'   This variant is tailored for pye with key modifications:
#'   1. Non-Monotone Condition: Uses the non-monotone line search
#'      condition, which generally allows for faster convergence by
#'      accepting temporary increases in the objective function.
#'   2. Non-Reversible Selection: Implements a non-reversible variable
#'      selection rule. Once a parameter is set to zero (after a warm-up
#'      phase), it is permanently excluded from the active set to enforce
#'      sparsity.
#'
#' @param x0 A named numeric vector representing the starting point for the
#'   optimization. It is highly recommended to use the zero vector
#'   to encourage a sparse solution.
#' @param c_pos The index position of the constant term (intercept) in the
#'   \eqn{x0} vector. Use \eqn{NULL} if no constant is included.
#'   Default is \eqn{NULL}.
#' @param delta_fx A function that computes the gradient of the smooth
#'   (loss) component, \eqn{f(x)}, of the objective function.
#' @param proxx A function that computes the proximal operator related
#'   to the non-smooth (penalty) component, \eqn{g(x)}.
#' @param Fx A function that computes the full objective function,
#'   \eqn{F(x) = f(x) + g(x)}.
#' @param lambda The penalization parameter (\eqn{\lambda}) related to \eqn{g(x)}.
#'   Used primarily for tracing and reporting. Default is \eqn{NULL}.
#' @param penalty The type of penalty (e.g., "L1", "SCAD"). Used for tracing
#'   and reporting. Default is \eqn{NULL}.
#' @param fold An optional numeric fold number, typically used when the
#'   function is called within a cross-validation loop. Default is \eqn{NULL}.
#' @param stepsizeShrink The shrinking factor for the step-size \eqn{\alpha}
#'   in the backtracking line search. Must be in $(0, 1)$. A value closer
#'   to 1 increases accuracy but slows convergence. Default is \eqn{0.8}.
#' @param max_alpha The maximum value allowed for the step-size \eqn{\alpha}.
#'   Default is \eqn{10000}.
#' @param min_alpha The minimum value allowed for the step-size \eqn{\alpha}.
#'   If \eqn{\alpha_x + \alpha_y} falls below this threshold, the algorithm
#'   is considered converged. Default is \eqn{1e-10}.
#' @param delta The convergence criterion parameter used in the line-search
#'   condition (see \eqn{mmAPG}). Default is \eqn{1e-5}.
#' @param trace An integer to control output verbosity:
#'   \eqn{2} = print details for every step;
#'   \eqn{1} = print only final result and convergence message;
#'   \eqn{0} = no output. Default is \eqn{2}.
#' @param seed Numeric seed for reproducibility. Default is \eqn{1}.
#' @param max_iter The maximum number of iterations. Default is \eqn{10000}.
#' @param convergence_error The absolute tolerance used for the 'minimum
#'   change' convergence criterion: convergence is assumed if the objective
#'   function changes by less than this value for 5 consecutive iterations.
#'   Default is \eqn{1e-7}.
#' @param eta The non-monotonicity control parameter \eqn{\eta \in [0, 1]}.
#'   It dictates how much the current objective function value can deviate
#'   from the best previous value. Li and Lin (2015) suggest $0.8$.
#'   Default is \eqn{0.8}.
#' @param zeros_stay_zeros_from_iteration The iteration number after which
#'   any coefficient that becomes zero is permanently fixed to zero.
#'   This enforces sparsity. Default is \eqn{5}.
#' @param max.print The number of non-zero elements to display in the trace
#'   output. Default is \eqn{10}.
#'
#' @return A \eqn{list} containing the following elements:
#' \item{x1}{The final estimated vector of parameters, including the constant.}
#' \item{tot_iters}{The total number of main algorithm iterations performed.}
#' \item{backtrack_iters}{The total number of backtracking steps executed.}
#' \item{estimation_time}{The total time taken for the estimation (in minutes).}
#' \item{reason_for_exit}{The convergence reason ("max_iter", "min_change",
#'   or "min_alpha").}
#'
#' @references
#' Li, C., & Lin, Z. (2015). Accelerated Proximal Gradient Methods for
#' Nonconvex Programming. *Advances in Neural Information Processing
#' Systems*, *28*.
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' ID <- rownames(df)
#' df1 <- cbind(ID, df[, c(y, X)], drop = FALSE)
#' penalty <- "L12"
#' c_zero_fixed <- TRUE
#' lambda <- 0.1
#' kernel <- "gaussian"
#' if (penalty == "SCAD") {a=3.7} else {a=3.0}
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' #wrappers
#' delta_fx <- function(x) {
#'  if (c_zero_fixed == TRUE) {
#'    x[names(x) == "c"] <- 0
#'   }
#'   result <- -pye_KS(df = df1[, names(df1) != "ID", drop = FALSE],
#'     X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")],
#'      lambda = lambda, c = x[(names(x) == "c")], kernel = kernel,
#'      alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty)$gr_yi
#'   if (c_zero_fixed == TRUE) {
#'     result[names(result) == "c"] <- 0
#'   }
#'   return(result)
#' }
#'
#' proxx <- function(x, eta) {
#'  if (c_zero_fixed == TRUE) {
#'     x[names(x) == "c"] <- 0
#'   }
#'   result <- c(prox_penalty(betas = x[!(names(x) == "c")], lambda = eta * lambda,
#'     alpha = 0.5, a = a), x[(names(x) == "c")])
#'   return(result)
#' }
#'
#' Fx <- function(x) {
#'   if (c_zero_fixed == TRUE) {
#'     x[(names(x) == "c")] <- 0
#'   }
#'   result <- -getElement(pye_KS(df = df1[, names(df1) != "ID", drop = FALSE],
#'     X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")],
#'     lambda = lambda, c = x[(names(x) == "c")], kernel = kernel,
#'     alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty), paste0("pye_", penalty))
#'   return(result)
#' }
#'
#' #starting point:
#' x0 <- c(rep(0, length(X)), 0) #the first zero is the constant term
#' names(x0) <- c(X, "c")
#' position_c <- which(names(x0) == "c")
#'
#' estim <- mnmAPG(x0 = x0, c_pos = position_c, delta_fx = delta_fx, proxx = proxx, Fx = Fx,
#'      lambda = lambda, penalty = penalty, max_iter = 10, trace = 2)
#' print(estim)
#'
#' @export

#Nonmonotone APG with line search
mnmAPG <- function(x0, c_pos = NULL, delta_fx, proxx, Fx,
                   lambda = NULL, penalty = NULL, fold = NULL,
                   stepsizeShrink = 0.8, eta = 0.8, max_alpha = 10000,
                   min_alpha = 1e-10, delta = 1e-5, trace = 2,
                   seed = 1, max_iter = 10000, convergence_error = 1e-7,
                   zeros_stay_zeros_from_iteration = 20, max.print = 10) {

  start_time <- Sys.time()
  
	if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
  # Set max.print temporarily
  old_options <- options(max.print = max.print)
  on.exit(options(old_options))

  if (length(names(x0)) == 0) {stop("x0 needs to have a name vector corresponding to the variable names.")}

  preparing_the_solution <- rep(0, length(x0))
  names(preparing_the_solution) <- names(x0)

  #initializations
  t1 <- 1
  t0 <- 0
  z0 <- x1 <- v1 <- x0
  c1 <- Fx(x1)
  q1 <- 1
  alpha_x <- 0
	reason_for_exit <- ""

  #counter
  i <- 1

  #set seed
  set.seed(seed)

  #counting the loops of the line search algorithm
  totalBacktracks <- 0
  backtrackCount <- 0

  while ((i < max_iter)) {

    y1 <- x1 + (t0 / t1) * (z0 - x1) + ((t0 - 1) / t1) * (x1 - x0)
    #keep not considering the lastly deleted betas
    zeros <- FALSE
    if (i > zeros_stay_zeros_from_iteration) {
      c_proxy <- c_pos
      if (length(c_pos) == 0) {
        c_proxy <- 100000000000000000000000000
      }
      if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros) {
      if (length(c_pos) == 0) {
        y1[x1 == 0] <- 0
      } else {
        y1[-c_pos][x1[-c_pos] == 0] <- 0 #we leave c the possibility vary
      }
    }

    if (i == 1) {
      s1 <- y1
    } else {
      s1 <- y1 - y0
    }
    dfx_y1 <- delta_fx(y1)
    if (i == 1) {
      r1 <- dfx_y1
      #to try not to stuck to the initial zero point when the derivative is too low
      if (sum(abs(dfx_y1)) < 0.00000001) {
        dfx_y1 <- dfx_y1 * (10^(min(round(abs(log10(abs(dfx_y1)) + 1)))))
      }
    } else {
      dfx_y0 <- delta_fx(y0)
      r1 <- dfx_y1 - dfx_y0
    }
    if (sum(s1) == 0) {
      alpha_y <- 1
    } else {
      alpha_y <- min(max_alpha, abs((t(s1) %*% t(t(s1))) / (t(s1) %*% t(t(r1)))), na.rm = TRUE) #abs since sometimes might be negative
    }
    # alternatively: alpha_y = (t(s1) %*% t(t(r1))) / (t(r1) %*% t(t(r1)))
    if (is.nan(alpha_y)) {stop("alpha_y is NaN -> Go and discover why!")}

    Fx_y1 <- Fx(y1)
    #line search 1
    while (TRUE) {

      z1 <- proxx(y1 - alpha_y * dfx_y1, alpha_y)
      if (i > 3) {
        if (length(c_pos) == 0) {
          z1[y1 == 0] <- 0
        } else {
          z1[-c_pos][y1[-c_pos] == 0] <- 0 #we leave c the possibility vary
        }
      }

      backtrackCount <- backtrackCount + 1

      #condition
      Fx_z1 <- Fx(z1)
      barrier <- Fx_y1 - delta * norm(z1 - y1, type = "2")^2
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
    if (length(c_pos) != 0) {
      if (sum(z1[-c_pos]) == 0) {z1[c_pos] <- 0}
    }

    #second condition
    barrier <- c1 - delta * norm(z1 - y1, type = "2")^2
    condition <- (round(Fx_z1, 10) <= round(barrier, 10))

    #cat("Fx_z1", Fx_z1, "\n")
    #cat("barrier", barrier, "\n")
    #cat("condition", condition, "\n")

    if (condition) {

      x0 <- x1
      x1 <- z1
      y0 <- y1
      Fx_x1 <- Fx_z1

      if (i == 1) {
        best_x1 <- z1
        best_Fx_x1 <- Fx_z1
      } else if (round(Fx_z1, 10) <= round(best_Fx_x1, 10)) {
        best_x1 <- z1
        best_Fx_x1 <- Fx_z1
				if (trace == 2) {
				 cat("\n I am updating the best x1: \n")
         print(best_x1)
         cat("\n")
				}
      }

    } else {

      #set the stepsize alpha_x
      if (i == 1) {
        s1 <- x1
      } else {
        s1 <- x1 - y0
      }
      dfx_x1 <- delta_fx(x1)
      if (i == 1) {
        r1 <- dfx_x1
        #to try not to stuck to the initial zero point when the derivative is too low
        if (sum(abs(dfx_x1)) < 0.00000001) {
          dfx_x1 <- dfx_x1 * (10^(min(round(abs(log10(abs(dfx_x1)) + 1)))))
        }
      } else {
        dfx_y0 <- delta_fx(y0)
        r1 <- dfx_x1 - dfx_y0
      }
      if (sum(s1) == 0) {
        alpha_x <- 1
      } else {
        alpha_x <- min(max_alpha, abs((t(s1) %*% t(t(s1))) / (t(s1) %*% t(t(r1)))), na.rm = TRUE) #abs since sometimes might be negative
      }
      # alternatively: alpha_x = (t(s1) %*% t(t(r1))) / (t(r1) %*% t(t(r1)))
      if (is.nan(alpha_x)) {stop("alpha_x is NaN -> Go and discover why!")}

      #Fx_x1 <- Fx(x1)
      #line search 2
      while (TRUE) {

        v1 <- proxx(x1 - alpha_x * dfx_x1, alpha_x)
        zeros <- FALSE
        if (i > zeros_stay_zeros_from_iteration) {
          c_proxy <- c_pos
          if (length(c_pos) == 0) {
            c_proxy <- 100000000000000000000000000
          }
          if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
            zeros <- TRUE
          }
        }
        if (zeros) {
          if (length(c_pos) == 0) {
            v1[x1 == 0] <- 0
          } else {
            v1[-c_pos][x1[-c_pos] == 0] <- 0 #we leave c the possibility vary
          }
        }

        backtrackCount <- backtrackCount + 1

        #condition
        Fx_v1 <- Fx(v1)
        barrier <- c1 - delta * norm(v1 - x1, type = "2")^2
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
      if (length(c_pos) != 0) {
        if (sum(v1[-c_pos]) == 0) {v1[c_pos] <- 0}
      }

      #sum the loops of the line search
      totalBacktracks <- totalBacktracks + backtrackCount

      #results and one step forward
      t0 <- t1
      t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
      z0 <- z1
      y0 <- y1
      q0 <- q1
      q1 <- eta * q0 + 1
      c1 <- (eta * q0 * c1 + Fx_x1) / q1
      x0 <- x1
      Fx_z1 <- Fx(z1)
      Fx_v1 <- Fx(v1)

      #cat("Fx_z1",Fx_z1, "\n")
      #cat("Fx_v1",Fx_v1, "\n")

      if (Fx_z1 <= Fx_v1) {
        x1 <- z1
        Fx_x1 <- Fx_z1
      } else {
        x1 <- v1
        Fx_x1 <- Fx_v1
      }

      if (i == 1) {
        best_x1 <- x1
        best_Fx_x1 <- Fx_x1
      } else if (round(Fx_x1, 10) <= round(best_Fx_x1, 10)) {
        best_x1 <- x1
        best_Fx_x1 <- Fx_x1
        if (trace == 2) {
				  cat("\n I am updating the best x1: \n")
          print(best_x1)
          cat("\n")
				}
      }
    }

    zeros <- FALSE
    if (i > zeros_stay_zeros_from_iteration) {
      c_proxy <- c_pos
      if (length(c_pos) == 0) {
        c_proxy <- 100000000000000000000000000
      }
      if (!((sum(x1[-c_proxy]) == 0) && (i < 4))) { #if we are at the beginning (i < 4) and with all zeros we allow to vary
        zeros <- TRUE
      }
    }
    if (zeros) {
      if (length(c_pos) == 0) {
        if (sum(x1) == 0) {
          x1 <- x1[1]
        } else {
          x1 <- x1[x1 != 0]
        }
        if (length(x0) != length(x1)) {exit_if_5_times_the_same <- 0} #if the number of selected var changes, we restart the exit_if_5_times_the_same counter
        x0 <- x0[names(x0) %in% names(x1)]
        v1 <- v1[names(v1) %in% names(x1)]
        z0 <- z0[names(z0) %in% names(x1)]
        z1 <- z1[names(z1) %in% names(x1)]
        y0 <- y0[names(y0) %in% names(x1)]
      }
    }
    #end-test

    #print the partial results
    if (trace == 2) {
      cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the non-monotone version. \n")
      if (!is.null(fold)) {
        cat("Fold:", fold, "; \n")
      }
      if (length(c_pos) == 0) {
        visualize_x1 <- x1[which(x1 != 0)]
      } else {
        visualize_x1 <- c(x1[c_pos], x1[which(x1[-c_pos] != 0)]) #we leave c the possibility vary
      }
      cat("iter:", i, ifelse(length(penalty) != 0, paste0("; penalty:", penalty), ""), ifelse(length(lambda) != 0, paste0("; lambda:", lambda), ""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
      print(visualize_x1)
      cat("\n")
    }

    #min_change convergence criteria
    if (i == 1) {
      Fx_x1_best <- Fx_x1 #save the best
      exit_if_5_times_the_same <- 0
    } else {
      if (Fx_x1_best - Fx_x1 < convergence_error) {

        exit_if_5_times_the_same <- exit_if_5_times_the_same + 1
        if (exit_if_5_times_the_same > 4) {
				  reason_for_exit <- "min_change"
					if (trace %in% c(1, 2)) {
            cat("mnmAPG converged because the convergence error is below the threshold \n")
          }
					break
        }
      } else {
        Fx_x1_best <- Fx_x1
        exit_if_5_times_the_same <- 0
      }
    }

    #if both the two alphas are very small, the new increment will be very poor, so we consider the algorithm converged
    if ((sum(alpha_x, alpha_y) < min_alpha) || ((0 < alpha_x) && (alpha_x < min_alpha))) {
		  reason_for_exit <- "min_alpha"
      if (trace %in% c(1, 2)) {
        cat("mnmAPG converged because alpha is below the threshold \n")
      }
      break
    }

    #counters
    i <- i + 1

    if (i == max_iter) {
		  reason_for_exit <- "max_iter"
      if (trace %in% c(1, 2)) {
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
  if (trace == 2) {
    cat("-> Final result: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the non-monotone version. \n")
    if (!is.null(fold)) {
      cat("Fold:", fold, "; \n")
    }
    if (length(c_pos) == 0) {
      visualize_x1 <- x1[which(x1 != 0)]
    } else {
      visualize_x1 <- c(x1[c_pos], x1[which(x1[-c_pos] != 0)]) #we leave c the possibility vary
    }
    cat("iters:", i, "; backtraking iters:", totalBacktracks, ifelse(length(penalty) != 0, paste0("; penalty:", penalty), ""), ifelse(length(lambda) != 0, paste0("; lambda:", lambda), ""), "; alpha_y:", alpha_y, "; alpha_x:", alpha_x, "; F(x1):", Fx_x1, "; \n")
    print(visualize_x1)
    cat("\n")
  }

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  return(list(x1 = x1, 
	            tot_iters = i, 
							backtrack_iters = totalBacktracks, 
							estimation_time = estimation_time,
							reason_for_exit = reason_for_exit))
}







#The next part of the code is internal and aims to test the differences
#in term of time and precision of using
#diffent values of stepsizeShrink and min_alpha, to better understand the
#cost-benefit of the choice if this params
#' @noRd
#' @keywords internal
test_parameters <- function(df, X = NULL, y = "y", stepsizeShrink = c(0.5, 0.6, 0.7, 0.8),
                            min_alpha = c(1e-7, 1e-10, 1e-12), lambda = 0.5, kernel = "gaussian", alpha = 0.5,
                            a1 = 3.7, a2 = 3, penalty = "L1") {

   # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

  # Handle y parameter
  if (inherits(y, "data.frame")) y <- names(y)[1]
  if (!is.character(y) || length(y) != 1) stop("'y' must be a single column name.")
  if (!(y %in% names(df))) stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")

  # Handle X and y parameters
  if (is.null(X)) X <- setdiff(names(df), y)
  if (inherits(X, "data.frame")) X <- names(X)
  if (!is.character(X) || length(X) == 0) stop("'X' must be a character vector of column names or a data.frame.")
  if (!all(X %in% names(df))) stop("Not all regressors in 'X' are found in 'df'.")

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  #df1 <- cbind(ID, df[, (names(df) %in% c(y, X))]) #OLD
  df1 <- cbind(ID, df[, c(y, X)], drop = FALSE)

  betas1_initial_zeros <- rep(0, ncol(df1) - 1)
  names <- c(colnames(df1[, 3:length(df1)]), "c")

  betas_start <- betas1_initial_zeros
  names(betas_start) <- names

  delta_fx <- function(x) {
    result <- -pye_KS(df = df1, X = X, y = y, betas = x[-length(x)], lambda = lambda, c = x[length(x)], kernel = kernel, alpha = alpha,
                      a1 = a1, a2 = a2, penalty = penalty)$gr_yi
    return(result)
  }

  a <- if (penalty == "SCAD") a1 else a2
  prox_penalty <- get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  proxx <- function(x, eta) {
    result <- c(prox_penalty(betas = x[-length(x)], lambda = eta * lambda, alpha = alpha, a = a), x[length(x)])
    return(result)
  }

  Fx <- function(x) {
    result <- -getElement(pye_KS(df = df1, X = X, y = y, betas = x[-length(x)], lambda = lambda, c = x[length(x)], kernel = kernel,
                                 alpha = alpha, a1 = a1, a2 = a2, penalty = penalty), paste0("pye_", penalty))
    return(result)
  }

  # lunch the test
  for (iii in stepsizeShrink) {
    for (iiii in min_alpha) {
      #monotone
      assign(paste0("mmAPG_test_", iii, "_", iiii), mmAPG(x0 = betas_start, c_pos = length(betas_start), delta_fx = delta_fx, proxx = proxx, Fx = Fx, lambda = lambda,
                                                   penalty = "L1", stepsizeShrink = iii, min_alpha = iiii, trace = 1, seed = 1, max_iter = 100))
      assign(paste0("pye_KS_value_mm_", iii, "_", iiii), pye_KS(df = df1, X = X, y = y, betas = get(paste0("mmAPG_test_", iii, "_", iiii))$x1[-length(get(paste0("mmAPG_test_", iii, "_", iiii))$x1)], lambda = lambda,
                                                                c = get(paste0("mmAPG_test_", iii, "_", iiii))$x1[length(get(paste0("mmAPG_test_", iii, "_", iiii))$x1)], kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty))

      #nonmonotone
      assign(paste0("mnmAPG_test_", iii, "_", iiii), mnmAPG(x0 = betas_start, c_pos = length(betas_start), delta_fx = delta_fx, proxx = proxx, Fx = Fx, lambda = lambda,
                                                          penalty = "L1", stepsizeShrink = iii, min_alpha = iiii, trace = 1, seed = 1, max_iter = 100))
      assign(paste0("pye_KS_value_mnm_", iii, "_", iiii), pye_KS(df = df1, X = X, y = y, betas = get(paste0("mnmAPG_test_", iii, "_", iiii))$x1[-length(get(paste0("mnmAPG_test_", iii, "_", iiii))$x1)], lambda = lambda,
                                                             c = get(paste0("mnmAPG_test_", iii, "_", iiii))$x1[length(get(paste0("mnmAPG_test_", iii, "_", iiii))$x1)], kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty))

    }
  }

  result <- lapply(stepsizeShrink, function(x) lapply(min_alpha, function(y) list(get(paste0("mmAPG_test_", x, "_", y)), get(paste0("pye_KS_value_mm_", iii, "_", iiii)), get(paste0("mnmAPG_test_", x, "_", y)), get(paste0("pye_KS_value_mnm_", iii, "_", iiii)))))
  return(result)
}
