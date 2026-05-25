#' @title Proximal Operator for the Non-Convex \eqn{\textrm{L}_{1 / 2}} Penalty
#'
#' @description Computes the proximal operator for the non-convex \eqn{\textrm{L}_{1 / 2}}
#'   penalty. This operator serves as the backward step in the
#'   Forward-Backward Splitting method (or Proximal Gradient method) used to
#'   optimize the Penalized Youden Index (pye) or other objectives with this
#'   penalty. The solution is found using a fixed-point iterative method.
#'
#' @param betas \code{numeric} vector to which the proximal operator should be applied.
#' @param lambda \code{numeric}. The non-negative penalization parameter \eqn{\lambda}.
#' @param q \code{numeric}. The exponent for the \eqn{\textrm{L}_{1 / 2}} penalty, typically
#'   \eqn{q=0.5}. Must be \eqn{0 < q < 1}. Default is \eqn{0.5}.
#' @param max_iter \code{integer}. Maximum number of iterations for the fixed-point
#'   optimization loop for \eqn{x_{\text{bar}}}. Default is \eqn{10000}.
#' @param tol \code{numeric}. Tolerance for the convergence of the fixed-point
#'   optimization loop. Default is \code{1e-10}.
#' @param ... Additional arguments passed to methods
#'
#' @return A \code{numeric} vector of the same length as \code{betas},
#'   containing the result of applying the proximal operator.
#'
#' @details
#' The \eqn{\textrm{L}_{1 / 2}} penalty is non-convex and promotes sparser solutions than
#' the Lasso (\eqn{\textrm{L}_{1}}). The proximal operator for a non-convex penalty
#' often lacks a closed-form solution.
#'
#' This implementation computes the solution by:
#' \enumerate{
#' \item Identifying a threshold \eqn{h_{\lambda}} such that if \eqn{|\beta_j| \le h_{\lambda}},
#'   the solution is \eqn{\text{prox}_{\lambda \textrm{L}_{1 / 2}}(\beta_j) = 0}.
#' \item If \eqn{|\beta_j| > h_{\lambda}}, the non-zero solution 
#'   \eqn{\text{prox}_{\lambda \textrm{L}_{1 / 2}}(\beta_j)} is found using a 
#'   fixed-point iteration method to solve for \eqn{x_{\text{bar}}},
#'   where \eqn{|\beta_j| = x_{\text{bar}} + \lambda q x_{\text{bar}}^{q-1}}.
#' }
#' The final proximal value is \eqn{\text{sign}(\beta_j) \cdot x_{\text{bar}}}.
#'
#' @examples
#' # Example usage:
#' betas <- seq(-2, 2, by = 0.2)
#' lambda <- 0.5
#' q <- 0.5
#'
#' # Apply the proximal operator
#' prox_betas <- proximal_operator_L12(betas = betas, lambda = lambda, q = q)
#' print(prox_betas)
#'
#' # Another example with different lambda
#' lambda2 <- 0.1
#' prox_betas2 <- proximal_operator_L12(betas = betas, lambda = lambda2)
#' print(prox_betas2)
#'
#' @export
proximal_operator_L12 <- function(betas, lambda, q = 0.5, max_iter = 10000, tol = 1e-10, ...) {

  # Internal function to optimize to find x_bar
  # This is a fixed-point iteration for the L1 / 2 proximal operator.
  # It solves for x_bar in the equation |beta| = x_bar + lambda * q * x_bar^(q-1)
  # When x_bar is positive.
  x_bar_step <- function(x_bar, lambda, q, beta) {
    # Ensure x_bar is positive for x_bar^(q-1)
    if (x_bar <= 0) return(0) # or handle appropriately for iterative method
    return(abs(beta) - lambda * q * x_bar^(q - 1))
  }

  # Inner function to compute single beta's proximal operator
  single_L12_internal <- function(beta_j, lambda, q, max_iter, tol) {

    # This threshold h_lambda is specific to the L1/2 penalty
    # and determines if the solution is zero or non-zero.
    x_lambda <- (2 * lambda * (1 - q))^(1 / (2 - q))
    h_lambda <- x_lambda + q * lambda * x_lambda^(q - 1)

    # Handle potential NaN from power of negative or zero if lambda is tiny and q=0.5.
    # For robust calculation, check for non-finite values rather than just NaN.
    if (!is.finite(h_lambda)) { # covers NaN, Inf, -Inf
      h_lambda <- 0 # Assuming threshold is 0 in such pathological cases
    }

    prox_L12 <- 0 # Default to 0 (shrinkage to zero)

    if (abs(beta_j) > h_lambda) {
      # Initialize x_bar_0. A common choice is |beta_j|.
      # (x_lambda + abs(beta_j)) / 2 is also a reasonable heuristic.
      #x_bar_0 <- abs(beta_j) # other option, not used
      x_bar_0 <- (x_lambda + abs(beta_j)) / 2

      x_bar_before <- x_bar_0
      x_bar <- x_bar_step(x_bar = x_bar_before, lambda = lambda, q = q, beta = beta_j)

      # Optimization loop (fixed-point iteration)
      iter <- 0
      while (abs(x_bar - x_bar_before) > tol && iter < max_iter) {
        x_bar_before <- x_bar
        x_bar <- x_bar_step(x_bar = x_bar_before, lambda = lambda, q = q, beta = beta_j)
        iter <- iter + 1

        # Safeguard against non-positive x_bar during iteration for x_bar^(q-1)
        # If x_bar becomes non-positive, it means the iteration has failed or
        # converged to a trivial solution (zero). This can happen if the
        # initial guess or fixed-point function is not suitable.
        if (!is.finite(x_bar) || x_bar <= 0) {
            x_bar <- 0 # Force convergence to 0 if problematic
            break
        }
      }

      # If it did not converge within max_iter, consider it converged to 0
      # For publication, it's safer to have a fallback.
      if (iter >= max_iter && abs(x_bar - x_bar_before) > tol) {
          warning("Fixed-point iteration for x_bar did not converge within ",
                  max_iter, " iterations for beta_j = ", beta_j)
          # A common robust fallback might be the soft-thresholding solution
          # or simply zero if it's very small.
          x_bar <- 0
      }

      prox_L12 <- sign(beta_j) * x_bar
    } else {
      # if abs(beta_j) <= h_lambda, prox_L12 remains 0
    }
    return(prox_L12)
  }

  # Input validation
  if (!is.numeric(betas)) {
    stop("`betas` must be a numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.")
  }
  if (!is.numeric(q) || length(q) != 1 || q <= 0 || q >= 1) {
    stop("`q` must be a single numeric value between 0 (exclusive) and 1 (exclusive).")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter <= 0 ||
      max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a single positive integer.")
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("`tol` must be a single positive numeric value.")
  }

  # Apply the internal function to each element of betas
  # using `vapply` for type safety and efficiency.
  output <- vapply(betas, single_L12_internal, FUN.VALUE = numeric(1),
                   lambda = lambda, q = q,
                   max_iter = max_iter, tol = tol)

  return(output)
}






#' @title Proximal Operator for the Convex \eqn{\textrm{L}_{1}} (Lasso) Penalty
#'
#' @description Computes the proximal operator for the convex \eqn{\textrm{L}_{1}} (Lasso)
#'   penalty. This operator is the soft-thresholding operator and serves
#'   as the backward step in Proximal Gradient optimization.
#'
#' @param betas \code{numeric} vector to which the proximal operator should be applied.
#' @param lambda \code{numeric}. The non-negative penalization parameter \eqn{\lambda}.
#' @param ... Additional arguments passed to methods
#'
#' @return A \code{numeric} vector of the same length as \code{betas},
#'   containing the result of applying the soft-thresholding operator.
#'
#' @details
#' The \eqn{\textrm{L}_{1}} penalty, also known as the Lasso penalty, is a convex
#' penalty commonly used to enforce sparsity (set coefficients to zero).
#' Its proximal operator has a closed-form solution known as the
#' soft-thresholding operator:
#'
#' \deqn{\text{prox}_{\lambda \textrm{L}_1}(\beta) = \text{sign}(\beta) \cdot \max(0, |\beta| - \lambda)}
#'
#' This function applies this operator element-wise to the input \code{betas} vector.
#'
#' @examples
#' # Example usage:
#' betas <- seq(-2, 2, by = 0.2)
#' lambda <- 0.5
#' proximal_operator_L1(betas = betas, lambda = lambda)
#'
#' # Another example with different lambda
#' lambda2 <- 0.1
#' proximal_operator_L1(betas = betas, lambda = lambda2)
#'
#' @export
proximal_operator_L1 <- function(betas, lambda, ...) {

  # Input validation
  if (!is.numeric(betas)) {
    stop("`betas` must be a numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.")
  }

  # Apply the soft-thresholding operator directly (vectorized)
  # sign(betas) gives the sign for each element
  # abs(betas) - lambda subtracts lambda from absolute values
  # pmax(..., 0) ensures that values less than 0 become 0 (thresholding)
  # The product combines the sign with the thresholded magnitude.
  output <- sign(betas) * pmax(abs(betas) - lambda, 0)

  return(output)
}





#' @title Proximal Operator for the Elastic-Net Penalty
#'
#' @description Computes the proximal operator for the Elastic-Net penalty,
#'   which is a convex combination of \eqn{\textrm{L}_{1}} (Lasso) and
#'   \eqn{\textrm{L}_2} (Ridge). This operator serves as the backward step in
#'   Proximal Gradient optimization.
#'
#' @param betas \code{numeric} vector to which the proximal operator should
#'   be applied.
#' @param lambda \code{numeric}. The overall non-negative penalization parameter
#'   \eqn{\lambda}.
#' @param alpha \code{numeric}. The mixing parameter \eqn{\alpha \in [0, 1]}.
#'   \eqn{\alpha = 1} is pure \eqn{\textrm{L}_{1}} (Lasso); 
#'   \eqn{\alpha = 0} is pure \eqn{\textrm{L}_2} (Ridge).
#' @param ... Additional arguments passed to methods
#'
#' @return A \code{numeric} vector of the same length as \code{betas},
#'   containing the result of applying the Elastic-Net proximal operator.
#'
#' @details
#' The Elastic-Net penalty is a convex combination of the \eqn{\textrm{L}_{1}}
#' and \eqn{\textrm{L}_2} penalties. The combined penalty \eqn{P(\beta)} is:
#' \deqn{P(\beta) = \lambda \left( \alpha |\beta| + (1-\alpha) \frac{1}{2} \beta^2 \right)}
#' The proximal operator for the Elastic-Net penalty has a closed-form solution:
#' \deqn{\text{prox}_{P}(\beta) = \frac{1}{1 + \lambda (1-\alpha)} \cdot \text{sign}(\beta) \cdot \max(0, |\beta| - \lambda \alpha)}
#' This function applies this operator element-wise to the input \code{betas} vector.
#'
#' @examples
#' library(pye)
#' betas <- seq(-2, 2, by = 0.2)
#' lambda <- 0.5
#' alpha <- 0.5
#' proximal_operator_EN(betas = betas, lambda = lambda, alpha = alpha)
#'
#' # Equivalent to Lasso
#' prox_lasso <- proximal_operator_EN(betas = betas, lambda = 0.5, alpha = 1)
#'
#' # Equivalent to a scaled Ridge
#' prox_ridge <- proximal_operator_EN(betas = betas, lambda = 0.5, alpha = 0)
#'
#' @export
proximal_operator_EN <- function(betas, lambda, alpha, ...) {

  # Input validation
  if (!is.numeric(betas)) {
    stop("`betas` must be a numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single numeric value between 0 and 1 (inclusive).")
  }

  # Calculate the scaling factor
  scaling_factor <- 1 / (1 + lambda * (1 - alpha))

  # Apply the Elastic-Net proximal operator directly (vectorized)
  # This combines soft-thresholding (for L1 part) with scaling (for L2 part)
  output <- scaling_factor * sign(betas) * pmax(abs(betas) - lambda * alpha, 0)

  return(output)
}





#' @title Proximal Operator for the Non-Convex SCAD Penalty
#'
#' @description Computes the proximal operator for the non-convex Smoothly
#'   Clipped Absolute Deviation (SCAD) penalty. This operator serves as the
#'   backward step in Proximal Gradient optimization.
#'
#' @param betas \code{numeric} vector to which the proximal operator should be applied.
#' @param lambda \code{numeric}. The non-negative penalization parameter \eqn{\lambda}.
#' @param a \code{numeric}. The hyperparameter \eqn{a > 2} for the SCAD penalty.
#'   It controls the concavity of the penalty function. \eqn{a=3.7} is a common choice.
#' @param ... Additional arguments passed to methods
#'
#' @return A \code{numeric} vector of the same length as \code{betas},
#'   containing the result of applying the SCAD proximal operator.
#'
#' @details
#' The SCAD penalty is a non-convex penalty designed for sparse and unbiased
#' estimation. Its proximal operator is defined piecewise based on \eqn{\lambda} and $a$:
#'
#' \enumerate{
#' \item If \eqn{|\beta| \le 2\lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = \text{sign}(\beta) \cdot \max(0, |\beta| - \lambda)}
#' \item If \eqn{2\lambda < |\beta| \le a\lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = \frac{(a - 1)\beta - \text{sign}(\beta)a\lambda}{a-2}}
#' \item If \eqn{|\beta| > a\lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = \beta}
#' }
#'
#' This function applies this operator element-wise to the input \code{betas} vector.
#'
#' @examples
#' # Example usage:
#' betas <- seq(-2, 2, by = 0.2)
#' lambda <- 0.5
#' a <- 3.7 # Common choice for 'a'
#' prox_betas <- proximal_operator_SCAD(betas = betas, lambda = lambda, a = a)
#' print(prox_betas)
#'
#' # Example with a different 'a' value
#' prox_betas2 <- proximal_operator_SCAD(betas = betas, lambda = 0.3, a = 3)
#' print(prox_betas2)
#'
#' @export
proximal_operator_SCAD <- function(betas, lambda, a, ...) {

  # Input validation
  if (!is.numeric(betas)) {
    stop("`betas` must be a numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.")
  }
  if (!is.numeric(a) || length(a) != 1 || a <= 2) {
    stop("`a` must be a single numeric value greater than 2.")
  }

  abs_betas <- abs(betas)
  sgn_betas <- sign(betas)

  # Initialize output vector
  prox_SCAD_output <- numeric(length(betas))

  # Case 1: |betas| <= 2*lambda
  # Equivalent to L1 soft-thresholding
  cond1 <- abs_betas <= 2 * lambda
  prox_SCAD_output[cond1] <- sgn_betas[cond1] * pmax(abs_betas[cond1] - lambda, 0)

  # Case 2: 2*lambda < |betas| <= a * lambda
  cond2 <- (abs_betas > 2 * lambda) & (abs_betas <= a * lambda)
  prox_SCAD_output[cond2] <- ((a - 1) * betas[cond2] - sgn_betas[cond2] * a * lambda) / (a - 2)

  # Case 3: |betas| > a * lambda
  cond3 <- abs_betas > a * lambda
  prox_SCAD_output[cond3] <- betas[cond3]

  # Preserve names
  names(prox_SCAD_output) <- names(betas)

  return(prox_SCAD_output)
}

# Other way to write SCAD for a single beta_j
#single_SCAD <- function(beta_j, lambda, a) {
#  if (abs(beta_j) <= 2*lambda) {prox_SCAD = sign(beta_j)*max(abs(beta_j)-lambda, 0)
#  } else if ((2*lambda < abs(beta_j)) & (abs(beta_j) <= a * lambda)) {prox_SCAD = ((a - 1)*beta_j-sign(beta_j)*a * lambda) / (a-2)
#  } else if (abs(beta_j) > a * lambda) {prox_SCAD = beta_j}
#
#  return(prox_SCAD)
#}








#' @title Proximal Operator for the Non-Convex MCP Penalty
#'
#' @description Computes the proximal operator for the non-convex Minimax
#'   Concave Penalty (MCP). This operator serves as the backward step
#'   in Proximal Gradient optimization.
#'
#' @param betas \code{numeric} vector to which the proximal operator should be applied.
#' @param lambda \code{numeric}. The non-negative penalization parameter \eqn{\lambda}.
#' @param a \code{numeric}. The hyperparameter \eqn{a > 1} for the MCP penalty.
#'   It controls the threshold after which the penalty is zero. $a=3.0$ is a common choice.
#' @param ... Additional arguments passed to methods
#'
#' @return A \code{numeric} vector of the same length as \code{betas},
#'   containing the result of applying the MCP proximal operator.
#'
#' @details
#' The MCP penalty is a non-convex penalty that smoothly transitions from the
#' \eqn{\textrm{L}_{1}} penalty (for small coefficients) to zero (for large coefficients),
#' thereby reducing estimation bias. Its proximal operator is defined piecewise:
#'
#' \enumerate{
#' \item If \eqn{|\beta| \le \lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = 0}
#' \item If \eqn{\lambda < |\beta| \le a\lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = \frac{\beta - \text{sign}(\beta)\lambda}{1 - 1/a}}
#' \item If \eqn{|\beta| > a\lambda}:
#'   \deqn{\text{prox}_{\lambda P}(\beta) = \beta}
#' }
#'
#' This function applies this operator element-wise to the input \code{betas} vector.
#'
#' @examples
#' # Example usage:
#' betas <- seq(-2, 2, by = 0.2)
#' lambda <- 0.5
#' a <- 3.0 # Common choice for 'a'
#' prox_betas <- proximal_operator_MCP(betas = betas, lambda = lambda, a = a)
#' print(prox_betas)
#'
#' # Example with a different 'a' value
#' prox_betas2 <- proximal_operator_MCP(betas = betas, lambda = 0.3, a = 2)
#' print(prox_betas2)
#'
#' @export
proximal_operator_MCP <- function(betas, lambda, a, ...) {

  # Input validation
  if (!is.numeric(betas)) {
    stop("`betas` must be a numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.")
  }
  if (!is.numeric(a) || length(a) != 1 || a <= 1) {
    stop("`a` must be a single numeric value greater than 1.")
  }

  abs_betas <- abs(betas)
  sgn_betas <- sign(betas)

  # Initialize output vector
  prox_MCP_output <- numeric(length(betas))

  # Case 1: |betas| <= lambda
  cond1 <- abs_betas <= lambda
  prox_MCP_output[cond1] <- 0

  # Case 2: lambda < |betas| <= a * lambda
  cond2 <- (abs_betas > lambda) & (abs_betas <= a * lambda)
  prox_MCP_output[cond2] <- (betas[cond2] - lambda * sgn_betas[cond2]) / (1 - (1 / a))

  # Case 3: |betas| > a * lambda
  cond3 <- abs_betas > a * lambda
  prox_MCP_output[cond3] <- betas[cond3]

  # Preserve names
  names(prox_MCP_output) <- names(betas)

  return(prox_MCP_output)
}

# Other way to write MCP for a single beta_j
#single_MCP <- function(beta_j, lambda, a) {
#  if (abs(beta_j) <= lambda) {prox_MCP = 0}
#  else if ((lambda < abs(beta_j)) & (abs(beta_j) <= a * lambda)) {prox_MCP = (beta_j-lambda * sign(beta_j)) / (1-(1/a))}
#  else if (abs(beta_j) > a * lambda) {prox_MCP = beta_j}
#
#  return(prox_MCP)
#}
