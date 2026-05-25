# single MCP (Minimax Concave Penalty) penalizzation
single_MCP_function <- function(beta, lambda, a) {

  #for every single element of beta
  if (abs(beta) <= a * lambda) {
   value <- lambda * abs(beta) - (beta^2) / (2 * a)
  } else {
   value <- (a * lambda^2) / 2
  }

  return(value)
}


#' @title Minimax Concave Penalty (MCP) Function Value
#'
#' @description Calculates the value of the Minimax Concave Penalty (MCP)
#'   function, a regularization penalty often used in sparse regression models.
#'
#' @param betas \code{numeric} vector. The vector of coefficients (e.g.,
#'   \eqn{\beta}) on which to apply the MCP penalization.
#' @param lambda \code{numeric}. The penalty parameter (\eqn{\lambda > 0}).
#' @param a \code{numeric}. The hyperparameter of the MCP penalization
#'   (\eqn{a > 1}). This controls the concavity and asymptote of the penalty.
#'
#' @return A single \code{numeric} value representing the sum of the
#'   MCP penalties applied to each element in the \code{betas} vector.
#'
#' @details
#' The MCP penalty aims to apply an unbiased penalty to large coefficients
#' while still performing continuous shrinkage, similar to SCAD.
#' The penalty for a single coefficient \eqn{|\beta|} is defined as:
#' \itemize{
#'   \item \eqn{\lambda |\beta| - \frac{\beta^2}{2a} \qquad \text{if } |\beta| \le a\lambda}
#'   \item \eqn{\frac{a\lambda^2}{2} \qquad \qquad \text{if } |\beta| > a\lambda}
#' }
#' This function computes the sum of this penalty over all elements in
#' \code{betas} using an assumed internal function, \code{single_MCP_function}.
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by = 0.2)
#' lambda <- 0.5
#' a <- 3.0
#' MCP_function(betas = betas, lambda = lambda, a = a)
#'
#' @export
MCP_function <- function(betas, lambda, a) {

  penalties <- sapply(betas, function (x) single_MCP_function(x, lambda, a))
  penalty <- sum(penalties)

  return(penalty)
}
