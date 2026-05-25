# single SCAD (Smoothly Clipped Absolute Deviation) penalization
single_SCAD_function <- function(beta, lambda, a) {

  #for every single element of beta
  if (abs(beta) <= lambda) {
    value <- lambda * (abs(beta))
  } else if ((lambda < abs(beta)) && (abs(beta) <= a * lambda)) {
    value <- (2 * a * lambda * abs(beta) - beta^2 - lambda^2) / (2 * (a - 1))
  } else {
    value <- (lambda^2) * (a + 1) / 2
  }

  return(value)
}

#' @title Smoothly Clipped Absolute Deviation (SCAD) Function Value
#'
#' @description Calculates the value of the Smoothly Clipped Absolute
#'   Deviation (SCAD) penalty function, a popular type of non-convex
#'   regularization used for feature selection and estimation.
#'
#' @param betas \code{numeric} vector. The vector of coefficients (e.g.,
#'   \eqn{\beta}) on which to apply the SCAD penalization.
#' @param lambda \code{numeric}. The penalty parameter (\eqn{\lambda > 0}).
#' @param a \code{numeric}. The hyperparameter of the SCAD penalization
#'   (\eqn{a > 2}). This parameter controls the clipping point and the rate
#'   at which the penalty tapers off. Common values include \eqn{a=3.7}.
#'
#' @return A single \code{numeric} value representing the sum of the
#'   SCAD penalties applied to each element in the \code{betas} vector.
#'
#' @details
#' The SCAD penalty applies hard-thresholding to large coefficients to
#' reduce estimation bias while retaining the sparsity property of Lasso.
#' The first derivative of the penalty function \eqn{p_\lambda(|\beta|)} is:
#' \itemize{
#'   \item \eqn{\lambda \qquad \qquad \qquad \qquad \qquad \text{if } |\beta| \le \lambda}
#'   \item \eqn{\frac{a\lambda - |\beta|}{a - 1} \qquad \qquad \text{if } \lambda < |\beta| \le a\lambda}
#'   \item \eqn{0 \qquad \qquad \qquad \qquad \qquad \text{if } |\beta| > a\lambda}
#' }
#' This function computes the sum of the integral of this derivative (the
#' penalty function value itself) over all elements in \code{betas} using
#' an assumed internal function, \code{single_SCAD_function}.
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by = 0.2)
#' lambda <- 0.5
#' a <- 3.7
#' SCAD_function(betas = betas, lambda = lambda, a = a)
#'
#' @export
SCAD_function <- function(betas, lambda, a) {

  penalties <- sapply(betas, function (x) single_SCAD_function(x, lambda, a))
  penalty <- sum(penalties)

  return(penalty)
}