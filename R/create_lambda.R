#' @title Generate a Sequence of Penalty Parameters (\eqn{\lambda})
#'
#' @description Creates a sequence of logarithmically or linearly
#' spaced values for the regularization penalty parameter, \eqn{\lambda},
#' typically used in penalized regression models like Lasso or Elastic Net.
#' The sequence runs from a specified maximum value (\code{lmax}) down to a
#' minimum value (\code{lmin}).
#'
#' @param n \code{integer}. The number of \eqn{\lambda} values to generate in
#'   the sequence. Must be \eqn{\ge 2}.
#' @param lmax \code{numeric}. The maximum (starting) value for \eqn{\lambda}.
#'   This is usually a value near where all coefficients are zero. Must be
#'   $> \code{lmin}$.
#' @param lmin \code{numeric}. The minimum (ending) value for \eqn{\lambda}.
#'   This is usually a small positive number to ensure model sparsity. Must be
#'   \eqn{\ge 0}.
#'
#' @return A \code{numeric} vector of length \code{n} containing the chosen
#'   \eqn{\lambda} values, sorted in descending order (from \code{lmax} to
#'   \code{lmin}). The sequence is logarithmically spaced unless \code{lmax}
#'   and \code{lmin} are equal, in which case it behaves like a linear sequence.
#'
#' @details
#' The function's primary purpose is to generate a sequence on a logarithmic
#' scale, which is standard practice for tuning regularization parameters, as
#' it allows for a finer grid search near zero while still covering large
#' values.
#'
#' Specifically:
#' \itemize{
#'   \item If \eqn{\log(\code{lmin} / \code{lmax}) \ne 0}, the function generates a
#'     geometric progression (logarithmic scale).
#'   \item If \eqn{\log(\code{lmin} / \code{lmax}) = 0} (i.e., \eqn{\code{lmax} = \code{lmin}}),
#'     it generates a sequence using an arithmetic progression (linear scale),
#'     though typically \code{lmax} and \code{lmin} should be different for
#'     cross-validation.
#' }
#' @examples
#' lambdas <- create_lambda(n = 100, lmax = 10, lmin = 1)
#' print(lambdas)
#'
#' @export

create_lambda <- function(n = 100, lmax = 10, lmin = 1) {
  lambda <- double(n)
  lambda[1] <- lmax
  if (log(lmin / lmax) == 0) {
    lstep <- (lmax - lmin) / (n - 1)
    l <- 2
    while (l <= n) {
      lambda[l] <- lambda[l - 1] - lstep
      l <- l + 1
    }
  } else {
    lstep <- log(lmin / lmax) / (n - 1)
    l <- 2
    while (l <= n) {
      lambda[l] <- lambda[l - 1] * exp(lstep)
      l <- l + 1
    }
  }
  return(lambda)
}
