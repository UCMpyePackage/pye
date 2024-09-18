#MCP (Minimax Concave Penalty ) penalizzation
single_MCP_function <- function(beta, lambda, a){

  #for every single element of beta
  if (abs(beta) <= a*lambda){
   value = lambda*abs(beta) - (beta^2)/(2*a)
  } else {
   value = (a*lambda^2)/2
  }

  return(value)
}


#' @title MCP_function
#'
#' @description MCP penalization function
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param a the hyperparameter of the penalization MCP
#'
#' @return the value of the penalty MCP function
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' a <- 3.0
#' MCP_function(betas=betas, lambda=lambda, a=a)
#' 
#' @export
MCP_function <- function(betas, lambda, a){

  penalties <- sapply(betas, function (x) single_MCP_function(x, lambda, a))
  penalty <- sum(penalties)

  return(penalty)
}


#MCP_derivative <- function(beta, lambda, a){
#
#  lambda=rep(lambda, length(beta))
#
#  ifelse(abs(beta) <= a*lambda,
#         value <- (lambda - abs(beta)/a)*sign(beta), value <- 0)
#
#  return(value)
#}
