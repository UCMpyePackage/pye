#SCAD (Smoothly Clipped Absolute Deviation) penalization
single_SCAD_function <- function(beta, lambda, a){
  
  #for every single element of beta
  if (abs(beta) <= lambda){ 
    value = lambda*(abs(beta))
  } else if ((lambda < abs(beta)) & (abs(beta) <= a*lambda)) {
    value= (2*a*lambda*abs(beta) - beta^2 - lambda^2)/(2*(a-1)) 
  } else {
    value= (lambda^2) * (a + 1)/2
  }

  return(value)
}

#' @title SCAD_function
#'
#' @description SCAD penalization function
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param a the hyperparameter of the penalization SCAD
#'
#' @return the value of the penalty SCAD function
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' a <- 3.7
#' SCAD_function(betas=betas, lambda=lambda, a=a)
#' 
#' @export
SCAD_function <- function(betas, lambda, a){
  
  penalties <- sapply(betas, function (x) single_SCAD_function(x, lambda, a))
  penalty <- sum(penalties)
  
  return(penalty)
}

#SCAD_derivative <- function(beta, lambda, a){
#  
#  lambda=rep(lambda, length(beta))
#  
#  ifelse(abs(beta) <= lambda, value <- lambda*sign(beta), 
#         ifelse((lambda < abs(beta)) & (abs(beta) <= a*lambda),
#                value <- (a*lambda - abs(beta))/(a-1), value <- 0))
#  
#  return(value)
#}
