#L1/2
#function to optimize to find x_bar
x_bar_step <- function (x_bar, lambda, q, beta){
  step <- abs(beta) - lambda*q*x_bar^(q-1)
  return(step)
}

#for a single beta_j
single_L12 <- function(beta_j, lambda, q=0.5){

  x_lambda=(2*lambda*(1-q))^(1/(2-q))
  h_lambda=x_lambda + q*lambda*x_lambda^(q-1)

  if (is.nan(h_lambda)){h_lambda=0}

  #find x_bar
  x_bar=0
  if (abs(beta_j) > h_lambda){
    #starting point
    x_bar_0=(x_lambda + abs(beta_j))/2
    x_bar_before = x_bar_step(x_bar_0, lambda=lambda, q=q, beta=beta_j)
    x_bar = x_bar_step(x_bar_before, lambda=lambda, q=q, beta=beta_j)

    #optimization
    while (x_bar_before != x_bar){
      x_bar_before = x_bar
      x_bar = x_bar_step(x_bar_before, lambda=lambda, q=q, beta=beta_j)
    }
  }

  #compute the proximal operator
  if(abs(beta_j) < h_lambda){prox_L12 <- 0
  } else if(abs(beta_j) == h_lambda) {prox_L12 <- 0 #Prox_L12=sign(beta_j)*x_lambda
  } else if(abs(beta_j) > h_lambda) {prox_L12 <- sign(beta_j)*x_bar}
  return(prox_L12)
}

#Proximal Operators for each of the selected penalties
#They need to identify the b(n+1) given b(n) in the estimation
#of the parameter of pye in its smoothed version


#' @title proximal_operator_L12
#'
#' @description Proximal operator related to the penalty L1/2.
#' The "backward" step of the "forward-backward splitting method",
#' also named Proximal methods, used to optimize pye.
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param ... possibility to insert additional variables if 
#' necessary, not affecting the result
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' proximal_operator_L12(betas=betas, lambda=lambda)
#' 
#' @export
proximal_operator_L12 <- function(betas, lambda,...){

  output=sapply(betas, function (x) (single_L12(beta_j=x, lambda=lambda)))

  return(output)
}



#L1
#for a single beta_j
single_L1 <- function(beta_j, lambda){

  prox_L1 = sign(beta_j)*max(abs(beta_j)- lambda, 0)

  return(prox_L1)
}

#' @title proximal_operator_L1
#'
#' @description Proximal operator related to the penalty L1.
#' The "backward" step of the "forward-backward splitting method",
#' also named Proximal methods, used to optimize pye.
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param ... possibility to insert additional variables if 
#' necessary, not affecting the result
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' proximal_operator_L1(betas=betas, lambda=lambda)
#' 
#' @export
#for all the betas
proximal_operator_L1 <- function(betas, lambda,...){

    output=sapply(betas, function (x) (single_L1(beta_j=x, lambda=lambda)))

    return(output)
}



#Elastic-Net
#for a single beta_j
single_EN <- function(beta_j, lambda, alpha){

  prox_EN = (1/(1+(lambda*(1-alpha))))*sign(beta_j)*max(abs(beta_j)-lambda*alpha, 0)

  return(prox_EN)
}

#' @title proximal_operator_EN
#'
#' @description Proximal operator related to the penalty Elastic-Net.
#' The "backward" step of the "forward-backward splitting method",
#' also named Proximal methods, used to optimize pye.
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param alpha the perameter of the Elastic-Net penalty to weight
#' the penalization L1 with respect of L2
#' @param ... possibility to insert additional variables if 
#' necessary, not affecting the result
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' alpha <- 0.5
#' proximal_operator_EN(betas=betas, lambda=lambda, alpha=alpha)
#' 
#' @export
#for all the betas
proximal_operator_EN <- function(betas, lambda, alpha,...){

  output=sapply(betas, function (x) (single_EN(beta_j=x, lambda=lambda, alpha=alpha)))

  return(output)
}



#SCAD
#for a single beta_j
single_SCAD <- function(beta_j, lambda, a){
  if(abs(beta_j) <= 2*lambda){prox_SCAD = sign(beta_j)*max(abs(beta_j)-lambda, 0)
  } else if ((2*lambda < abs(beta_j)) & (abs(beta_j) <= a*lambda)){prox_SCAD = ((a-1)*beta_j-sign(beta_j)*a*lambda)/(a-2)
  } else if (abs(beta_j) > a*lambda){prox_SCAD = beta_j}

  return(prox_SCAD)
}

#' @title proximal_operator_SCAD
#'
#' @description Proximal operator related to the penalty SCAD.
#' The "backward" step of the "forward-backward splitting method",
#' also named Proximal methods, used to optimize pye.
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param a the hyperparameter of the penalization SCAD
#' @param ... possibility to insert additional variables if 
#' necessary, not affecting the result
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' a <- 3.7
#' proximal_operator_SCAD(betas=betas, lambda=lambda, a=a)
#' 
#' @export
#for all the betas
proximal_operator_SCAD <- function(betas, lambda, a, ...){

  output=sapply(betas, function (x) (single_SCAD(beta_j=x, lambda=lambda, a=a)))

  return(output)
}



#MCP
#for a single beta_j
single_MCP <- function(beta_j, lambda, a){
  if(abs(beta_j) <= lambda){prox_MCP = 0}
  else if((lambda < abs(beta_j)) & (abs(beta_j) <= a*lambda)){prox_MCP = (beta_j-lambda*sign(beta_j))/(1-(1/a))}
  else if(abs(beta_j) > a*lambda){prox_MCP = beta_j}

  return(prox_MCP)
}

#' @title proximal_operator_MCP
#'
#' @description Proximal operator related to the penalty MCP.
#' The "backward" step of the "forward-backward splitting method",
#' also named Proximal methods, used to optimize pye.
#'
#' @param betas vertor to apply the proximal operator
#' @param lambda the penalization parameter
#' @param a the hyperparameter of the penalization MCP
#' @param ... possibility to insert additional variables if 
#' necessary, not affecting the result
#'
#' @examples
#' library(pye)
#'
#' betas <- seq(0, 2, by=0.2)
#' lambda <- 0.5
#' a <- 3.0
#' proximal_operator_MCP(betas=betas, lambda=lambda, a=a)
#' 
#' @export
#for all the betas
proximal_operator_MCP <- function(betas, lambda, a,...){

  output=sapply(betas, function (x) (single_MCP(beta_j=x, lambda=lambda, a=a)))

  return(output)
}



