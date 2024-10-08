% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pye_proximal_operators.R
\name{proximal_operator_SCAD}
\alias{proximal_operator_SCAD}
\title{proximal_operator_SCAD}
\usage{
proximal_operator_SCAD(betas, lambda, a, ...)
}
\arguments{
\item{betas}{vertor to apply the proximal operator}

\item{lambda}{the penalization parameter}

\item{a}{the hyperparameter of the penalization SCAD}

\item{...}{possibility to insert additional variables if 
necessary, not affecting the result}
}
\description{
Proximal operator related to the penalty SCAD.
The "backward" step of the "forward-backward splitting method",
also named Proximal methods, used to optimize pye.
}
\examples{
library(pye)

betas <- seq(0, 2, by=0.2)
lambda <- 0.5
a <- 3.7
proximal_operator_SCAD(betas=betas, lambda=lambda, a=a)

}
