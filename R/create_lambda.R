#' @title create_lambda
#'
#' @description Function to create a sequence of n lambdas, given a lambda-max and a lambda-min
#'
#' @param n Number of lambdas
#' @param lmax Maximum lambda value
#' @param lmin Minimum lambda value
#'
#' @return A a vector containing n chosen values of lambdas
#'
#' @examples
#' lambdas <- create_lambda(n=100, lmax=10, lmin=1)
#' @export

create_lambda <- function(n=100, lmax=10, lmin=1){
  lambda = double(n)
  lambda[1] = lmax
  if (log(lmin/lmax) == 0){
    lstep <- (lmax - lmin)/(n-1)
    l=2
    while (l <= n){
      lambda[l] = lambda[l-1]-lstep
      l=l+1
    }
  } else {
    lstep = log(lmin/lmax)/(n - 1)
    l=2
    while (l <= n){
      lambda[l] = lambda[l-1]*exp(lstep)
      l=l+1
    }
  }
  return(lambda)
}
