lam_names <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}

#' @importFrom stats model.matrix
#' @importFrom ncvreg std
#' @useDynLib pye
ncvreg_modified <- function (X, y,
                    family = c("binomial"),
                    penalty = c("MCP", "SCAD"),
                    gamma = switch(penalty, SCAD = 3.7, 3),
					alpha = 1,
                    lambda,
                    eps = 1e-07,
                    max.iter = 10000000,
                    convex = TRUE, dfmax = p + 1,
                    penalty.factor = rep(1, ncol(X)),
                    warn = TRUE, returnX, ...)
{
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- stats::model.matrix(~0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call. = FALSE)
  }
  if (typeof(X) == "integer") storage.mode(X) <- "double"
  if (typeof(X)=="character") stop("X must be a numeric matrix", call.=FALSE)
  if ((!is.null(ncol(y))) & (ncol(y) > 1)) stop("y should be a vector of responses, not a matrix", call.=FALSE)
  if (!is.double(y)) {
    op <- options(warn=2) #stop the execution if there is a warning
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }
  if (!is.double(penalty.factor)) penalty.factor <- as.double(penalty.factor)

 # Error checking
  if (gamma <= 1 & penalty == "MCP") stop("gamma must be greater than 1 for the MC penalty", call. = FALSE)
  if (gamma <= 2 & penalty == "SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call. = FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call. = FALSE)
  if (length(penalty.factor) != ncol(X)) stop("penalty.factor does not match up with X", call. = FALSE)
  if (family == "binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call. = FALSE)
  if (family == "binomial" & !identical(sort(unique(y)), 0:1)) y <- as.double(y == max(y))
  if (length(y) != nrow(X)) stop("X and y do not have the same number of observations", call. = FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call. = FALSE)

  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## Set up XX, yy, lambda
  XX <- ncvreg::std(X)
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  yy <- y
  n <- length(yy)

  nlambda <- length(lambda)
  user.lambda <- TRUE

  if (sys.nframe() > 1) {
    cl <- sys.call(-1)[[1]]
    if (is.name(cl)) {
      if (cl == "local_mfdr")
        return(list(X = XX, y = yy))
    }
  }

  ## Fit
  res <- .Call("cdfit_glm", XX, yy, family, penalty,
               lambda, eps, as.integer(max.iter),
               as.double(gamma), PACKAGE="pye",
               penalty.factor, alpha,
               as.integer(dfmax),
               as.integer(user.lambda | any(penalty.factor == 0)),
               as.integer(warn))

  a <- res[[1]]
  b <- matrix(res[[2]], p, nlambda)
  loss <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)
  iter <- res[[5]]

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  a <- a[ind]
  b <- b[, ind, drop = FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow = (ncol(X) + 1), ncol = length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns + 1, ] <- bb
  beta[1, ] <- a - crossprod(attr(XX, "center")[ns], bb)

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep = "") else colnames(X)
  varnames <- c("(Intercept)", varnames)
  dimnames(beta) <- list(varnames, lam_names(lambda))
  obsnames <- if (is.null(rownames(X))) 1:nrow(X) else rownames(X)

  Eta <- as.matrix(Eta)
  dimnames(Eta) <- list(obsnames, lam_names(lambda))

  ## Output
  val <- structure(list(beta = beta, iter = iter, lambda = lambda,
                        penalty = penalty, family = family, gamma = gamma, alpha = alpha,
                        convex.min = convex.min, loss = loss,
                        linear.predictors = Eta, penalty.factor = penalty.factor,
                        n = n, y = y), class = "ncvreg")

  if (missing(returnX)) {
    if (utils::object.size(XX) > 1e+08) {
      warning("Due to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    }
    else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- XX
  }
  val
}

