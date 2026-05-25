# function to generate the correlated regressors following different distributions
#' @importFrom Matrix bdiag
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm qnorm
create_variables <- function(rows, cols, cov, mu = rep(0, cols), seed = 1) {

  # Input validation
  if (!is.numeric(rows) || length(rows) != 1 || rows <= 0 || rows != as.integer(rows)) {
    stop("`rows` must be a single positive integer.")
  }
  if (!is.numeric(cols) || length(cols) != 1 || cols <= 0 || cols != as.integer(cols)) {
    stop("`cols` must be a single positive integer.")
  }
  if (!is.numeric(cov) || length(cov) != 1 || cov < -1 || cov > 1) {
    stop("`cov` must be a single numeric value between -1 and 1.")
  }
  if (!is.numeric(mu) || length(mu) != cols) {
    stop(paste0("`mu` must be a numeric vector of length `cols` (", cols, ")."))
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed)) {
    stop("`seed` must be a single integer.")
  }
  if (cols %% 4 != 0) {
    stop("`cols` must be a multiple of 4 to create the block-diagonal covariance matrix.")
  }
  set.seed(seed)

  #block-diagonal variance-covariance matrix
  little_sigma <- matrix(cov, nrow = cols / 4, ncol = cols / 4) + diag(cols / 4) * (1 - cov)
  Sigma <- as.matrix(Matrix::bdiag(little_sigma, little_sigma, little_sigma, little_sigma))
  #diagonal matrix
  # Sigma <- matrix(0, nrow = cols, ncol = cols)
  # Sigma[1:(cols / 4), 1:(cols / 4)] <- matrix(cov, nrow = cols / 4, ncol = cols / 4) + diag(cols / 4) * (1 - cov)
  # Sigma[((cols / 4) + 1):(cols / 2), ((cols / 4) + 1):(cols / 2)] <- matrix(cov, nrow = cols / 4, ncol = cols / 4) + diag(cols / 4) * (1 - cov)
  # Sigma[((cols / 2) + 1):((cols / 4)*3), ((cols / 2) + 1):((cols / 4)*3)] <- matrix(cov, nrow = cols / 4, ncol = cols / 4)+ diag(cols / 4) * (1 - cov)
  # Sigma[(((cols / 4)*3) + 1):cols,(((cols / 4)*3) + 1):cols] <- matrix(cov, nrow = cols / 4, ncol = cols / 4)+ diag(cols / 4) * (1 - cov)

  #Creation
  normal <- MASS::mvrnorm(n = rows, mu = mu, Sigma = Sigma)
  #CDF
  pvars <- stats::pnorm(normal)

  #normal variables
  #normvars <- as.data.frame(stats::qnorm(pvars[, 1:(cols / 4)], mean = 3, sd= 2))
  #colnames(normvars) <- paste("norm", 1:(cols / 4), sep="")
  #keep only the normal variables
  normvars <- as.data.frame(stats::qnorm(pvars[, 1:cols], mean = 0, sd = 1))
  colnames(normvars) <- paste0("norm", 1:(cols))
  #apply normal copula
  #install.packages("copula", repos = "http://R-Forge.R-project.org")
  #involved packages to put in Description file: gsl, pcaPP, pspline, copula
  # myCop <- normalCopula(param=c(rep(cov, cols)), dim = cols, dispstr = "un")
  # myMvd <- mvdc(copula = myCop, margins = c(rep("normal", cols)),
  #               paramMargins = list(list(shape = 2, scale = 1),
  #                                 list(shape1=2, shape2=2),
  #                                list(df=5)) )

  #student-t variables
  #stuvars <- as.data.frame(qt(pvars[, 1:(cols / 4)], df = 1))
  #colnames(stuvars) <- paste("stt", 1:(cols / 4), sep="")
  #chi-squared variables
  #chivars <- as.data.frame(qchisq(pvars[, 1:(cols / 4)], df=7))
  #colnames(chivars) <- paste("chisq", 1:(cols / 4), sep="")
  #gamma variables
  #gammavars <- as.data.frame(qgamma(pvars[, 1:(cols / 4)], shape = 2, rate = 2))
  #colnames(gammavars) <- paste("gamma", 1:(cols / 4), sep="")
  #Bernoulli variables
  #bervars <- as.data.frame(apply(pvars[, ((cols / 4) + 1):(cols / 2)], c(1, 2) , function(x) if (x>0.5) {1} else {0}))
  #colnames(bervars) <- paste("ber", 1:(cols / 4), sep="")
  #exponential variables
  #expvars <- as.data.frame(qexp(pvars[, ((cols / 2) + 1):((cols / 4)*3)], rate = 0.5))
  #colnames(expvars) <- paste("exp", 1:(cols / 4), sep="")
  #poisson variables
  #poisvars <- as.data.frame(qpois(pvars[, (((cols / 4)*3) + 1):cols], 5))
  #colnames(poisvars) <- paste("pois", 1:(cols / 4), sep="")

  #merge
  #df4dist <- cbind(normvars, bervars, expvars, poisvars)
  df4dist <- cbind(normvars)

  return(df4dist)
}

# function to generate the correlated regressors following different distributions
#' @importFrom Matrix bdiag
#' @importFrom MASS mvrnorm
create_covariates <- function(rows, cols_cov, cov, mu_cov = rep(0, cols_cov), seed = 1) {
  set.seed(seed)
  #Var-cov matrix
  #library(Matrix)

  #block-diagonal matrix
  little_sigma <- matrix(cov, nrow = cols_cov / 4, ncol = cols_cov / 4) + diag(cols_cov / 4) * (1 - cov)
  Sigma <- as.matrix(Matrix::bdiag(little_sigma, little_sigma, little_sigma, little_sigma))

  #Creation
  normal <- MASS::mvrnorm(n = rows, mu = mu_cov, Sigma = Sigma)
  #CDF
  pvars <- stats::pnorm(normal)

  #normal variables
  normvars <- as.data.frame(stats::qnorm(pvars[, 1:cols_cov], mean = 0, sd = 1))
  colnames(normvars) <- paste0("norm_cov", 1:(cols_cov))
  #normvars1 <- as.data.frame(stats::qnorm(pvars[, 1:(cols_cov / 4)], mean = 0, sd = 1))
  #colnames(normvars1) <- paste("norm_cov", 1:(cols_cov / 4), sep="")
  #normvars2 <- as.data.frame(stats::qnorm(pvars[, ((cols_cov / 2) + 1):((cols_cov / 4)*3)], mean = 0, sd = 1))
  #colnames(normvars2) <- paste("norm_cov", ((cols_cov / 2) + 1):((cols_cov / 4)*3), sep="")
  #apply normal copula
  #install.packages("copula", repos = "http://R-Forge.R-project.org")
  #involved packages to put in Description file: gsl, pcaPP, pspline, copula
  # myCop <- normalCopula(param=c(rep(cov, cols_cov)), dim = cols_cov, dispstr = "un")
  # myMvd <- mvdc(copula = myCop, margins = c(rep("normal", cols_cov)),
  #               paramMargins = list(list(shape = 2, scale = 1),
  #                                 list(shape1=2, shape2=2),
  #                                list(df=5)) )

  #student-t variables
  #stuvars <- as.data.frame(qt(pvars[, 1:(cols_cov / 4)], df = 1))
  #colnames(stuvars) <- paste("stt", 1:(cols_cov / 4), sep="")
  #chi-squared variables
  #chivars <- as.data.frame(qchisq(pvars[, 1:(cols_cov / 4)], df=7))
  #colnames(chivars) <- paste("chisq", 1:(cols_cov / 4), sep="")
  #gamma variables
  #gammavars <- as.data.frame(qgamma(pvars[, 1:(cols_cov / 4)], shape = 2, rate = 2))
  #colnames(gammavars) <- paste("gamma", 1:(cols_cov / 4), sep="")
  #Bernoulli variables
  #bervars1 <- as.data.frame(apply(pvars[, ((cols_cov / 4) + 1):(cols_cov / 2)], c(1, 2) , function(x) if (x>0.5) {1} else {0}))
  #colnames(bervars1) <- paste("ber_cov", ((cols_cov / 4) + 1):(cols_cov / 2), sep="")
  #bervars2 <- as.data.frame(apply(pvars[, (((cols_cov / 4)*3) + 1):cols_cov], c(1, 2) , function(x) if (x>0.5) {1} else {0}))
  #colnames(bervars2) <- paste("ber_cov", (((cols_cov / 4)*3) + 1):cols_cov, sep="")
  #exponential variables
  #expvars <- as.data.frame(qexp(pvars[, ((cols_cov / 2) + 1):((cols_cov / 4)*3)], rate = 0.5))
  #colnames(expvars) <- paste("exp", 1:(cols_cov / 4), sep="")
  #poisson variables
  #poisvars <- as.data.frame(qpois(pvars[, (((cols_cov / 4)*3) + 1):cols_cov], 5))
  #colnames(poisvars) <- paste("pois", 1:(cols_cov / 4), sep="")

  #merge
  #df4dist <- cbind(normvars, bervars, expvars, poisvars)
  #df4dist <- cbind(normvars1, bervars1, normvars2, bervars2)
  df4dist <- normvars

  return(df4dist)
}

#' @title Generate Correlated Regressors and Covariates
#'
#' @description This function generates synthetic regressors (biomarkers)
#'   and covariates using a joint multivariate normal distribution. The
#'   correlation structure is constructed via a latent factor model,
#'   guaranteeing that the output covariance matrix is positive semi-definite.
#'   All generated variables are subsequently transformed to standard normal
#'   scores (mean 0, variance 1).
#'
#' @param rows \code{integer}. Number of observations (rows) to generate.
#' @param cols_reg \code{integer}. Number of regressor variables. Must be a
#'   multiple of 4 (required by subsequent sampling functions).
#' @param cols_cov \code{integer}. Number of covariate variables. Must be a
#'   multiple of 4.
#' @param max_corr \code{numeric}. Maximum absolute correlation allowed among
#'   all generated variables. Actual correlations are randomly drawn from a
#'   normal distribution centered at 0, with a standard deviation scaled by
#'   \code{max_corr}. Must be strictly between 0 and 1.
#' @param mu_reg \code{numeric} vector. Mean values for regressor variables.
#'   Its length must match \code{cols_reg}. Defaults to a zero vector.
#' @param mu_cov \code{numeric} vector. Mean values for covariate variables.
#'   Its length must match \code{cols_cov}. Defaults to a zero vector.
#' @param n_latent_factors \code{integer}. Number of latent factors used in
#'   the covariance matrix construction. Default is 2. Higher values can lead
#'   to more complex correlation structures.
#' @param seed \code{integer}. Seed for reproducibility. Default is 1.
#'
#' @return A \code{list} containing two data frames:
#'   \item{regressors}{Data frame of standardized normal regressors (mean 0,
#'     variance 1).}
#'   \item{covariates}{Data frame of standardized normal covariates (mean 0,
#'     variance 1).}
#'
#' @details
#' The correlation structure for all \code{cols_reg} + \code{cols_cov}
#' variables is defined using a latent factor model. Variables' loadings
#' onto \code{n_latent_factors} latent components are randomly drawn, then
#' clipped to ensure correlations are within [-\code{max_corr}, \code{max_corr}].
#' The covariance matrix \eqn{\Sigma} is derived from these loadings
#' (\eqn{L \times L^T}), ensuring it is positive semi-definite. Finally, the
#' diagonal elements of \eqn{\Sigma} are explicitly set to 1 to enforce unit
#' variance across all generated variables.
#'
#' The generated multivariate normal variables undergo a rank-preserving
#' transformation to standard normal scores (using \code{pnorm} and
#' \code{qnorm}) to ensure stable distribution properties.
#'
#' @examples
#' # Simulate 100 observations with 8 regressors and 8 covariates
#' sim <- create_data_all(
#'   rows = 100,
#'   cols_reg = 8,
#'   cols_cov = 8,
#'   max_corr = 0.6,
#'   seed = 123
#' )
#'
#' # Inspect the structure of simulated outputs
#' str(sim$regressors)
#' str(sim$covariates)
#'
#' # Visualize the correlation distribution
#' cor_mat <- cor(cbind(sim$regressors, sim$covariates))
#' hist(cor_mat[lower.tri(cor_mat)],
#'      main = "Histogram of Pairwise Correlations",
#'      xlab = "Correlation", col = "steelblue")
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm qnorm rnorm
#' @export
create_data_all <- function(rows, cols_reg, cols_cov,
                            max_corr = 0.5, #covar_reg = 0.3, covar_cov = 0.3, covar_cross = 0.3,
                            mu_reg  = rep(0, cols_reg),
                            mu_cov = rep(0, cols_cov),
                            n_latent_factors = 2,
                            seed = 1) {

  # --- Input Validation ---
  # Check that both rg and cv dimensions are divisible by 4
  if (cols_reg %% 4 != 0 || cols_cov %% 4 != 0) {
    stop("`cols_reg` and `cols_cov` must both be divisible by 4 for 4 blocks each.")
  }
  # Validate rows
  if (!is.numeric(rows) || length(rows) != 1 || rows <= 0 || rows != floor(rows)) {
    stop("`rows` must be a single positive integer.")
  }
  # Validate cols_reg
  if (!is.numeric(cols_reg) || length(cols_reg) != 1 || cols_reg <= 0 || cols_reg != floor(cols_reg)) {
    stop("`cols_reg` must be a single positive integer.")
  }
  # Validate cols_cov
  if (!is.numeric(cols_cov) || length(cols_cov) != 1 || cols_cov <= 0 || cols_cov != floor(cols_cov)) {
    stop("`cols_cov` must be a single positive integer.")
  }
  # Validate divisibility constraint
  if (cols_reg %% 4 != 0 || cols_cov %% 4 != 0) {
    stop("`cols_reg` and `cols_cov` must be divisible by 4 for 4 blocks each.")
  }
  # Validate max_corr
  if (!is.numeric(max_corr) || length(max_corr) != 1 || max_corr <= 0 || max_corr >= 1) {
    stop("`max_corr` must be a numeric value strictly between 0 and 1.")
  }
  # Validate mu_reg
  if (!is.numeric(mu_reg) || length(mu_reg) != cols_reg) {
    stop("`mu_reg` must be a numeric vector of length `cols_reg`.")
  }
  # Validate mu_cov
  if (!is.numeric(mu_cov) || length(mu_cov) != cols_cov) {
    stop("`mu_cov` must be a numeric vector of length `cols_cov`.")
  }
  # Validate seed
  if (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed)) {
    stop("`seed` must be a single integer.")
  }

  set.seed(seed)

  # Create a random covariate matrix using a latent factor model:
  # Total number of variables
  cols <- cols_reg + cols_cov
  # Number of latent factors (can be tuned: e.g., 5, 10)
  k <- n_latent_factors
  # --- Construct Covariance Matrix using Latent Factor Model ---
  # Factor loadings: each row corresponds to a variable, each column to a latent factor.
  # Loadings are drawn from a normal distribution with mean 0 and standard deviation
  # scaled by `max_corr`. This allows for random correlation strengths.
  loadings <- matrix(rnorm(cols * k, mean = 0, sd = max_corr / 3), nrow = cols)
  # Clip values between −max_corr and max_corr
  loadings <- pmax(pmin(loadings, max_corr), -max_corr)
  # Construct Sigma = L × L^T → automatically positive semi-definite
  # (Since it's a product of a matrix and its transpose, Σ is guaranteed to be
  # symmetric and positive semi-definite, no matter what the loadings are.)
  Sigma <- loadings %*% t(loadings)
  # Set diagonal to 1 for unit variance
  diag(Sigma) <- 1
  #Sigma <- cov2cor(loadings %*% t(loadings))

  mu <- c(mu_reg, mu_cov)
  normal_all <- MASS::mvrnorm(n = rows, mu = mu, Sigma = Sigma)

  # check corr
  # max(abs(cor(normal_all)[lower.tri(cor(normal_all))]))

  pvars <- stats::pnorm(normal_all)
	eps <- 1e-12
  pvars <- pmin(pmax(pvars, eps), 1 - eps)

  normvars <- as.data.frame(stats::qnorm(pvars, mean = 0, sd = 1))

  #TEST TO DELETE
  #cor_matrix <- cor(normvars)
  #lower_values <- cor_matrix[lower.tri(cor_matrix)]
  #hist(lower_values)
  #cor_matrix[1998:2003, 1998:2003]

  #cor_matrix <- cor(df_breastCA$df_breastCA_all[,tail(df_breastCA$X, 200)], df_breastCA$df_breastCA_all[df_breastCA$C])
  #lower_values <- cor_matrix[lower.tri(cor_matrix)]
  #hist(lower_values)
  #cor_matrix

  #cor_matrix <- cor(df_breastCA$df_breastCA_all[df_breastCA$X])
  #lower_values <- cor_matrix[lower.tri(cor_matrix)]
  #hist(lower_values)
  #cor_matrix

  #cor_matrix <- cor(simMicroarrayData_cov03_dim200_covariates20$df)
  #lower_values <- cor_matrix[lower.tri(cor_matrix)]
  #hist(lower_values)
  #cor_matrix[1998:2003, 1998:2003]

  #cor_matrix <- cor(simMicroarrayData_cov03_dim200_covariates2000$df)
  #lower_values <- cor_matrix[lower.tri(cor_matrix)]
  #hist(lower_values)
  #cor_matrix[1998:2003, 1998:2003]

  rg <- normvars[, 1:cols_reg]
  cv <- normvars[, (cols_reg + 1):(cols_reg + cols_cov)]

  colnames(rg) <- paste0("norm", 1:cols_reg)
  colnames(cv) <- paste0("norm_cov", 1:cols_cov)

  #df <- cbind(rg, cv)

  return(list(regressors = rg, covariates = cv))
}



#' @title Create Synthetic High-Dimensional Sample with Binary Target
#'
#' @description Creates a synthetic dataset with a binary target variable
#'   (\code{y}) and a set of correlated normal regressors (\code{X}). This function
#'   is useful for testing models in high-dimensional settings, particularly
#'   those that do not incorporate covariates.
#'
#' @param rows_train \code{numeric}. Number of rows for the training sample.
#'   Default is 50.
#' @param cols \code{numeric}. Total number of regressor variables in both
#'   training and test samples. Must be a multiple of 4. Default is 2000.
#' @param cov \code{numeric}. The correlation coefficient used to define the
#'   internal correlation within each of the four block-diagonal structures
#'   for generating regressors. Higher values result in more correlated
#'   features. Must be between -1 and 1. Default is 0.5.
#' @param mu \code{numeric} vector. Mean vector for the multivariate normal
#'   distribution used to generate the regressors. Its length must match
#'   \code{cols}. Default is a vector of zeros of length \code{cols}.
#' @param rows_test \code{numeric}. Number of rows for the test sample.
#'   Default is 1000. It is recommended to create a test sample larger than
#'   the training sample.
#' @param seed \code{numeric}. An integer to fix the random seed for
#'   reproducibility. Default is 1.
#' @param varsN \code{numeric} vector or \code{NULL}. A vector of two distinct
#'   integers specifying the column indices of the "real regressors" (used in
#'   the latent variable \code{z}) from the first \code{cols}/4 block.
#'   If \code{NULL}, two are randomly selected. Indices are relative to the
#'   full set of \code{cols} variables.
#' @param varsB \code{numeric} vector or \code{NULL}. Indices of "real
#'   regressors" from the second \code{cols}/4 block.
#' @param varsE \code{numeric} vector or \code{NULL}. Indices of "real
#'   regressors" from the third \code{cols}/4 block.
#' @param varsP \code{numeric} vector or \code{NULL}. Indices of "real
#'   regressors" from the fourth \code{cols}/4 block.
#'
#' @return A \code{list} containing the following elements:
#'   \item{df}{A data frame with the generated target variable (\code{y}) and all
#'     regressors, before scaling.}
#'   \item{df_scaled}{A list containing the full dataset, with regressors
#'     scaled (mean 0, variance 1) and the target \code{y} unscaled.}
#'   \item{train_df_scaled}{The data frame for the training set, with scaled
#'     regressors and unscaled target.}
#'   \item{test_df_scaled}{The data frame for the test set, with scaled
#'     regressors and unscaled target.}
#'   \item{nregressors}{A numeric vector indicating the \code{cols} column
#'     indices of the 8 selected "real regressors".}
#'   \item{regressors}{A character vector of the names of the 8 selected
#'     "real regressors".}
#'   \item{coefficients}{A numeric vector of the fixed coefficients used for the
#'     "real regressors" in the linear combination to generate \code{z} (fixed as
#'     c(-16, -4, 12, 8, -8, -12, 4, 16)).}
#'   \item{linearformula}{A character string representing the linear formula
#'     used to generate the latent variable \code{z}.}
#'   \item{z}{A numeric vector of the latent linear predictor values for all
#'     observations.}
#'   \item{cutoff}{A numeric value, the \eqn{\mathbf{0.7}} quantile of \eqn{z},
#'     used as the threshold to binarize \eqn{y}.}
#'
#' @details
#' The function generates \code{cols} regressors, all following a standard
#' normal distribution, structured into four equal blocks with internal
#' correlation defined by \code{cov} using a block-diagonal covariance matrix.
#'
#' A latent variable \code{z} is created as a linear combination of exactly eight
#' of these regressors (two from each block). The binary target variable \eqn{y}
#' is then generated by comparing \eqn{z} to its 0.7 quantile: \eqn{y=1} if
#' \eqn{z > \text{quantile}(z, 0.7)}
#' and \eqn{y=0} otherwise. Data is finally split into training and test sets.
#'
#' @examples
#' library(pye)
#' # Create a small sample with 20 regressors and specified 'real' variables
#' df <- create_sample(rows_train = 200, cols = 20, rows_test = 20,
#'   varsN = c(2, 4), varsB = c(6, 8), varsE = c(12, 14), varsP = c(16, 18))
#' head(df$df)
#'
#' @importFrom stats quantile
#' @export
create_sample <- function (rows_train = 50, cols = 2000, cov = 0.5, mu = rep(0, cols),
                           rows_test = 1000, seed = 1,
                           varsN = NULL, varsB = NULL, varsE = NULL, varsP = NULL) {

  # Input validation (carried over from internal create_variables and added
  # for other params)
  if (!is.numeric(rows_train) || length(rows_train) != 1 || rows_train <= 0 ||
      rows_train != as.integer(rows_train)) {
    stop("`rows_train` must be a single positive integer.")
  }
  if (!is.numeric(cols) || length(cols) != 1 || cols <= 0 ||
      cols != as.integer(cols)) {
    stop("`cols` must be a single positive integer.")
  }
  if (cols %% 4 != 0) {
    stop("`cols` must be a multiple of 4.")
  }
  if (cols < 8 ) {
    stop("`cols` must be at least 8 to randomly select 2 'real regressors' ", "from each block.")
  }
  if (!is.numeric(cov) || length(cov) != 1 || cov < -1 || cov > 1) {
    stop("`cov` must be a single numeric value between -1 and 1.")
  }
  if (!is.numeric(mu) || length(mu) != cols) {
    stop(paste0("`mu` must be a numeric vector of length `cols` (",
                cols, ")."))
  }
  if (!is.numeric(rows_test) || length(rows_test) != 1 || rows_test <= 0 ||
      rows_test != as.integer(rows_test)) {
    stop("`rows_test` must be a single positive integer.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed)) {
    stop("`seed` must be a single integer.")
  }

  # --- Handle selection of real regressors (varsN, varsB, varsE, varsP) ---
  validate_vars_param <- function(vars_param, block_start_idx, block_end_idx) {
    if (length(vars_param) != 0) {
      if (!is.numeric(vars_param) || length(vars_param) != 2 ||
          any(vars_param <= 0) || any(vars_param != as.integer(vars_param))) {
        stop(paste0("`", substitute(vars_param), "` must be a numeric vector of two positive integers."))
      }
      if (any(vars_param < block_start_idx) ||
          any(vars_param > block_end_idx)) {
        stop(paste0("`", substitute(vars_param), "` values must be within the range [",
                    block_start_idx, ", ", block_end_idx, "]."))
      }
      if (vars_param[1] == vars_param[2]) {
        stop(paste0("`", substitute(vars_param), "` values must be distinct."))
      }
    }
  }

  validate_vars_param(varsN, 1, (cols / 4))
  validate_vars_param(varsB, (cols / 4 + 1), (cols / 2))
  validate_vars_param(varsE, (cols / 2 + 1), (cols / 4 * 3))
  validate_vars_param(varsP, (cols / 4 * 3 + 1), cols)

  #set the seed
  set.seed(seed)
  if (length(varsN) == 0) { varsN <- sort(sample(1:(cols / 4), 2, replace = FALSE))}
  if (length(varsB) == 0) { varsB <- sort(sample((cols / 4 + 1):(cols / 2), 2, replace = FALSE))}
  if (length(varsE) == 0) { varsE <- sort(sample((cols / 2 + 1):(cols / 4 * 3), 2, replace = FALSE))}
  if (length(varsP) == 0) { varsP <- sort(sample((cols / 4 * 3 + 1): cols, 2, replace = FALSE))}
  # With default values the result is: 167, 324, 629, 918, 1299, 1471, 1770, 1966

  rows <- rows_train + rows_test

  df <- create_variables(rows = rows, cols = cols, cov = cov, mu = mu, seed = seed)

  nregressors <- sort(c(varsN, varsB, varsE, varsP))
  regressors <- colnames(df)[nregressors]

  # Coefficients for the real regressors (fixed as per original code)
  betas <- c(-16, -4, 12, 8, -8, -12, 4, 16)

  # Create the latent variable 'z'
  # Directly select columns by index and multiply by betas
  df["z"] <- as.matrix(df[, nregressors]) %*% betas
  cutoff <- stats::quantile(df$z, 0.7)
  df["y"] <- ifelse(df$z > cutoff, 1, 0)
  z <- df$z

  # Reorder columns: target 'y' and delete 'z'
  df <- df[, c("y", setdiff(names(df), c("y", "z")))] #z is removed from he final dataset

  df_scaled <- scaling_df_for_pye(df = df, X = colnames(df[-1]), y = "y")

  # Split data into training and test sets
  split <- sample(rep(1:rows), size = rows_train, replace = FALSE)

  train_df_scaled <- df_scaled$df_scaled[split, ]
  test_df_scaled <- df_scaled$df_scaled[-split, ]

  linearformula <- "z <- as.matrix(df[, nregressors])%*%coefficients"

  return(list(df = df, df_scaled = df_scaled, train_df_scaled = train_df_scaled,
              test_df_scaled = test_df_scaled, nregressors = nregressors, regressors = regressors,
              coefficients = betas, linearformula = linearformula, z = z, cutoff = cutoff))
}


#' @title Create a Synthetic Data Set with Covariates and Binary Target
#'
#' @description Creates a synthetic data set featuring a binary target
#'   variable, a large set of correlated regressors (biomarkers), and
#'   a set of correlated covariates. This synthetic environment is
#'   designed for testing methods that account for covariates in binary
#'   classification, such as the Penalized Youden Index (pye).
#'   Regressors and covariates are generated with controlled correlation
#'   structures based on a latent factor model.
#'
#' @param rows_train \code{integer}. Number of observations for the
#'   training sample. Default is 50.
#' @param cols \code{integer}. Total number of regressor variables
#'   (biomarkers). Must be a multiple of 4. Default is 2000.
#' @param cols_cov \code{integer}. Total number of covariate variables.
#'   Must be a multiple of 4. Default is 20. Increase this for high-
#'   dimensional covariate settings.
#' @param max_rho \code{numeric}. The maximum correlation coefficient used
#'   in the latent factor model for generating variables. Higher values lead
#'   to more correlated features. Must be between 0 and 1. Default is 0.3.
#' @param mu \code{numeric} vector. Mean vector for the multivariate
#'   normal distribution used to generate the regressors. Its length
#'   must match \code{cols}. Default is a vector of zeros.
#' @param mu_cov \code{numeric} vector. Mean vector for the multivariate
#'   normal distribution used to generate the covariates. Its length
#'   must match \code{cols_cov}. Default is a vector of zeros.
#' @param rows_test \code{integer}. Number of observations for the test
#'   sample. Default is 1000.
#' @param seed \code{integer}. An integer to fix the random seed for
#'   reproducibility. Default is 1.
#' @param varsN \code{numeric} vector or \code{NULL}. Indices of the two
#'   "real regressors" chosen from the first \code{cols}/4 block
#'   of variables. If \code{NULL}, two are randomly selected.
#' @param varsB \code{numeric} vector or \code{NULL}. Indices of the two
#'   "real regressors" chosen from the second \code{cols}/4 block
#'   of variables. If \code{NULL}, two are randomly selected.
#' @param varsE \code{numeric} vector or \code{NULL}. Indices of the two
#'   "real regressors" chosen from the third \code{cols}/4 block
#'   of variables. If \code{NULL}, two are randomly selected.
#' @param varsP \code{numeric} vector or \code{NULL}. Indices of the two
#'   "real regressors" chosen from the fourth \code{cols}/4 block
#'   of variables. If \code{NULL}, two are randomly selected.
#' @param varN_cov \code{numeric} or \code{NULL}. Index of the single
#'   "real covariate" chosen from the first \code{cols_cov}/4 block.
#'   If \code{NULL}, one is randomly selected.
#' @param varB_cov \code{numeric} or \code{NULL}. Index of the single
#'   "real covariate" chosen from the second \code{cols_cov}/4 block.
#'   If \code{NULL}, one is randomly selected.
#' @param varE_cov \code{numeric} or \code{NULL}. Index of the single
#'   "real covariate" chosen from the third \code{cols_cov}/4 block.
#'   If \code{NULL}, one is randomly selected.
#' @param varP_cov \code{numeric} or \code{NULL}. Index of the single
#'   "real covariate" chosen from the fourth \code{cols_cov}/4 block.
#'   If \code{NULL}, one is randomly selected.
#' @param n_latent_factors \code{integer}. Number of latent factors used in
#'   the covariance matrix construction. Default is 2. Higher values create
#'   more complex correlation structures.
#'
#' @return A \code{list} containing the following data elements:
#'   \item{df}{A data frame of the full sample (\code{rows} = \code{rows_train} + \code{rows_test})
#'     containing the generated target variable (\code{y}), all standardized
#'     regressors (\code{X}), and all standardized covariates (\code{C}).}
#'   \item{train_df_scaled}{The training subset of \code{df}.}
#'   \item{test_df_scaled}{The test subset of \code{df}.}
#'   \item{nregressors}{\code{numeric} vector of the column indices of the
#'     selected "real regressors" from the full \code{cols} set.}
#'   \item{regressors}{\code{character} vector of the names of the selected
#'     "real regressors."}
#'   \item{coefficients_regressors}{\code{numeric} vector of the coefficients
#'     used for the "real regressors" in the latent variable \code{z}.}
#'   \item{ncovariates}{\code{numeric} vector of the column indices of the
#'     selected "real covariates" from the full \code{cols_cov} set.}
#'   \item{covariates}{\code{character} vector of the names of the selected
#'     "real covariates."}
#'   \item{coefficients_covariates}{\code{numeric} vector of the coefficients
#'     used for the "real covariates" in the dynamic \code{cutoff}.}
#'   \item{target}{The name of the target variable (always "y").}
#'   \item{linearformula_z}{Character string of the formula used to generate
#'     the latent predictor \code{z}.}
#'   \item{linearformula_cutoff}{Character string of the formula used to
#'     generate the dynamic cutoff.}
#'   \item{X}{\code{character} vector of all regressor variable names.}
#'   \item{y}{\code{character} string, the name of the target variable.}
#'   \item{C}{\code{character} vector of all covariate variable names.}
#'   \item{z}{\code{numeric} vector of the latent linear predictor values for
#'     all observations.}
#'   \item{cutoff}{\code{numeric} vector of the dynamic cutoff values for all
#'     observations.}
#'
#' @details
#' The binary target variable \code{y} is generated using a Probit-like
#' model: \eqn{y = 1} if \eqn{z \ge \text{cutoff}}, and \eqn{y = 0} otherwise.
#'
#' The latent linear predictor \eqn{\mathbf{z}} is a linear combination of
#' the few designated "real regressors":
#' \deqn{\mathbf{z} = \mathbf{X}_{real} \mathbf{\beta}_{real}}
#'
#' The dynamic cutoff for classification is a linear combination of
#' the few designated "real covariates":
#' \deqn{text{cutoff} = \mathbf{C}_{real} \mathbf{\gamma}_{real}}
#'
#' All generated regressors and covariates are standardized (mean 0, variance 1)
#' and generated with a controlled correlation structure via a latent factor
#' model. The function includes input validation to enforce the block-based
#' structure (e.g., \code{cols} and \code{cols_cov} must be multiples of 4).
#'
#' @examples
#' library(pye)
#' sample_data <- create_sample_with_covariates(cols = 20, cols_cov = 20,
#'   rows_test = 10,
#'   varsN = c(2, 4), varsB = c(6, 8), varsE = c(12, 14), varsP = c(16, 18),
#'   varN_cov = 2, varB_cov = 6, varE_cov = 15, varP_cov = 19)
#' names(sample_data)
#' head(sample_data$df)
#'
#' @export
create_sample_with_covariates <- function(rows_train = 50, cols = 2000, cols_cov = 20, max_rho = 0.3, #covar_cov = 0.5,
                                          mu = rep(0, cols), mu_cov = rep(0, cols_cov),
                                          rows_test = 1000, seed = 1,
                                          varsN = NULL, varsB = NULL, varsE = NULL, varsP = NULL,
                                          varN_cov = NULL, varB_cov = NULL, varE_cov = NULL, varP_cov = NULL,
                                          n_latent_factors = 2) {

  # Input validation
  if (!is.numeric(rows_train) || length(rows_train) != 1 || rows_train <= 0 ||
      rows_train != as.integer(rows_train)) {
    stop("`rows_train` must be a single positive integer.")
  }
  if (!is.numeric(cols) || length(cols) != 1 || cols <= 0 ||
      cols != as.integer(cols)) {
    stop("`cols` must be a single positive integer.")
  }
  if (cols %% 4 != 0) {
    stop("`cols` must be a multiple of 4.")
  }
  if (!is.numeric(cols_cov) || length(cols_cov) != 1 || cols_cov <= 0 ||
      cols_cov != as.integer(cols_cov)) {
    stop("`cols_cov` must be a single positive integer.")
  }
  if (cols_cov %% 4 != 0) {
    stop("`cols_cov` must be a multiple of 4.")
  }
  if (!is.numeric(max_rho) || length(max_rho) != 1 || max_rho < 0 || max_rho > 1) {
    stop("`max_rho` (maximum variable correlation) must be a single numeric value between 0 and 1.")
  }
  #if (!is.numeric(covar_cov) || length(covar_cov) != 1 || covar_cov < -1 ||
  #    covar_cov > 1) {
  #  stop("`covar_cov` (covariate covariance) must be a single numeric value ",
  #       "between -1 and 1.")
  #}
  if (!is.numeric(mu) || length(mu) != cols) {
    stop(paste0("`mu` must be a numeric vector of length `cols` (",
                cols, ")."))
  }
  if (!is.numeric(mu_cov) || length(mu_cov) != cols_cov) {
    stop(paste0("`mu_cov` must be a numeric vector of length `cols_cov` (",
                cols_cov, ")."))
  }
  if (!is.numeric(rows_test) || length(rows_test) != 1 || rows_test < 0 ||
      rows_test != as.integer(rows_test)) {
    stop("`rows_test` must be a single positive integer or zero.")
  }
  if (!is.numeric(seed) || length(seed) != 1 ||
      seed != as.integer(seed)) {
    stop("`seed` must be a single integer.")
  }


  #set the seed
  set.seed(seed)
  if (is.null(varsN)) {varsN <- c(1, 2)}
  if (is.null(varsB)) {varsB <- c((cols / 4 + 1), (cols / 4 + 2))}
  if (is.null(varsE)) {varsE <- c((cols / 2 + 1), (cols / 2 + 2))}
  if (is.null(varsP)) {varsP <- c((cols / 4 * 3 + 1), (cols / 4 * 3 + 2))}
  #if (length(varsN) == 0) { varsN <- sort(sample(1:(cols / 4), 2, replace = FALSE))} #var v167 (norm167) and v324 (norm324)
  #if (length(varsB) == 0) { varsB <- sort(sample((cols / 4 + 1):(cols / 2), 2, replace = FALSE))} #var v629 (ber129) and v918 (ber418)
  #if (length(varsE) == 0) { varsE <- sort(sample((cols / 2 + 1):(cols / 4 * 3), 2, replace = FALSE))} #var v1299 (exp299) and v1471 (exp471)
  #if (length(varsP) == 0) { varsP <- sort(sample((cols / 4 * 3 + 1): cols, 2, replace = FALSE))} #var v1770 (pois270) and v1966 (pois466)

  if (is.null(varN_cov)) {varN_cov <- 1}
  if (is.null(varB_cov)) {varB_cov <- (cols_cov / 4 + 1)}
  if (is.null(varE_cov)) {varE_cov <- (cols_cov / 2 + 1)}
  if (is.null(varP_cov)) {varP_cov <- (cols_cov / 4 * 3 + 1)}
  #if (length(varN_cov) == 0) { varN_cov <- sort(sample(1:(cols_cov / 4), 1, replace = FALSE))}
  #if (length(varB_cov) == 0) { varB_cov <- sort(sample((cols_cov / 4 + 1):(cols_cov / 2), 1, replace = FALSE))}
  #if (length(varE_cov) == 0) { varE_cov <- sort(sample((cols_cov / 2 + 1):(cols_cov / 4 * 3), 1, replace = FALSE))}
  #if (length(varP_cov) == 0) { varP_cov <- sort(sample((cols_cov / 4 * 3 + 1): cols_cov, 1, replace = FALSE))}


  # --- Handle selection of real regressors (varsN, varsB, varsE, varsP) ---
  validate_vars_param <- function(vars_param, block_start_idx, block_end_idx) {
    if (length(vars_param) != 0) {
      if (!is.numeric(vars_param) || length(vars_param) != 2 ||
          any(vars_param <= 0) || any(vars_param != as.integer(vars_param))) {
        stop(paste0("`", substitute(vars_param), "` must be a numeric vector of two positive integers."))
      }
      if (any(vars_param < block_start_idx) ||
          any(vars_param > block_end_idx)) {
        stop(paste0("`", substitute(vars_param), "` values must be within the range [",
                    block_start_idx, ", ", block_end_idx, "]."))
      }
      if (vars_param[1] == vars_param[2]) {
        stop(paste0("`", substitute(vars_param), "` values must be distinct."))
      }
    }
  }

  validate_vars_param(varsN, 1, (cols / 4))
  validate_vars_param(varsB, (cols / 4 + 1), (cols / 2))
  validate_vars_param(varsE, (cols / 2 + 1), (cols / 4 * 3))
  validate_vars_param(varsP, (cols / 4 * 3 + 1), cols)

  # --- Handle selection of real covariates (varN_cov, varB_cov, varE_cov, varP_cov) ---
  validate_cov_vars <- function(vars_param, block_start_idx, block_end_idx) {
    if (!is.null(vars_param)) {
      if (!is.numeric(vars_param) || length(vars_param) != 1 ||
          vars_param <= 0 || vars_param != as.integer(vars_param)) {
        stop(paste0("`", substitute(vars_param), "` must be a single positive integer."))
      }
      if (vars_param < block_start_idx || vars_param > block_end_idx) {
        stop(paste0("`", substitute(vars_param), "` value must be within the range [",
                    block_start_idx, ", ", block_end_idx, "]."))
      }
    }
  }

  validate_cov_vars(varN_cov, 1, (cols_cov / 4))
  validate_cov_vars(varB_cov, (cols_cov / 4 + 1), (cols_cov / 2))
  validate_cov_vars(varE_cov, (cols_cov / 2 + 1), (cols_cov / 4 * 3))
  validate_cov_vars(varP_cov, (cols_cov / 4 * 3 + 1), cols_cov)

  #increase the sample of other rows_test var (to test the model)
  rows <- rows_train + rows_test

  # Generate covariates
  #cv <- create_covariates(rows = rows, cols_cov = cols_cov, cov = covar_cov, mu_cov = mu_cov, seed = seed)
  # Compute additive signal from covariates
  # Normalize addOn_for_mu to have mean 0 and standard deviation 1.
  #addOn_for_mu <- (rowMeans(cv + sin(cv)) - mean(rowMeans(cv + sin(cv))))/ stats::sd(rowMeans(cv + sin(cv)))
  #addOn_for_mu <- scale(rowMeans(cv + sin(cv)))
  # We add a noise to not correlate addOn_for_mu

  # Regressors
  #rg <- as.matrix(create_variables(rows = rows, cols = cols, cov = max_rho, mu = mu, seed = seed+1))

  # Generate regressors, incorporating the additive effect on the mean
  #df <- 0.6*addOn_for_mu + rg
  # With this approach there was too much collinearity in the generated data

  #Let's use a different approach:
  df_all <- create_data_all(rows = rows, cols_reg = cols, cols_cov = cols_cov,
                            max_corr = max_rho,
                            mu_reg  = mu,
                            mu_cov = mu_cov,
                            seed = seed,
                            n_latent_factors = n_latent_factors)
  df <- df_all$regressors
  cv <- df_all$covariates


  # X and C (all variable names for regressors and covariates)
  X <- names(df)
  C <- names(cv)

  #create the response variable as
  #select 10 random number between 1 and 2000
  #set.seed(seed)

  # Define "real regressors" indices and names
  #old
  nregressors <- rep(0, cols)
  for (n in c(varsN, varsB, varsE, varsP)) {
    nregressors[n] <- n
  }
  regressors <- names(df[, which(nregressors != 0)])
  #new
  #nregressors <- sort(c(varsN, varsB, varsE, varsP))
  #regressors <- colnames(df)[nregressors]

  # Define "real covariates" indices and names
  #old
  ncovariates <- rep(0, cols_cov)
  for (n in c(varN_cov, varB_cov, varE_cov, varP_cov)) {
    ncovariates[n] <- n
  }
  covariates <- names(cv[, which(ncovariates != 0)])
  #new
  #ncovariates <- sort(c(varN_cov, varB_cov, varE_cov, varP_cov))
  #covariates <- colnames(cv)[ncovariates]

  #intercept
  #b0=2

  # Coefficients for real regressors
  betas <- c(1, -2, 3, -4, 4, -3, 2, -1) #chosen ad hoc
  #betas <- c(5, -6, 7, -8, 8, -7, 6, -5) #chosen ad hoc
  coefficients <- rep(0, cols)
  for (i in seq_along(betas)) {
    coefficients[c(which(nregressors != 0))][i] <- betas[i]
  }

  #coefficients covariates
  gammas_cv <- c(1, -2, 2, -1) #chosen ad hoc
  #gammas_cv <- c(3, -4, 4, -3) #chosen ad hoc
  coefficients_cv <- rep(0, cols_cov)
  for (i in seq_along(gammas_cv)) {
    coefficients_cv[c(which(ncovariates != 0))][i] <- gammas_cv[i]
  }

  # Create the latent variable 'z' from real regressors
  df["z"] <- as.matrix(df[, nregressors]) %*% betas
  # Define the dynamic cutoff 'cutoff_cv' from real covariates
  cutoff_cv <- as.matrix(cv[, ncovariates]) %*% gammas_cv
  # Create the binary target variable 'y'
  df['y'] <- ifelse(df$z >= cutoff_cv, 1, 0)
  z <- data.frame(ID = rownames(df), z = df$z)
  target <- "y"
  y <- target

  #check correlation
  #corr.xy <- data.frame(matrix(ncol = cols+1, nrow = 0))#-4 since the const here is included. I want it always be the first var to be computed!
  #names(corr.xy) <- names(df[, 1:(cols+1)])
  #corr.xy[1:3, ] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(lapply(seq(1, cols+1, 1), function (x) stats::cor(df[,'y'], df[, x],  method = h)))))
  #row.names(corr.xy) <- c("pearson", "kendall", "spearman")
  #mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
  #mean.corr.xy <- mean.corr.xy[,order(abs(mean.corr.xy),decreasing = TRUE)]

  #cor(df[,'z'], df[, 1:(ncol(df)-2)])

  #Create an ID
  #df['ID']=as.numeric(rownames(df))
  # Combine regressors (df) and covariates (cv) into a single data frame, reordering 'y'
  df1 <- cbind(df[, c("y", setdiff(names(df), c("y", "z")))], cv) #z is removed from he final dataset
  # I do not need to scale df1 since the generative process created variables from correlated N(0, 1).
  #df_scaled <- scaling_df_for_pye (df = df1, X = colnames(df1[-1]), y = "y")
  df_scaled <- df1

  # Split data into training and test sets
  # Add retry logic to ensure enough observations per class in training set
  split_attempts <- 0
  max_split_attempts <- 100 # Maximum attempts to get a balanced split
  min_obs_per_class <- 10  # Minimum desired observations for each class
  repeat {
    split <- sample(rep(1:rows), size = rows_train, replace = FALSE)
    sum_ones <- sum(df_scaled[split, ][y])
    sum_zeros <- sum(df_scaled[split, ][y] == 0)
    if (sum_ones >= min_obs_per_class && sum_zeros >= min_obs_per_class) {
      break # Condition met, exit loop
    }
    split_attempts <- split_attempts + 1
    if (split_attempts >= max_split_attempts) {
      warning(
        "Could not achieve at least ", min_obs_per_class,
        " observations per class in training set after ",
        max_split_attempts, " attempts. Proceeding with current split, ",
        "consider increasing `rows_train`."
      )
      #break # Give up after too many attempts
    }
  }

  train_df_scaled <- df_scaled[split, ]
  test_df_scaled <- df_scaled[-split, ]

  # Prepare output formulas
  linearformula_z <- "z <- as.matrix(df[, nregressors]) %*% betas"
  linearformula_cutoff <- "as.matrix(cv[, ncovariates]) %*% gammas_cv"

  return(list(df = df1,
              train_df_scaled = train_df_scaled,
              test_df_scaled = test_df_scaled,
              nregressors = nregressors,
              regressors = regressors,
              coefficients_regressors = coefficients, #coefficients_regressors = betas,
              ncovariates = ncovariates,
              covariates = covariates,
              coefficients_covariates = coefficients_cv, #coefficients_covariates = gammas_cv,
              target = target,
              linearformula_z = linearformula_z,
              linearformula_cutoff = linearformula_cutoff,
              X = X, y = y, C = C, z = z,
              cutoff = cutoff_cv))
}



#' @title Load and Preprocess Longitudinal PBC Data for Time-Dependent Analysis
#'
#' @description Loads and preprocesses the Mayo Clinic Primary Biliary
#'   Cholangitis (PBC) Longitudinal Data (`pbcseq` from the `survival`
#'   package), making it suitable for use in time-dependent analysis
#'   (e.g., with the \code{longpye} package). The function performs NA
#'   imputation, factor conversion, data filtering, and a subject-level
#'   train-test split to prevent data leakage.
#'
#' @param T_min \code{integer}. The minimum number of complete visits
#'   a subject must have to be included in the dataset. Default is 4,
#'   as subjects with fewer visits are often excluded for robustness.
#'   Must be a positive integer \eqn{\ge 1}.
#' @param T_max \code{integer}. The latest visit number to consider.
#'   Data points for visits beyond \code{T_max} for any subject will
#'   be excluded. The maximum possible visit number is 16. Default is 6.
#' @param seed \code{integer}. An integer to fix the random seed for
#'   reproducible splitting. Default is 1.
#'
#' @return A \code{list} containing the prepared PBC dataset(s):
#'   \item{df_PBC_final_unscaled}{A \code{data.frame} with all included
#'     subjects and visits (up to \code{T_max}). Variables are
#'     preprocessed (NAs imputed, factors converted) but not scaled.}
#'   \item{df_PBC_scaled_all}{The full preprocessed dataset returned by
#'     \code{scaling_df_for_pye}, with numerical features scaled (mean 0,
#'     variance 1) for all subjects and visits.}
#'   \item{train_df_PBC_scaled}{A \code{data.frame} for the training set
#'     (70\% of subjects), with scaled features.}
#'   \item{test_df_PBC_scaled}{A \code{data.frame} for the test set
#'     (30\% of subjects), with scaled features.}
#'   \item{X}{\code{character} vector of names of the regressor (feature)
#'     variables.}
#'   \item{y}{\code{character} string, the name of the target variable
#'     (always "y").}
#'   \item{t}{\code{character} string, the name of the time variable
#'     (always "count").}
#'   \item{id}{\code{character} string, the name of the subject ID variable
#'     (always "id").}
#'
#' @details
#' This function prepares the \code{survival::pbcseq} dataset for use in
#' longitudinal classification models. The dataset contains repeated
#' measures from the Mayo Clinic PBC trial.
#'
#' \strong{Data Preprocessing Steps:}
#' \itemize{
#'   \item \strong{Exclusion}: Patients who underwent a liver transplant
#'     (\code{status} == 1) are excluded.
#'   \item \strong{Target Creation}: A binary target variable \code{y} is
#'     created: 0 for alive/censored (\code{status} < 2), 1 for dead
#'     (\code{status} = 2).
#'   \item \strong{Time Variable}: A \code{count} variable is added to
#'     represent the visit number (\eqn{1, 2, 3, \ldots}).
#'   \item \strong{Filtering}: Only subjects with \eqn{\ge} \code{T_min} visits
#'     and visits up to \eqn{\le} \code{T_max} are retained.
#'   \item \strong{Categorical Variables}: Factors (\code{sex}, \code{ascites}, \ldots)
#'     are converted into dummy variables using \code{model.matrix}.
#'   \item \strong{Missing Data}: Missing numerical values in \code{chol},
#'     \code{alk.phos}, and \code{platelet} are imputed with the column mean.
#'     A binary flag variable (e.g., \code{cholNA}) is created to indicate
#'     original missingness for each imputed feature.
#'   \item \strong{Scaling}: Numerical features are scaled (mean 0, variance 1).
#'   \item \strong{Split}: A 70/30 train-test split is performed at the
#'     subject (ID) level to ensure no patient's visits appear in both
#'     sets, thereby preventing data leakage.
#' }
#'
#' @examples
#' \donttest{
#' # Requires the 'survival' package to be attached or available
#' pbc_data <- pye:::PBC_Mayo_Clinic_data_for_longpye(T_min = 4, T_max = 6)
#'
#' # Check dimensions and the separation of unique IDs
#' cat("Train visits:", nrow(pbc_data$train_df_PBC_scaled), "\\n")
#' cat("Test visits:", nrow(pbc_data$test_df_PBC_scaled), "\\n")
#' cat("Unique train IDs:", length(unique(pbc_data$train_df_PBC_scaled$id)), "\\n")
#' cat("Unique test IDs:", length(unique(pbc_data$test_df_PBC_scaled$id)), "\\n")
#' }
#'
#' @import survival
#' @importFrom stats aggregate ave
#' @noRd
#' @keywords internal
PBC_Mayo_Clinic_data_for_longpye <- function(T_min = 4, T_max = 6, seed = 1) {

  # Input validation
  T_min_int <- as.integer(T_min)
  if (T_min != T_min_int) {
    warning("`T_min` is not an integer; only the integer part has been ",
            "considered, i.e., T_min = ", T_min_int)
  }
  T_min <- T_min_int

  T_max_int <- as.integer(T_max)
  if (T_max != T_max_int) {
    warning("`T_max` is not an integer; only the integer part has been ",
            "considered, i.e., T_max = ", T_max_int)
  }
  T_max <- T_max_int

  if (T_min < 1) {
    stop("`T_min` must be a positive integer (>= 1).")
  }
  if (T_min > T_max) {
    stop("`T_min` cannot be greater than `T_max`.")
  }

  if (T_max > 16) {
    warning("The dataset contains at most 16 visits for some patients. ",
            "`T_max` has been set to 16.")
    T_max <- 16
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed)) {
    stop("`seed` must be a single integer.")
  }

  # Load the dataset
  #utils::data(pbc, package = "survival")
  df_PBC <- survival::pbcseq #this is the dataset to be used, not pbc
  #df_PBC2 <- pbc #we could use pbc just tu add copper and trig variables at the starting time point
  df_PBC_ordered <- df_PBC[order(df_PBC$id, df_PBC$futime, df_PBC$day), ]
  #unique(df_PBC_ordered$id)
  #19 variables, 312 subjects (IDs) with 1.945 clinical visits recorded


  #DETAILS:
  #Primary sclerosing cholangitis is an autoimmune disease leading to destruction
  #of the small bile ducts in the liver. Progression is slow but inexhortable,
  #eventually leading to cirrhosis and liver decompensation. The condition has been
  #recognized since at least 1851 and was named "primary biliary cirrhosis" in 1949.
  #Because cirrhosis is a feature only of advanced disease, a change of its name
  #to "primary biliary cholangitis" was proposed by patient advocacy groups in 2014.

  #This data is from the Mayo Clinic trial in PBC conducted between 1974 and 1984.
  #A total of 424 PBC patients, referred to Mayo Clinic during that ten-year
  #interval, met eligibility criteria for the randomized placebo controlled trial
  #of the drug D-penicillamine. The first 312 cases in the data set participated
  #in the randomized trial and contain largely complete data. The additional 112
  #cases did not participate in the clinical trial, but consented to have basic
  #measurements recorded and to be followed for survival. Six of those cases were
  #lost to follow-up shortly after diagnosis, so the data here are on an
  #additional 106 cases as well as the 312 randomized participants.

  # --- Variable Preprocessing and NA handling ---

  # 1. Binary variables from factors (sex, ascites, hepato, spiders, edema)
  #    `model.matrix` handles factor conversion and NA-to-dummy.
  # 2. Impute missing numerical variables with mean and add NA flag
  #VARIABLE DESCRIPTION and DATA PREPARATION:
  #1) id: case number
  #2) futime: number of days between registration and the earlier of death, transplantion, or study analysis in July, 1986
  #3) status: status at endpoint - 0=alive, 1=transplanted, 2=dead
  #4) trt: drug - 1 = D-penicillamine, 0 = placebo
  #5) age in years, at registration
  #6) sex: m/f
  #df_PBC_ordered$sex
  sex <- data.frame(sex = ifelse(df_PBC_ordered$sex == "m", 1, 0))
  #7) day: number of days between enrollment and this visit date, remaining values on the line of data refer to this visit date.
  #8) ascites: presence of ascites: 0=no 1=yes and NAs
  #df_PBC_ordered$ascites
  ascites <- data.frame(ascites = factor(df_PBC_ordered$ascites, exclude = NULL))
  ascites <- model.matrix(~. - 1, data = ascites, contrasts.arg = list(ascites = stats::contrasts(ascites$ascites, contrasts = FALSE)))
  #9) hepato: presence of hepatomegaly or enlarged liver 0=no 1=yes
  #df_PBC_ordered$hepato
  hepato <- data.frame(hepato = factor(df_PBC_ordered$hepato, exclude = NULL))
  hepato <- model.matrix(~. - 1, data = hepato, contrasts.arg = list(hepato = stats::contrasts(hepato$hepato, contrasts = FALSE)))
  #10) spiders: presence of spiders 0=no 1=yes (blood vessel malformations in the skin)
  #df_PBC_ordered$spiders
  spiders <- data.frame(spiders = factor(df_PBC_ordered$spiders, exclude = NULL))
  spiders <- model.matrix(~. - 1, data = spiders, contrasts.arg = list(spiders = stats::contrasts(spiders$spiders, contrasts = FALSE)))
  #11) edema: presence of edema 0=no edema and no diuretic therapy for edema; 0.5 = edema present without diuretics (untreated), or edema resolved by diuretics (successfully treated); 1 = edema despite diuretic therapy
  #df_PBC_ordered$edema
  edema <- data.frame(edema = factor(df_PBC_ordered$edema, exclude = NULL))
  edema <- model.matrix(~. - 1, data = edema, contrasts.arg = list(edema = stats::contrasts(edema$edema, contrasts = FALSE)))
  #12) bili: serum bilirubin in mg/dl
  #13) chol: serum cholesterol in mg/dl
  cholNA <- data.frame(cholNA = ifelse(is.na(df_PBC_ordered$chol), 1, 0)) #create a flag variable when chol is NA
  chol <- df_PBC_ordered["chol"]
  chol[is.na(chol)] <- mean(chol[!is.na(chol)]) #input the mean value to chol
  #14) albumin: serum albumin in gm/dl
  #15) alk.phos: alkaline phosphatase in U/liter
  #df_PBC_ordered$alk.phos
  alk.phosNA <- data.frame(alk.phosNA = ifelse(is.na(df_PBC_ordered$alk.phos), 1, 0)) #create a flag variable when alk.phos is NA
  alk.phos <- df_PBC_ordered["alk.phos"]
  alk.phos[is.na(alk.phos)] <- mean(alk.phos[!is.na(alk.phos)]) #input the mean value to alk.phos
  #16) ast: aspartate aminotransferase, once called SGOT in U/ml (serum glutamic-oxaloacetic transaminase, the enzyme name has subsequently changed to "ALT" in the medical literature)
  #17) platelet: platelets per cubic ml / 1000
  #df_PBC_ordered$platelet
  plateletNA <- data.frame(plateletNA = ifelse(is.na(df_PBC_ordered$platelet), 1, 0)) #create a flag variable when platelet is NA
  platelet <- df_PBC_ordered["platelet"]
  platelet[is.na(platelet)] <- mean(platelet[!is.na(platelet)]) #input the mean value to platelet
  #18) protime: prothrombin time in seconds (standardised blood clotting time)
  #19) stage: histologic stage of disease (needs biopsy)
  #df_PBC_ordered$stage
  stage <- data.frame(stage = factor(df_PBC_ordered$stage, exclude = NULL))
  stage <- model.matrix(~. - 1, data = stage, contrasts.arg = list(stage = stats::contrasts(stage$stage, contrasts = FALSE)))

  # Combine all processed variables
  df_PBC_cleaned <- cbind(df_PBC_ordered[, c("id", "futime", "status", "trt", "age")], sex, df_PBC_ordered["day"],
                              ascites, hepato, spiders, edema, df_PBC_ordered["bili"],
                              chol, cholNA, df_PBC_ordered["albumin"], alk.phos, alk.phosNA,
                              df_PBC_ordered["ast"], platelet, plateletNA, df_PBC_ordered["protime"],
                              stage)

  #check
  #sum(!complete.cases(df_PBC_cleaned))
  #no more NAs

  #TARGET VARIABLE "status": status at endpoint, 0/1 / 2 for censored, transplant, dead
  #table(df_PBC_cleaned$status)
  #  0    1    2
  #1073  147  725

  # Filter out transplanted patients (status == 1) (they can bias the outcome)
  df_PBC_ordered2 <- df_PBC_cleaned[!(df_PBC_cleaned$status == 1), ]
  #unique(df_PBC_ordered2$id)
  #283 subjects (IDs) with 1798 clinical visits recorded

  #TARGET VARIABLE "status": status at endpoint, 0/1 / 2 for censored, transplant, dead
  #table(df_PBC_ordered2$status)
  #  0    2
  #1073  725

  # Define the target variable 'y' (0 = alive/censored, 1 = dead)
  df_PBC_ordered2["y"] <- ifelse(df_PBC_ordered2$status < 2, 0, 1)

  #let's consider only the ones with no left censored
  #unique(df_PBC_ordered2[df_PBC_ordered2$day == 0, ][, "id"])
  #all the subjects had the starting visit at t = 0

  # Add visit 'count' for each patient
  df_PBC_ordered2 <- transform(df_PBC_ordered2, count = stats::ave(id, id, FUN = seq_along))

  # --- Filter by T_min and T_max ---

  # 1. Keep only subjects with at least T_min visits
  # Get max count for each subject
  max_visits_per_id <- stats::aggregate(count ~ id, data = df_PBC_ordered2, FUN = max)
  ids_to_keep_Tmin <- max_visits_per_id$id[max_visits_per_id$count >= T_min]
  df_PBC_filtered_by_Tmin <- df_PBC_ordered2[df_PBC_ordered2$id %in% ids_to_keep_Tmin, ]

  #nrow(df_PBC_filtered_by_Tmin)
  #length(unique(df_PBC_filtered_by_Tmin$id))
  #table(df_PBC_filtered_by_Tmin[[y]])
  #distribution of count and y:
  #t1 <- df_PBC_filtered_by_Tmin[order(df_PBC_filtered_by_Tmin$id, -df_PBC_filtered_by_Tmin$count), ]
  #t1 <- t1[!duplicated(t1[, c('id')]), ][, ]
  #cbind(table(subset(t1, select = c('count', 'y'))), tot = table(t1$count))
  #times 14,15 and 16 has zero diseased.

  # 2. Keep only visits up to T_max for the remaining subjects
  df_PBC_final_unscaled <- df_PBC_filtered_by_Tmin[df_PBC_filtered_by_Tmin$count <= T_max, ]

  #nrow(df_PBC_final_unscaled)
  #table(df_PBC_final_unscaled[[y]])

  # Define variable sets for the returned data
  # Define variable sets for the returned data
  id <- "id"
  t <- "count"
  y <- "y"
  # All columns except ID, futime, status, y, and count are features (X)
  X <- colnames(df_PBC_final_unscaled)[!colnames(df_PBC_final_unscaled) %in% c(id, "futime", "status", y, t)]

  # --- Scale features ---
  df_PBC_scaled_all <- scaling_df_for_pye (df = df_PBC_final_unscaled, X = X, y = "y")

  # --- Subject-level Train-Test Split (70/30) ---
  unique_ids <- unique(df_PBC_scaled_all$df_scaled$id)
  num_train_ids <- round(length(unique_ids) * 0.7)
  # Ensure reproducible split
  set.seed(seed)
  train_ids <- sample(unique_ids, num_train_ids, replace = FALSE)

  train_df_PBC_scaled <- df_PBC_scaled_all$df_scaled[df_PBC_scaled_all$df_scaled$id %in% train_ids, ]
  test_df_PBC_scaled <- df_PBC_scaled_all$df_scaled[!(df_PBC_scaled_all$df_scaled$id %in% train_ids), ]

  return(list(df_PBC_final_unscaled = df_PBC_final_unscaled,
              df_PBC_scaled_all = df_PBC_scaled_all,
              train_df_PBC_scaled = train_df_PBC_scaled,
              test_df_PBC_scaled = test_df_PBC_scaled,
              X = X, y = y, t = t, id = id))
}
