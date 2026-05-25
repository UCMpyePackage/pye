#' @title Simulation Study for Penalized Models on Real Data
#'
#' @description This function performs a simulation study by repeatedly
#'   splitting a dataset into training and testing sets, estimating a
#'   competitor method (e.g., penalized logistic regression, penalized SVM,
#'   or penalized AUC-based method) on the training data, and evaluating its
#'   classification performance on the test set. Optionally, it can
#'   incorporate covariate information (covYI) for the cut-off point.
#'
#' @param n Integer. The number of simulation experiments to run.
#'   Default is 1000.
#' @param df A data frame containing the complete dataset, including the
#'   target variable, regressors, and optional covariates.
#' @param X Character vector. Column names from `df` to be used as
#'   regressors in the model estimation. It can be a data frame (whose
#'   column names will be extracted) or a character vector. Defaults to
#'   all columns in `df` not specified as `y` or `C`.
#' @param y Character string. The column name in `df` representing the
#'   binary target variable (0 or 1). It can be a data frame (whose
#'   first column name will be extracted) or a character string.
#'   Default is "y".
#' @param C Character vector or `NULL`. Column names from `df` to be used as
#'   covariate variables in covYI. It can be a data frame or a character
#'   vector. Default is `NULL`.
#' @param model_estimation_function Function. The function to use for model
#'   estimation. Expected values are "plr_estimation", "psvm_estimation",
#'   "AucPR_estimation".
#' @param model_prediction_function Function. The function to use for model
#'   prediction. Expected values are "plr_predict", "psvm_predict",
#'   "AucPR_predict".
#' @param model_type Character string. The specific model to use in the
#'   estimation. Examples include "logLasso", "logElasticNet", "logSCAD",
#'   "logMCP", "SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM", "AucPR_L1",
#'   "AucPR_EN".
#' @param lambda Numeric. The penalization parameter for the regressors `X`.
#'   Must be a single non-negative value.
#' @param tau Numeric. The penalization parameter for the covariates `C`
#'   in covYI. Default is 0, indicating no penalization. Must be a single
#'   non-negative value.
#' @param w_g Numeric. A value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param train_data_proportion Numeric. A numeric value between 0 and 1. The proportion
#'   of the total observations in `df` to be used for the training set in each
#'   simulation split. The split is performed using stratified sampling on `y`
#'   to maintain class balance. Default is 0.7 (70\% training, 30\% testing).
#' @param trace Integer (0, 1, or 2). Controls the verbosity of the output.
#'   2 for all steps, 1 for main results, 0 for no output. Default is 1.
#' @param alpha Numeric. A value between 0 and 1, the elastic-net mixing
#'   parameter for the primary model's penalty. Default is 0.5.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param used_cores Integer. Number of CPU cores to use for parallel
#'   processing. If 0, it automatically uses 70\% of available physical cores.
#'   If 1, no parallelization is adopted. Default is 1.
#' @param scaling Logical. If `TRUE`, the dataset `df` is scaled
#'   (mean 0, variance 1) before processing. Default is `FALSE`.
#' @param c_function_of_covariates Logical. If `TRUE`, `covYI` is used to
#'   estimate the cut-off point as a function of the covariate information `C`.
#'   If `FALSE`, the covariate information is ignored for the cutoff.
#'   Default is `FALSE`.
#' @param alpha_g Numeric. A value between 0 and 1, the elastic-net mixing
#'   parameter for the `covYI` penalty. Default is 0.5.
#' @param penalty_g Character string. The penalty type for `covYI`. To be
#'   chosen among "L12", "L1", "EN", "SCAD", "MCP". Default is "L1".
#' @param kernel_g Character string. The kernel type to use for the estimation
#'   of the density function in `covYI`. (Note: tested primarily for "gaussian").
#'   Default is "gaussian".
#' @param a1_g Numeric. The `a` parameter for SCAD and MCP penalties in `covYI`.
#'   Default is 3.7.
#' @param a2_g Numeric. The `b` parameter for the MCP penalty in `covYI`.
#'   Default is 3.0.
#' @param trend_g Character string. For `covYI` optimization. If "monotone",
#'   mmAPG is used; if "nonmonotone", mnmAPG is used. Default is "monotone".
#' @param gamma_start_input Numeric vector or `NULL`. A specific starting point
#'   for gamma coefficients in `covYI`. Default is `NULL`, meaning no custom
#'   starting point is provided.
#' @param gamma_start_default Character string. Sets the default starting point
#'   of gamma coefficients if `gamma_start_input` is `NULL`. If "zeros", it
#'   starts with all zero values; if "corr", it starts with the correlation of
#'   every covariate with the target variable. Default is "zeros".
#' @param regressors_gammas Numeric vector or `NULL`. A vector containing the
#'   indices of the true non-zero gamma coefficients (if known), used for
#'   performance evaluation of variable selection. Default is `NULL`.
#' @param max_iter_g Integer. Maximum number of iterations for the `covYI`
#'   optimization algorithm (mmAPG and mnmAPG). Default is 10000.
#' @param delta_g Numeric. Parameter for the convergence condition of the
#'   `covYI` optimization algorithm. Default is 1e-5.
#' @param max_alpha_g Numeric. Maximum value of the step-size parameter alpha
#'   in `covYI`'s backtracking line-search. Default is 10000.
#' @param stepsizeShrink_g Numeric. Parameter to adjust the step-size in the
#'   backtracking line-search for `covYI` optimization. Takes values between
#'   0 and 1. Closer to 1 implies more accurate estimation but longer
#'   computation time. Default is 0.8.
#' @param min_alpha_g Numeric. Minimum value of the step-size parameter alpha
#'   in `covYI`'s backtracking line-search. Default is 1e-12.
#' @param convergence_error_g Numeric. The error threshold to accept for
#'   considering the `covYI` algorithm converged. Default is 1e-7.
#' @param run_aauc Logical. If `FALSE`, the aAUC and aYI measures are not
#'   computed, which can save estimation time if not requested. Default is `FALSE`.
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If `NULL`, output goes to the console.
#'   Default is "log_sim_Other_real.txt".
#'
#' @return A list containing the aggregated classification measures and other
#'   results related to every simulation experiment.
#'   \item{model_type}{Character string. The type of model used for estimation.}
#'   \item{simulation_time}{Numeric. Total time taken for the entire simulation
#'     study in minutes.}
#'   \item{estimation_time_original_method}{Numeric. Mean estimation time for
#'     the primary model across all simulations.}
#'   \item{estimation_time_covYI}{Numeric. Mean estimation time for the covYI
#'     method across all simulations (0 if `c_function_of_covariates` is `FALSE`).}
#'   \item{used_cores}{Integer. The number of CPU cores used for parallelization.}
#'   \item{n}{Integer. The number of simulation experiments performed.}
#'   \item{lambda}{Numeric. The penalization parameter for the primary model.}
#'   \item{tau}{Numeric. The penalization parameter for the covYI method.}
#'   \item{auc}{Matrix (n x 2). AUC values for training and testing sets
#'     across all simulations for the primary model.}
#'   \item{youden_index}{Matrix (n x 2). Youden's Index values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{sensitivity}{Matrix (n x 2). Sensitivity values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{specificity}{Matrix (n x 2). Specificity values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{geometric_mean}{Matrix (n x 2). Geometric Mean values for training
#'     and testing sets across all simulations for the primary model.}
#'   \item{fdr}{Matrix (n x 2). False Discovery Rate values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{mcc}{Matrix (n x 2). Matthews Correlation Coefficient values for
#'     training and testing sets across all simulations for the primary model.}
#'   \item{corrclass}{Matrix (n x 2). Correct Classification Rate values for
#'     training and testing sets across all simulations for the primary model.}
#'   \item{auc_covYI, aauc_covYI, aYI_covYI, youden_index_covYI,
#'     sensitivity_covYI, specificity_covYI, geometric_mean_covYI,
#'     fdr_covYI, mcc_covYI, corrclass_covYI}{Matrices (n x 2).
#'     Classification measures (AUC, aAUC, aYI, Youden's Index, Sensitivity,
#'     Specificity, Geometric Mean, FDR, MCC, Correct Classification Rate)
#'     for training and testing sets across all simulations for the covYI
#'     method (only if `c_function_of_covariates` is `TRUE`).}
#'   \item{n_total_var_betas}{Matrix (n x 1). Total number of regressors for
#'     betas in each simulation.}
#'   \item{n_predicted_zeros_betas}{Matrix (n x 1). Number of predicted zero
#'     beta coefficients in each simulation.}
#'   \item{n_predicted_non_zeros_betas}{Matrix (n x 1). Number of predicted
#'     non-zero beta coefficients in each simulation.}
#'   \item{n_caught_betas}{Matrix (n x 1). Number of true non-zero beta
#'     coefficients correctly identified as non-zero in each simulation.}
#'   \item{n_non_caught_betas}{Matrix (n x 1). Number of true non-zero beta
#'     coefficients incorrectly identified as zero in each simulation.}
#'   \item{n_caught_zero_betas}{Matrix (n x 1). Number of true zero beta
#'     coefficients correctly identified as zero in each simulation.}
#'   \item{n_zero_not_caught_betas}{Matrix (n x 1). Number of true zero beta
#'     coefficients incorrectly identified as non-zero in each simulation.}
#'   \item{n_total_var_gammas}{Matrix (n x 1). Total number of covariates for
#'     gammas in each simulation (only if `c_function_of_covariates` is `TRUE`).}
#'   \item{n_predicted_zeros_gammas}{Matrix (n x 1). Number of predicted zero
#'     gamma coefficients in each simulation (only if `c_function_of_covariates`
#'     is `TRUE`).}
#'   \item{n_predicted_non_zeros_gammas}{Matrix (n x 1). Number of predicted
#'     non-zero gamma coefficients in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_caught_gammas}{Matrix (n x 1). Number of true non-zero gamma
#'     coefficients correctly identified as non-zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_non_caught_gammas}{Matrix (n x 1). Number of true non-zero gamma
#'     coefficients incorrectly identified as zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_caught_zero_gammas}{Matrix (n x 1). Number of true zero gamma
#'     coefficients correctly identified as zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_zero_not_caught_gammas}{Matrix (n x 1). Number of true zero gamma
#'     coefficients incorrectly identified as non-zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{betas_times_selected}{Matrix (p x 1). For each regressor, the number
#'     of times its coefficient was estimated as non-zero across all simulations.
#'     `p` is the number of regressors in `X`.}
#'   \item{gammas_times_selected}{Matrix (q x 1). For each covariate, the number
#'     of times its coefficient was estimated as non-zero across all simulations
#'     (only if `c_function_of_covariates` is `TRUE`). `q` is the number of
#'     covariates in `C`.}
#'   \item{roc_spec_points}{Numeric vector. The fixed specificity points used
#'     for extracting ROC curve sensitivity values.}
#'   \item{roc_sens_points_train}{List of numeric vectors. Sensitivity points
#'     for ROC curves on training data for each simulation, corresponding to
#'     `roc_spec_points`.}
#'   \item{roc_sens_points_test}{List of numeric vectors. Sensitivity points
#'     for ROC curves on testing data for each simulation, corresponding to
#'     `roc_spec_points`.}
#'   \item{regressors_betas}{Numeric vector or `NULL`. The input `regressors_betas`
#'     parameter, indicating true non-zero beta indices.}
#'   \item{regressors_gammas}{Numeric vector or `NULL`. The input `regressors_gammas`
#'     parameter, indicating true non-zero gamma indices.}
#'   \item{input_parameters}{\code{character vector}. A list containing the
#'     input parameters.}
#'   \item{c_function_of_covariates}{Logical. The input `c_function_of_covariates`
#'     parameter.}
#'   \item{run_aauc}{Logical. The input `run_aauc` parameter.}
#'   \item{betas}{List of numeric vectors. The estimated beta coefficients from
#'     each simulation run. Each element of the list is a named numeric vector
#'     of coefficients for that simulation.}
#'   \item{gammas}{List of numeric vectors. The estimated gamma coefficients from
#'     each simulation run (only if `c_function_of_covariates` is `TRUE`). Each
#'     element of the list is a named numeric vector of coefficients for that
#'     simulation, including the 'const' term.}
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_betas <- sim_data$nregressors
#' regressors_gammas <- sim_data$ncovariates
#' n <- 4 #number of simulations
#' c_function_of_covariates <- TRUE
#'
#' #Possible choices for "method":
#' #"logLasso", "logElasticNet", "logSCAD", "logMCP",  "SCADSVM", "ElasticSCADSVM",
#' #"l1SVM", "enSVM", "AucPR_L1", "AucPR_EN"
#' method <- "logLasso"
#' model_estimation_function <- plr_estimation
#' model_prediction_function <- plr_predict
#'
#' # Run a small simulation study
#' results_study <- model_simulation_study(
#'   n = n, # Small number of simulations for example
#'   df = df,
#'   X = X,
#'   y = y,
#'   C = C,
#'   lambda = 0.1,
#'   tau = 0.05,
#'   model_estimation_function = model_estimation_function,
#'   model_prediction_function = model_prediction_function,
#'   model_type = method,
#'   trace = 1,
#'   gamma_start_default = "zeros",
#'   alpha_g = 0.5,
#'   penalty_g = "L12", # Options: "L12", "L1", "EN", "SCAD", "MCP"
#'   used_cores = 1,
#'   c_function_of_covariates = c_function_of_covariates,
#'   regressors_betas = regressors_betas,
#'   regressors_gammas = regressors_gammas,
#'   max_iter_g = 8 #<---- Reduced to speed up estimation
#' )
#'
#' # Print some of the results
#' cat("Simulation Time: ", format(results_study$simulation_time, digits = 4))
#' cat("\nCCR first-step method:\n")
#' print(results_study$corrclass)
#' cat("Betas Times Selected:\n")
#' print(results_study$betas_times_selected[results_study$betas_times_selected > 1])
#' if (results_study$input_parameters$c_function_of_covariates) {
#'   cat("CCR covYI:\n")
#'   print(results_study$corrclass_covYI)
#'   cat("\nGammas Times Selected:\n")
#'   print(results_study$gammas_times_selected[results_study$gammas_times_selected > 1])
#' }
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom pROC roc coords
#' @importFrom tools file_path_sans_ext file_ext
#' @export
model_simulation_study <- function(n = 1000, df, X = NULL, y = "y", C = NULL,
                                             model_estimation_function,
																						 model_prediction_function,
                                             model_type, lambda, tau = 0, w_g = 0.5,
																						 train_data_proportion = 0.7,
																						 trace = 1, alpha = 0.5,
																						 regressors_betas = NULL,
                                             used_cores = 1, scaling = FALSE,
                                             c_function_of_covariates = FALSE,
                                             alpha_g = 0.5, penalty_g = "L1",
																						 kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                                             trend_g = "monotone",
																						 gamma_start_input = NULL,
																						 gamma_start_default = "zeros",
                                             regressors_gammas = NULL,
																						 max_iter_g = 10000,
                                             delta_g = 1e-5, max_alpha_g = 10000,
                                             stepsizeShrink_g = 0.8, min_alpha_g = 1e-12,
																						 convergence_error_g = 1e-7,
                                             run_aauc = FALSE,
																						 log_file = "log_sim_Other_real.txt") {
  # Start calculation of estimation time
  start_time <- Sys.time()

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

  # Handle X, y, C parameters
  if (inherits(y, "data.frame")) y <- names(y)[1]
  if (!is.character(y) || length(y) != 1) stop("'y' must be a single column name.")
  if (!(y %in% names(df))) stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")

  if (!is.null(C)) {
    if (inherits(C, "data.frame")) C <- names(C)
    if (!is.character(C)) stop("'C' must be a character vector or NULL.")
    if (length(C) == 0) C <- NULL
    if (!is.null(C) && !all(C[C != "const"] %in% names(df))) {
      stop("Not all covariates in 'C' are found in 'df'.")
    }
  }

  if (is.null(X)) X <- setdiff(names(df), c(y, C))
  if (inherits(X, "data.frame")) X <- names(X)
  if (!is.character(X) || length(X) == 0) stop("'X' must be a character vector of column names or a data.frame.")
  if (!all(X %in% names(df))) stop("Not all regressors in 'X' are found in 'df'.")


  # Check c_function_of_covariates
  if (!is.logical(c_function_of_covariates)) {
    stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")
  }

  # Check if tau is provided and C is valid when c_function_of_covariates is TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
    if (is.null(C) || length(C) == 0) { stop("Parameter 'C' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

  # Check for "ID" column conflict: 'ID' is used internally
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  # Create an 'ID' column from row names and select relevant columns
  ID <- rownames(df)
  #df1 <- cbind(ID, df[, (names(df) %in% c(y,X,C))]) #OLD
  df1 <- cbind(ID, df[, c(y, X, C), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Further input validation for numeric parameters
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty_g %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty_g'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(max_iter_g) || length(max_iter_g) != 1 || max_iter_g < 2 || !is.integer(as.integer(max_iter_g))) {stop("Parameter 'max_iter_g' needs to be an integer and at least 2.")}
  if (!is.numeric(trace) || length(trace) != 1 || !(trace %in% c(0, 1, 2))) {stop("The parameter 'trace' has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!is.numeric(n) || length(n) != 1 || n < 2 || n != floor(n)) {stop("Parameter 'n' must be a single integer value and at least 2.")}
  if (!is.logical(scaling)) {stop("Parameter 'scaling' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_function_of_covariates)) {stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores < 0 || !is.integer(as.integer(used_cores))) {stop("Parameter 'used_cores' must be a single non-negative integer value.")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  if (!(kernel_g %in% valid_kernels)) {stop(paste0("Parameter 'kernel_g' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}


  # Standardize df1 if scaling is TRUE
  if (scaling == TRUE) {
    df1 <- scaling_df_for_pye (df = df1, X = colnames(df1[, names(df1) %in% c(X, C)]), y = "y")$df_scaled
  }

  # Generate seeds for each simulation run
  seeds <- 1:n

  # Names for results matrices
  names <- paste("seed", seeds, sep = "=")

  # Initialize matrices to store results for various performance metrics
  auc <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  youden_index <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  sensitivity <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  specificity <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  geometric_mean <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  fdr <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  mcc <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  corrclass <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  auc_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  aauc_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  aYI_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  youden_index_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  sensitivity_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  specificity_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  geometric_mean_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  fdr_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  mcc_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  corrclass_covYI <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))

  n_total_var_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_total_var_betas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_zeros"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_non_zeros"))
  n_caught_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_caught_betas"))
  n_non_caught_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_non_caught_betas"))
  n_caught_zero_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_caught_zero_betas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_zero_not_caught_betas"))
  n_total_var_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_total_var_gammas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_non_zeros_gammas"))
  n_caught_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_caught_gammas"))
  n_non_caught_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_non_caught_gammas"))
  n_caught_zero_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_zero_not_caught_gammas"))
  # betas_times_selected: how many times each beta regressor was selected (non-zero)
  betas_times_selected <- matrix(0, nrow = length(X), ncol = 1, dimnames = list(X, "n_times_beta_diff_zero"))
  # gammas_times_selected: how many times each gamma covariate was selected (non-zero)
  gammas_times_selected <- matrix(0, nrow = length(C), ncol = 1, dimnames = list(C, "n_times_gamma_diff_zero"))


  # Function to be computed in each simulation iteration (either sequentially or in parallel)
  func <- function(seed, n, df, X, y, C,
	                 lambda, tau, w_g,
									 train_data_proportion,
                   model_estimation_function,
                   model_prediction_function,
                   model_type,
                   trace, alpha,
                   a1_g, a2_g,
                   alpha_g,
                   penalty_g,
                   regressors_betas,
                   regressors_gammas,
                   c_function_of_covariates,
                   run_aauc,
                   max_iter_g,
                   delta_g,
                   max_alpha_g,
                   stepsizeShrink_g,
                   min_alpha_g,
                   convergence_error_g,
                   gamma_start_input,
                   gamma_start_default,
                   trend_g,
                   kernel_g) {

    if (trace %in% c(1, 2)) {
      cat("\n--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda = ", lambda, "; model_type = ", model_type, "\n")
      if (c_function_of_covariates == TRUE) {
        cat("tau = ", tau, "; weight_g:", w_g, "; penalty_g = ", penalty_g, "\n")
      }
    }

    # Set seed for reproducibility of the current simulation iteration
    set.seed(seed)

    # Stratified sampling for train/test split to ensure class balance
    df_sort <- df[order(df[[y]]), c("ID", y)] # Sort by target variable
    rows_y0 <- which(df_sort[[y]] == 0)
    rows_y1 <- which(df_sort[[y]] == 1)
    # Calculate sizes for stratified split (70% train, 30% test)
    size_0 <- as.integer(round(length(rows_y0) * train_data_proportion, 0))
    size_1 <- as.integer(round(length(rows_y1) * train_data_proportion, 0))
    split_0 <- sample(rows_y0, size = size_0, replace = FALSE)
    split_1 <- sample(rows_y1, size = size_1, replace = FALSE)
    split <- sort(c(split_0, split_1))
    train_df <- df[df$ID %in% df_sort[split, "ID"], ]
    test_df <- df[df$ID %in% df_sort[-split, "ID"], ]

    # --- Model Estimation on Training Data ---
    # Remove 'ID' column before passing to estimation function
    train_solution <- model_estimation_function(df = train_df[, names(train_df) != "ID", drop = FALSE],
                                                X = X, y = y,
                                                lambda = lambda,
                                                alpha = alpha, # Specific to Elastic-Net
                                                regressors_betas = regressors_betas,
                                                model_type = model_type,
                                                trace = trace)


    estimation_time_original_method <- train_solution$estimation_time
    z_hat <- train_solution$z_hat

    train_covYI_solution <- NULL
    estimation_time_covYI <- 0

    # --- covYI Estimation (if c_function_of_covariates is TRUE) ---
    if (c_function_of_covariates == TRUE) {
      gamma_start_input1 <- gamma_start_input
      if (!is.null(gamma_start_input) && (length(gamma_start_input) != (1 + length(C)))) {
        warning("gamma_start_input length does not match 'const' + C length. Defaulting to 'NULL'", call. = FALSE)
        gamma_start_input1 <- NULL
      }

			if (length(gamma_start_input1) == 0) {
			#if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
				gamma_start_input1 <- c(train_solution$c_hat, rep(0, length(C)))
				names(gamma_start_input1) <- c("const", C)
			} else {
				names(gamma_start_input1) <- c("const", C)
			}

      # Run covYI_KS_estimation. z_hat from primary model is passed as the 'z' variable.
      train_covYI_solution <- covYI_KS_estimation(df = cbind(train_df[, names(train_df) != "ID", drop = FALSE], z_hat = train_solution$z_hat[, "z_hat"]),
                                                  z = "z_hat",
																									y = y,
																									C = C,
																									tau = tau,
																									w = w_g,
                                                  gamma_start_input = gamma_start_input1,
                                                  gamma_start_default = gamma_start_default,
																									trace = trace,
                                                  alpha = alpha_g,
																									a1 = a1_g, a2 = a2_g,
																									penalty = penalty_g,
                                                  max_iter = max_iter_g,
                                                  min_alpha = min_alpha_g,
                                                  convergence_error = convergence_error_g,
                                                  regressors_gammas = regressors_gammas,
                                                  trend = trend_g,
                                                  stepsizeShrink = stepsizeShrink_g,
                                                  delta = delta_g, max_alpha = max_alpha_g, kernel = kernel_g,
                                                  run_aauc = run_aauc)

      estimation_time_covYI <- train_covYI_solution$estimation_time
			niter_covYI <- train_covYI_solution$niter
      z_hat <- train_covYI_solution$z_hat

    }

    # ROC curve @ certain levels on TRAIN data
    roc_spec_points <- seq(0, 1, by = 0.05)
    est_roc <- pROC::roc(as.numeric(train_df[[y]]), z_hat[, "z_hat"], levels = c(0, 1), direction = "<", quiet = TRUE)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)

    # --- Model Prediction on TEST data ---
    test_solution <- model_prediction_function(df = test_df[, names(test_df) != "ID", drop = FALSE],
                                               y = y,
                                               model_to_use = train_solution,
                                               trace = trace)

    z_hat <- test_solution$z_hat

    test_covYI_solution <- NULL

    # --- covYI Prediction on TEST data (if c_function_of_covariates is TRUE) ---
    if (c_function_of_covariates == TRUE) {
      #put "const" in C
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df = cbind(test_df[, names(test_df) != "ID", drop = FALSE], z_hat = test_solution$z_hat[, "z_hat"]),
                                       z = "z_hat", y = y, C = C1,
																			 gammas = train_covYI_solution$gammas_hat,
                                       tau = tau, w = w_g,
																			 kernel = kernel_g, alpha = alpha_g,
																			 a1 = a1_g, a2 = a2_g,
                                       penalty = penalty_g,
																			 prediction = TRUE, run_aauc = run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1, 2)) {
        cat("-> Results on the TEST SET on seed:", seed, "\n")
        cat("-> algorithm: covYI_KS_proximal_gradient_method ; ")
        visualize_gammas <- train_covYI_solution$gammas_hat[which(train_covYI_solution$gammas_hat != 0)]
        cat("tau:", tau, "; weight:", w_g, "; penalty:", penalty_g, "; covYI_KS:", getElement(test_covYI_solution,  paste0("covYI_KS_", penalty_g)), "; youden_index:", test_covYI_solution$youden_index, "; aYI:", test_covYI_solution$aYI, "; sensitivity:", test_covYI_solution$sensitivity, "; specificity:", test_covYI_solution$specificity, "; geometric_mean:", test_covYI_solution$geometric_mean, "; fdr:", test_covYI_solution$fdr, "; mcc:", test_covYI_solution$mcc, "; auc:", test_covYI_solution$auc, "; aauc:", test_covYI_solution$aauc, "; corrclass:", test_covYI_solution$corrclass, " \n")
        cat("TP:", test_covYI_solution$TP, "; TN:", test_covYI_solution$TN, "; FP:", test_covYI_solution$FP, "; FN:", test_covYI_solution$FN, "; gammas: \n")
        print(visualize_gammas)
				cat("Estimation time:", estimation_time_covYI, "; Number of iterations:", niter_covYI, "\n\n\n")
      }

    }

    #ROC curve @ certain levels
    est_roc <-  pROC::roc(as.numeric(getElement(test_df, y)), z_hat[, "z_hat"], levels = c(0, 1), direction = "<", quiet = TRUE)
    roc_sens_points_test <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)

    return(list(estimation_time_original_method = estimation_time_original_method,
                estimation_time_covYI = estimation_time_covYI,
                seed = seed,
                train_solution = train_solution,
                train_covYI_solution = train_covYI_solution,
                test_solution = test_solution,
                test_covYI_solution = test_covYI_solution,
                roc_spec_points = roc_spec_points,
                roc_sens_points_train = roc_sens_points_train,
                roc_sens_points_test = roc_sens_points_test))
  }

  if (trace %in% c(1, 2)) {
    cat("------------------------------------------------------------------\n")
    cat("|         Starting simulation study with", n, "simulations         |\n")
    cat("------------------------------------------------------------------\n")
  }

  # --- Parallel / Sequential Execution ---
  cl <- NULL # Initialize cluster object to NULL
  if (used_cores > 1) {
    max.cores <- parallel::detectCores(logical = FALSE)
    if (used_cores > max.cores) {
      warning("The number of specified cores (", used_cores, ") is larger than the number of physical cores available (", max.cores, ")!")
    }

    setup_strategy <- ifelse(.Platform$OS.type == "windows", "sequential", "parallel")

    if (!is.null(log_file)) {
      if (!is.character(log_file) || length(log_file) != 1) {stop("log_file must be a character string specifying the path to the log file.")}
      # Validate that the directory exists or can be created
      log_dir <- dirname(log_file)
      # Check if directory exists, allowing current directory
      if (!dir.exists(log_dir) && log_dir != ".") stop("Directory for log_file does not exist: ", log_dir)

      # Add timestamp to log file to prevent overwriting
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      log_base <- tools::file_path_sans_ext(basename(log_file))
      log_ext  <- tools::file_ext(log_file)
      final_log_file <- file.path(log_dir, paste0(log_base, "_", timestamp, ifelse(log_ext != "", paste0(".", log_ext), "")))

      cl <- parallel::makeCluster(used_cores, outfile = final_log_file, setup_strategy = setup_strategy)
      if (trace > 0) message("Parallel computing output is being written to: ", final_log_file)
    } else {
      cl <- parallel::makeCluster(used_cores, setup_strategy = setup_strategy) # No outfile, output goes to console
    }
    if (!("cluster" %in% class(cl))) stop("cl is not of class 'cl'; see ?makeCluster")

    # Ensure cluster is stopped on function exit, even if errors occur
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (trace %in% c(1, 2)) cat("Parallel computing for cross-validation started on", length(cl), "cores.\n")
    # Load the package on all workers
    parallel::clusterCall(cl, function() library(pye))

    # --- Fit the simulation ---
    parallel::clusterExport(cl, c("n", "df1", "X", "y", "C", "lambda", "seeds", "tau", "w_g",
		                              "model_estimation_function", "model_prediction_function",
																	"trace", "model_type", "alpha", "a1_g", "a2_g",
                                  "regressors_betas", "regressors_gammas", "alpha_g", "penalty_g",
                                  "c_function_of_covariates", "run_aauc", "trend_g", "kernel_g", "max_iter_g",
                                  "delta_g", "max_alpha_g", "stepsizeShrink_g", "min_alpha_g",
                                  "convergence_error_g", "run_aauc", "n", "func", "covYI_KS_estimation",
                                  "covYI_KS"), envir = environment())

    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed = x, n = n, df = df1,
		                                                                    X = X, y = y, C = C,
                                                                        lambda = lambda,
																																				tau = tau, w_g = w_g,
																																				train_data_proportion = train_data_proportion,
                                                                        model_estimation_function = model_estimation_function,
                                                                        model_prediction_function = model_prediction_function,
                                                                        model_type = model_type,
                                                                        trace = trace, alpha = alpha,
                                                                        a1_g = a1_g, a2_g = a2_g,
                                                                        alpha_g = alpha_g,
                                                                        penalty_g = penalty_g,
                                                                        regressors_betas = regressors_betas,
                                                                        regressors_gammas = regressors_gammas,
                                                                        c_function_of_covariates = c_function_of_covariates,
                                                                        run_aauc = run_aauc,
                                                                        trend_g = trend_g,
                                                                        kernel_g = kernel_g,
                                                                        max_iter_g = max_iter_g, delta_g = delta_g,
                                                                        max_alpha_g = max_alpha_g,
                                                                        stepsizeShrink_g = stepsizeShrink_g,
                                                                        min_alpha_g = min_alpha_g,
                                                                        convergence_error_g = convergence_error_g,
                                                                        gamma_start_input = gamma_start_input,
                                                                        gamma_start_default = gamma_start_default))
  } else {
    # Sequential execution
    cat("Running simulation in sequential mode (used_cores = 1).\n")
    simulation_study <- lapply(seeds, function(x) func(seed = x, n = n, df = df1, X = X, y = y, C = C,
                                                       lambda = lambda,
																											 tau = tau, w_g = w_g,
																											 train_data_proportion = train_data_proportion,
                                                       model_estimation_function = model_estimation_function,
                                                       model_prediction_function = model_prediction_function,
                                                       model_type = model_type,
                                                       trace = trace, alpha = alpha,
                                                       a1_g = a1_g, a2_g = a2_g,
                                                       alpha_g = alpha_g,
                                                       penalty_g = penalty_g,
                                                       regressors_betas = regressors_betas,
                                                       regressors_gammas = regressors_gammas,
                                                       c_function_of_covariates = c_function_of_covariates,
                                                       run_aauc = run_aauc,
                                                       trend_g = trend_g,
                                                       kernel_g = kernel_g,
                                                       max_iter_g = max_iter_g, delta_g = delta_g,
                                                       max_alpha_g = max_alpha_g,
                                                       stepsizeShrink_g = stepsizeShrink_g,
                                                       min_alpha_g = min_alpha_g,
                                                       convergence_error_g = convergence_error_g,
                                                       gamma_start_input = gamma_start_input,
                                                       gamma_start_default = gamma_start_default))

  }

  estimation_time_original_method <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  #fill the matrices
  #measures on the train set
  list_of_measures <- c("auc", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][, 1], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i]))))
    mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }

  if (c_function_of_covariates == TRUE) {
    list_of_measures_covYI <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, list_of_measures_covYI[i])))
      #assign(list_of_measures[i][, 1], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      assign(paste0(list_of_measures_covYI[i], "_covYI"), mes)
    }

    n_total_var_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_total_var_gammas")))
    n_predicted_zeros_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_zeros_gammas")))
    n_predicted_non_zeros_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_non_zeros_gammas")))
    n_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_gammas")))
    n_non_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_non_caught_gammas")))
    n_caught_zero_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_zero_gammas")))
    n_zero_not_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_zero_not_caught_gammas")))
    gammas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, paste0("gammas_hat_", penalty_g))) != 0)

  }

  n_total_var_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var")))
  n_predicted_zeros_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros")))
  n_predicted_non_zeros_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros")))
  n_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_non_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_caught_zero_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero")))
  n_zero_not_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught")))
  betas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "betas_hat")) != 0)

  # Measures on test set
  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][, 2], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i]))))
    mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if (c_function_of_covariates == TRUE) {
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      #assign(list_of_measures[i][, 2], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures_covYI[i])))
      assign(paste0(list_of_measures_covYI[i], "_covYI"), mes)
    }
  }

  # ROC curve data
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_test)

  betas_start <- getElement(simulation_study[[1]]$train_solution, "betas_start")
  betas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_solution, "betas_hat"))
  if (c_function_of_covariates == TRUE) {
    gammas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))
  } else {
    gammas <- NULL
  }

  # End computing estimation time
  simulation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(1, 2)) {
    cat("Total Simulation Time: ", format(simulation_time, digits = 4), " mins.\n")
  }

  results <- list(model_type = model_type,
                  simulation_time = simulation_time,
                  estimation_time_original_method = estimation_time_original_method,
                  estimation_time_covYI = estimation_time_covYI,
                  used_cores = used_cores,
                  n = n, lambda = lambda,
                  tau = tau,
                  auc = auc,
                  youden_index = youden_index,
                  sensitivity = sensitivity, specificity = specificity,
                  geometric_mean = geometric_mean, fdr = fdr,
                  mcc = mcc, corrclass = corrclass,
                  auc_covYI = auc_covYI,
                  aauc_covYI = aauc_covYI,
                  aYI_covYI = aYI_covYI,
                  youden_index_covYI = youden_index_covYI,
                  sensitivity_covYI = sensitivity_covYI,
                  specificity_covYI = specificity_covYI,
                  geometric_mean_covYI = geometric_mean_covYI,
                  fdr_covYI = fdr_covYI,
                  mcc_covYI = mcc_covYI,
                  corrclass_covYI = corrclass_covYI,
                  n_total_var_betas = n_total_var_betas,
                  n_predicted_zeros_betas = n_predicted_zeros_betas,
                  n_predicted_non_zeros_betas = n_predicted_non_zeros_betas,
                  n_caught_betas = n_caught_betas,
                  n_non_caught_betas = n_non_caught_betas,
                  n_caught_zero_betas = n_caught_zero_betas,
                  n_zero_not_caught_betas = n_zero_not_caught_betas,
                  n_total_var_gammas = n_total_var_gammas,
                  n_predicted_zeros_gammas = n_predicted_zeros_gammas,
                  n_predicted_non_zeros_gammas = n_predicted_non_zeros_gammas,
                  n_caught_gammas = n_caught_gammas,
                  n_non_caught_gammas = n_non_caught_gammas,
                  n_caught_zero_gammas = n_caught_zero_gammas,
                  n_zero_not_caught_gammas = n_zero_not_caught_gammas,
                  betas_times_selected = betas_times_selected,
                  gammas_times_selected = gammas_times_selected,
                  betas_start = betas_start,
                  roc_spec_points = roc_spec_points,
                  roc_sens_points_train = roc_sens_points_train,
                  roc_sens_points_test = roc_sens_points_test,
                  regressors_betas = regressors_betas,
                  regressors_gammas = regressors_gammas,
									input_parameters = list(
                    n = n,
                    df = df,
                    X = X,
                    y = y,
                    C = C,
                    model_estimation_function = model_estimation_function,
                    model_prediction_function = model_prediction_function,
                    model_type = model_type,
                    lambda = lambda,
                    tau = tau,
                    w_g = w_g,
                    trace = trace,
                    alpha = alpha,
                    regressors_betas = regressors_betas,
                    used_cores = used_cores,
                    scaling = scaling,
                    c_function_of_covariates = c_function_of_covariates,
                    alpha_g = alpha_g,
                    penalty_g = penalty_g,
                    kernel_g = kernel_g,
                    a1_g = a1_g,
                    a2_g = a2_g,
                    trend_g = trend_g,
                    gamma_start_input = gamma_start_input,
                    gamma_start_default = gamma_start_default,
                    regressors_gammas = regressors_gammas,
                    max_iter_g = max_iter_g,
                    delta_g = delta_g,
                    max_alpha_g = max_alpha_g,
                    stepsizeShrink_g = stepsizeShrink_g,
                    min_alpha_g = min_alpha_g,
                    convergence_error_g = convergence_error_g,
                    run_aauc = run_aauc,
                    log_file = log_file
                  ),
                  betas = betas, gammas = gammas)

  return(results)
}













#' @title Simulation Study for Penalized Models on Synthetic Data
#'
#' @description This function performs a simulation study by repeatedly
#'   generating a synthetic dataset, splitting it into training and testing
#'   sets, estimating a competitor method (e.g., penalized logistic regression,
#'   penalized SVM, or penalized AUC-based method) on the training data, and
#'   evaluating its classification performance on the test set. Optionally, it
#'   can incorporate covariate information (covYI) for the cut-off point.
#'
#' @param n Integer. The number of simulation experiments to run.
#'   Default is 1000.
#' @param rows_train Integer. Number of rows for the training sample.
#'   Default is 50, suitable for a high-dimensional setting.
#' @param rows_test Integer. Number of rows for the test sample. Default is 1000.
#'   It's recommended to create a test sample much larger than the training
#'   sample to evaluate the method/model on more data.
#' @param cols Integer. Number of regressor variables for both training and
#'   test samples. Default is 2000, suitable for a high-dimensional setting.
#' @param cols_cov Integer. Number of covariate variables for both training and
#'   test samples. Default is 20. Increase for a higher-dimensional setting.
#' @param max_rho Numeric. The maximum correlation coefficient used
#'   in the latent factor model for generating variables. Higher values lead
#'   to more correlated features. Must be between 0 and 1. Default is 0.5.
#' @param mu Numeric vector. Mean of the (multivariate) normal distribution
#'   of the regressors. Default is 0.
#' @param mu_cov Numeric vector. Mean of the (multivariate) normal distribution
#'   of the covariates. Default is 0.
#' @param model_estimation_function Function. The function to use for model
#'   estimation. Expected values are "plr_estimation", "psvm_estimation",
#'   "AucPR_estimation".
#' @param model_prediction_function Function. The function to use for model
#'   prediction. Expected values are "plr_predict", "psvm_predict",
#'   "AucPR_predict".
#' @param model_type Character string. The specific model to use in the
#'   estimation. Examples include "logLasso", "logElasticNet", "logSCAD",
#'   "logMCP", "SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM", "AucPR_L1",
#'   "AucPR_EN".
#' @param lambda Numeric. The penalization parameter for the regressors `X`.
#'   Must be a single non-negative value.
#' @param tau Numeric. The penalization parameter for the covariates `C`
#'   in covYI. Default is 0, indicating no penalization. Must be a single
#'   non-negative value.
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param trace Integer (0, 1, or 2). Controls the verbosity of the output.
#'   2 for all steps, 1 for main results, 0 for no output. Default is 1.
#' @param alpha Numeric. A value between 0 and 1, the elastic-net mixing
#'   parameter for the primary model's penalty. Default is 0.5.
#' @param used_cores Integer. Number of CPU cores to use for parallel
#'   processing. If 0, it automatically uses 70\% of available physical cores.
#'   If 1, no parallelization is adopted. Default is 1.
#' @param c_function_of_covariates Logical. If `TRUE`, `covYI` is used to
#'   estimate the cut-off point as a function of the covariate information `C`.
#'   If `FALSE`, the covariate information is ignored for the cutoff.
#'   Default is `FALSE`.
#' @param alpha_g Numeric. A value between 0 and 1, the elastic-net mixing
#'   parameter for the `covYI` penalty. Default is 0.5.
#' @param penalty_g Character string. The penalty type for `covYI`. To be
#'   chosen among "L12", "L1", "EN", "SCAD", "MCP". Default is "L1".
#' @param kernel_g Character string. The kernel type to use for the estimation
#'   of the density function in `covYI`. (Note: tested primarily for "gaussian").
#'   Default is "gaussian".
#' @param a1_g Numeric. The `a` parameter for SCAD and MCP penalties in `covYI`.
#'   Default is 3.7.
#' @param a2_g Numeric. The `b` parameter for the MCP penalty in `covYI`.
#'   Default is 3.0.
#' @param trend_g Character string. For `covYI` optimization. If "monotone",
#'   mmAPG is used; if "nonmonotone", mnmAPG is used. Default is "monotone".
#' @param gamma_start_input Numeric vector or `NULL`. A specific starting point
#'   for gamma coefficients in `covYI`. Default is `NULL`, meaning no custom
#'   starting point is provided.
#' @param gamma_start_default Character string. Sets the default starting point
#'   of gamma coefficients if `gamma_start_input` is `NULL`. If "zeros", it
#'   starts with all zero values; if "corr", it starts with the correlation of
#'   every covariate with the target variable. Default is "zeros".
#' @param max_iter_g Integer. Maximum number of iterations for the `covYI`
#'   optimization algorithm (mmAPG and mnmAPG). Default is 10000.
#' @param delta_g Numeric. Parameter for the convergence condition of the
#'   `covYI` optimization algorithm. Default is 1e-5.
#' @param max_alpha_g Numeric. Maximum value of the step-size parameter alpha
#'   in `covYI`'s backtracking line-search. Default is 10000.
#' @param stepsizeShrink_g Numeric. Parameter to adjust the step-size in the
#'   backtracking line-search for `covYI` optimization. Takes values between
#'   0 and 1. Closer to 1 implies more accurate estimation but longer
#'   computation time. Default is 0.8.
#' @param min_alpha_g Numeric. Minimum value of the step-size parameter alpha
#'   in `covYI`'s backtracking line-search. Default is 1e-12.
#' @param convergence_error_g Numeric. The error threshold to accept for
#'   considering the `covYI` algorithm converged. Default is 1e-7.
#' @param run_aauc Logical. If `FALSE`, the aAUC and aYI measures are not
#'   computed, which can save estimation time if not requested. Default is `FALSE`.
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If `NULL`, output goes to the console.
#'   Default is "log_sim_Other_synthetic.txt".
#'
#' @return A list containing the aggregated classification measures and other
#'   results related to every simulation experiment.
#'   \item{model_type}{Character string. The type of model used for estimation.}
#'   \item{simulation_time}{Numeric. Total time taken for the entire simulation
#'     study in minutes.}
#'   \item{estimation_time_original_method}{Numeric. Mean estimation time for
#'     the primary model across all simulations.}
#'   \item{estimation_time_covYI}{Numeric. Mean estimation time for the covYI
#'     method across all simulations (0 if `c_function_of_covariates` is `FALSE`).}
#'   \item{used_cores}{Integer. The number of CPU cores used for parallelization.}
#'   \item{n}{Integer. The number of simulation experiments performed.}
#'   \item{lambda}{Numeric. The penalization parameter for the primary model.}
#'   \item{tau}{Numeric. The penalization parameter for the covYI method.}
#'   \item{auc}{Matrix (n x 2). AUC values for training and testing sets
#'     across all simulations for the primary model.}
#'   \item{youden_index}{Matrix (n x 2). Youden's Index values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{sensitivity}{Matrix (n x 2). Sensitivity values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{specificity}{Matrix (n x 2). Specificity values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{geometric_mean}{Matrix (n x 2). Geometric Mean values for training
#'     and testing sets across all simulations for the primary model.}
#'   \item{fdr}{Matrix (n x 2). False Discovery Rate values for training and
#'     testing sets across all simulations for the primary model.}
#'   \item{mcc}{Matrix (n x 2). Matthews Correlation Coefficient values for
#'     training and testing sets across all simulations for the primary model.}
#'   \item{corrclass}{Matrix (n x 2). Correct Classification Rate values for
#'     training and testing sets across all simulations for the primary model.}
#'   \item{auc_covYI, aauc_covYI, aYI_covYI, youden_index_covYI,
#'     sensitivity_covYI, specificity_covYI, geometric_mean_covYI,
#'     fdr_covYI, mcc_covYI, corrclass_covYI}{Matrices (n x 2).
#'     Classification measures (AUC, aAUC, aYI, Youden's Index, Sensitivity,
#'     Specificity, Geometric Mean, FDR, MCC, Correct Classification Rate)
#'     for training and testing sets across all simulations for the covYI
#'     method (only if `c_function_of_covariates` is `TRUE`).}
#'   \item{n_total_var_betas}{Matrix (n x 1). Total number of regressors for
#'     betas in each simulation.}
#'   \item{n_predicted_zeros_betas}{Matrix (n x 1). Number of predicted zero
#'     beta coefficients in each simulation.}
#'   \item{n_predicted_non_zeros_betas}{Matrix (n x 1). Number of predicted
#'     non-zero beta coefficients in each simulation.}
#'   \item{n_caught_betas}{Matrix (n x 1). Number of true non-zero beta
#'     coefficients correctly identified as non-zero in each simulation.}
#'   \item{n_non_caught_betas}{Matrix (n x 1). Number of true non-zero beta
#'     coefficients incorrectly identified as zero in each simulation.}
#'   \item{n_caught_zero_betas}{Matrix (n x 1). Number of true zero beta
#'     coefficients correctly identified as zero in each simulation.}
#'   \item{n_zero_not_caught_betas}{Matrix (n x 1). Number of true zero beta
#'     coefficients incorrectly identified as non-zero in each simulation.}
#'   \item{n_total_var_gammas}{Matrix (n x 1). Total number of covariates for
#'     gammas in each simulation (only if `c_function_of_covariates` is `TRUE`).}
#'   \item{n_predicted_zeros_gammas}{Matrix (n x 1). Number of predicted zero
#'     gamma coefficients in each simulation (only if `c_function_of_covariates`
#'     is `TRUE`).}
#'   \item{n_predicted_non_zeros_gammas}{Matrix (n x 1). Number of predicted
#'     non-zero gamma coefficients in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_caught_gammas}{Matrix (n x 1). Number of true non-zero gamma
#'     coefficients correctly identified as non-zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_non_caught_gammas}{Matrix (n x 1). Number of true non-zero gamma
#'     coefficients incorrectly identified as zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_caught_zero_gammas}{Matrix (n x 1). Number of true zero gamma
#'     coefficients correctly identified as zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{n_zero_not_caught_gammas}{Matrix (n x 1). Number of true zero gamma
#'     coefficients incorrectly identified as non-zero in each simulation (only if
#'     `c_function_of_covariates` is `TRUE`).}
#'   \item{betas_times_selected}{Matrix (p x 1). For each regressor, the number
#'     of times its coefficient was estimated as non-zero across all simulations.
#'     `p` is the number of regressors in `X`.}
#'   \item{gammas_times_selected}{Matrix (q x 1). For each covariate, the number
#'     of times its coefficient was estimated as non-zero across all simulations
#'     (only if `c_function_of_covariates` is `TRUE`). `q` is the number of
#'     covariates in `C`.}
#'   \item{roc_spec_points}{Numeric vector. The fixed specificity points used
#'     for extracting ROC curve sensitivity values.}
#'   \item{roc_sens_points_train}{List of numeric vectors. Sensitivity points
#'     for ROC curves on training data for each simulation, corresponding to
#'     `roc_spec_points`.}
#'   \item{roc_sens_points_test}{List of numeric vectors. Sensitivity points
#'     for ROC curves on testing data for each simulation, corresponding to
#'     `roc_spec_points`.}
#'   \item{regressors_betas}{List of numeric vectors. The indices of the true
#'     non-zero beta coefficients for each simulation, as determined by the
#'     data generation process.}
#'   \item{regressors_gammas}{List of numeric vectors. The indices of the true
#'     non-zero gamma coefficients for each simulation, as determined by the
#'     data generation process (only if `c_function_of_covariates` is `TRUE`).}
#'   \item{input_parameters}{\code{character vector}. A list containing the
#'     input parameters.}
#'   \item{c_function_of_covariates}{Logical. The input `c_function_of_covariates`
#'     parameter.}
#'   \item{run_aauc}{Logical. The input `run_aauc` parameter.}
#'   \item{betas}{List of numeric vectors. The estimated beta coefficients from
#'     each simulation run. Each element of the list is a named numeric vector
#'     of coefficients for that simulation.}
#'   \item{gammas}{List of numeric vectors. The estimated gamma coefficients from
#'     each simulation run (only if `c_function_of_covariates` is `TRUE`). Each
#'     element of the list is a named numeric vector of coefficients for that
#'     simulation, including the 'const' term.}
#'
#' @examples
#' library(pye)
#'
#' # Define simulation parameters
#' cols <- 200 # Number of X variables
#' cols_cov <- 20 # Number of C variables
#' mu <- rep(0, cols) # Mean of each X var
#' mu_cov <- rep(0, cols_cov) # Mean of each C var
#'
#' # Number of simulation runs for this example
#' n <- 5 #number of simulations
#' lambda_to_use <- 0.15
#' tau_to_use <- 0.1
#' c_function_of_covariates <- TRUE
#'
#' #Possible choices for "method":
#' #"logLasso", "logElasticNet", "logSCAD", "logMCP",  "SCADSVM", "ElasticSCADSVM",
#' #"l1SVM", "enSVM", "AucPR_L1", "AucPR_EN"
#' method <- "logLasso"
#' model_estimation_function <- plr_estimation
#' model_prediction_function <- plr_predict
#'
#' sim_result <- model_simulation_study_synthetic_data(n = n,
#'   rows_train = 50,
#'   rows_test = 1000,
#'   cols = cols,
#'   cols_cov = cols_cov,
#'   max_rho = 0.2,
#'   mu = mu ,
#'   mu_cov = mu_cov,
#'   lambda = lambda_to_use,
#'   model_estimation_function = model_estimation_function,
#'   model_prediction_function = model_prediction_function,
#'   model_type = method,
#'   tau = tau_to_use,
#'   gamma_start_default = "zeros",
#'   a1_g = 3.7, a2_g = 3, # SCAD/MCP param for covYI
#'   penalty_g = "L12",
#'   trace = 1,
#'   alpha_g = 0.5,
#'   used_cores = 1,
#'   c_function_of_covariates = c_function_of_covariates,
#'   run_aauc = FALSE,
#'   max_iter_g = 10 #<---- reduced to speed up estimation
#' )
#'
#' # You can now access the results, e.g.:
#' print(sim_result$auc)
#' print(sim_result$betas_times_selected[sim_result$betas_times_selected > 1])
#' print(sim_result$aauc)
#' print(sim_result$gammas_times_selected[sim_result$gammas_times_selected > 1])
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom pROC coords roc
#' @importFrom tools file_path_sans_ext file_ext

#' @noRd
#' @keywords internal
model_simulation_study_synthetic_data <- function(n = 1000, rows_train = 50, rows_test = 1000, cols = 2000, cols_cov = 20,
                                                  max_rho = 0.5, mu = rep(0, cols), mu_cov = rep(0, cols_cov),
                                                  model_estimation_function, model_prediction_function, model_type,
                                                  lambda, tau = 0, w_g = 0.5, trace = 1, alpha = 0.5, used_cores = 1,
                                                  c_function_of_covariates = FALSE,
                                                  alpha_g = 0.5, penalty_g = "L1", kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                                                  trend_g = "monotone", gamma_start_input = NULL, gamma_start_default = "zeros",
                                                  max_iter_g = 10000, delta_g = 1e-5, max_alpha_g = 10000,
                                                  stepsizeShrink_g = 0.8, min_alpha_g = 1e-12, convergence_error_g = 1e-7,
                                                  run_aauc = FALSE, log_file = "log_sim_Other_synthetic.txt") {

  #start calc sim. time
  start_time <- Sys.time()

  #starting controls
  # Validate tau
  if (length(tau) != 1 || !is.numeric(tau) || tau < 0) {stop("Parameter 'tau' must be a single non-negative numeric value.")}
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (n < 2) {stop("n needs to be at least 2.")}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}
  if (!is.numeric(n) || length(n) != 1 || n < 2 || n != floor(n)) {stop("Parameter 'n' must be a single integer value and at least 2.")}
	if (!is.numeric(max_rho) || length(max_rho) != 1 || max_rho < 0 || max_rho > 1) {stop("Parameter 'max_rho' must be a single numeric value between 0 and 1.")}
	if (!is.numeric(alpha_g) || length(alpha_g) != 1 || alpha_g < 0 || alpha_g > 1) {stop("Parameter 'alpha_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(a1_g) || length(a1_g) != 1 || a1_g < 0) {stop("Parameter 'a1_g' must be a single non-negative numeric value.")}
  if (!is.numeric(a2_g) || length(a2_g) != 1 || a2_g < 0) {stop("Parameter 'a2_g' must be a single non-negative numeric value.")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel_g %in% valid_kernels)) {stop(paste0("Parameter 'kernel_g' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
	valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty_g %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty_g'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
  if (!(trend_g %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend_g' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  if (!(gamma_start_default %in% c("zeros", "corr"))) {stop("Parameter 'gamma_start_default' must be 'zeros' or 'corr'.")}
  if (!is.numeric(max_iter_g) || length(max_iter_g) != 1 || max_iter_g < 2 || !is.integer(as.integer(max_iter_g))) {stop("Parameter 'max_iter_g' needs to be an integer and at least 2.")}
  if (!is.numeric(delta_g) || length(delta_g) != 1 || delta_g <= 0) {stop("Parameter 'delta_g' must be a single positive numeric value.")}
  if (!is.numeric(max_alpha_g) || length(max_alpha_g) != 1 || max_alpha_g <= 0) {stop("Parameter 'max_alpha_g' must be a single positive numeric value.")}
  if (!is.numeric(stepsizeShrink_g) || length(stepsizeShrink_g) != 1 || stepsizeShrink_g <= 0 || stepsizeShrink_g >= 1) {stop("Parameter 'stepsizeShrink_g' must be a single numeric value between 0 and 1 (exclusive).")}
  if (!is.numeric(min_alpha_g) || length(min_alpha_g) != 1 || min_alpha_g <= 0) {stop("Parameter 'min_alpha_g' must be a single positive numeric value.")}
  if (!is.numeric(convergence_error_g) || length(convergence_error_g) != 1 || convergence_error_g <= 0) {stop("Parameter 'convergence_error_g' must be a single positive numeric value.")}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}

  # Check c_function_of_covariates
  if (!is.logical(c_function_of_covariates)) {
    stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")
  }

  # Check if tau exists when c_function_of_covariates = TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}

  # Generate seeds
  seeds <- 1:n

  # Names of the columns
  names <- paste("seed", seeds, sep = "=")

  # Function to be computed
  func <- function(seed, n, rows_train, rows_test,
                   mu = mu, mu_cov = mu_cov,
                   cols, cols_cov, max_rho,
									 lambda, tau, w_g,
                   model_estimation_function,
                   model_prediction_function,
                   model_type, trace, alpha,
                   gamma_start_input,
                   gamma_start_default,
                   max_iter_g,
                   delta_g,
                   min_alpha_g,
                   max_alpha_g,
                   stepsizeShrink_g,
                   trend_g, alpha_g,
                   penalty_g,
                   a1_g, a2_g,
                   kernel_g,
                   convergence_error_g,
                   c_function_of_covariates,
                   run_aauc) {

    if (trace %in% c(1, 2)) {
      cat("\n--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda = ", lambda, "; model_type = ", model_type, "\n")
      if (c_function_of_covariates == TRUE) {
        cat("tau = ", tau, "; weight_g:", w_g, "; penalty_g = ", penalty_g, "\n")
      }
    }

    # Create the dataframe (the seed makes different all the simulations)
    df <- create_sample_with_covariates(rows_train = rows_train, cols = cols, cols_cov = cols_cov, max_rho = max_rho, mu = mu,
                                        mu_cov = mu_cov, rows_test = rows_test, seed = seed)
                                        #varsN=varsN, varsB=varsB, varsE=varsE, varsP=varsP,
                                        #varN_cov=varN_cov, varB_cov=varB_cov, varE_cov=varE_cov, varP_cov=varP_cov)
    train_df <- df$train_df_scaled
    test_df <- df$test_df_scaled

    X <- df[[X]]
    y <- df[[y]]
    C <- df[[C]]
    regressors_betas <- df$nregressors
    regressors_gammas <- df$ncovariates # Do not include constant here, covYI handles it

    # Train the model
    train_solution <- model_estimation_function(df = train_df, X = X, y = y,
                                                lambda = lambda,
                                                alpha = alpha,
                                                regressors_betas = regressors_betas,
                                                model_type = model_type,
                                                trace = trace)

    estimation_time_original_method <- train_solution$estimation_time
    z_hat <- train_solution$z_hat

    train_covYI_solution <- NULL
    estimation_time_covYI <- 0

    # covYI estimation
    # Computing c
    if (c_function_of_covariates == TRUE) {
      gamma_start_input1 <- gamma_start_input
      if (!is.null(gamma_start_input) && (length(gamma_start_input) != (1 + length(C)))) {
        warning("gamma_start_input length does not match 'const' + C length. Defaulting to 'NULL'", call. = FALSE)
        gamma_start_input1 <- NULL
      }

			if (length(gamma_start_input1) == 0) {
			#if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
				gamma_start_input1 <- c(train_solution$c_hat, rep(0, length(C)))
				names(gamma_start_input1) <- c("const", C)
			} else {
				names(gamma_start_input1) <- c("const", C)
			}

      # Check C for null/empty here, as it's generated dynamically
      if (is.null(C) || length(C) == 0) {
        stop("C cannot be NULL or empty when c_function_of_covariates is TRUE, even in synthetic data generation.")
      }

      train_covYI_solution <- covYI_KS_estimation(df = cbind(train_df[, names(train_df) != "ID", drop = FALSE],
                                                  z_hat = train_solution$z_hat[, "z_hat"]),
                                                  z = "z_hat", y = y, C = C, tau = tau, w = w_g,
                                                  gamma_start_input = gamma_start_input1,
                                                  gamma_start_default = gamma_start_default,
                                                  trace = trace,
                                                  alpha = alpha_g,
                                                  a1 = a1_g, a2 = a2_g,
                                                  penalty = penalty_g,
                                                  max_iter = max_iter_g,
                                                  min_alpha = min_alpha_g,
                                                  convergence_error = convergence_error_g,
                                                  regressors_gammas = regressors_gammas,
                                                  trend = trend_g,
                                                  stepsizeShrink = stepsizeShrink_g,
                                                  delta = delta_g, max_alpha = max_alpha_g,
                                                  kernel = kernel_g,
                                                  run_aauc = run_aauc)

      estimation_time_covYI <- train_covYI_solution$estimation_time
			niter_covYI <- train_covYI_solution$niter
      z_hat <- train_covYI_solution$z_hat

    }

    #ROC curve @ certain levels
    est_roc <-  pROC::roc(as.numeric(getElement(train_df, y)), z_hat[, "z_hat"], levels = c(0, 1), direction = "<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)

    #test the results
    test_solution <- model_prediction_function(df = test_df,
                                               y = y,
                                               model_to_use = train_solution,
                                               trace = trace)

    z_hat <- test_solution$z_hat

    test_covYI_solution <- NULL

    if (c_function_of_covariates == TRUE) {
      #put "const" in C
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df = cbind(test_df[, names(test_df) != "ID", drop = FALSE], z_hat = test_solution$z_hat[, "z_hat"]),
                                      z = "z_hat", y = y, C = C1,
                                      gammas = train_covYI_solution$gammas_hat,
                                      tau = tau, w = w_g, kernel = kernel_g,
                                      alpha = alpha_g, a1 = a1_g, a2 = a2_g,
                                      penalty = penalty_g,
                                      prediction = TRUE, run_aauc = run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1, 2)) {
        cat("-> Results on the TEST SET on seed:", seed, "\n")
        cat("-> algorithm: covYI_KS_proximal_gradient_method ; ")
        visualize_gammas <- train_covYI_solution$gammas_hat[which(train_covYI_solution$gammas_hat != 0)]
        cat("tau:", tau, "; weight:", w_g, "; penalty:", penalty_g, "; covYI_KS:", getElement(test_covYI_solution,  paste0("covYI_KS_", penalty_g)), "; youden_index:", test_covYI_solution$youden_index, "; aYI:", test_covYI_solution$aYI, "; sensitivity:", test_covYI_solution$sensitivity, "; specificity:", test_covYI_solution$specificity, "; geometric_mean:", test_covYI_solution$geometric_mean, "; fdr:", test_covYI_solution$fdr, "; mcc:", test_covYI_solution$mcc, "; auc:", test_covYI_solution$auc, "; aauc:", test_covYI_solution$aauc, "; corrclass:", test_covYI_solution$corrclass, " \n")
        cat("TP:", test_covYI_solution$TP, "; TN:", test_covYI_solution$TN, "; FP:", test_covYI_solution$FP, "; FN:", test_covYI_solution$FN, "; gammas: \n")
        print(visualize_gammas)
				cat("Estimation time:", estimation_time_covYI, "; Number of iterations:", niter_covYI, "\n\n\n")
      }

    }

    #ROC curve @ certain levels
    est_roc <-  pROC::roc(as.numeric(getElement(test_df, y)), z_hat[, "z_hat"], levels = c(0, 1), direction = "<")
    roc_sens_points_test <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)


    return(list(estimation_time_original_method = estimation_time_original_method,
                estimation_time_covYI = estimation_time_covYI,
                seed = seed,
                train_solution = train_solution,
                train_covYI_solution = train_covYI_solution,
                test_solution = test_solution,
                test_covYI_solution = test_covYI_solution,
                roc_spec_points = roc_spec_points,
                roc_sens_points_train = roc_sens_points_train,
                roc_sens_points_test = roc_sens_points_test,
                regressors_betas = regressors_betas,
                regressors_gammas = regressors_gammas))
  }

  #simulate one dataset to take the names of the variables:
  #create the dataframe
  #set.seed(seed)
  #choose the regressors. I keep it commented since I keep default ones (below the default code)
  #varsN <- c(1, 2)
  #varsB <- c((cols / 4 + 1), (cols / 4 + 2))
  #varsE <- c((cols / 2 + 1), (cols / 2 + 2))
  #varsP <- c((cols / 4 * 3 + 1), (cols / 4 * 3 + 2))
  #varsN_cov <- 1
  #varsB_cov <- (cols_cov / 4 + 1)
  #varsE_cov <- (cols_cov / 2 + 1)
  #varsP_cov <- (cols_cov / 4 * 3 + 1)
  #compute the df (the seed makes different all the simulations)
  df_for_the_names <- create_sample_with_covariates(rows_train = rows_train, cols = cols, cols_cov = cols_cov, max_rho = max_rho, mu = mu,
                                                    mu_cov = mu_cov, rows_test = rows_test, seed = 1)
                                                    #varsN=varsN, varsB=varsB, varsE=varsE, varsP=varsP,
                                                    #varN_cov=varN_cov, varB_cov=varB_cov, varE_cov=varE_cov, varP_cov=varP_cov)

  if (trace %in% c(1, 2)) {
    cat("Real regressors of betas are: ", df_for_the_names$regressors, "\n")
    if (c_function_of_covariates == TRUE) {
      cat("Real regressors of gammas are: ", df_for_the_names$covariates, "\n\n\n")
    }
  }

  X <- df_for_the_names[[X]]
  #y <- df_for_the_names[[y]]
  C <- df_for_the_names[[C]]

  # Measures (train and test)
  auc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  youden_index <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  sensitivity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  specificity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  geometric_mean <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  fdr <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  mcc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  corrclass <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  auc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  aauc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  aYI_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  youden_index_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  sensitivity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  specificity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  geometric_mean_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  fdr_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  mcc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  corrclass_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames = list(names, c("train", "test")))
  n_total_var_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_total_var_betas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_predicted_zeros"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_predicted_non_zeros"))
  n_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_caught_betas"))
  n_non_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_non_caught_betas"))
  n_caught_zero_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_caught_zero_betas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_zero_not_caught_betas"))
  n_total_var_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_total_var_gammas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_predicted_non_zeros_gammas"))
  n_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_caught_gammas"))
  n_non_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_non_caught_gammas"))
  n_caught_zero_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames = list(names, "n_zero_not_caught_gammas"))
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))])), nrow =  ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))]), ncol = 1, dimnames = list(colnames(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))]), "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))])), nrow =  ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))]), ncol = 1, dimnames = list(colnames(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))]), "n_times_gamma_diff_zero"))
  #betas <- vector(mode = "list", length=length(seeds))
  #names(betas) <- names

  if (trace %in% c(1, 2)) {
    cat("------------------------------------------------------------------\n")
    cat("|         Starting simulation study with", n, "simulations         |\n")
    cat("------------------------------------------------------------------\n")
  }

  # --- Parallel / Sequential Execution ---
  cl <- NULL # Initialize cluster object to NULL
  if (used_cores > 1) {
    max.cores <- parallel::detectCores(logical = FALSE)
    if (used_cores > max.cores) {
      warning("The number of specified cores (", used_cores, ") is larger than the number of physical cores available (", max.cores, ")!")
    }

    setup_strategy <- ifelse(.Platform$OS.type == "windows", "sequential", "parallel")

    if (!is.null(log_file)) {
      if (!is.character(log_file) || length(log_file) != 1) {
        stop("log_file must be a character string specifying the path to the log file.")
      }
      # Validate that the directory exists or can be created
      log_dir <- dirname(log_file)
      # Check if directory exists, allowing current directory
      if (!dir.exists(log_dir) && log_dir != ".") stop("Directory for log_file does not exist: ", log_dir)

      # Add timestamp to log file to prevent overwriting
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      log_base <- tools::file_path_sans_ext(basename(log_file))
      log_ext  <- tools::file_ext(log_file)
      final_log_file <- file.path(log_dir, paste0(log_base, "_", timestamp, ifelse(log_ext != "", paste0(".", log_ext), "")))

      cl <- parallel::makeCluster(used_cores, outfile = final_log_file, setup_strategy = setup_strategy)
      if (trace > 0) message("Parallel computing output is being written to: ", final_log_file)
    } else {
      cl <- parallel::makeCluster(used_cores, setup_strategy = setup_strategy) # No outfile, output goes to console
    }
    if (!("cluster" %in% class(cl))) stop("cl is not of class 'cl'; see ?makeCluster")

    # Ensure cluster is stopped on function exit, even if errors occur
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (trace %in% c(1, 2)) cat("Parallel computing for cross-validation started on", length(cl), "cores.\n")

    # Load the package on all workers
    parallel::clusterCall(cl, function() library(pye))

    # --- Fit the simulation ---
    parallel::clusterExport(cl, c("n", "rows_train", "rows_test", "cols", "cols_cov", "max_rho", "lambda", "seeds", "tau", "w_g",
                                  "model_estimation_function", "model_prediction_function", "trace", "model_type",
                                  "alpha", "alpha_g", "penalty_g", "a1_g", "a2_g", "trend_g",
                                  "gamma_start_input", "gamma_start_default", "kernel_g",
                                  "max_iter_g", "delta_g", "max_alpha_g", "stepsizeShrink_g", "min_alpha_g", "convergence_error_g",
                                  "c_function_of_covariates", "run_aauc", "mu", "mu_cov", "func"), envir = environment())
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed = x, n = n, rows_train = rows_train, rows_test = rows_test,
                                                                        cols = cols, cols_cov = cols_cov, max_rho = max_rho,
                                                                        mu = mu, mu_cov = mu_cov,
                                                                        lambda = lambda, tau = tau,
																																				w_g = w_g,
                                                                        model_estimation_function = model_estimation_function,
                                                                        model_prediction_function = model_prediction_function,
                                                                        model_type = model_type,
                                                                        trace = trace, alpha = alpha,
                                                                        a1_g = a1_g, a2_g = a2_g,
                                                                        alpha_g = alpha_g,
                                                                        penalty_g = penalty_g,
                                                                        trend_g = trend_g,
                                                                        kernel_g = kernel_g,
                                                                        gamma_start_input = gamma_start_input,
                                                                        gamma_start_default = gamma_start_default,
                                                                        max_iter_g = max_iter_g,
                                                                        delta_g = delta_g, max_alpha_g = max_alpha_g,
                                                                        stepsizeShrink_g = stepsizeShrink_g,
                                                                        min_alpha_g = min_alpha_g,
                                                                        convergence_error_g = convergence_error_g,
                                                                        c_function_of_covariates = c_function_of_covariates,
                                                                        run_aauc = run_aauc))
  } else {
		# Sequential execution
		cat("Running simulation in sequential mode (used_cores = 1).\n")
    simulation_study <- lapply(seeds, function(x) func(seed = x, n = n, rows_train = rows_train, rows_test = rows_test,
                                                       cols = cols, cols_cov = cols_cov, max_rho = max_rho,
                                                       mu = mu, mu_cov = mu_cov,
                                                       lambda = lambda, tau = tau,
																											 w_g = w_g,
                                                       model_estimation_function = model_estimation_function,
                                                       model_prediction_function = model_prediction_function,
                                                       model_type = model_type,
                                                       trace = trace, alpha = alpha,
                                                       a1_g = a1_g, a2_g = a2_g,
                                                       alpha_g = alpha_g,
                                                       penalty_g = penalty_g,
                                                       trend_g = trend_g,
                                                       kernel_g = kernel_g,
                                                       gamma_start_input = gamma_start_input,
                                                       gamma_start_default = gamma_start_default,
                                                       max_iter_g = max_iter_g,
                                                       delta_g = delta_g, max_alpha_g = max_alpha_g,
                                                       stepsizeShrink_g = stepsizeShrink_g,
                                                       min_alpha_g = min_alpha_g,
                                                       convergence_error_g = convergence_error_g,
                                                       c_function_of_covariates = c_function_of_covariates,
                                                       run_aauc = run_aauc))

  }

  estimation_time_original_method <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  #fill the matrices
  #measures on the train set
  list_of_measures <- c("auc", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][, 1], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i]))))
    mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if (c_function_of_covariates == TRUE) {
    list_of_measures_covYI <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, list_of_measures_covYI[i])))
      #assign(list_of_measures[i][, 1], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      assign(paste0(list_of_measures_covYI[i], "_covYI"), mes)
    }

    n_total_var_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_total_var_gammas")))
    n_predicted_zeros_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_zeros_gammas")))
    n_predicted_non_zeros_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_non_zeros_gammas")))
    n_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_gammas")))
    n_non_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_non_caught_gammas")))
    n_caught_zero_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_zero_gammas")))
    n_zero_not_caught_gammas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_zero_not_caught_gammas")))
    gammas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, paste0("gammas_hat_", penalty_g))) != 0)

  }

  n_total_var_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var")))
  n_predicted_zeros_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros")))
  n_predicted_non_zeros_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros")))
  n_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_non_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_caught_zero_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero")))
  n_zero_not_caught_betas[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught")))
  betas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "betas_hat")) != 0)

  # Measures on test set
  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][, 2], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i]))))
    mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if (c_function_of_covariates == TRUE) {
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      #auc[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      #assign(list_of_measures[i][, 2], unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures_covYI[i])))
      assign(paste0(list_of_measures_covYI[i], "_covYI"), mes)
    }
  }

  #roc curve
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_test)

  betas_start <- getElement(simulation_study[[1]]$train_solution, "betas_start")
  betas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_solution, "betas_hat"))
  if (c_function_of_covariates == TRUE) {
    gammas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))
  } else {
    gammas <- NULL
  }

  regressors_betas <- simulation_study[[1]]$regressors_betas
  regressors_gammas <- simulation_study[[1]]$regressors_gammas

  #end computing sim. time
  simulation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(1, 2)) {
    cat("Total Simulation Time: ", format(simulation_time, digits = 4), " mins.\n")
  }

  results <- list(model_type = model_type,
                  simulation_time = simulation_time,
                  estimation_time_original_method = estimation_time_original_method,
                  estimation_time_covYI = estimation_time_covYI,
                  used_cores = used_cores,
                  n = n,
                  lambda = lambda,
                  tau = tau,
                  auc = auc,
                  youden_index = youden_index,
                  sensitivity = sensitivity,
                  specificity = specificity,
                  geometric_mean = geometric_mean,
                  fdr = fdr,
                  mcc = mcc,
                  corrclass = corrclass,
                  auc_covYI = auc_covYI,
                  aauc_covYI = aauc_covYI,
                  aYI_covYI = aYI_covYI,
                  youden_index_covYI = youden_index_covYI,
                  sensitivity_covYI = sensitivity_covYI,
                  specificity_covYI = specificity_covYI,
                  geometric_mean_covYI = geometric_mean_covYI,
                  fdr_covYI = fdr_covYI, mcc_covYI = mcc_covYI,
                  corrclass_covYI = corrclass_covYI,
                  n_total_var_betas = n_total_var_betas,
                  n_predicted_zeros_betas = n_predicted_zeros_betas,
                  n_predicted_non_zeros_betas = n_predicted_non_zeros_betas,
                  n_caught_betas = n_caught_betas,
                  n_non_caught_betas = n_non_caught_betas,
                  n_caught_zero_betas = n_caught_zero_betas,
                  n_zero_not_caught_betas = n_zero_not_caught_betas,
                  n_total_var_gammas = n_total_var_gammas,
                  n_predicted_zeros_gammas = n_predicted_zeros_gammas,
                  n_predicted_non_zeros_gammas = n_predicted_non_zeros_gammas,
                  n_caught_gammas = n_caught_gammas,
                  n_non_caught_gammas = n_non_caught_gammas,
                  n_caught_zero_gammas = n_caught_zero_gammas,
                  n_zero_not_caught_gammas = n_zero_not_caught_gammas,
                  betas_times_selected = betas_times_selected,
                  gammas_times_selected = gammas_times_selected,
									betas_start = betas_start,
                  roc_spec_points = roc_spec_points,
                  roc_sens_points_train = roc_sens_points_train,
                  roc_sens_points_test = roc_sens_points_test,
                  regressors_betas = regressors_betas,
                  regressors_gammas = regressors_gammas,
									input_parameters = list(
                    n = n,
                    rows_train = rows_train,
                    rows_test = rows_test,
                    cols = cols,
                    cols_cov = cols_cov,
                    max_rho = max_rho,
                    mu = mu,
                    mu_cov = mu_cov,
                    model_estimation_function = model_estimation_function,
                    model_prediction_function = model_prediction_function,
                    model_type = model_type,
                    lambda = lambda,
                    tau = tau,
                    w_g = w_g,
                    trace = trace,
                    alpha = alpha,
                    used_cores = used_cores,
                    c_function_of_covariates = c_function_of_covariates,
                    alpha_g = alpha_g,
                    penalty_g = penalty_g,
                    kernel_g = kernel_g,
                    a1_g = a1_g,
                    a2_g = a2_g,
                    trend_g = trend_g,
                    gamma_start_input = gamma_start_input,
                    gamma_start_default = gamma_start_default,
                    max_iter_g = max_iter_g,
                    delta_g = delta_g,
                    max_alpha_g = max_alpha_g,
                    stepsizeShrink_g = stepsizeShrink_g,
                    min_alpha_g = min_alpha_g,
                    convergence_error_g = convergence_error_g,
                    run_aauc = run_aauc,
                    log_file = log_file
                  ),
                  betas = betas,
                  gammas = gammas)

  return(results)
}
