#' @title pye_KS_with_print
#'
#' @description Internal helper function to print and return classification
#'   measures for the pye-KS method. This function is designed for internal
#'   use within simulation studies to provide detailed output during the
#'   testing phase.
#'
#' @param df The input data frame, excluding the ID column, containing the
#'   covariates (X) and the target variable (y).
#' @param X A character vector of column names representing the regressors.
#' @param y A character string, the column name of the target variable
#'   (must be binary: 0 or 1).
#' @param betas A numeric vector of estimated beta coefficients from the
#'   pye KS estimation.
#' @param lambda A numeric value representing the penalization parameter
#'   for betas.
#' @param c A numeric value representing the estimated cut-off point.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in pye. Sensitivity is weighted by `w` and
#'   specificity by `1 - w`. Default is 0.5 (which corresponds to the
#'   standard Youden Index).
#' @param sim_n An integer indicating the current simulation number.
#' @param n An integer indicating the total number of simulations.
#' @param alpha A numeric value for the Elastic-Net mixing parameter
#'   (0 for L1, 1 for L2).
#' @param a1 A numeric value for the SCAD penalty parameter.
#' @param a2 A numeric value for the MCP penalty parameter.
#' @param penalty A character string specifying the penalty type used
#'   ("L12", "L1", "EN", "SCAD", "MCP").
#' @param est_time A numeric value indicating the estimation time for the
#'   pye KS method.
#' @param niter An integer indicating the number of iterations taken by
#'   the optimization algorithm.
#' @param kernel A character string specifying the kernel type used for
#'   density estimation (e.g., "gaussian").
#' @param trace An integer (0, 1, or 2) controlling the verbosity of output.
#'   2 for full details, 1 for partial, 0 for none.
#'
#' @return A list containing:
#'   \item{pye_KS_result}{The results from pye_KS including classification metrics}
#'   \item{est_time}{Estimation time}
#'   \item{niter}{Number of iterations}
#'
#' @noRd
#' @keywords internal
pye_KS_with_print <- function(df, X, y,
                              betas, lambda,
                              c, w, sim_n, n,
                              alpha, a1, a2,
                              penalty, est_time,
                              niter, kernel, trace) {

  pye_KS_result <- pye_KS(df = df,
                          X = X,
                          y = y,
                          betas = betas,
                          lambda = lambda,
                          c = c,
													w = w,
                          alpha = alpha,
                          a1 = a1,
                          a2 = a2,
                          penalty = penalty,
                          prediction = TRUE,
                          kernel = kernel)

  # Validate inputs
  if (!is.data.frame(df)) stop("df must be a data frame")
  if (!is.numeric(betas)) stop("betas must be numeric")
  if (!is.numeric(lambda)) stop("lambda must be numeric")
  if (!is.numeric(c)) stop("c must be numeric")
	if (!is.numeric(w) || length(w) != 1 || w < 0 || w > 1) {stop("Parameter 'w' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(sim_n) || !is.numeric(n)) stop("sim_n and n must be numeric")
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if (!is.numeric(a1) || !is.numeric(a2)) stop("a1 and a2 must be numeric")
  if (!(penalty %in% c("L12", "L1", "EN", "SCAD", "MCP"))) stop("penalty must be one of: L12, L1, EN, SCAD, MCP")
  if (!is.numeric(est_time)) stop("est_time must be numeric")
  if (!is.numeric(niter)) stop("niter must be numeric")
  if (!(trace %in% c(0, 1, 2))) stop("trace must be 0, 1 or 2")

  if (trace %in% c(1, 2)) {
    cat("-> Results on the TEST SET \n")
    cat("-> algorithm: pye_KS_proximal_gradient_method ; ")
    cat("simulation n.", sim_n, "of", n , "; ")
    visualize_betas <- c(betas[which(betas != 0)], c)
    cat("lambda:", lambda, "; weight:", w, "; penalty:", penalty, "; pye_KS:", getElement(pye_KS_result,  paste0("pye_", penalty)), "; youden_index:", pye_KS_result$youden_index, "; sensitivity:", pye_KS_result$sensitivity, "; fdr:", pye_KS_result$fdr, "; mcc:", pye_KS_result$mcc, "; auc:", pye_KS_result$auc, "; corrclass:", pye_KS_result$corrclass, " \n")
    cat("TP:", pye_KS_result$TP, "; TN:", pye_KS_result$TN, "; FP:", pye_KS_result$FP, "; FN:", pye_KS_result$FN, ";  betas: \n")
    print(visualize_betas)
    cat("Estimation time:", est_time, "; Number of iterations:", niter, "\n\n\n")
  }

  return(pye_KS_result)
}

#' @title Simulation Study for Penalized Youden Index (pye) on Real Data
#'
#' @description This function performs a simulation study (repeated train/test
#'   splits) to estimate a penalized Youden Index (pye KS) model, and
#'   evaluates its performance on hold-out test data. It supports various
#'   penalties and optional cut-off point adjustment using covariates (covYI).
#'
#' @param n An integer, the number of simulation experiments to run.
#'   Default is 1000.
#' @param df A data frame containing the complete dataset (y, X, C).
#' @param X A character vector of column names from `df` for regressors
#'   in the pye KS model. Defaults to all columns not `y` or `C`.
#' @param y A character string, the column name for the binary target
#'   variable (0 or 1). Default is "y".
#' @param C A character vector of column names from `df` for covariates
#'   in the covYI model. Default is `NULL`.
#' @param lambda A numeric value, the penalization parameter \eqn{\lambda}
#'   for regressors \eqn{X} (pye KS model).
#' @param tau A numeric value, the penalization parameter \eqn{\tau} for
#'   covariates \eqn{C} (covYI model). Default is 0 (no penalization).
#' @param w A numeric value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in pye. Sensitivity is weighted by `w` and specificity by
#'   `1 - w`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param w_g A numeric value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param train_data_proportion A numeric value between 0 and 1. The proportion
#'   of the total observations in `df` to be used for the training set in each
#'   simulation split. The split is performed using stratified sampling on `y`
#'   to maintain class balance. Default is 0.7 (70\% training, 30\% testing).
#' @param beta_start_input A numeric vector for custom starting \eqn{\beta}
#'   coefficients. Default is `NULL`.
#' @param beta_start_default A character string for default starting \eqn{\beta}s
#'   if `beta_start_input` is `NULL`. Options: "zeros", "corr".
#'   Default is "zeros".
#' @param trace An integer (0, 1, or 2) controlling output verbosity.
#'   0=none, 1=main, 2=all steps. Default is 1.
#' @param alpha A numeric value [0, 1], the elastic-net mixing parameter
#'   for pye KS. Default is 0.5.
#' @param a1 A numeric value, the $a$ parameter for SCAD/MCP in pye KS.
#'   Default is 3.7.
#' @param a2 A numeric value, the $b$ parameter for MCP in pye KS.
#'   Default is 3.0.
#' @param penalty \code{character}. The penalty type for regressors (\code{X})
#'   in \code{pye}. Options: \code{"L12"}, \code{"L1"} (Lasso),
#'   \code{"EN"} (Elastic-Net), \code{"SCAD"}, and \code{"MCP"}. Default is \code{"L1"}.
#' @param regressors_betas numeric vector. The "true" \eqn{\beta} coefficients,
#'   if known, for variable selection accuracy. Default is `NULL`.
#' @param used_cores An integer specifying CPU cores for parallelization.
#'   0 uses \eqn{70\%} of available cores. Default is 1 (no parallel).
#' @param scaling A logical value. If `TRUE`, the dataset `df` is scaled
#'   before processing. Default is `FALSE`.
#' @param c_zero_fixed A logical value. If `TRUE`, the cut-off point \code{c}
#'   is fixed at zero. Default is `FALSE`.
#' @param max_iter An integer, max iterations for pye KS optimization.
#'   Default is 10000.
#' @param trend A character string specifying the pye KS optimization
#'   algorithm trend. Options: "monotone" (mmAPG), "nonmonotone" (mnmAPG).
#'   Default is "monotone".
#' @param delta A numeric value, convergence tolerance for pye KS.
#'   Default is 1e-5.
#' @param max_alpha A numeric value, max step-size for pye KS
#'   backtracking. Default is 10000.
#' @param stepsizeShrink A numeric value [0, 1], shrinkage factor for
#'   pye KS backtracking step-size. Default is 0.8.
#' @param min_alpha A numeric value, min step-size for pye KS
#'   backtracking. Default is 1e-10.
#' @param convergence_error A numeric value, the convergence threshold
#'   in pye KS. Default is 1e-7.
#' @param kernel A character string specifying the kernel type for pye KS
#'   density estimation (e.g., "gaussian"). Default is "gaussian".
#' @param c_function_of_covariates A logical value. If `TRUE`, covYI is
#'   used to estimate cut-off \code{c} as a function of \code{C}. Default is `FALSE`.
#' @param alpha_g A numeric value [0, 1], the elastic-net mixing parameter
#'   for covYI. Default is 0.5.
#' @param penalty_g A character string specifying the penalty type for
#'   covYI. Must be "L12", "L1", "EN", "SCAD", or "MCP". Default is "L1".
#' @param kernel_g A character string specifying the kernel type for covYI
#'   density estimation. Default is "gaussian".
#' @param a1_g A numeric value, the $a$ parameter for SCAD/MCP in covYI.
#'   Default is 3.7.
#' @param a2_g A numeric value, the $b$ parameter for MCP in covYI.
#'   Default is 3.0.
#' @param trend_g A character string specifying the covYI optimization
#'   algorithm trend. Default is "monotone".
#' @param gamma_start_input A numeric vector for custom starting \eqn{\gamma}
#'   coefficients. Default is `NULL`.
#' @param gamma_start_default A character string for default starting \eqn{\gamma}s
#'   if `gamma_start_input` is `NULL`. Options: "zeros", "corr".
#'   Default is "zeros".
#' @param regressors_gammas A numeric vector representing the true \eqn{\gamma}
#'   coefficients (if known), for evaluation. Default is `NULL`.
#' @param max_iter_g An integer, max iterations for covYI optimization.
#'   Default is 10000.
#' @param delta_g A numeric value, convergence tolerance for covYI.
#'   Default is 1e-5.
#' @param max_alpha_g A numeric value, max step-size for covYI
#'   backtracking. Default is 10000.
#' @param stepsizeShrink_g A numeric value [0, 1], shrinkage factor for
#'   covYI backtracking step-size. Default is 0.8.
#' @param min_alpha_g A numeric value, min step-size for covYI
#'   backtracking. Default is 1e-12.
#' @param convergence_error_g A numeric value, the convergence threshold
#'   in covYI. Default is 1e-7.
#' @param run_aauc A logical value. If `FALSE`, aAUC and aYI measures
#'   are not computed to save time. Default is `FALSE`.
#' @param log_file Character. Path to a file for logging parallel
#'   output. Default is "log_sim_pye_real.txt".
#'
#' @return A \code{list} containing the aggregated results of the simulation
#'   study, including:
#' \item{simulation_time}{Total time taken for the entire simulation (time unit
#'   is in the output).}
#' \item{estimation_time_original_method}{Mean estimation time for the pye KS
#'   model across all simulations.}
#' \item{estimation_time_covYI}{Mean estimation time for the covYI model, if
#'   \code{c_function_of_covariates} is \code{TRUE}.}
#' \item{used_cores}{Number of CPU cores used for parallel processing.}
#' \item{n}{Number of simulation experiments performed.}
#' \item{lambda, tau}{The penalization parameters used for \code{X} and \code{C}
#'   respectively.}
#' \item{pye_L12, pye_L1, pye_EN, pye_SCAD, pye_MCP}{Matrices (n x 2, for
#'   train/test) containing the pye KS objective function values, based on the
#'   selected \code{penalty} and the fitted betas.}
#' \item{auc, youden_index, sensitivity, specificity, geometric_mean, fdr,
#'   mcc, corrclass}{Matrices (n x 2, for train/test) of classification
#'   performance measures for the pye KS model.}
#' \item{auc_covYI, aauc_covYI, aYI_covYI, youden_index_covYI,
#'   sensitivity_covYI, specificity_covYI, geometric_mean_covYI, fdr_covYI,
#'   mcc_covYI, corrclass_covYI}{Matrices (n x 2, for train/test) of
#'   classification performance measures for the covYI model, if
#'   \code{c_function_of_covariates} is \code{TRUE}.}
#' \item{n_total_var_betas, n_predicted_zeros_betas, n_predicted_non_zeros_betas,
#'   n_caught_betas, n_non_caught_betas, n_caught_zero_betas,
#'   n_zero_not_caught_betas}{Matrices related to variable selection performance
#'   (comparison with \code{regressors_betas}), where relevant.}
#' \item{n_total_var_gammas, n_predicted_zeros_gammas,
#'   n_predicted_non_zeros_gammas, n_caught_gammas, n_non_caught_gammas,
#'   n_caught_zero_gammas, n_zero_not_caught_gammas}{Matrices related to
#'   variable selection performance for gammas (comparison with
#'   \code{regressors_gammas}), if \code{c_function_of_covariates} is
#'   \code{TRUE}.}
#' \item{betas_times_selected}{A vector counting the number of times (out of
#'   \code{n}) each \code{X} regressor was estimated as non-zero.}
#' \item{gammas_times_selected}{A vector counting the number of times (out of
#'   \code{n}) each \code{C} covariate was estimated as non-zero, if
#'   \code{c_function_of_covariates} is \code{TRUE}.}
#' \item{betas_start}{The starting beta coefficients used.}
#' \item{c_zero_fixed}{A logical value indicating if the cut-off point \code{c}
#'   was fixed at zero.}
#' \item{roc_spec_points}{A vector of specificity points (e.g.,
#'   \code{seq(0, 1, by = 0.05)}) used to compute the ROC curves.}
#' \item{roc_sens_points_train}{A \code{list} where each element contains the
#'   sensitivity points corresponding to \code{roc_spec_points} on the training
#'   data for each simulation.}
#' \item{roc_sens_points_test}{A \code{list} where each element contains the
#'   sensitivity points corresponding to \code{roc_spec_points} on the testing
#'   data for each simulation.}
#' \item{regressors_betas, regressors_gammas}{The true regressor vectors used
#'   for evaluation (if provided).}
#' \item{input_parameters}{\code{character vector}. A list containing the
#'   input parameters.}
#' \item{c_function_of_covariates, run_aauc}{The input parameters indicating if
#'   \code{covYI} was used and if advanced AUC metrics were run.}
#' \item{betas}{A \code{list} containing the vector of estimated beta
#'   coefficients from each simulation run.}
#' \item{gammas}{A \code{list} containing the vector of estimated gamma
#'   coefficients from each simulation run, if \code{c_function_of_covariates}
#'   is \code{TRUE}.}
#'
#' @examples
#' library(pye)
#'
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_betas <- sim_data$nregressors
#' regressors_gammas <- sim_data$ncovariates
#'
#' # Run a small simulation study
#' results_study <- pye_KS_simulation_study(
#'   n = 2, # Small number of simulations for example
#'   df = df,
#'   X = X,
#'   y = y,
#'   C = C,
#'   lambda = 1.3,
#'   tau = 0.8,
#'   beta_start_default = "zeros",
#'   gamma_start_default = "zeros",
#'   trace = 1,
#'   alpha = 0.5,
#'   alpha_g = 0.5,
#'   penalty = "L1",
#'   penalty_g = "L1",
#'   kernel = "gaussian",
#'   used_cores = 1, # Use 1 core for example
#'   c_function_of_covariates = TRUE,
#'   c_zero_fixed = FALSE,
#'   run_aauc = FALSE,
#'   max_iter = 5, # Reduced iterations for example
#'   max_iter_g = 5 # Reduced iterations for example
#' )
#'
#' # Print some of the results
#' cat("Simulation Time: ", format(results_study$simulation_time, digits = 4))
#' cat("CCR:\n")
#' print(results_study$corrclass)
#' cat("Betas Times Selected:\n")
#' print(results_study$betas_times_selected[results_study$betas_times_selected > 0])
#' if (results_study$input_parameters$c_function_of_covariates) {
#'   cat("CCR with covYI:\n")
#'   print(results_study$corrclass_covYI)
#'   cat("\nGammas Times Selected:\n")
#'   print(results_study$gammas_times_selected[results_study$gammas_times_selected > 0])
#' }
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom pROC roc coords
#' @importFrom tools file_path_sans_ext file_ext

#' @export
pye_KS_simulation_study <- function(n = 1000, df, X = NULL, y = "y", C = NULL,
                                           lambda, tau = 0, w = 0.5, w_g = 0.5,
																					 train_data_proportion = 0.7,
                                           beta_start_input = NULL,
                                           beta_start_default = "zeros", trace = 1,
                                           alpha = 0.5, a1 = 3.7, a2 = 3, penalty = "L1",
                                           regressors_betas = NULL,
                                           used_cores = 1, scaling = FALSE,
                                           c_zero_fixed = FALSE,
                                           max_iter = 10000, trend = "monotone", delta = 1e-5,
                                           max_alpha = 10000, stepsizeShrink = 0.8,
                                           min_alpha = 1e-10, convergence_error = 1e-7, kernel = "gaussian",
                                           c_function_of_covariates = FALSE,
                                           alpha_g = 0.5, penalty_g = "L1",
                                           kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                                           trend_g = "monotone", gamma_start_input = NULL,
                                           gamma_start_default = "zeros",
                                           regressors_gammas = NULL, max_iter_g = 10000,
                                           delta_g = 1e-5, max_alpha_g = 10000,
                                           stepsizeShrink_g = 0.8, min_alpha_g = 1e-12,
                                           convergence_error_g = 1e-7,
                                           run_aauc = FALSE, log_file = "log_sim_pye_real.txt") {

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

  # Check if tau exists when c_function_of_covariates = TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
    if (is.null(C) || length(C) == 0) { stop("Parameter 'C' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[, c(y, X, C), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Further input validation
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
  if (!(penalty_g %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty_g'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
	if (!is.numeric(w) || length(w) != 1 || w < 0 || w > 1) {stop("Parameter 'w' must be a single numeric value between 0 and 1.")}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 2 || !is.integer(as.integer(max_iter))) {stop("Parameter 'max_iter' needs to be an integer and at least 2.")}
  if (!is.numeric(max_iter_g) || length(max_iter_g) != 1 || max_iter_g < 2 || !is.integer(as.integer(max_iter_g))) {stop("Parameter 'max_iter_g' needs to be an integer and at least 2.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter 'trace' has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  if (!(trend_g %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend_g' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  if (!is.numeric(n) || length(n) != 1 || n < 2 || n != floor(n)) {stop("Parameter 'n' must be a single integer value and at least 2.")}
  if (!is.logical(scaling)) {stop("Parameter 'scaling' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_zero_fixed)) {stop("Parameter 'c_zero_fixed' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_function_of_covariates)) {stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel %in% valid_kernels)) {stop(paste0("Parameter 'kernel' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  if (!(kernel_g %in% valid_kernels)) {stop(paste0("Parameter 'kernel_g' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}

  #standardize df1
  if (scaling == TRUE) {
    df1 <- scaling_df_for_pye (df = df1, X = colnames(df1[, names(df1) %in% c(X, C)]), y = "y")$df_scaled
  }

  # Generate seeds
  seeds <- 1:n

  # Initialize matrices to store results
  names <- paste("seed", seeds, sep = "=")
  pye_L12 <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_L1 <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_EN <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_SCAD <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_MCP <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
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
  n_predicted_zeros_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_zeros_betas"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = n, ncol = 1, dimnames = list(names, "n_predicted_non_zeros_betas"))
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
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(0, nrow = length(X), ncol = 1, dimnames = list(X, "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(0, nrow = length(C), ncol = 1, dimnames = list(C, "n_times_gamma_diff_zero"))

  #function to be computed
  func <- function(seed, n, df, X, y, C, lambda, tau,
	                 w, w_g,
									 train_data_proportion,
                   beta_start_input,
                   beta_start_default,
                   gamma_start_input,
                   gamma_start_default,
                   trace, alpha,
                   alpha_g, penalty_g,
                   a1, a2, a1_g, a2_g,
                   penalty,
                   regressors_betas,
                   regressors_gammas,
                   trend, trend_g,
                   kernel,
                   c_zero_fixed,
                   c_function_of_covariates,
                   run_aauc,
                   max_iter, max_iter_g,
                   delta, delta_g,
                   max_alpha, max_alpha_g,
                   stepsizeShrink,
                   stepsizeShrink_g,
                   min_alpha,
                   min_alpha_g,
                   convergence_error,
                   convergence_error_g) {

    if (trace %in% c(1, 2)) {
      cat("\n--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda = ", lambda, "; weight:", w, "; penalty = ", penalty, "\n")
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

    # Train pye KS
    # Remove 'ID' column before passing to pye_KS_estimation
    train_solution <- pye_KS_estimation(df = train_df[, names(train_df) != "ID", drop = FALSE],
                                        X = X, y = y,
                                        lambda = lambda,
																				w = w,
                                        beta_start_input = beta_start_input,
                                        beta_start_default = beta_start_default,
                                        trace = trace,
                                        alpha = alpha,
                                        a1 = a1, a2 = a2,
                                        max_iter = max_iter,
                                        penalty = penalty,
                                        regressors_betas = regressors_betas,
                                        trend = trend,
                                        stepsizeShrink = stepsizeShrink,
                                        delta = delta,
                                        max_alpha = max_alpha,
                                        min_alpha = min_alpha,
                                        convergence_error = convergence_error,
                                        kernel = kernel,
                                        c_zero_fixed = c_zero_fixed)

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
                                                  delta = delta_g,
                                                  max_alpha = max_alpha_g,
                                                  kernel = kernel_g,
                                                  run_aauc = run_aauc)

      estimation_time_covYI <- train_covYI_solution$estimation_time
			niter_covYI <- train_covYI_solution$niter
      z_hat <- train_covYI_solution$z_hat

    }

    # ROC curve @ certain levels on TRAIN data
    roc_spec_points <- seq(0, 1, by = 0.05)
    est_roc <- pROC::roc(as.numeric(train_df[[y]]), z_hat[, "z_hat"], levels = c(0, 1), direction = "<", quiet = TRUE)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)

    # Estimated betas and c
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    c <- getElement(train_solution, "c_hat")
    niter <- train_solution$niter

    # --- Model Prediction on TEST data ---
    test_solution <- pye_KS_with_print(df = test_df[, names(test_df) != "ID", drop = FALSE],
                                       X = X, y = y, betas = betas,
                                       lambda = lambda, c = c,
																			 w = w,
                                       sim_n = seed, n = n,
                                       alpha = alpha,
                                       a1 = a1, a2 = a2,
                                       penalty = penalty,
                                       est_time = estimation_time_original_method,
                                       niter = niter,
                                       kernel = kernel,
                                       trace = trace)

    z_hat <- test_solution$z_hat

    test_covYI_solution <- NULL

    # --- covYI Prediction on TEST data (if c_function_of_covariates is TRUE) ---
    if (c_function_of_covariates == TRUE) {
      #put "const" in C
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df = cbind(test_df[, names(test_df) != "ID", drop = FALSE],
                                      z_hat = test_solution$z_hat[, "z_hat"]),
                                      z = "z_hat", y = y, C = C1,
																			w = w_g,
                                      gammas = train_covYI_solution$gammas_hat,
                                      tau = tau, kernel = kernel_g,
                                      alpha = alpha_g, a1 = a1_g, a2 = a2_g,
                                      penalty = penalty_g,
                                      prediction = TRUE,
                                      run_aauc = run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1, 2)) {
        cat("-> Results on the TEST SET\n")
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
    parallel::clusterExport(cl, c("n", "df1", "X", "y", "C", "lambda", "tau", "w", "w_g",
                                  "beta_start_input", "beta_start_default", "gamma_start_input", "gamma_start_default",
                                  "trace", "a1", "a2", "a1_g", "a2_g",
                                  "alpha", "alpha_g", "penalty", "penalty_g", "trend", "kernel", "c_zero_fixed",
                                  "regressors_betas", "regressors_gammas", "c_function_of_covariates", "run_aauc",
                                  "trend_g", "max_iter", "max_iter_g", "delta", "delta_g", "max_alpha", "max_alpha_g",
                                  "stepsizeShrink", "stepsizeShrink_g", "min_alpha", "min_alpha_g",
                                  "convergence_error", "convergence_error_g", "func", "pye_KS_with_print"
                                   ), envir = environment())

    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed = x, n = n, df = df1, X = X, y = y, C = C,
                                                                        lambda = lambda, tau = tau,
                                                                        train_data_proportion = train_data_proportion,
																																				w = w, w_g = w_g,
                                                                        beta_start_input = beta_start_input,
                                                                        beta_start_default = beta_start_default,
                                                                        gamma_start_input = gamma_start_input,
                                                                        gamma_start_default = gamma_start_default,
                                                                        c_zero_fixed = c_zero_fixed,
                                                                        trace = trace, alpha = alpha,
                                                                        alpha_g = alpha_g,
                                                                        a1 = a1, a2 = a2,
                                                                        a1_g = a1_g, a2_g = a2_g,
                                                                        penalty = penalty,
                                                                        penalty_g = penalty_g,
                                                                        regressors_betas = regressors_betas,
                                                                        regressors_gammas = regressors_gammas,
                                                                        trend = trend, trend_g = trend_g,
                                                                        kernel = kernel,
                                                                        c_function_of_covariates = c_function_of_covariates,
                                                                        run_aauc = run_aauc,
                                                                        max_iter = max_iter,
                                                                        max_iter_g = max_iter_g,
                                                                        delta = delta,
                                                                        delta_g = delta_g,
                                                                        max_alpha = max_alpha,
                                                                        max_alpha_g = max_alpha_g,
                                                                        stepsizeShrink = stepsizeShrink,
                                                                        stepsizeShrink_g = stepsizeShrink_g,
                                                                        min_alpha = min_alpha,
                                                                        min_alpha_g = min_alpha_g,
                                                                        convergence_error = convergence_error,
                                                                        convergence_error_g = convergence_error_g))
  } else {
    # Sequential execution
    cat("Running simulation in sequential mode (used_cores = 1).\n")
    simulation_study <- lapply(seeds, function(x) func(seed = x, n = n, df = df1, X = X, y = y, C = C,
                                                       lambda = lambda, tau = tau,
                                                       train_data_proportion = train_data_proportion,
																											 w = w, w_g = w_g,
                                                       beta_start_input = beta_start_input,
                                                       beta_start_default = beta_start_default,
                                                       gamma_start_input = gamma_start_input,
                                                       gamma_start_default = gamma_start_default,
                                                       c_zero_fixed = c_zero_fixed,
                                                       trace = trace, alpha = alpha,
                                                       alpha_g = alpha_g,
                                                       a1 = a1, a2 = a2,
                                                       a1_g = a1_g, a2_g = a2_g,
                                                       penalty = penalty,
                                                       penalty_g = penalty_g,
                                                       regressors_betas = regressors_betas,
                                                       regressors_gammas = regressors_gammas,
                                                       trend = trend, trend_g = trend_g,
                                                       kernel = kernel,
                                                       c_function_of_covariates = c_function_of_covariates,
                                                       run_aauc = run_aauc, max_iter = max_iter,
                                                       max_iter_g = max_iter_g, delta = delta,
                                                       delta_g = delta_g, max_alpha = max_alpha,
                                                       max_alpha_g = max_alpha_g,
                                                       stepsizeShrink = stepsizeShrink,
                                                       stepsizeShrink_g = stepsizeShrink_g,
                                                       min_alpha = min_alpha, min_alpha_g = min_alpha_g,
                                                       convergence_error = convergence_error,
                                                       convergence_error_g = convergence_error_g))
  }

  estimation_time_original_method <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  #fill the matrices
  #measures on the train set
  temp_pye <- get(paste0("pye_", penalty))
  temp_pye[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, paste0("pye_KS_", penalty))))

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
  betas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty))) != 0)

  # Measures on test set
  temp_pye[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, paste0("pye_", penalty))))
  assign(paste0("pye_", penalty), temp_pye)

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
  betas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))
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

  results <- list(simulation_time = simulation_time,
                  estimation_time_original_method = estimation_time_original_method,
                  estimation_time_covYI = estimation_time_covYI,
                  used_cores = used_cores, n = n, lambda = lambda, tau = tau,
                  pye_L12 = pye_L12, pye_L1 = pye_L1, pye_EN = pye_EN, pye_SCAD = pye_SCAD, pye_MCP = pye_MCP,
                  auc = auc, youden_index = youden_index, sensitivity = sensitivity,
                  specificity = specificity, geometric_mean = geometric_mean,
                  fdr = fdr, mcc = mcc, corrclass = corrclass,
                  auc_covYI = auc_covYI, aauc_covYI = aauc_covYI,
									aYI_covYI = aYI_covYI, youden_index_covYI = youden_index_covYI,
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
									c_zero_fixed = c_zero_fixed,
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
                    lambda = lambda,
                    tau = tau,
                    w = w,
                    w_g = w_g,
                    beta_start_input = beta_start_input,
                    beta_start_default = beta_start_default,
                    trace = trace,
                    alpha = alpha,
                    a1 = a1,
                    a2 = a2,
                    penalty = penalty,
                    regressors_betas = regressors_betas,
                    used_cores = used_cores,
                    scaling = scaling,
                    c_zero_fixed = c_zero_fixed,
                    max_iter = max_iter,
                    trend = trend,
                    delta = delta,
                    max_alpha = max_alpha,
                    stepsizeShrink = stepsizeShrink,
                    min_alpha = min_alpha,
                    convergence_error = convergence_error,
                    kernel = kernel,
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
									betas = betas,
									gammas = gammas)

  return(results)
}
















#' @title Simulation Study for Penalized Youden Index (pye) on Synthetic Data
#'
#' @description This function conducts a simulation study for the \link{pye_KS}
#' (Penalized Youden Index Estimation) method using synthetic data.
#' It generates multiple synthetic datasets, splits each into training
#' and testing sets, estimates pye models on the training data
#' (optionally incorporating covariates via \link{covYI_KS}), and
#' evaluates classification performance measures on the test sets.
#' Parallel computation is supported.
#'
#' @param n \code{integer}. The number of simulation experiments to run. Must
#'   be a single integer \eqn{\ge 2}. Default is 1000.
#' @param rows_train \code{integer}. The number of observations for the
#'   training set in each simulation. Must be a single integer \eqn{\ge 1}.
#'   Default is \eqn{50}.
#' @param rows_test \code{integer}. The number of observations for the
#'   test set. Recommended to be much larger than \code{rows_train} to
#'   adequately evaluate model performance. Must be a single integer \eqn{\ge 1}.
#'   Default is \eqn{1000}.
#' @param cols \code{integer}. The number of regressor variables (\code{X} features)
#'   for both samples. Must be a single integer \eqn{\ge 1}. Default is 2000.
#' @param cols_cov \code{integer}. The number of covariate variables (\code{C} features)
#'   for both samples. Must be a single integer \eqn{\ge 0}. Default is \eqn{20}.
#' @param max_rho \code{numeric}. The maximum correlation coefficient used
#'   in the latent factor model for generating variables. Higher values lead
#'   to more correlated features. Must be between 0 and 1. Default is 0.5.
#' @param mu \code{numeric} vector. The mean vector for the multivariate normal
#'   distribution of the regressors (\code{X}). If \code{NULL}, defaults to a
#'   vector of zeros with length \code{cols}.
#' @param mu_cov \code{numeric} vector. The mean vector for the covariates (\code{C}).
#'   If \code{NULL}, defaults to a vector of zeros with length \code{cols_cov}.
#' @param lambda \code{numeric}. The penalization parameter \eqn{\lambda} for the
#'   regressors (\eqn{X}) in the \code{pye} estimation. Must be a single non-negative
#'   numeric value.
#' @param tau \code{numeric}. The penalization parameter \eqn{\tau} for the covariates
#'   (\code{C}) when using \code{covYI}. Must be a single non-negative numeric
#'   value. Ignored if \code{c_function_of_covariates} is \code{FALSE}
#'   Default is \eqn{0}.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in pye. Sensitivity is weighted by `w` and specificity by
#'   `1 - w`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param trace \code{integer}. Controls the level of output messages during the
#'   simulation. \eqn{0}: No output. 1: Partial output (progress and key results).
#'   2: Full output (detailed progress for each iteration). Default is \eqn{1}.
#' @param used_cores \code{integer}. The number of CPU cores for parallel
#'   computation. If 1, execution is sequential. If \eqn{\le 0}, it attempts
#'   to detect available cores and use \eqn{70\%} of physical cores. Default is \eqn{1}.
#' @param c_function_of_covariates \code{logical}. If \code{TRUE}, the cut-off
#'   point \code{c} is estimated as a function of the covariates using \code{\link{covYI_KS}}.
#'   If \code{FALSE}, covariates are ignored. Default is \code{FALSE}.
#' @param c_zero_fixed \code{logical}. If \code{TRUE}, the cut-off point \code{c} is
#'   fixed at zero, simplifying the optimization. Default is \code{FALSE}.
#' @param beta_start_input \code{numeric} vector. Optional initial starting
#'   point for the \eqn{\beta} coefficients. If \code{NULL}, \code{beta_start_default}
#'   is used. Default is \code{NULL}.
#' @param beta_start_default \code{character}. Specifies the default starting point
#'   for \eqn{\beta} coefficients if \code{beta_start_input} is \code{NULL}.
#'   \code{"zeros"}: Starts with a vector of all zeros.
#'   \code{"corr"}: Starts with values based on the correlation of each regressor
#'   with the target variable \code{y}. Default is \code{"zeros"}.
#' @param max_iter \code{integer}. The maximum number of iterations for the
#'   \code{pye} optimization algorithms (mmAPG/mnmAPG). Must be an integer \eqn{\ge 2}.
#'   Default is 10000.
#' @param trend \code{character}. Specifies the type of optimization algorithm
#'   for \code{pye}. \code{"monotone"} (mmAPG) or \code{"nonmonotone"} (mnmAPG).
#'   Default is \code{"monotone"}.
#' @param delta \code{numeric}. Parameter for the convergence condition of the
#'   \code{pye} optimization algorithm. Must be a single positive numeric value.
#'   Default is \code{1e-5}.
#' @param max_alpha \code{numeric}. The maximum value for the step-size \eqn{\alpha}
#'   in the backtracking line-search within \code{pye} optimization.
#'   Must be a single positive numeric value. Default is \eqn{10000}.
#' @param stepsizeShrink \code{numeric}. Parameter (\eqn{0 < \cdot < 1}) to adjust the
#'   step-size in the backtracking line-search for \code{pye}. Default is 0.8.
#' @param min_alpha \code{numeric}. The minimum value for the step-size \eqn{\alpha}
#'   in the backtracking line-search within \code{pye}. Must be a single positive
#'   numeric value. Default is \code{1e-10}.
#' @param convergence_error \code{numeric}. The tolerance for the \code{pye}
#'   optimization convergence check. Must be a single positive numeric value.
#'   Default is \code{1e-7}.
#' @param alpha \code{numeric}. The \eqn{\alpha} mixing parameter (\eqn{0 \le \alpha \le 1})
#'   for the Elastic-Net penalty (if \code{penalty} is "EN") in \code{pye}.
#'   Default is \eqn{0.5}.
#' @param alpha_g \code{numeric}. The \eqn{\alpha} mixing parameter (\eqn{0 \le \alpha \le 1})
#'   for the Elastic-Net penalty (if \code{penalty_g} is "EN") in \code{covYI}.
#'   Default is \eqn{0.5}.
#' @param a1 \code{numeric}. Parameter for the SCAD and MCP penalties in \code{pye}.
#'   Must be a single non-negative numeric value. Default is 3.7.
#' @param a2 \code{numeric}. Parameter for the MCP penalty in \code{pye}. Must
#'   be a single non-negative numeric value. Default is 3.0.
#' @param kernel \code{character}. The kernel type for density estimation in \code{pye}.
#'   \code{"gaussian"} is recommended. Other options: \code{"normal"}, \code{"uniform"},
#'   \code{"epanechnikov"}, etc. Default is \code{"gaussian"}.
#' @param penalty \code{character}. The penalty type for regressors (\code{X}) in \code{pye}.
#'   Options: \code{"L12"}, \code{"L1"} (Lasso), \code{"EN"} (Elastic-Net),
#'   \code{"SCAD"}, and \code{"MCP"}. Default is \code{"L1"}.
#' @param penalty_g \code{character}. The penalty type for covariates (\code{C}) in
#'   \code{covYI} estimation. Same options as \code{penalty}. Default is \code{"L1"}.
#' @param kernel_g \code{character}. The kernel type for density estimation in \code{covYI}.
#'   Same options as \code{kernel}. Default is \code{"gaussian"}.
#' @param a1_g \code{numeric}. Parameter for the SCAD and MCP penalties in \code{covYI}.
#'   Must be a single non-negative numeric value. Default is 3.7.
#' @param a2_g \code{numeric}. Parameter for the MCP penalty in \code{covYI}.
#'   Must be a single non-negative numeric value. Default is 3.0.
#' @param trend_g \code{character}. Specifies the optimization algorithm for \code{covYI}.
#'   \code{"monotone"} (mmAPG) or \code{"nonmonotone"} (mnmAPG).
#'   Default is \code{"monotone"}.
#' @param gamma_start_input \code{numeric} vector. Optional initial starting
#'   point for the \eqn{\gamma} coefficients in \code{covYI}. If \code{NULL},
#'   \code{gamma_start_default} is used. Default is \code{NULL}.
#' @param gamma_start_default \code{character}. Specifies the default starting point
#'   for \eqn{\gamma} coefficients if \code{gamma_start_input} is \code{NULL}.
#'   \code{"zeros"} or \code{"corr"} (correlation with \code{y}). Default is \code{"zeros"}.
#' @param max_iter_g \code{integer}. The maximum number of iterations for the
#'   optimization algorithms used in \code{covYI}. Must be an integer \eqn{\ge 2}.
#'   Default is 10000.
#' @param delta_g \code{numeric}. Parameter for the convergence condition of the
#'   \code{covYI} optimization algorithm. Must be a single positive numeric value.
#'   Default is \code{1e-5}.
#' @param max_alpha_g \code{numeric}. The maximum step-size \eqn{\alpha} in the
#'   backtracking line-search within \code{covYI}. Default is \eqn{10000}.
#' @param stepsizeShrink_g \code{numeric}. Parameter (\eqn{0 < \cdot < 1}) to adjust the
#'   step-size in the backtracking line-search for \code{covYI}. Default is \eqn{0.8}.
#' @param min_alpha_g \code{numeric}. The minimum step-size \eqn{\alpha} in the
#'   backtracking line-search within \code{covYI}. Default is \code{1e-12}.
#' @param convergence_error_g \code{numeric}. The tolerance for the \code{covYI}
#'   optimization convergence check. Default is \code{1e-7}.
#' @param run_aauc \code{logical}. If \code{FALSE}, the aAUC and aYI (adjusted
#'   measures) for \code{covYI} are not computed, saving estimation time.
#'   Default is \code{FALSE}.
#' @param log_file \code{character}. Path to a file for logging output from parallel
#'   workers. Default is \code{"log_sim_pye_synthetic.txt"}.
#'
#' @return A \code{list} containing the aggregated results from the simulation study
#'   across all \code{n} experiments. Key components include:
#'   \item{simulation_time}{Total time taken for the entire simulation study.}
#'   \item{estimation_time_original_method}{Average estimation time for the
#'     \code{pye} method.}
#'   \item{estimation_time_covYI}{Average estimation time for the \code{covYI}
#'     method (if covariates are enabled).}
#'   \item{used_cores}{Number of cores utilized for parallel processing.}
#'   \item{n}{Input parameter \code{n} of the simulation.}
#'   \item{lambda}{Input parameter \code{lambda} of the simulation.}
#'   \item{tau}{Input parameter \code{tau} of the simulation.}
#'   \item{pye_L12, pye_L1, pye_EN, pye_SCAD, pye_MCP}{Matrices (rows = seeds, 2
#'     columns for train/test) of pye values for each penalty type. Only the matrix
#'     corresponding to the \code{penalty} parameter will contain data.}
#'   \item{auc, youden_index, sensitivity, specificity, geometric_mean, fdr,
#'     mcc, corrclass}{Matrices (rows = seeds, 2 columns for train/test) of
#'     performance measures for the \code{pye} method.}
#'   \item{auc_covYI, aauc_covYI, aYI_covYI, youden_index_covYI, sensitivity_covYI,
#'     specificity_covYI, geometric_mean_covYI, fdr_covYI, mcc_covYI,
#'     corrclass_covYI}{Matrices (rows = seeds, 2 columns for train/test) of
#'     performance measures for the \code{covYI} method (if covariates are enabled).}
#'   \item{n_total_var_betas, n_predicted_zeros_betas, ..., n_zero_not_caught_betas}{
#'     Matrices providing counts related to \eqn{\beta} coefficient selection success.}
#'   \item{n_total_var_gammas, n_predicted_zeros_gammas, ..., n_zero_not_caught_gammas}{
#'     Matrices providing counts related to \eqn{\gamma} coefficient selection success
#'     (if covariates are enabled).}
#'   \item{betas_times_selected}{Vector indicating the frequency each \eqn{\beta}
#'     regressor was selected across all simulations.}
#'   \item{gammas_times_selected}{Vector indicating the frequency each \eqn{\gamma}
#'     regressor was selected (if covariates are enabled).}
#'   \item{roc_spec_points}{Vector of specificity values at which ROC sensitivity
#'     points are recorded.}
#'   \item{roc_sens_points_train, roc_sens_points_test}{Lists of sensitivity
#'     values for the ROC curve on the training/test set for each simulation.}
#'   \item{regressors_betas}{The indices of the true non-zero \eqn{\beta} coefficients.}
#'   \item{regressors_gammas}{The indices of the true non-zero \eqn{\gamma} coefficients.}
#'   \item{input_parameters}{\code{character vector}. A list containing the
#'     input parameters.}
#'   \item{c_function_of_covariates, run_aauc}{Logical flags indicating feature use.}
#'   \item{betas, gammas}{Lists of all estimated coefficient vectors from each
#'     simulation (\eqn{\gamma} only if covariates are enabled).}
#'
#' @examples
#' # A small-scale example for demonstration purposes.
#' # For meaningful results, 'n' should be much larger (e.g., 500).
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
#' lambda_to_use <- 1
#' tau_to_use <- 0.4
#' c_function_of_covariates <- TRUE
#'
#' sim_result <- pye_KS_simulation_study_synthetic_data(
#'   n = n,
#'   rows_train = 200,
#'   rows_test = 1000,
#'   cols = cols,
#'   cols_cov = cols_cov,
#'   max_rho = 0.2,
#'   mu = mu ,
#'   mu_cov = mu_cov,
#'   lambda = lambda_to_use,
#'   tau = tau_to_use,
#'   trend = "monotone",
#'   beta_start_default = "zeros",
#'   gamma_start_default = "zeros",
#'   trace = 1,
#'   a1 = 3.7, a2 = 3, # SCAD/MCP param for pye
#'   a1_g = 3.7, a2_g = 3, # SCAD/MCP param for covYI
#'   penalty = "L1", # Options: "L12", "L1", "EN", "SCAD", "MCP"
#'   penalty_g = "SCAD", # Options: "L12", "L1", "EN", "SCAD", "MCP"
#'   kernel = "gaussian",
#'   used_cores = 1,
#'   c_function_of_covariates = c_function_of_covariates,
#'   run_aauc = TRUE,
#'   max_iter = 10, # Reduced for a quick example
#'   max_iter_g = 10 # Reduced for a quick example
#' )
#'
#' # You can now access the results, e.g.:
#' print(sim_result$auc)
#' print(sim_result$betas_times_selected[sim_result$betas_times_selected > 1])
#' print(sim_result$aauc)
#' print(sim_result$gammas_times_selected[sim_result$gammas_times_selected > 1])
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom pROC roc coords
#' @importFrom tools file_path_sans_ext file_ext

#' @noRd
#' @keywords internal
pye_KS_simulation_study_synthetic_data <- function(n = 1000, rows_train = 50, rows_test = 1000, cols = 2000, cols_cov = 20,
                                                 max_rho = 0.5, mu = rep(0, cols), mu_cov = rep(0, cols_cov),
                                                 lambda, tau = 0, w = 0.5, w_g = 0.5, trace = 1, used_cores = 1,
                                                 c_function_of_covariates = FALSE, c_zero_fixed = FALSE,
                                                 beta_start_input = NULL, beta_start_default = "zeros",
                                                 max_iter = 10000, trend = "monotone", delta = 1e-5, max_alpha = 10000,
                                                 stepsizeShrink = 0.8, min_alpha = 1e-10, convergence_error = 1e-7,
                                                 alpha = 0.5, alpha_g = 0.5, a1 = 3.7, a2 = 3, kernel = "gaussian",
                                                 penalty = "L1", penalty_g = "L1", kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                                                 trend_g = "monotone", gamma_start_input = NULL, gamma_start_default = "zeros",
                                                 max_iter_g = 10000, delta_g = 1e-5, max_alpha_g = 10000,
                                                 stepsizeShrink_g = 0.8, min_alpha_g = 1e-12, convergence_error_g = 1e-7,
                                                 run_aauc = FALSE, log_file = "log_sim_pye_synthetic.txt") {

  # Start calculating simulation time
  start_time <- Sys.time()

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 2 || n != floor(n)) {stop("Parameter 'n' must be a single integer value and at least 2.")}
  if (!is.numeric(rows_train) || length(rows_train) != 1 || rows_train < 1 || !is.integer(as.integer(rows_train))) {stop("Parameter 'rows_train' must be a single integer value and at least 1.")}
  if (!is.numeric(rows_test) || length(rows_test) != 1 || rows_test < 1 || !is.integer(as.integer(rows_test))) {stop("Parameter 'rows_test' must be a single integer value and at least 1.")}
  if (!is.numeric(cols) || length(cols) != 1 || cols < 1 || !is.integer(as.integer(cols))) {stop("Parameter 'cols' must be a single integer value and at least 1.")}
  if (!is.numeric(cols_cov) || length(cols_cov) != 1 || cols_cov < 0 || !is.integer(as.integer(cols_cov))) {stop("Parameter 'cols_cov' must be a single integer value and at least 0.")}
  if (!is.numeric(max_rho) || length(max_rho) != 1 || max_rho < 0 || max_rho > 1) {stop("Parameter 'max_rho' must be a single numeric value between 0 and 1.")}

  if (is.null(mu)) {
    mu <- rep(0, cols)
  } else if (!is.numeric(mu) || length(mu) != cols) {
    stop("Parameter 'mu' must be a numeric vector of length 'cols'.")
  }
  if (is.null(mu_cov)) {
    mu_cov <- rep(0, cols_cov)
  } else if (!is.numeric(mu_cov) || length(mu_cov) != cols_cov) {
    stop("Parameter 'mu_cov' must be a numeric vector of length 'cols_cov'.")
  }

  # Check if tau exists when c_function_of_covariates = TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
    if (is.null(C) || length(C) == 0) { stop("Parameter 'C' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
	if (!is.numeric(w) || length(w) != 1 || w < 0 || w > 1) {stop("Parameter 'w' must be a single numeric value between 0 and 1.")}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter 'trace' has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}
  if (!is.logical(c_function_of_covariates)) {stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_zero_fixed)) {stop("Parameter 'c_zero_fixed' must be a logical (TRUE/FALSE).")}
  if (!(beta_start_default %in% c("zeros", "corr"))) {stop("Parameter 'beta_start_default' must be 'zeros' or 'corr'.")}
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 2 || !is.integer(as.integer(max_iter))) {stop("Parameter 'max_iter' needs to be an integer and at least 2.")}
  if (!(trend %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0) {stop("Parameter 'delta' must be a single positive numeric value.")}
  if (!is.numeric(max_alpha) || length(max_alpha) != 1 || max_alpha <= 0) {stop("Parameter 'max_alpha' must be a single positive numeric value.")}
  if (!is.numeric(stepsizeShrink) || length(stepsizeShrink) != 1 || stepsizeShrink <= 0 || stepsizeShrink >= 1) {stop("Parameter 'stepsizeShrink' must be a single numeric value between 0 and 1 (exclusive).")}
  if (!is.numeric(min_alpha) || length(min_alpha) != 1 || min_alpha <= 0) {stop("Parameter 'min_alpha' must be a single positive numeric value.")}
  if (!is.numeric(convergence_error) || length(convergence_error) != 1 || convergence_error <= 0) {stop("Parameter 'convergence_error' must be a single positive numeric value.")}
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {stop("Parameter 'alpha' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(alpha_g) || length(alpha_g) != 1 || alpha_g < 0 || alpha_g > 1) {stop("Parameter 'alpha_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(a1) || length(a1) != 1 || a1 < 0) {stop("Parameter 'a1' must be a single non-negative numeric value.")}
  if (!is.numeric(a2) || length(a2) != 1 || a2 < 0) {stop("Parameter 'a2' must be a single non-negative numeric value.")}
  if (!is.numeric(a1_g) || length(a1_g) != 1 || a1_g < 0) {stop("Parameter 'a1_g' must be a single non-negative numeric value.")}
  if (!is.numeric(a2_g) || length(a2_g) != 1 || a2_g < 0) {stop("Parameter 'a2_g' must be a single non-negative numeric value.")}

  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel %in% valid_kernels)) {stop(paste0("Parameter 'kernel' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  if (!(kernel_g %in% valid_kernels)) {stop(paste0("Parameter 'kernel_g' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}

  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
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

  # Check if tau exists when c_function_of_covariates = TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

  # If cols_cov is 0, c_function_of_covariates is not performed and set to FALSE.
  if (cols_cov == 0) {
    c_function_of_covariates <- FALSE
    cat("Setting 'c_function_of_covariates' equal to FALSE since 'cols_cov' is equal to zero")
  }

  # Generate seeds
  seeds <- 1:n

  # Simulate one dataset to take the names of the variables and check initial conditions
  df_for_the_names <- create_sample_with_covariates(rows_train = rows_train, cols = cols, cols_cov = cols_cov, max_rho = max_rho, mu = mu,
                                                    mu_cov = mu_cov, rows_test = rows_test, seed = 1)
  X <- df_for_the_names[[X]]
  #y <- df_for_the_names[[y]]
  C <- df_for_the_names[[C]]

  if (trace %in% c(1, 2)) {
    cat("Real regressors of betas are: ", df_for_the_names$regressors, "\n")
    if (c_function_of_covariates == TRUE) {
      cat("Real regressors of gammas are: ", df_for_the_names$covariates, "\n\n\n")
    }
  }

  # Initialize matrices to store results
  names <- paste("seed", seeds, sep = "=")
  pye_L12 <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_L1 <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_EN <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_SCAD <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
  pye_MCP <- matrix(NA, nrow = n, ncol = 2, dimnames = list(names, c("train", "test")))
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
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(0, nrow = length(X), ncol = 1, dimnames = list(X, "n_times_beta_diff_zero"))
  #betas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))])), nrow =  ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))]), ncol = 1, dimnames = list(colnames(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(X))]), "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(0, nrow = length(C), ncol = 1, dimnames = list(C, "n_times_gamma_diff_zero"))
  #gammas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))])), nrow =  ncol(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))]), ncol = 1, dimnames = list(colnames(df_for_the_names$train_df_scaled[, (names(df_for_the_names$train_df_scaled) %in% c(C))]), "n_times_gamma_diff_zero"))


  func <- function(seed, n, rows_train, rows_test,
                   cols, cols_cov, max_rho, mu, mu_cov,
                   lambda, tau, w, w_g,
                   beta_start_input, beta_start_default,
                   gamma_start_input, gamma_start_default,
                   max_iter, max_iter_g,
                   trace,
                   alpha, alpha_g,
                   a1, a2, a1_g, a2_g,
                   penalty, penalty_g,
                   trend, trend_g,
                   kernel, kernel_g,
                   delta, delta_g,
                   max_alpha, max_alpha_g,
                   stepsizeShrink, stepsizeShrink_g,
                   min_alpha, min_alpha_g,
                   convergence_error, convergence_error_g,
                   c_zero_fixed, c_function_of_covariates, run_aauc) {

    if (trace %in% c(1, 2)) {
      cat("\n--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda = ", lambda, "; weight:", w, "; penalty = ", penalty, "\n")
      if (c_function_of_covariates == TRUE) {
        cat("tau = ", tau, "; weight_g:", w_g, "; penalty_g = ", penalty_g, "\n")
      }
    }

    # Create the dataframe
    set.seed(seed) # Use seed_val for reproducibility of each simulation
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
    regressors_gammas <- df$ncovariates

    # Train pye KS
    train_solution <- pye_KS_estimation(df = train_df, X = X, y = y, lambda = lambda, w = w,
                                        beta_start_input = beta_start_input, beta_start_default = beta_start_default,
                                        trace = trace, alpha = alpha, a1 = a1, a2 = a2, max_iter = max_iter, penalty = penalty,
                                        regressors_betas = regressors_betas, trend = trend,
                                        stepsizeShrink = stepsizeShrink, delta = delta,
                                        max_alpha = max_alpha, min_alpha = min_alpha,
                                        convergence_error = convergence_error,
                                        kernel = kernel, c_zero_fixed = c_zero_fixed)

    estimation_time_original_method <- train_solution$estimation_time
    z_hat <- train_solution$z_hat

    train_covYI_solution <- NULL
    estimation_time_covYI <- 0

    # Computing c with covYI estimation
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

      train_covYI_solution <- covYI_KS_estimation(df = cbind(train_df[, names(train_df) != "ID", drop = FALSE], z_hat = train_solution$z_hat[, "z_hat"]),
                                                  z = "z_hat", y = y, C = C, tau = tau, w = w_g,
                                                  gamma_start_input = gamma_start_input1,
                                                  gamma_start_default = gamma_start_default, trace = trace,
                                                  alpha = alpha_g, a1 = a1_g, a2 = a2_g, penalty = penalty_g,
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

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(train_df[[y]]), z_hat[, "z_hat"], levels = c(0, 1), direction = "<", quiet = TRUE)
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input = "specificity", ret = "sensitivity", transpose = FALSE)

    #identify the estimated betas
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    c <- getElement(train_solution, "c_hat")
    niter <- train_solution$niter

    #test the results
    test_solution <- pye_KS_with_print(df = test_df, X = X, y = y, betas = betas,
                                       lambda = lambda, c = c, w = w,
																			 sim_n = seed, n = n, alpha = alpha, a1 = a1, a2 = a2,
                                       penalty = penalty, est_time = estimation_time_original_method, niter = niter, kernel = kernel, trace = trace)

    z_hat <- test_solution$z_hat

    test_covYI_solution <- NULL

    if (c_function_of_covariates == TRUE) {
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df = cbind(test_df[, names(test_df) != "ID", drop = FALSE], z_hat = test_solution$z_hat[, "z_hat"]),
                                      z = "z_hat", y = y, C = C1, gammas = train_covYI_solution$gammas_hat,
                                      tau = tau, w = w_g, kernel = kernel_g, alpha = alpha_g, a1 = a1_g, a2 = a2_g,
                                      penalty = penalty_g, prediction = TRUE, run_aauc = run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1, 2)) {
        cat("-> Results on the TEST SET \n")
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
    roc_spec_points <- seq(0, 1, by = 0.05)
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
    parallel::clusterExport(cl, c("n", "rows_train", "rows_test", "cols", "cols_cov", "max_rho", "mu", "mu_cov",
                                  "lambda", "tau", "w", "w_g", "beta_start_input", "beta_start_default",
                                  "gamma_start_input", "gamma_start_default", "trace", "max_iter", "max_iter_g",
                                  "a1", "a2", "a1_g", "a2_g", "alpha", "alpha_g",
                                  "penalty", "penalty_g", "trend", "trend_g", "kernel", "kernel_g", "c_zero_fixed",
                                  "c_function_of_covariates", "run_aauc",
                                  "delta", "delta_g", "max_alpha", "max_alpha_g",
                                  "stepsizeShrink", "stepsizeShrink_g", "min_alpha", "min_alpha_g",
                                  "convergence_error", "convergence_error_g"), envir = environment())
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed = x, n = n, rows_train = rows_train, rows_test = rows_test,
                                                                        cols = cols, cols_cov = cols_cov, max_rho = max_rho,
                                                                        mu = mu, mu_cov = mu_cov,
                                                                        lambda = lambda, tau = tau,
																																				w = w, w_g = w_g,
                                                                        beta_start_input = beta_start_input,
                                                                        beta_start_default = beta_start_default,
                                                                        gamma_start_input = gamma_start_input,
                                                                        gamma_start_default = gamma_start_default,
                                                                        max_iter = max_iter,
                                                                        max_iter_g = max_iter_g,
                                                                        trace = trace,
                                                                        alpha = alpha,
                                                                        alpha_g = alpha_g,
                                                                        a1 = a1, a2 = a2, a1_g = a1_g, a2_g = a2_g,
                                                                        penalty = penalty, penalty_g = penalty_g,
                                                                        trend = trend, trend_g = trend_g,
                                                                        kernel = kernel,
                                                                        kernel_g = kernel_g, delta = delta, delta_g = delta_g,
                                                                        max_alpha = max_alpha, max_alpha_g = max_alpha_g,
                                                                        stepsizeShrink = stepsizeShrink,
                                                                        stepsizeShrink_g = stepsizeShrink_g,
                                                                        min_alpha = min_alpha, min_alpha_g = min_alpha_g,
                                                                        convergence_error = convergence_error,
                                                                        convergence_error_g = convergence_error_g,
                                                                        c_zero_fixed = c_zero_fixed,
                                                                        c_function_of_covariates = c_function_of_covariates,
                                                                        run_aauc = run_aauc))
  } else {
    # Sequential execution
    cat("Running simulation in sequential mode (used_cores = 1).\n")
    simulation_study <- lapply(seeds, function(x) func(seed = x, n = n, rows_train = rows_train, rows_test = rows_test,
                                                       cols = cols, cols_cov = cols_cov, max_rho = max_rho,
                                                       mu = mu, mu_cov = mu_cov,
                                                       lambda = lambda, tau = tau,
																											 w = w, w_g = w_g,
                                                       beta_start_input = beta_start_input,
                                                       beta_start_default = beta_start_default,
                                                       gamma_start_input = gamma_start_input,
                                                       gamma_start_default = gamma_start_default,
                                                       max_iter = max_iter, max_iter_g = max_iter_g,
                                                       trace = trace, alpha = alpha, alpha_g = alpha_g,
                                                       a1 = a1, a2 = a2, a1_g = a1_g, a2_g = a2_g,
                                                       penalty = penalty, penalty_g = penalty_g,
                                                       trend = trend, trend_g = trend_g,
                                                       kernel = kernel, kernel_g = kernel_g,
                                                       delta = delta, delta_g = delta_g,
                                                       max_alpha = max_alpha, max_alpha_g = max_alpha_g,
                                                       stepsizeShrink = stepsizeShrink,
                                                       stepsizeShrink_g = stepsizeShrink_g,
                                                       min_alpha = min_alpha, min_alpha_g = min_alpha_g,
                                                       convergence_error = convergence_error,
                                                       convergence_error_g = convergence_error_g,
                                                       c_zero_fixed = c_zero_fixed,
                                                       c_function_of_covariates = c_function_of_covariates,
                                                       run_aauc = run_aauc))
  }

  estimation_time_original_method <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(seeds, function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  # Fill the matrices with results
  # Measures on the train set (pye)
  temp_pye <- get(paste0("pye_", penalty))
  temp_pye[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, paste0("pye_KS_", penalty))))

  list_of_measures <- c("auc", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }

  # Measures on the train set (covYI)
  if (c_function_of_covariates == TRUE) {
    #measures on the train set
    list_of_measures_covYI <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      mes[, 1] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$train_covYI_solution, list_of_measures_covYI[i])))
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
  betas_times_selected <- rowSums(sapply(seeds, function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty))) != 0)

  # Measures on test set (pye)
  temp_pye[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, paste0("pye_", penalty))))
  assign(paste0("pye_", penalty), temp_pye)

  for (i in seq_along(list_of_measures)) {
    mes <- get(list_of_measures[i])
    mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if (c_function_of_covariates == TRUE) {
    for (i in seq_along(list_of_measures_covYI)) {
      mes <- get(paste0(list_of_measures_covYI[i], "_covYI"))
      mes[, 2] <- unlist(lapply(seeds, function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures_covYI[i])))
      assign(paste0(list_of_measures_covYI[i], "_covYI"), mes)
    }
  }

  # ROC curve data
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(seeds, function(x) simulation_study[[x]]$roc_sens_points_test)

  betas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))
  if (c_function_of_covariates == TRUE) {
    gammas <- lapply(seeds, function(k) getElement(simulation_study[[k]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))
  } else {
    gammas <- NULL
  }

  regressors_betas <- simulation_study[[1]]$regressors_betas
  regressors_gammas <- simulation_study[[1]]$regressors_gammas

  # End computing sim. time
  simulation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(1, 2)) {
    cat("Total Simulation Time: ", format(simulation_time, digits = 4), " mins.\n")
  }

  results <- list(simulation_time = simulation_time,
                  estimation_time_original_method = estimation_time_original_method,
                  estimation_time_covYI = estimation_time_covYI,
                  used_cores = used_cores, n = n, lambda = lambda, tau = tau,
                  pye_L12 = pye_L12, pye_L1 = pye_L1, pye_EN = pye_EN,
									pye_SCAD = pye_SCAD, pye_MCP = pye_MCP,
                  auc = auc, youden_index = youden_index,
                  sensitivity = sensitivity,
									specificity = specificity,
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
                    lambda = lambda,
                    tau = tau,
                    w = w,
                    w_g = w_g,
                    trace = trace,
                    used_cores = used_cores,
                    c_function_of_covariates = c_function_of_covariates,
                    c_zero_fixed = c_zero_fixed,
                    beta_start_input = beta_start_input,
                    beta_start_default = beta_start_default,
                    max_iter = max_iter,
                    trend = trend,
                    delta = delta,
                    max_alpha = max_alpha,
                    stepsizeShrink = stepsizeShrink,
                    min_alpha = min_alpha,
                    convergence_error = convergence_error,
                    alpha = alpha,
                    alpha_g = alpha_g,
                    a1 = a1,
                    a2 = a2,
                    kernel = kernel,
                    penalty = penalty,
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
