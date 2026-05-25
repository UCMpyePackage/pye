#Bioinformatics, 2014, titled "AucPR: An AUC-based approach
#through penalized regression for disease prediction using
#high-dimensional omics data".


#' @title Calculate Elastic Net Coefficients
#' @description Calculates coefficients for elastic net regression
#' @param cvm Cross-validation mean
#' @param mu Mean values
#' @param lambda Regularization parameter (default: 0.1)
#' @param alpha Elastic net mixing parameter (default: 1)
#' @return Numeric vector of coefficients
#' @noRd
#' @keywords internal
Coefs_Elas <- function(cvm, mu,
                       lambda = 0.1, alpha = 1,
                       family = "gaussian",
                       intercept = FALSE) {
  if (!is.numeric(cvm) || !is.numeric(mu))  stop("cvm and mu must be numeric")
  if (lambda < 0 || alpha < 0 || alpha > 1) stop("Invalid lambda or alpha parameters")

  fitm <- glmnet::glmnet(x = cvm, y = mu, alpha = alpha, intercept = intercept,
                         lambda = lambda, family = family)
  coefs <- fitm$beta
  return(as.numeric(coefs))
}

#' @title AucPR Estimation for Coefficient and Feature Selection
#'
#' @description function to estimate the optimal value of betas using the
#' AucPR method.
#'
#' @description This function estimates the optimal regression coefficients
#' (betas) for a binary classification problem by maximizing the penalized
#' Area Under the Precision-Recall Curve (AUC-PR). The penalization is either
#' L1 (Lasso) or Elastic-Net, which helps in feature selection and regularization.
#' The function returns the estimated coefficients along with various
#' performance metrics.
#'
#' @param df A data frame containing the complete dataset.
#' @param X A character vector of column names from `df` to be used as
#'   regressor variables. Alternatively, a data frame containing only the
#'   regressors. If not specified, all columns not in `y` are used.
#' @param y A character string specifying the column name of the target
#'   variable in `df`. It must be a binomial variable with values 0 or 1.
#'   Alternatively, a data frame containing only the target variable. Default
#'   is "y".
#' @param lambda A single numeric value for the regularization parameter. A
#'   larger value results in stronger penalization.
#' @param alpha The elastic-net mixing parameter. A value of 1 corresponds
#'   to L1 penalization (Lasso), while a value between 0 and 1 corresponds
#'   to Elastic-Net. Default is 1.
#' @param c_to_use character or numeric value specifying the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If set
#'   to `"Youden"`, the function computes the optimal cut-point using Youden's Index.
#'   If a numeric value is provided, it is used directly as the threshold.
#'   Default is "Youden".
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param trace An integer controlling the amount of output. `0` for no
#'   output, `1` for a summary of the results, and `2` for a detailed
#'   breakdown of all steps. Default is 1.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing the optimal betas and a collection of performance
#'   metrics.
#'     \item{model_type}{A character string indicating the final model type used,
#'       either "AucPR_L1" or "AucPR_EN".}
#'     \item{betas_hat}{The estimated vector of optimal beta coefficients.}
#'     \item{X_model}{A character vector containing the names of all regressor
#'     variables (\code{X}) used in the model estimation. This vector defines
#'     the canonical order of the regressors used to align the coefficients
#'     in \code{betas_hat}.}
#'     \item{youden_index}{The Youden Index value at the optimal cut-point,
#'       which is a measure of the overall effectiveness of a diagnostic test.}
#'     \item{sensitivity}{The sensitivity (True Positive Rate) at the optimal
#'       cut-point.}
#'     \item{specificity}{The specificity (True Negative Rate) at the optimal
#'       cut-point.}
#'     \item{geometric_mean}{The geometric mean of sensitivity and
#'       specificity.}
#'     \item{fdr}{The False Discovery Rate at the optimal cut-point.}
#'     \item{mcc}{The Matthews Correlation Coefficient, a measure of the
#'       quality of binary classification.}
#'     \item{corrclass}{The overall correct classification rate.}
#'     \item{auc}{The Area Under the ROC Curve (AUC).}
#'     \item{lambda}{The penalization parameter used in the estimation.}
#'     \item{alpha}{The alpha parameter used in the estimation.}
#'     \item{c_hat}{The optimal cut-point for the `z_hat` score.}
#'     \item{z_hat}{A data frame with the estimated linear combination of
#'       regressors for each observation.}
#'     \item{y_hat}{A data frame with the predicted binary outcomes (0 or 1)
#'       for each observation.}
#'     \item{n_betas}{The number of non-zero beta coefficients estimated.}
#'     \item{n_total_var}{The total number of regressors considered.}
#'     \item{n_predicted_zeros}{The number of estimated betas that are zero.}
#'     \item{n_predicted_non_zeros}{The number of estimated betas that are
#'       non-zero.}
#'     \item{n_caught_betas}{The number of true non-zero betas correctly
#'       identified as non-zero (if `regressors_betas` is provided).}
#'     \item{n_non_caught_betas}{The number of true non-zero betas incorrectly
#'       identified as zero.}
#'     \item{n_caught_zero}{The number of true zero betas correctly identified
#'       as zero.}
#'     \item{n_zero_not_caught}{The number of true zero betas incorrectly
#'       identified as non-zero.}
#'     \item{estimation_time}{The time taken for the estimation, in minutes.}
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' regressors_betas <- sim_data$nregressors
#' alpha <- 1
#' lambda <- 0.1
#'
#' AucPR_result <- AucPR_estimation(df = df, X = X, y = y, lambda = lambda, alpha = alpha,
#'         regressors_betas = regressors_betas, trace = 1)
#'
#' print(AucPR_result)
#'
#' @importFrom glmnet glmnet
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom stats cov
#' @importFrom stats setNames
#' @export
AucPR_estimation <- function (df,
                              X = NULL,
                              y = "y",
                              lambda,
                              alpha = 1,
                              c_to_use = "Youden",
                              fold = NULL,
                              regressors_betas = NULL,
                              trace = 1,
                              max.print = 10,
															...) {


   if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
  # Set max.print temporarily
  old_options <- options(max.print = max.print)
  on.exit(options(old_options))

  # Start calculation of estimation time
  start_time <- Sys.time()

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

   # Handle y parameter
  if (inherits(y, "data.frame")) y <- names(y)[1]
  if (!is.character(y) || length(y) != 1) stop("'y' must be a single column name.")
  if (!(y %in% names(df))) stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")

  # Handle X and y parameters
  if (is.null(X)) X <- setdiff(names(df), y)
  if (inherits(X, "data.frame")) X <- names(X)
  if (!is.character(X) || length(X) == 0) stop("'X' must be a character vector of column names or a data.frame.")
  if (!all(X %in% names(df))) stop("Not all regressors in 'X' are found in 'df'.")


  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Further input validation
  # Validate alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha > 1) {stop("The parameter 'alpha' must be a single numeric value between 0 (exclusive) and 1 (inclusive).")}
  model_type <- if (alpha == 1) "AucPR_L1" else "AucPR_EN"
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  if (!is.null(fold) && (!is.numeric(fold) || length(fold) != 1)) {stop("The parameter 'fold' must be a single numeric value or NULL.")}
  if (!is.null(regressors_betas) && !is.numeric(regressors_betas)) {stop("The parameter 'regressors_betas' must be a numeric vector.")}
  if (!is.numeric(trace) || length(trace) != 1 || !(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(is.character(c_to_use) && c_to_use == "Youden") && !(is.numeric(c_to_use) && length(c_to_use) == 1 && !is.na(c_to_use))) stop("Parameter 'c_to_use' must be either the string 'Youden' or a single numeric value.")

  # Separate cases and controls
  cases <- as.matrix(df1[which(df1[[y]] == 1), X, drop = FALSE])
  controls <- as.matrix(df1[which(df1[[y]] == 0), X, drop = FALSE])

  if (nrow(cases) < 2 || nrow(controls) < 2) stop("Need at least two observations in each class to compute covariance.")

  # Calculate betas and z_hat score
  #n = nrow(cases)
  #m = nrow(controls)
  cvmd <- stats::cov(cases)
  cvmn <- stats::cov(controls)
  dd <- cvmd + cvmn
  #len = length(lambda)

  mu1 <- colMeans(cases) - colMeans(controls)
  maxabs <- max(abs(mu1))
  if (maxabs == 0) mu <- mu1 else mu <- mu1 / maxabs

  betas_hat <- Coefs_Elas(dd, mu, lambda, alpha = alpha)
  names(betas_hat) <- X

  # Standardize coefficients
  maxcoef <- max(abs(betas_hat))
  if (maxcoef > 0) betas_hat <- betas_hat / maxcoef

  # Predict
  z_hat <- as.matrix(df1[, X, drop = FALSE]) %*% betas_hat
  names(z_hat) <- "z_hat"
  df_hat <- cbind(df1, z_hat)

  # Find optimal measures, some of them based on the Youden's Index
  opt <- OptimalCutpoints::optimal.cutpoints(data = df_hat, X = "z_hat", status = y, methods = "Youden", tag.healthy = 0)

  # Which cut-point to use?
  if (c_to_use == "Youden") {
    # Select the optimal cut-point (using the mean if multiple are found)
    c <- mean(opt$Youden$Global$optimal.cutoff$cutoff)
  } else {
    c <- c_to_use
  }

  # Predict binary outcomes based on the optimal cutpoint
  y_hat <- ifelse(df_hat["z_hat"] < c, 0, 1)
  colnames(y_hat) <- "y_hat"
  df_hat <- cbind(df_hat, y_hat)

  # Calculate performance measures
  auc <- mean(opt$Youden$Global$measures.acc$AUC)
  youden_index <- mean(opt$Youden$Global$optimal.criterion)

  # Confusion matrix components
  TP <- sum(ifelse(df_hat[y] == 1 & df_hat["z_hat"] >=  c, 1, 0))
  TN <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] < c, 1, 0))
  FP <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] >=  c, 1, 0))
  FN <- sum(ifelse(df_hat[y] == 1 & df_hat["z_hat"] < c, 1, 0))

  # Sensitivity, Specificity, and other metrics
  spec <- TN / (TN + FP)
  fnr <- FN / (FN + TP)
  # Sensitivity
  sensitivity <- 1 - fnr
  # Geometric mean
  gm <- sqrt(spec * sensitivity)
  # FDR
  fdr <- FP / (FP + TP)
  # MCC
  mcc <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  # CCR
  corrclass <- (TP + TN) / nrow(df_hat)

  # If all betas are zero, measures are set to zero to avoid NaN
  if (sum(betas_hat) == 0) {
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
    gm <- 0
  }

  # Prepare the output
  betas_hat <- betas_hat
  names(betas_hat) <- X
  auc <- auc
  corrclass <- corrclass
  lambda <- lambda
  c_hat <- c
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  if (trace %in% c(1, 2)) {
    cat("Estimation done using the Penalized AUC method type:", model_type, "; \n")
    cat("-> ")
    if (!is.null(fold)) {cat("fold =", fold, "; ")}
    if (!is.na(lambda)) {cat("lambda:", lambda, ";")}
    visualize_betas <- c(betas_hat[which(betas_hat != 0)], stats::setNames(c_hat, "c"))
    cat("youden_index:", youden_index, "; sensitivity:", sensitivity, "; specificity:", spec, "; geometric_mean:", gm, "; fdr:", fdr, "; mcc:", mcc, "; auc:", auc, "; corrclass:", corrclass, "; \n")
    cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, "; betas_hat: \n")
    print(visualize_betas)
  }

  # Number of selected betas
  n_betas <- sum(betas_hat != 0)
  # Total number of variables
  n_total_var <- length(betas_hat)
  # n of predicted zeros
  n_predicted_zeros <- sum(betas_hat == 0)
  # n of predicted betas different from 0
  n_predicted_non_zeros <- sum(betas_hat != 0)

  # compute other measures
  if (!is.null(regressors_betas)) {
    if (length(betas_hat) != length(regressors_betas)) {stop("The length of 'betas_hat' does not match the true betas in regressors_betas.")}
    # n of beta caught
    n_caught_betas <- sum(betas_hat[regressors_betas != 0] != 0)
    # n of beta not caught
    n_non_caught_betas <- sum(betas_hat[regressors_betas != 0] == 0)
    # n of zeros caught
    n_caught_zero <- sum(betas_hat[regressors_betas == 0] == 0)
    # n of zeros not caught
    n_zero_not_caught <- sum(betas_hat[regressors_betas == 0] != 0)
  } else {
    n_caught_betas <- NA
    n_non_caught_betas <- NA
    n_caught_zero <- NA
    n_zero_not_caught <- NA
  }

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(1, 2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  return(list(model_type = model_type,
              c_to_use = c_to_use,
              betas_hat = betas_hat,
              X_model = X,
              youden_index = youden_index,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr, mcc = mcc,
              corrclass = corrclass,
              auc = auc, lambda = lambda, alpha = alpha,
              c_hat = c_hat, z_hat = z_hat, y_hat = y_hat,
              n_betas = n_betas,
              n_total_var = n_total_var,
              n_predicted_zeros = n_predicted_zeros,
              n_predicted_non_zeros = n_predicted_non_zeros,
              n_caught_betas = n_caught_betas,
              n_non_caught_betas = n_non_caught_betas,
              n_caught_zero = n_caught_zero,
              n_zero_not_caught = n_zero_not_caught,
              estimation_time = estimation_time))
}






#' @title Prediction and Performance Evaluation for Penalized AucPR Models
#'
#' @description This function applies a previously trained penalized AUC-PR
#' model to new data to generate predictions and evaluate performance metrics.
#' The model to be used is the result of the `AucPR_estimation` function.
#'
#' @param df A data frame containing the new dataset for prediction.
#' @param y A character string specifying the column name of the target
#'   variable in `df`. It must be a binomial variable with values 0 or 1.
#'   Default is "y".
#' @param model_to_use A list containing the parameters from a model trained
#'   by the `AucPR_estimation` function, specifically the `betas_hat` and
#'   `c_hat` components.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param trace An integer controlling the amount of output. `0` for no
#'   output, `1` for a summary of the results, and `2` for a detailed
#'   breakdown of all steps. Default is 1.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#' @param c_function_of_covariates A logical flag to suppress printing of
#'   results. If `TRUE`, the function will not print the summary to the
#'   console, as the output is expected to be handled by an external function
#'   (e.g., `covYI`). Default is `FALSE`.
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing the predictions and a collection of performance
#'   metrics, including:
#'     \item{model_type}{The model type ("AucPR_L1" or "AucPR_EN").}
#'     \item{betas_hat}{The estimated vector of optimal beta coefficients.}
#'     \item{youden_index}{The Youden Index value at the optimal cut-point.}
#'     \item{sensitivity}{The sensitivity at the optimal cut-point.}
#'     \item{specificity}{The specificity at the optimal cut-point.}
#'     \item{geometric_mean}{The geometric mean of sensitivity and specificity.}
#'     \item{fdr}{The False Discovery Rate at the optimal cut-point.}
#'     \item{mcc}{The Matthews Correlation Coefficient.}
#'     \item{corrclass}{The overall correct classification rate.}
#'     \item{auc}{The Area Under the ROC Curve (AUC).}
#'     \item{lambda}{The penalization parameter used in the estimation.}
#'     \item{alpha}{The alpha parameter used in the estimation.}
#'     \item{c_hat}{The optimal cut-point for the `z_hat` score.}
#'     \item{z_hat}{A data frame with the estimated linear combination of
#'       regressors for each observation.}
#'     \item{y_hat}{A data frame with the predicted binary outcomes (0 or 1).}
#'     \item{n_total_var}{The total number of regressors considered.}
#'     \item{n_predicted_zeros}{The number of estimated betas that are zero.}
#'     \item{n_predicted_non_zeros}{The number of estimated betas that are
#'       non-zero.}
#'     \item{n_caught_betas}{The number of true non-zero betas correctly
#'       identified as non-zero (if `regressors_betas` is provided).}
#'     \item{n_non_caught_betas}{The number of true non-zero betas incorrectly
#'       identified as zero.}
#'     \item{n_caught_zero}{The number of true zero betas correctly identified
#'       as zero.}
#'     \item{n_zero_not_caught}{The number of true zero betas incorrectly
#'       identified as non-zero.}
#'     \item{TP, TN, FP, FN}{The counts from the confusion matrix.}
#'     \item{estimation_time}{The time taken for the prediction.}
#'
#' @examples
#' library(pye)
#'
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' train_df <- sim_data$train_df_scaled
#' test_df <- sim_data$test_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' regressors_betas <- sim_data$nregressors
#' alpha <- 0.5
#' lambda <- 0.1
#' c_function_of_covariates <- FALSE
#'
#' # 1. Train the model using AucPR_estimation
#' model <- AucPR_estimation(df = train_df, y = y, X = X, alpha = alpha,
#'    lambda = lambda, regressors_betas = regressors_betas, trace = 1)
#'
#' # 2. Apply the trained model to new data using AucPR_predict
#' predictions <- AucPR_predict(df = test_df, model_to_use = model,
#'    trace = 1, regressors_betas = regressors_betas,
#'    c_function_of_covariates = c_function_of_covariates)
#'
#' print(predictions)
#'
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @export
AucPR_predict <- function (df,
                           y = "y",
                           model_to_use,
                           fold = NULL,
                           regressors_betas = NULL,
                           trace = 1,
                           max.print = 10,
                           c_function_of_covariates = FALSE,
													 ...) {

	if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
  old_options <- options(max.print = max.print)
  on.exit(options(old_options))

  # Start calculation of estimation time
  start_time <- Sys.time()

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

  # Handle X parameter
  #if (inherits(X, "data.frame")) {
  #  X <- names(X)
  #} else if (!inherits(X, "character")) {
  #  stop("Parameter 'X' must be a character vector of column names or a data.frame.")
  #}
  #if (length(X) == 0) {
  #  stop("No regressors (X) specified or found in the data frame.")
  #}
  #if (!all(X %in% names(df))) {
  #  stop("Not all specified regressors (X) are found in the input data frame 'df'.")
  #}

  # Handle y parameter
  if (inherits(y, "data.frame")) y <- names(y)[1]
  if (!is.character(y) || length(y) != 1) stop("'y' must be a single column name.")
  if (!(y %in% names(df))) stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")

  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}
  if (!is.list(model_to_use) || !all(c("betas_hat", "c_hat", "model_type", "lambda", "alpha") %in% names(model_to_use))) {
    stop("The parameter 'model_to_use' must be a valid list result from the 'AucPR_estimation' function.")
  }
  if (!is.null(fold) && (!is.numeric(fold) || length(fold) != 1)) {stop("The parameter 'fold' must be a single numeric value or NULL.")}
  if (!is.null(regressors_betas) && !is.numeric(regressors_betas)) {stop("The parameter 'regressors_betas' must be a numeric vector.")}
  if (!is.numeric(trace) || length(trace) != 1 || !(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!is.logical(c_function_of_covariates) || length(c_function_of_covariates) != 1) {stop("The parameter 'c_function_of_covariates' must be a single logical value (TRUE/FALSE).")}

  # Extract model components
  betas_hat <- model_to_use$betas_hat
  c_hat <- model_to_use$c_hat
  X <- model_to_use$X_model

  # Ensure all variables required by the model (X) are present in the new data (df)
  if (length(X) == 0) {stop("The model contains no regressors (X). Cannot perform prediction.")}
  if (!all(X %in% names(df))) {
    missing_vars <- X[!X %in% names(df)]
    stop("The following required regressors from the trained model are missing in the new data: ", paste(missing_vars, collapse = ", "))
  }

  # Prepare data for prediction
  ID <- rownames(df)
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Predict z_hat score
  z_hat <- as.matrix(df1[, names(betas_hat), drop = FALSE]) %*% betas_hat
  names(z_hat) <- "z_hat"
  df_hat <- cbind(df1, z_hat)

  # Predict binary outcomes
  y_hat <- ifelse(df_hat["z_hat"] < c_hat, 0, 1)
  colnames(y_hat) <- "y_hat"
  df_hat <- cbind(df_hat, y_hat)

  # Calculate performance measures
  # AUC e YI
  opt <- OptimalCutpoints::optimal.cutpoints(data = df_hat, X = "z_hat", status = y, methods = "Youden", tag.healthy = 0)
  auc <-  mean(opt$Youden$Global$measures.acc$AUC)
  youden_index <- mean(opt$Youden$Global$optimal.criterion)

  # Confusion matrix components
  TP <- sum(ifelse(df_hat[y] == 1 & df_hat["z_hat"] >=  c_hat, 1, 0))
  TN <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] < c_hat, 1, 0))
  FP <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] >=  c_hat, 1, 0))
  FN <- sum(ifelse(df_hat[y] == 1 & df_hat["z_hat"] < c_hat, 1, 0))

  # Sensitivity, Specificity, and other metrics, handling division by zero
  spec <- TN / (TN + FP)
  fnr <- FN / (FN + TP)
  #youden_index <- spec - fnr
  #sensitivity
  sensitivity <- 1 - fnr
  #geometric mean
  gm <- sqrt(spec * sensitivity)
  #FDR
  fdr <- FP / (FP + TP)
  #MCC
  mcc <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  #compute the correct classiication:
  corrclass <- (TP + TN) / nrow(df_hat)

  #if all the betas are zero, the measures are zeros
  if (sum(betas_hat) == 0) {
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    spec <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  # Prepare the output
  lambda <- model_to_use$lambda
  alpha <- model_to_use$alpha
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  if (trace %in% c(1, 2)) {
    #print the results only if c_function_of_covariates = FALSE
    if (c_function_of_covariates == FALSE) {
      cat("Prediction executed using the AucPR model type:", model_to_use$model_type, "; \n")
      cat("-> ")
      if (!is.null(fold)) {cat("fold =", fold, "; ")}
      if (!is.na(lambda)) {cat("lambda:", lambda, ";")}
      visualize_betas <- betas_hat[which(betas_hat != 0)]
      cat("youden_index:", youden_index, "; sensitivity:", sensitivity, "; specificity:", spec, "; geometric_mean:", gm, "; fdr:", fdr, "; mcc:", mcc, "; auc:", auc, "; corrclass:", corrclass, "; \n")
      cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, "; betas_hat: \n")
      print(visualize_betas)
    }
  }

  # Number of selected betas
  n_betas <- sum(betas_hat != 0)
  # Total number of variables
  n_total_var <- length(betas_hat)
  # n of predicted zeros
  n_predicted_zeros <- sum(betas_hat == 0)
  # n of predicted betas different from 0
  n_predicted_non_zeros <- sum(betas_hat != 0)

  # Comparison with true betas (if provided)
  n_caught_betas <- NA
  n_non_caught_betas <- NA
  n_caught_zero <- NA
  n_zero_not_caught <- NA
  if (!is.null(regressors_betas)) {
    if (length(betas_hat) != length(regressors_betas)) {
      stop("The length of 'betas_hat' does not match the true betas in regressors_betas.")
    }
    # n of true non-zero betas caught as non-zero
    n_caught_betas <- sum(betas_hat[regressors_betas != 0] != 0)
    # n of true non-zero betas missed (False Negatives)
    n_non_caught_betas <- sum(betas_hat[regressors_betas != 0] == 0)
    # n of true zero betas caught as zero (True Negatives)
    n_caught_zero <- sum(betas_hat[regressors_betas == 0] == 0)
    # n of true zero betas incorrectly caught as non-zero (False Positives)
    n_zero_not_caught <- sum(betas_hat[regressors_betas == 0] != 0)
  }

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
  #print the results only if c_function_of_covariates = FALSE
  if (c_function_of_covariates == FALSE) {
    if (trace %in% c(1, 2)) {
      cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
    }
  }

  return(list(model_type = model_to_use$model_type,
              betas_hat = model_to_use$betas_hat,
              youden_index = youden_index,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr,
              mcc = mcc,
              corrclass = corrclass,
              auc = auc,
              lambda = lambda,
              alpha = alpha,
              c_hat = c_hat,
              z_hat = z_hat,
              y_hat = y_hat,
              n_betas = n_betas,
              n_total_var = n_total_var,
              n_predicted_zeros = n_predicted_zeros,
              n_predicted_non_zeros = n_predicted_non_zeros,
              n_caught_betas = n_caught_betas,
              n_non_caught_betas = n_non_caught_betas,
              n_caught_zero = n_caught_zero,
              n_zero_not_caught = n_zero_not_caught,
              TP = TP, TN = TN, FP = FP, FN = FN,
              estimation_time = estimation_time))
}


#create the output class of the AucPR.cv function
setClass(Class = "AucPR_cross_validation_output",
         representation(
           cv_time = "ANY",
           auc = "ANY",
           aauc = "ANY",
           aYI = "ANY",
           youden_index = "ANY",
           sensitivity = "ANY",
           specificity = "ANY",
           geometric_mean = "ANY",
           fdr = "ANY",
           mcc = "ANY",
           corrclass = "ANY",
           n_betas = "ANY",
           n_gammas = "ANY",
           betas = "ANY",
           gammas = "ANY"
         )
)

#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom methods new
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext file_ext
#' @keywords internal
#' @noRd
AucPR.cv <- function (df, X, y, C,
                      alpha, lambda, tau,
											w_g,
                      c_to_use, alpha_g,
                      penalty_g, folds_i, k,
                      regressors_betas,
                      trace, model_type, kernel_g,
                      a1_g, a2_g, trend_g,
                      gamma_start_input, gamma_start_default,
                      regressors_gammas = NULL,
                      max_iter_g, delta_g,
                      max_alpha_g, stepsizeShrink_g,
                      min_alpha_g,
                      convergence_error_g,
                      auc, aauc, aYI,
                      youden_index, sensitivity,
                      specificity,
                      geometric_mean, fdr, mcc,
                      corrclass, n_betas,
                      n_gammas, used_cores,
                      c_function_of_covariates,
                      simultaneous, run_aauc,
                      log_file) {

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  greek <- if (c_function_of_covariates) "tau" else "lambda"
  if (trace %in% c(1, 2)) {
    cat("----------------------------------------------------------------\n")
    cat("|       starting with the", k, "-th fold for the CV of ", greek, "      |\n")
    cat("----------------------------------------------------------------\n")
  }

  #if ((c_function_of_covariates == TRUE) & (length(lambda) == 1)) {
  #  #in this case we have to estimate the combination z_hat over all the dataset
  #  train_df1 <- df
  #} else {
  #  train_df1 <- train_df
  #}

  # --- Beta Estimation ---
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

    # --- Fit primary models ---
    parallel::clusterExport(cl, c("train_df", "X", "y", "regressors_betas", "k", "trace", "alpha", "lambda", "c_to_use"), envir = environment())
    fitted_models <- parallel::parLapply(cl, lambda, function(x) AucPR_estimation(df = train_df, X = X, y = y,
                                                                                   alpha = alpha,
                                                                                   lambda = x,
                                                                                   c_to_use = c_to_use,
                                                                                   regressors_betas = regressors_betas,
                                                                                   fold = k,
                                                                                   trace = trace))

  } else {

    fitted_models <- lapply(lambda, function(x) AucPR_estimation(df = train_df, X = X, y = y,
                                                                 alpha = alpha,
                                                                 lambda = x,
                                                                 c_to_use = c_to_use,
                                                                 regressors_betas = regressors_betas,
                                                                 fold = k,
                                                                 trace = trace))

  }

  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat"))

  if (length(gamma_start_input) == 0) {
  #if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
    gamma_start_input1 <- lapply(seq_along(lambda), function(x) {
                            gamma_start_input1 <- c(getElement(fitted_models[[x]], "c_hat"), rep(0, length(C)))
                            names(gamma_start_input1) <- c("const", C)
                            gamma_start_input1
                            })
  } else {
    gamma_start_input1 <- lapply(seq_along(lambda), function(x) {
                            gamma_start_input1 <- gamma_start_input
                            names(gamma_start_input1) <- c("const", C)
                            gamma_start_input1
                            })
  }

  # --- Fit covariate-adjusted cut-point models (covYI) or collect results ---
  gammas_hat <- lapply(seq_along(lambda), function(y) NA)
  names(gammas_hat) <- lambda
  if (c_function_of_covariates == FALSE) {
    # If c_function_of_covariates is FALSE, collect performance measures on the train set
    for (i in seq_along(lambda)) {
      auc[[i]]$train[k, ] <- getElement(fitted_models[[i]], "auc")
      aauc[[i]]$train[k, ]  <- NA
      aYI[[i]]$train[k, ]  <- NA
      youden_index[[i]]$train[k, ] <- getElement(fitted_models[[i]], "youden_index")
      sensitivity[[i]]$train[k, ] <- getElement(fitted_models[[i]], "sensitivity")
      specificity[[i]]$train[k, ]  <- getElement(fitted_models[[i]], "specificity")
      geometric_mean[[i]]$train[k, ]  <- getElement(fitted_models[[i]], "geometric_mean")
      fdr[[i]]$train[k, ]  <- getElement(fitted_models[[i]], "fdr")
      mcc[[i]]$train[k, ]  <- getElement(fitted_models[[i]], "mcc")
      corrclass[[i]]$train[k, ] <- getElement(fitted_models[[i]], "corrclass")
    }
  } else {
    # If c_function_of_covariates is TRUE, fit covariate-adjusted cut-point models (covYI)

    # --- Parallel / Sequential Execution ---
    if (!is.null(cl)) { # Parallel execution (reusing the cluster)
      if (trace > 0) cat("Start parallel computing for cross-validation of covYI, I am using:", length(cl), "cores.\n")
      if (used_cores > max.cores) {
        warning("The number of specified cores (", used_cores, ") is larger than the number of physical cores available (", max.cores, ")!")
      }
      parallel::clusterExport(cl, c("train_df", "tau", "y", "C", "trace", "penalty_g", "w_g",
                                    "alpha_g", "a1_g", "a2_g", "trend_g", "kernel_g", "k",
                                    "gamma_start_input1", "gamma_start_default", "z_hat",
                                    "regressors_gammas", "max_iter_g", "delta_g", "max_alpha_g", "stepsizeShrink_g",
                                    "min_alpha_g", "convergence_error_g", "run_aauc"), envir = environment())

      if (length(lambda) >= length(tau)) {
        if (trace > 0) cat("Parallel loop is on lambda for covYI.\n")
        fitted_gammas <- parallel::parLapply(cl, seq_along(z_hat), function (x) lapply(tau, function(t) covYI_KS_estimation(df = cbind(train_df, z_hat = z_hat[[x]]$z_hat),
                                                                                                                 z = "z_hat", y = y, C = C, tau = t, w = w_g,
                                                                                                                 gamma_start_input = gamma_start_input1[[x]],
                                                                                                                 gamma_start_default = gamma_start_default, trace = trace,
                                                                                                                 alpha = alpha_g, a1 = a1_g, a2 = a2_g, penalty = penalty_g,
                                                                                                                 max_iter = max_iter_g,
                                                                                                                 min_alpha = min_alpha_g,
                                                                                                                 convergence_error = convergence_error_g,
                                                                                                                 regressors_gammas = regressors_gammas, fold = k,
                                                                                                                 trend = trend_g,
                                                                                                                 stepsizeShrink = stepsizeShrink_g,
                                                                                                                 delta = delta_g, max_alpha = max_alpha_g, kernel = kernel_g,
                                                                                                                 run_aauc = run_aauc)))
      } else {
        if (trace > 0) cat("Parallel loop is on tau for covYI.\n")
        fitted_gammas <- lapply(seq_along(z_hat), function (x) parallel::parLapply(cl, tau, function(t) covYI_KS_estimation(df = cbind(train_df, z_hat = z_hat[[x]]$z_hat),
                                                                                                                 z = "z_hat", y = y, C = C, tau = t, w = w_g,
                                                                                                                 gamma_start_input = gamma_start_input1[[x]],
                                                                                                                 gamma_start_default = gamma_start_default, trace = trace,
                                                                                                                 alpha = alpha_g, a1 = a1_g, a2 = a2_g, penalty = penalty_g,
                                                                                                                 max_iter = max_iter_g,
                                                                                                                 min_alpha = min_alpha_g,
                                                                                                                 convergence_error = convergence_error_g,
                                                                                                                 regressors_gammas = regressors_gammas, fold = k,
                                                                                                                 trend = trend_g,
                                                                                                                 stepsizeShrink = stepsizeShrink_g,
                                                                                                                 delta = delta_g, max_alpha = max_alpha_g, kernel = kernel_g,
                                                                                                                 run_aauc = run_aauc)))
      }

    } else {

      fitted_gammas <- lapply(seq_along(z_hat), function (x) lapply(tau, function(t) covYI_KS_estimation(df = cbind(train_df, z_hat = z_hat[[x]]$z_hat),
                                                z = "z_hat", y = y, C = C, tau = t, w = w_g,
                                                gamma_start_input = gamma_start_input1[[x]],
                                                gamma_start_default = gamma_start_default, trace = trace,
                                                alpha = alpha_g, a1 = a1_g, a2 = a2_g, penalty = penalty_g,
                                                max_iter = max_iter_g,
                                                min_alpha = min_alpha_g,
                                                convergence_error = convergence_error_g,
                                                regressors_gammas = regressors_gammas, fold = k,
                                                trend = trend_g,
                                                stepsizeShrink = stepsizeShrink_g,
                                                delta = delta_g, max_alpha = max_alpha_g, kernel = kernel_g,
                                                run_aauc = run_aauc)))
    }

    names(fitted_gammas) <- paste("lambda", lambda, sep = "=")
    for (zy in seq_along(lambda)) {
      names(fitted_gammas[[zy]]) <- paste("tau", tau, sep = "=")
    }

    # Organize results
    names <- paste("tau", tau, sep = "=")
    gammas_hat <- lapply(lambda, function(x) sapply(names, function(xx) NULL))
    names(gammas_hat) <- paste("lambda", lambda, sep = "=")

    for (i in seq_along(lambda)) {
      # Measures on the train set - these are multiple tables based on the number of considered TAUs
      auc[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "auc")))
      aauc[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "aauc")))
      aYI[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "aYI")))
      youden_index[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "sensitivity")))
      specificity[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "fdr")))
      mcc[[i]]$train[k, ]  <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "mcc")))
      corrclass[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "corrclass")))

      n_gammas[[i]][k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "n_gammas")))
      gammas_hat[[i]] <- lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("gammas_hat_", penalty_g)))
      names(gammas_hat[[i]]) <- paste("tau", tau, sep = "=")
    }

  }

  # n_betas is only based on lambda, not tau!
  n_betas[k, ] <- unlist(lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "n_betas")))

  # create a list of the results - it is only based on lambda, not tau!
  betas_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "betas_hat"))
  names(betas_hat) <- paste("lambda", lambda, sep = "=")

  # c_hat is the fixed value of c in case we don't use the covariates to estimate a patient's specific cut-off value
  c_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "c_hat"))
  names(c_hat) <- paste("lambda", lambda, sep = "=")

  betas <- betas_hat
  gammas <- gammas_hat

  # --- Test Set Evaluation ---
  all_measures_test <- mapply(function(z) AucPR_predict(df = test_df,
                                                        y = y,
                                                        model_to_use = fitted_models[[z]],
                                                        fold = k, trace = trace,
                                                        c_function_of_covariates = c_function_of_covariates),
                              seq_along(lambda))

  # --- Organize test set results ---
  if (c_function_of_covariates == FALSE) {
    # Measures on the test set - if we don't compute c with covariates, we are only dependent of lambda:
    for (i in seq_along(lambda)) {
      auc[[i]]$test[k, ] <- all_measures_test["auc", ][[i]]
      aauc[[i]]$test[k, ] <- NA
      aYI[[i]]$test[k, ] <- NA
      youden_index[[i]]$test[k, ] <- all_measures_test["youden_index", ][[i]]
      sensitivity[[i]]$test[k, ] <- all_measures_test["sensitivity", ][[i]]
      specificity[[i]]$test[k, ] <- all_measures_test["specificity", ][[i]]
      geometric_mean[[i]]$test[k, ] <- all_measures_test["geometric_mean", ][[i]]
      fdr[[i]]$test[k, ] <- all_measures_test["fdr", ][[i]]
      mcc[[i]]$test[k, ] <- all_measures_test["mcc", ][[i]]
      corrclass[[i]]$test[k, ] <- all_measures_test["corrclass", ][[i]]
    }
  } else {
    # Put "const" in C if not already there - it is needed in covYI function
    # Ensure 'const' is the first element in the covariate vector C1 for correct coefficient alignment.
    # First, remove 'const' from C if it exists to avoid duplication.
    C_without_const <- C[C != "const"]
    # Then, prepend 'const' to the vector.
    C1 <- c("const", C_without_const)

    if ("const" %in% names(test_df)) {
      # If 'const' column already exists, just select the required columns in the correct order.
      test_df1 <- test_df[, c(y, C1), drop = FALSE]
    } else {
      # If 'const' column does not exist, create it and bind it to the dataframe.
      const <- rep(1, nrow(test_df))
      test_df1 <- cbind(test_df[, y, drop = FALSE], const, test_df[, C_without_const, drop = FALSE])
      # Ensure column order matches C1
      test_df1 <- test_df1[, c(y, C1)]
    }

    # use covYI to compute the value of c
    listing <- lapply(lambda, function(xx) lapply(tau, function(x) NULL))
    for (l in seq_along(lambda)) {for (t in seq_along(tau)) { listing[[l]][[t]] <- list(gammas_hat[[l]][[t]], tau[[t]])}}

    z_hat <- lapply(seq_along(lambda), function(x) getElement(all_measures_test["z_hat", ][[x]], "z_hat"))

    cov_results <- lapply(seq_along(lambda), function(xx) lapply(listing[[xx]], function(x) covYI_KS(df = cbind(test_df1, z_hat = z_hat[[xx]]),
                                                                                                    z = "z_hat",
                                                                                                    y = y, C = C1,
                                                                                                    gammas = x[[1]],
                                                                                                    tau = x[[2]],
																																																		w = w_g,
                                                                                                    kernel = kernel_g,
                                                                                                    alpha = alpha_g,
                                                                                                    a1 = a1_g, a2 = a2_g,
                                                                                                    penalty = penalty_g,
                                                                                                    prediction = TRUE,
                                                                                                    run_aauc = run_aauc)))

    # Name the lambdas
    names(cov_results) <- paste("lambda", lambda, sep = "=")
    for (i in seq_along(lambda)) {
      # Name the taus
      names(cov_results[[i]]) <- paste("tau", tau, sep = "=")

      if (trace %in% c(1, 2)) {
        # Print the results on the Test Set
        cat("\n")
        cat("Final results on TEST SET of covYI for the FOLD:", k, ": \n")
        cat("With lambda:", lambda[i], " and method:", model_type, "\n")
        visualize_betas <- c(betas_hat[[i]][which(betas_hat[[i]] != 0)], stats::setNames(c_hat[[i]], "c"))
        cat(" youden_index:", all_measures_test["youden_index", ][[i]], "; sensitivity:", all_measures_test["sensitivity", ][[i]], "; specificity:", all_measures_test["specificity", ][[i]], "; geometric_mean:", all_measures_test["geometric_mean", ][[i]], "; fdr:", all_measures_test["fdr", ][[i]], "; mcc:", all_measures_test["mcc", ][[i]], "; auc:", all_measures_test["auc", ][[i]], "; corrclass:", all_measures_test["corrclass", ][[i]], "; \n")
        cat("TP:", all_measures_test["TP", ][[i]], "; TN:", all_measures_test["TN", ][[i]], "; FP:", all_measures_test["FP", ][[i]], "; FN:", all_measures_test["FN", ][[i]], ";  betas: \n")
        print(visualize_betas)
        cat("Cross-validation time:", fitted_models[[i]]$estimation_time, "\n\n")

        for (ii in seq_along(cov_results[[i]])) {
          cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend_g, "version ; \n")
          visualize_gammas <- c(gammas_hat[[i]][[ii]][which(gammas_hat[[i]][[ii]] != 0)])
          cat("tau:", tau[ii], "; weight:", w_g, "; penalty:", penalty_g, "; ", paste0("covYI_KS_", penalty_g), ":" , getElement(cov_results[[i]][[ii]], paste0("covYI_KS_", penalty_g)), "; youden_index:", cov_results[[i]][[ii]]$youden_index, "; aYI:", cov_results[[i]][[ii]]$aYI, "; sensitivity:", cov_results[[i]][[ii]]$sensitivity, "; specificity:", cov_results[[i]][[ii]]$specificity, "; geometric_mean:", cov_results[[i]][[ii]]$geometric_mean, "; fdr:", cov_results[[i]][[ii]]$fdr, "; mcc:", cov_results[[i]][[ii]]$mcc, "; auc:", cov_results[[i]][[ii]]$auc, "; aauc:", cov_results[[i]][[ii]]$aauc, "; corrclass:", cov_results[[i]][[ii]]$corrclass, "; \n")
          cat("TP:", cov_results[[i]][[ii]]$TP, "; TN:", cov_results[[i]][[ii]]$TN, "; FP:", cov_results[[i]][[ii]]$FP, "; FN:", cov_results[[i]][[ii]]$FN, "; gammas: \n")
          print(visualize_gammas)
          cat("Cross-validation time:", fitted_gammas[[i]][[ii]]$estimation_time, "; Number of iterations:", fitted_gammas[[i]][[ii]]$niter, "\n\n\n")
        }
				cat("\n")
      }

      # Store test measures
      auc[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "auc")))
      aauc[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "aauc")))
      aYI[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "aYI")))
      youden_index[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "sensitivity")))
      specificity[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "fdr")))
      mcc[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "mcc")))
      corrclass[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], "corrclass")))
    }
  }

  return(methods::new("AucPR_cross_validation_output", auc = auc,
                                                       aauc = aauc,
                                                       aYI = aYI,
                                                       youden_index = youden_index,
                                                       sensitivity = sensitivity,
                                                       specificity = specificity,
                                                       geometric_mean = geometric_mean,
                                                       fdr = fdr,
                                                       mcc = mcc,
                                                       corrclass = corrclass,
                                                       n_betas = n_betas,
                                                       n_gammas = n_gammas,
                                                       betas = betas,
                                                       gammas = gammas))
}







#' @title Cross-Validation for Optimal AucPR Regularization Parameter Selection
#'
#' @description This function performs k-fold cross-validation to determine the optimal
#' values for the penalization parameters `lambda` and `tau` (if applicable) for the `AucPR` model.
#' The primary goal is to select the hyperparameters that yield the best model
#' performance based on a variety of accuracy measures. The function first
#' divides the dataset into stratified folds to maintain the class distribution of the
#' target variable. It then iterates through the specified `lambda` and `tau`
#' values, training and evaluating the model on each fold. The model's performance
#' is assessed using metrics such as AUC, aAUC, aYI, Youden's Index,
#' sensitivity, specificity, and others. The function returns the results of
#' this process, including the performance metrics for each fold and the
#' hyperparameter values that produced the best overall test performance for
#' each measure.
#'
#'
#' @param n_folds Integer. Number of folds for cross-validation.
#' @param df Data frame. Input dataset containing predictors and target.
#' @param X Character vector or data frame. Names of predictor variables.
#'   Defaults to all columns not in `y` or `C`.
#' @param y Character or data frame. Name of the binary target variable (0/1).
#'   Defaults to `"y"`.
#' @param C Character vector or data frame. Names of covariate variables for
#'   `covYI`. Default is `NULL`.
#' @param alpha Numeric in (0, 1]. Elastic-net mixing parameter for AucPR.
#'   `1` = Lasso, `(0, 1)` = Elastic-Net.
#' @param lambda Numeric vector. Penalization parameter(s) for predictors `X`.
#' @param tau Numeric vector. Penalization parameter(s) for covariates `C` in
#'   `covYI`. Ignored if `c_function_of_covariates = FALSE`.
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the 
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by 
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param c_to_use character or numeric value specifying the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If set
#'   to `"Youden"`, the function computes the optimal cut-point using Youden's Index.
#'   If a numeric value is provided, it is used directly as the threshold.
#'   Default is "Youden".
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param trace Integer (0, 1, 2). Level of printed output.
#' @param seed Integer. Random seed for reproducibility.
#' @param used_cores Integer. Number of CPU cores for parallelization.
#' @param scaling Logical. If `TRUE`, scale predictors and covariates.
#' @param c_function_of_covariates Logical. If `TRUE`, estimate cut-off as a
#'   function of covariates via `covYI`.
#' @param simultaneous Logical. If `TRUE`, estimate betas and gammas jointly.
#' @param measure_to_select_lambda Character. Metric to select `lambda` when
#'   `simultaneous = FALSE`. One of `"auc"`, `"aauc"`, `"aYI"`, `"yi"`, `"sen"`,
#'   `"spc"`, `"gm"`, `"fdr"`, `"mcc"`, `"ccr"`.
#' @param alpha_g Numeric. Elastic-net mixing parameter for `covYI`.
#' @param penalty_g Character. Penalty type for `covYI`: `"L12"`, `"L1"`, `"EN"`,
#'   `"SCAD"`, `"MCP"`.
#' @param kernel_g Character. Kernel type for density estimation in `covYI`.
#' @param a1_g,a2_g Numeric. Parameters for SCAD/MCP penalties in `covYI`.
#' @param trend_g Character. `"monotone"` uses mmAPG, `"nonmonotone"` uses mnmAPG.
#' @param gamma_start_input Numeric vector. Starting point for gammas.
#' @param gamma_start_default Character. `"zeros"` or `"corr"`.
#' @param regressors_gammas Numeric vector. True gamma coefficients (optional).
#' @param max_iter_g Integer. Max iterations for `covYI` optimization.
#' @param delta_g,max_alpha_g,stepsizeShrink_g,min_alpha_g,convergence_error_g
#'   Numeric. Optimization parameters for `covYI`.
#' @param run_aauc Logical. If `FALSE`, skip aAUC/aYI computation.
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If `NULL`, output goes to the console.
#'   Default is `"log_AUC_models.txt"`.
#'
#' @return A list containing the results of the cross-validation process. The list
#' includes:
#'   \item{model_type}{The type of model used (e.g., "AucPR_EN").}
#'   \item{alpha}{The alpha value used for the AucPR model.}
#'   \item{penalty_g}{The penalty type used for the `covYI` model.}
#'   \item{cv_time}{The total time taken for the cross-validation,
#'     in minutes}.
#'   \item{auc_first_step, aauc_first_step, aYI_first_step, 
#'     youden_index_first_step, sensitivity_first_step,
#'     specificity_first_step, geometric_mean_first_step, fdr_first_step, 
#'     mcc_first_step, corrclass_first_step}{A series of nested lists,
#'      each containing performance matrices for the respective measure
#'      of the first-step estimated method.}
#'   \item{auc, aauc, aYI, youden_index, sensitivity,
#'     specificity, geometric_mean, fdr, mcc, corrclass}{A series of
#'     nested lists, each containing performance matrices for the 
#'     respective measure of the final estimated method.}
#'   \item{n_betas}{A matrix showing the number of non-zero beta
#'     coefficients for each fold and lambda value.}
#'   \item{n_gammas}{A list of matrices showing the number of
#'     non-zero gamma coefficients for each fold, lambda, and tau value.}
#'   \item{betas}{A list of the estimated beta coefficients for
#'     each fold.}
#'   \item{betas_star}{A list of the omptimal estimated beta
#'     coefficients used for the estimation of the gammas, for each fold.}
#'   \item{gammas}{A list of the estimated gamma coefficients for
#'     each fold.}
#'   \item{lambda_hat_*}{A vector containing the optimal lambda
#'     value selected by each performance measure.}
#'   \item{tau_hat_*}{A vector containing the optimal tau value
#'     selected by each performance measure (only if `c_function_of_covariates`
#'     is TRUE).}
#'
#' @examples
#' # 1. Simulate a small dataset
#' set.seed(123)
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_betas <- sim_data$nregressors
#' regressors_gammas <- sim_data$ncovariates
#'
#' # 2. Define a small grid for hyperparameters
#' # In a real scenario, these would be calibrated (e.g., using calibrate_lambda_max)
#' # and the sequence would be longer.
#' lambda_seq <- c(0.1, 0.05)
#' tau_seq <- c(0.1, 0.05)
#'
#' # 3. Run cross-validation
#' cv_results <- AucPR_compute_cv(alpha = 0.5,
#'                               trace = 1, n_folds = 2,
#'                               df = df, X = X, y = y, C = C, 
#'                               lambda = lambda_seq, tau = tau_seq,
#'                               regressors_betas = regressors_betas,
#'                               regressors_gammas = regressors_gammas,
#'                               c_function_of_covariates = TRUE,
#'                               simultaneous = TRUE,
#'                               measure_to_select_lambda = "ccr",
#'                               penalty_g = "L1", max_iter_g = 5)
#'
#' # 4. Inspect the results
#' cat("Cross-validation finished in:", round(cv_results$cv_time, 2), "minutes\n")
#' cat("Optimal lambda based on AUC:", cv_results$lambda_hat_auc, "\n")
#' cat("Optimal tau based on CCR:", cv_results$tau_hat_ccr, "\n")
#'
#' @export
AucPR_compute_cv <- function (n_folds, df,
                              X = NULL,
                              y = "y",
                              C = NULL,
                              alpha = 0.5,
                              lambda, tau = 0,
															w_g = 0.5,
                              c_to_use = "Youden",
                              regressors_betas = NULL, trace = 1,
                              seed = 1, used_cores = 1, scaling = FALSE,
                              c_function_of_covariates = FALSE,
                              simultaneous = FALSE,
                              measure_to_select_lambda = "ccr",
                              alpha_g = 0.5, penalty_g = "L1",
                              kernel_g = "gaussian",
                              a1_g = 3.7, a2_g = 3,
                              trend_g = "monotone",
                              gamma_start_input = NULL,
                              gamma_start_default = "zeros",
                              regressors_gammas = NULL,
                              max_iter_g = 10000,
                              delta_g = 1e-5, max_alpha_g = 10000,
                              stepsizeShrink_g = 0.8,
                              min_alpha_g = 1e-12,
                              convergence_error_g = 1e-7,
                              run_aauc = FALSE,
                              log_file = "log_AUC_models.txt") {

  # Start calculation of estimation time
  start_time <- Sys.time()

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")
  if (!is.numeric(n_folds) || length(n_folds) != 1 || n_folds < 2 || floor(n_folds) != n_folds) {stop("'n_folds' must be an single integer >= 2.")}

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

  tau_initial <- 0

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
    if (simultaneous == FALSE) {
      c_function_of_covariates <- FALSE
      tau_initial <- tau #if tau_initial != 0 then c_function_of_covariates = TRUE, and we have to re-change it later
      tau <- 0
    }
  } else { # If c_function_of_covariates is FALSE, tau is irrelevant
    tau <- 0
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
  # Validate alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha > 1) {stop("The parameter 'alpha' must be a single numeric value between 0 (exclusive) and 1 (inclusive).")}
  model_type <- if (alpha == 1) "AucPR_L1" else "AucPR_EN"
  if (!is.numeric(lambda) || length(lambda) < 1 || any(lambda < 0)) stop("Parameter 'lambda' must be a numeric vector of non-negative values.")

  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty_g %in% valid_penalties)) {stop("Parameter 'penalty_g' must be one of: ", paste(valid_penalties, collapse = ", "))}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  valid_measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
  if (!(measure_to_select_lambda %in% valid_measures)) {stop("Parameter 'measure_to_select_lambda' must be one of: ", paste(valid_measures, collapse = ", "))}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}
  if (!is.numeric(max_iter_g) || length(max_iter_g) != 1 || max_iter_g < 2 || floor(max_iter_g) != max_iter_g) {stop("Parameter 'max_iter_g' must be an integer greater than or equal to 2.")}
  if (!is.numeric(delta_g) || length(delta_g) != 1 || delta_g <= 0) {stop("Parameter 'delta_g' must be a single positive numeric value.")}
  if (!is.numeric(max_alpha_g) || length(max_alpha_g) != 1 || max_alpha_g <= 0) {stop("Parameter 'max_alpha_g' must be a single positive numeric value.")}
  if (!is.numeric(stepsizeShrink_g) || length(stepsizeShrink_g) != 1 || stepsizeShrink_g <= 0 || stepsizeShrink_g >= 1) {stop("Parameter 'stepsizeShrink_g' must be a single numeric value between 0 (exclusive) and 1 (exclusive).")}
  if (!is.numeric(min_alpha_g) || length(min_alpha_g) != 1 || min_alpha_g <= 0) {stop("Parameter 'min_alpha_g' must be a single positive numeric value.")}
  if (!is.numeric(convergence_error_g) || length(convergence_error_g) != 1 || convergence_error_g <= 0) {stop("Parameter 'convergence_error_g' must be a single positive numeric value.")}
  if (!is.logical(scaling)) {stop("Parameter 'scaling' must be a logical (TRUE/FALSE).")}
  if (!is.logical(simultaneous)) {stop("Parameter 'simultaneous' must be a logical (TRUE/FALSE).")}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}
  if (!(is.character(c_to_use) && c_to_use == "Youden") && !(is.numeric(c_to_use) && length(c_to_use) == 1 && !is.na(c_to_use))) stop("Parameter 'c_to_use' must be either the string 'Youden' or a single numeric value.")

  # standardize df1
  if (scaling == TRUE) {
    df1 <- scaling_df_for_pye (df = df1, X = colnames(df1[, names(df1) %in% c(X, C)]), y = "y")$df_scaled
  }

  set.seed(seed)
  # check if df1 is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1[[y]][df1[[y]] == 1]) < 2 * n_folds)) {stop("df1 contains too few 1s for this number of folds")
  } else if ((length(df1[[y]][df1[[y]] == 0]) < 2 * n_folds)) {stop("df1 contains too few 0s for this number of folds")}

  # Divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1, y)), 1:2]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 0, ])), replace = FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 1, ])), replace = FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0, fold_i_1))
  folds_i <- merge(df1[, 1:2], df_sort[, 1:3], by = 'ID', all = FALSE, sort = FALSE)[, 4]

  # Names of the columns
  lambdanames <- paste("lambda", lambda, sep = "=")
  taunames <- paste("tau", tau, sep = "=")
  foldnames <- paste("fold", 1:n_folds, sep = "=")
  # Accuracy measures (train and test)
  list_of_measures <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                        "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

  auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures) {
    assign(mes, lapply(lambda, function(x) list (
      train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)),
      test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames))
      )))
    eval(substitute(names(x) <- lambdanames, list(x = as.symbol(mes))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(foldnames, lambdanames))
  n_gammas <- lapply(lambda, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
  names(n_gammas) <- lambdanames

  if (trace > 0) {
    cat("Starting CV with model_type:", model_type, "\n")
  }

  # fill the matrices
  results <- mapply(function(k) AucPR.cv(df = df1[, names(df1) != "ID", drop = FALSE],
                                          X = X, y = y, C = C, alpha = alpha,
                                          lambda = lambda, tau = tau, w_g = w_g,
																					c_to_use = c_to_use, alpha_g = alpha_g,
                                          penalty_g = penalty_g, folds_i = folds_i, k = k,
                                          regressors_betas = regressors_betas,
                                          trace = trace, model_type = model_type,
                                          auc = auc, aauc = aauc, aYI = aYI,
                                          youden_index = youden_index,
                                          sensitivity = sensitivity,
                                          specificity = specificity,
                                          geometric_mean = geometric_mean, fdr = fdr,
                                          mcc = mcc, corrclass = corrclass,
                                          used_cores = used_cores,
                                          n_betas = n_betas,
                                          n_gammas = n_gammas,
                                          c_function_of_covariates = c_function_of_covariates,
                                          simultaneous = simultaneous,
                                          kernel_g = kernel_g, a1_g = a1_g, a2_g = a2_g,
                                          trend_g = trend_g,
                                          gamma_start_input = gamma_start_input,
                                          gamma_start_default = gamma_start_default,
                                          regressors_gammas = regressors_gammas,
                                          max_iter_g = max_iter_g, delta_g = delta_g,
                                          max_alpha_g = max_alpha_g,
                                          stepsizeShrink_g = stepsizeShrink_g,
                                          min_alpha_g = min_alpha_g,
                                          convergence_error_g = convergence_error_g,
                                          run_aauc = run_aauc, log_file = log_file),
                                        seq(1:n_folds))

  # Prepare the result
  wrapper <- function(results, mes, i, dataset, tau) {
    m <- as.matrix(sapply(c(1:n_folds), function(k) getElement(getElement(results[[k]], mes)[[i]], dataset)[k, , drop = FALSE]))
    if (length(tau) > 1) {
      m2 <- t(m)
    } else{
      m2 <- m
    }
    colnames(m2) <- taunames
    rownames(m2) <- foldnames
    return(m2)
  }

  for (mes in list_of_measures) {
    assign(mes[], lapply(seq_along(lambda), function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
    eval(substitute(names(x) <- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep = "="))), list(x = as.symbol(mes))))
  }

  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ])))
  n_gammas[] <- lapply(seq_along(lambda), function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ]))))

  betas <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas)))
  gammas <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

  if (length(lambda) == 1) {
    n_betas[] <- t(t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ]))))
    n_gammas[] <- lapply(seq_along(lambda), function(i) t(t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ])))))

    betas <- t(t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas))))
    gammas <- t(t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas))))
  }

  rownames(betas) <- rownames(gammas) <- foldnames

  measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
  list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                         "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

  tau_hat_auc <- tau_hat_aauc <- tau_hat_aYI <- tau_hat_yi <- tau_hat_sen <- NULL
  tau_hat_spc <- tau_hat_gm <- tau_hat_fdr <- tau_hat_mcc <- tau_hat_ccr <- NULL
  lambda_hat_auc <- lambda_hat_aauc <- lambda_hat_aYI <- lambda_hat_yi <- lambda_hat_sen <- NULL
  lambda_hat_spc <- lambda_hat_gm <- lambda_hat_fdr <- lambda_hat_mcc <- lambda_hat_ccr <- NULL
  auc_first_step <- aauc_first_step <- aYI_first_step <- youden_index_first_step <- sensitivity_first_step <- NULL
  specificity_first_step <- geometric_mean_first_step <- fdr_first_step <- mcc_first_step <- corrclass_first_step <- NULL
  for (i in seq_along(measures)) {
    #in general, create tau_hat equal to NA
    assign(paste0("tau_hat_", measures[i]), NA)
		#save the current table measures with another name
		assign(paste0(list_of_measures2[i], "_first_step"), get(list_of_measures2[i]))

    if (length(tau) != 1) {
      measures_matrix <- t(sapply(seq_along(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    } else {
      measures_matrix <- as.matrix(sapply(seq_along(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    }

    rownames(measures_matrix) <- lambdanames
    colnames(measures_matrix) <- taunames

    if (!all(is.na(measures_matrix))) {
      max_measures <- which(measures_matrix == max(measures_matrix, na.rm = TRUE), arr.ind = TRUE)[1, ]
      assign(paste0("lambda_hat_", measures[i]), lambda[max_measures[1]])

      if (c_function_of_covariates == TRUE) {
        if (length(tau) > 1) {
          assign(paste0("tau_hat_", measures[i]), tau[max_measures[2]])
        } else {
          assign(paste0("tau_hat_", measures[i]), tau)
        }
      }
    } else {
      assign(paste0("lambda_hat_", measures[i]), NA)
    }
  }

  n_betas_star <- NULL
  betas_star <- NULL
  lambda_star <- NA

  if (any(tau_initial != 0)) {
    c_function_of_covariates <- TRUE
    tau <- tau_initial

    # now we execute AucPR.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
    lambda_star <- get(paste0("lambda_hat_", measure_to_select_lambda))
    if (is.na(lambda_star)) {stop("Optimal lambda for sequential search was NA. Skipping tau search.")}

    # Re-initialize result structures for the second run
    taunames <- paste("tau", tau, sep = "=")
    lambdanames_star <- paste("lambda", lambda_star, sep = "=")
    # Accuracy measures (train and test)
    list_of_measures <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                          "fdr", "mcc", "corrclass")
    auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL

    for (mes in list_of_measures) {
      #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(foldnames, taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(foldnames, taunames)))
      assign(mes, lapply(lambda_star, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)),
                                                         test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))))
      eval(substitute(names(x) <- lambdanames_star, list(x = as.symbol(mes))))
    }

    n_betas_star <- matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(foldnames, lambdanames_star))
    n_gammas <- lapply(lambda_star, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
    names(n_gammas) <- lambdanames_star

    if (trace > 0) {
      cat(" \n Starting the CV of tau using as lambda:", lambda_star, ", that is the best value of lambda as per:", measure_to_select_lambda, "\n")
    }

    # Re-run CV with optimal lambda and full tau grid
    results <- mapply(function(k) AucPR.cv(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, C = C,
                                           alpha = alpha, 
																					 lambda = lambda_star, tau = tau, w_g = w_g,
                                           c_to_use = c_to_use, alpha_g = alpha_g,
                                           penalty_g = penalty_g, folds_i = folds_i, k = k,
                                           regressors_betas = regressors_betas,
                                           trace = trace, model_type = model_type,
                                           auc = auc, aauc = aauc, aYI = aYI,
                                           youden_index = youden_index,
                                           sensitivity = sensitivity,
                                           specificity = specificity,
                                           geometric_mean = geometric_mean, fdr = fdr,
                                           mcc = mcc, corrclass = corrclass,
                                           used_cores = used_cores,
                                           n_betas = n_betas_star, n_gammas = n_gammas,
                                           c_function_of_covariates = c_function_of_covariates,
                                           simultaneous = simultaneous,
                                           kernel_g = kernel_g, a1_g = a1_g, a2_g = a2_g,
                                           trend_g = trend_g,
                                           gamma_start_input = gamma_start_input,
                                           gamma_start_default = gamma_start_default,
                                           regressors_gammas = regressors_gammas,
                                           max_iter_g = max_iter_g, delta_g = delta_g,
                                           max_alpha_g = max_alpha_g,
                                           stepsizeShrink_g = stepsizeShrink_g,
                                           min_alpha_g = min_alpha_g,
                                           convergence_error_g = convergence_error_g,
                                           run_aauc = run_aauc, log_file = log_file),
                      seq(1:n_folds))

    # Prepare the result
    wrapper <- function(results, mes, i, dataset, tau) {
      m <- as.matrix(sapply(c(1:n_folds), function(k) getElement(getElement(results[[k]], mes)[[i]], dataset)[k, , drop = FALSE]))
      if (length(tau) > 1) {
        m2 <- t(m)
      } else{
        m2 <- m
      }
      colnames(m2) <- taunames
    rownames(m2) <- foldnames
      return(m2)
    }

    for (mes in list_of_measures) {
      assign(mes[], lapply(1, function(i) list(train = wrapper(results, mes, i, "train", tau),
                                               test = wrapper(results, mes, i, "test", tau))))
      eval(substitute(names(x) <- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep = "="))), list(x = as.symbol(mes))))
    }

    n_betas_star[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ])))
    n_gammas[] <- lapply(1, function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ]))))

    betas_star <- matrix(sapply(c(1:n_folds), function(k) results[[k]]@betas), byrow = TRUE, dimnames = list(foldnames, lambdanames_star))
    rownames(betas_star) <- foldnames
    gammas <- lapply(1, function(i) { x <- t(sapply(c(1:n_folds), function(k) results[[k]]@gammas[[i]])); rownames(x) <- foldnames; x})
    names(gammas) <- lambdanames_star

    measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
    list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                           "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

    for (i in seq_along(measures)) {
      #in general, create tau_hat equal to NA
      assign(paste0("tau_hat_", measures[i]), NA)

      if (length(tau) != 1) {
        measures_matrix <- t(sapply(1, function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
      } else {
        measures_matrix <- as.matrix(sapply(1, function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
      }

      rownames(measures_matrix) <- paste0("lambda=", lambda_star)
      colnames(measures_matrix) <- taunames

      if (!all(is.na(measures_matrix))) {
        #if the number of gammas is zero, we set the measure to zero
        mean_n_gammas <- t(sapply(n_gammas, function(i) colMeans(i)))
        measures_matrix[mean_n_gammas == 0] <- 0

        max_measures <- which(measures_matrix == max(measures_matrix, na.rm = TRUE), arr.ind = TRUE)[1, ]

        if (c_function_of_covariates == TRUE) {
          if (length(tau) > 1) {
            assign(paste0("tau_hat_", measures[i]), tau[max_measures[2]])
          } else {
            assign(paste0("tau_hat_", measures[i]), tau)
          }
        }
      }
    }
  }

  cv_time <- difftime(Sys.time(), start_time, units = "mins")

  if (trace %in% c(1, 2)) {
    cat("-------------------> END OF THE CROSS-VALIDATION OF THE AucPR METHOD <------------------- \n")
    cat("-------------> For the whoole Cross-validation it took:", cv_time, "minutes <------------- \n")
  }

  return(list(model_type = model_type,
              c_to_use = c_to_use,
              alpha = alpha,
              penalty_g = penalty_g,
							w_g = w_g,
              cv_time = cv_time,
              auc = auc,
              aauc = aauc,
              aYI = aYI,
              youden_index = youden_index,
              sensitivity = sensitivity,
              specificity = specificity,
              geometric_mean = geometric_mean,
              fdr = fdr,
              mcc = mcc,
              corrclass = corrclass,
              auc_first_step = auc_first_step,
              aauc_first_step = aauc_first_step,
              aYI_first_step = aYI_first_step,
              youden_index_first_step = youden_index_first_step,
              sensitivity_first_step = sensitivity_first_step,
              specificity_first_step = specificity_first_step,
              geometric_mean_first_step = geometric_mean_first_step,
              fdr_first_step = fdr_first_step,
              mcc_first_step = mcc_first_step,
              corrclass_first_step = corrclass_first_step,
              lambda_hat_yi = lambda_hat_yi,
              lambda_hat_auc = lambda_hat_auc,
              lambda_hat_aauc = lambda_hat_aauc,
              lambda_hat_aYI = lambda_hat_aYI,
              lambda_hat_ccr = lambda_hat_ccr,
              lambda_hat_sen = lambda_hat_sen,
              lambda_hat_spc = lambda_hat_spc,
              lambda_hat_gm = lambda_hat_gm,
              tau_hat_yi = tau_hat_yi,
              tau_hat_auc = tau_hat_auc,
              tau_hat_aauc = tau_hat_aauc,
              tau_hat_aYI = tau_hat_aYI,
              tau_hat_ccr = tau_hat_ccr,
              tau_hat_sen = tau_hat_sen,
              tau_hat_spc = tau_hat_spc,
              tau_hat_gm = tau_hat_gm,
              c_function_of_covariates = c_function_of_covariates,
              simultaneous = simultaneous,
              measure_to_select_lambda = measure_to_select_lambda,
              lambda_star = lambda_star,
              n_betas = n_betas,
              n_betas_star = n_betas_star,
              n_gammas = n_gammas,
              betas = betas,
              betas_star = betas_star,
              gammas = gammas))
}
