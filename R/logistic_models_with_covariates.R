#' @title glmnet Estimation for Coefficient and Feature Selection
#'
#' @description function to estimate the optimal regression coefficients
#' (betas) using penalized logistic regression methods such as Lasso,
#' Elastic-Net, SCAD, and MCP. This function fits a penalized logistic
#' regression model for a binary outcome.
#' The user chooses the type of penalization (L1, Elastic-Net,
#' SCAD, or MCP) and provides the fixed penalization parameter, \code{lambda}.
#' It returns the estimated coefficients and a set of comprehensive
#' classification performance metrics based on the optimal cut-point derived
#' from Youden's Index.
#'
#' @param df A data frame containing the complete dataset.
#' @param X A character vector of column names from \code{df} to be used as
#'   regressor variables. Alternatively, a data frame containing only the
#'   regressors. If not specified, all columns not in \code{y} are used.
#' @param y A character string specifying the column name of the target
#'   variable in \code{df}. It must be a binomial variable with values 0 or 1.
#'   Alternatively, a data frame containing only the target variable. Default
#'   is "y".
#' @param alpha The elastic-net mixing parameter. A value of 1 corresponds
#'   to L1 penalization (Lasso). A value between 0 and 1 is used for
#'   Elastic-Net. Note: This parameter is only used when
#'   \code{model_type = "logElasticNet"}; for \code{logLasso}, \code{alpha} is
#'   forced to 1, and for \code{logSCAD}/\code{logMCP}, it is ignored as they
#'   do not typically use an \code{alpha} mixing parameter. Default is 0.5.
#' @param lambda A single non-negative numeric value for the regularization
#'   parameter. A larger value results in stronger penalization.
#' @param model_type A character string specifying the penalized logistic
#'   model to use. Options are: "logLasso" (L1 penalty, via \code{glmnet}),
#'   "logElasticNet" (Elastic-Net penalty, via \code{glmnet}), "logSCAD" (SCAD
#'   penalty, via \code{ncvreg_modified} in \code{pye}), and "logMCP" (MCP
#'   penalty, via \code{ncvreg_modified} in \code{pye}). Default is "logLasso".
#' @param c_to_use character or numeric value specifying the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If set
#'   to `"Youden"`, the function computes the optimal cut-point using Youden's Index.
#'   If a numeric value is provided, it is used directly as the threshold.
#'   Default is "Youden".
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param trace An integer controlling the amount of output printed to the
#'   console. \code{0} for no output, \code{1} for a result summary, and
#'   \code{2} for a detailed breakdown. Default is 1.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#'
#' @return A list containing the model fit, estimated coefficients, and a
#'   collection of performance metrics.
#'   \item{model}{The fitted model object (from \code{glmnet} or
#'     \code{ncvreg_modified}).}
#'   \item{model_type}{A character string indicating the model type used.}
#'   \item{betas_hat}{The estimated vector of optimal beta coefficients
#'     (excluding the intercept).}
#'   \item{X_model}{A character vector containing the names of all regressor
#'   variables (\code{X}) used in the model estimation. This vector defines
#'   the canonical order of the regressors used to align the coefficients
#'   in \code{betas_hat}.}
#'   \item{youden_index}{The Youden Index value at the optimal cut-point.}
#'   \item{sensitivity}{The sensitivity (True Positive Rate) at the optimal
#'     cut-point.}
#'   \item{specificity}{The specificity (True Negative Rate) at the optimal
#'     cut-point (\code{spec}).}
#'   \item{geometric_mean}{The geometric mean of sensitivity and specificity
#'     (\code{gm}).}
#'   \item{fdr}{The False Discovery Rate at the optimal cut-point.}
#'   \item{mcc}{The Matthews Correlation Coefficient.}
#'   \item{corrclass}{The overall correct classification rate (CCR).}
#'   \item{auc}{The Area Under the ROC Curve (AUC).}
#'   \item{lambda}{The penalization parameter used.}
#'   \item{c_hat}{The optimal cut-point for the \code{z_hat} score.}
#'   \item{z_hat}{A data frame with the estimated linear combination of
#'     regressors (or probability) for each observation.}
#'   \item{y_hat}{A data frame with the predicted binary outcomes (0 or 1).}
#'   \item{n_betas}{The number of non-zero beta coefficients estimated.}
#'   \item{n_total_var}{The total number of regressors considered.}
#'   \item{n_predicted_zeros}{The number of estimated betas that are zero.}
#'   \item{n_predicted_non_zeros}{The number of estimated betas that are
#'     non-zero.}
#'   \item{n_caught_betas}{The number of true non-zero betas correctly
#'     identified as non-zero (if \code{regressors_betas} is provided).}
#'   \item{n_non_caught_betas}{The number of true non-zero betas incorrectly
#'     identified as zero (False Negatives).}
#'   \item{n_caught_zero}{The number of true zero betas correctly identified
#'     as zero (True Negatives).}
#'   \item{n_zero_not_caught}{The number of true zero betas incorrectly
#'     identified as non-zero (False Positives).}
#'   \item{TP, TN, FP, FN}{The counts of True Positives, True Negatives,
#'     False Positives, and False Negatives at the optimal cut-point.}
#'   \item{estimation_time}{The time taken for the estimation, in minutes.}
#'
#' @examples
#' library(pye)
#' # --- 1. Simulate data ---
#' # Simulate a small dataset with 50 observations, 200 total features,
#' # and 20 covariates (not used here, but good for context).
#' set.seed(42) # for reproducibility
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled # scaled training data
#' X <- sim_data$X                  # names of regressors
#' y <- sim_data$y                 # name of target variable
#' regressors_betas <- sim_data$nregressors # true non-zero betas count
#'
#' # --- 2. Estimate using logLasso ---
#' model_type_lasso <- "logLasso"
#' lambda_lasso <- 0.05
#'
#' lasso_result <- plr_estimation(df = df, X = X, y = y,
#'   model_type = model_type_lasso, lambda = lambda_lasso, alpha = 1,
#'   regressors_betas = regressors_betas, trace = 1)
#'
#' print(lasso_result$model_type)
#' print(lasso_result$corrclass)
#' print(lasso_result$n_betas)
#'
#' # --- 3. Estimate using logElasticNet ---
#' en_result <- plr_estimation(df = df, X = X, y = y,
#'   model_type = "logElasticNet", lambda = 0.05, alpha = 0.5,
#'   regressors_betas = regressors_betas, trace=0)
#'
#' print(en_result$model_type)
#' print(en_result$auc)
#'
#' @importFrom glmnet glmnet
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom Rmpfr mpfr
#' @importFrom stats coef predict setNames
#' @export
plr_estimation <- function (df,
                               X = NULL,
                               y = "y",
                               alpha = 0.5,
                               lambda,
                               model_type = "logLasso",
                               c_to_use = "Youden",
                               regressors_betas = NULL,
                               fold = NULL,
                               trace = 1,
                               max.print = 10) {

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
  #df1 <- cbind(ID, df[, (names(df) %in% c(y, X))]) #OLD
  df1 <- cbind(ID, df[, c(y, X)], drop = FALSE)


  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Additional parameter validation
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!is.null(fold) && (!is.numeric(fold) || length(fold) != 1)) {stop("The parameter 'fold' must be a single numeric value or NULL.")}
  if (!is.null(regressors_betas) && !is.numeric(regressors_betas)) {stop("The parameter 'regressors_betas' must be a numeric vector or NULL.")}
  if (!(is.character(c_to_use) && c_to_use == "Youden") && !(is.numeric(c_to_use) && length(c_to_use) == 1 && !is.na(c_to_use))) stop("Parameter 'c_to_use' must be either the string 'Youden' or a single numeric value.")

  # --- Model Type Validation ---
  valid_models <- c("logLasso", "logElasticNet", "logSCAD", "logMCP")
  if (!(model_type %in% valid_models)) {
    stop("The parameter 'model_type' must be one of: ", paste(valid_models, collapse = ", "))
  }
  # --- Alpha Validation and Standardization based on Model Type ---
  # 1. Check basic alpha properties (must be a single number/NA)
  if (!is.numeric(alpha) || length(alpha) != 1) {
    if (!(model_type %in% c("logSCAD", "logMCP") && is.na(alpha))) {
      stop("The parameter 'alpha' must be a single numeric value (or NA for SCAD/MCP).")
    }
  }
  if (model_type == "logLasso") {
    # Lasso is strictly L1, so alpha must be 1.
    if (alpha != 1) {
      if (trace %in% c(2)) {
        warning("For 'logLasso', 'alpha' is forced to 1, overriding the input value.")
      }
      alpha <- 1
    }
  } else if (model_type == "logElasticNet") {
    # Elastic-Net requires 0 < alpha < 1. Check range.
    if (alpha <= 0 || alpha >= 1) stop("For 'logElasticNet', 'alpha' must be a single numeric value between 0 (exclusive) and 1 (exclusive).")
  } else if (model_type %in% c("logSCAD", "logMCP")) {
    # SCAD/MCP do not use alpha mixing parameter. Set to NA for clarity in output.
    if (trace %in% c(2) && !is.na(alpha) && alpha!=0) {
      warning(paste0("For '", model_type, "', 'alpha' is not used for penalty mixing and is set to 0."))
			alpha <- 0
    }
  }


  # Prepare data for modeling functions
  df_x_as_matrix <- as.matrix(df1[, X, drop = FALSE])
  df_y_as_matrix <- as.matrix(df1[[y]])

  # Suppress warnings during model fitting
  suppressWarnings({

    # Fitting logistic lasso
    if (model_type == "logLasso") {

      # For Lasso, alpha must be 1
      alpha <- 1

      # Estimate the Lasso model
      fit_lasso <- glmnet::glmnet(df_x_as_matrix, df_y_as_matrix, family = "binomial", alpha = alpha, lambda = lambda, intercept = TRUE)

      # Extract coefficients (excluding intercept)
      betas_hat <- as.numeric(stats::coef (fit_lasso))[-1]
      names(betas_hat) <- X

      # Compute z_hat (predicted probabilities)
      z_hat <- stats::predict(fit_lasso, newx = df_x_as_matrix, type = "response", s = lambda) #type = "response", "class", "coefficients", "nonzero")
      colnames(z_hat) <- "z_hat"
      df_hat <- cbind(df1, z_hat)
      model <- fit_lasso

    } else if (model_type == "logElasticNet") {

      # Estimate the Elastic-Net model
      fit_EN <- glmnet::glmnet(df_x_as_matrix, df_y_as_matrix, family = "binomial", alpha = alpha, lambda = lambda, intercept = TRUE)

      # Extract coefficients (excluding intercept)
      betas_hat <- as.numeric(stats::coef (fit_EN))[-1]
      names(betas_hat) <- X

      # Compute z_hat (predicted probabilities)
      z_hat <- stats::predict(fit_EN, newx = df_x_as_matrix, type = "response", s = lambda) #type = "response", "class", "coefficients", "nonzero")
      colnames(z_hat) <- "z_hat"
      df_hat <- cbind(df1, z_hat)
      model <- fit_EN

    } else if (model_type %in% c("logSCAD", "logMCP")) {

      penalty_type <- gsub("log", "", model_type)
      eps <- 1e-7 # Initial convergence tolerance

      # Estimation using ncvreg_modified
      tmp <- try({
        fit_SCAD_or_MCP <- ncvreg_modified(X = df_x_as_matrix, #matrix
                                    y = df_y_as_matrix, #vector
                                    family = "binomial", #"gaussian", "binomial", "poisson"
                                    penalty = penalty_type, #"MCP", "SCAD", "lasso"
                                    #gamma=switch(penalty, SCAD=3.7, 3),
                                    alpha = 1,
                                    max.iter = 1000000,
                                    lambda = lambda,
                                    eps = eps,
                                    #dfmax = p+1,
                                    warn = FALSE)
        }, silent = TRUE)
      #if it produces an error it is because the algorithm do not converge
      #the only solution is to decrease eps
      while ((inherits(tmp, "try-error")) && (eps < 1)) {
        eps <- eps * 5
        tmp <- try({
          fit_SCAD_or_MCP <- ncvreg_modified(X = df_x_as_matrix, #matrix
                                             y = df_y_as_matrix, #vector
                                             family = "binomial", #"gaussian", "binomial", "poisson"
                                             penalty = penalty_type, #"MCP", "SCAD", "lasso"
                                             #gamma=switch(penalty, SCAD=3.7, 3),
                                             alpha = 1,
                                             max.iter = 100000,
                                             lambda = lambda, eps = eps,
                                             #dfmax = p+1,
                                             warn = FALSE)
          }, silent = TRUE)
      }
      if (inherits(tmp, "try-error") || eps > 1) {
        stop(paste0("The algorithm for ", model_type, " did not converge."))
      }

      # Extract coefficients (excluding intercept)
      betas_hat <- as.numeric(stats::coef (fit_SCAD_or_MCP))[-1]
      names(betas_hat) <- X

      # Compute z_hat (predicted probabilities)
      # Linear predictor: X*beta + intercept
      eta <-  Rmpfr::mpfr(sweep(df_x_as_matrix %*% fit_SCAD_or_MCP$beta[-1, drop = FALSE], 2, fit_SCAD_or_MCP$beta[1], "+"), precBits = 120)
      # Predicted probability: exp(eta) / (1 + exp(eta))
      z_hat <- as.matrix(as.numeric(exp(eta) / (1 + exp(eta))))
      colnames(z_hat) <- "z_hat"
      df_hat <- cbind(df1, z_hat)
      model <- fit_SCAD_or_MCP
      # For SCAD/MCP, we set alpha back to NA/NULL as it's not a primary parameter
      alpha <- NA

    } else {
      stop("The parameter model_type is not valid!")
    }
  })

  # Find optimal measures, some of them based on the Youden's Index
  opt <- OptimalCutpoints::optimal.cutpoints(data = df_hat, X = "z_hat", status = y, methods = "Youden", tag.healthy = 0)

  # Which cut-point to use?
  if (c_to_use == "Youden") {
    # Select the optimal cut-point (using the mean if multiple are found)
    c <- mean(opt$Youden$Global$optimal.cutoff$cutoff)
  } else {
    c <- c_to_use
  }

  # Predict binary outcomes based on the optimal cut-point
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

  # Prepare the output
  c_hat <- c
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  #if all the betas are zero, the measures are zeros
  if (sum(betas_hat) == 0) {
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    gm <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  if (trace %in% c(1, 2)) {
    cat("Estimation done using the Penalized Logistic method type:", model_type, "; \n")
    cat("-> ")
    if (!is.null(fold)) {cat("fold =", fold, "; ")}
    if (!is.na(lambda)) {cat("lambda:", lambda, "; ")}
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
  if (trace %in% c(1, 2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  return(list(model = model, model_type = model_type,
              c_to_use = c_to_use,
              betas_hat = betas_hat,
              X_model = X,
              youden_index = youden_index,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr, mcc = mcc, corrclass = corrclass,
              auc = auc,
              lambda = lambda,
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














#' @title Prediction and Performance Evaluation for Penalized glmnet Models
#'
#' @description This function applies a previously trained penalized logistic
#' model to new data to generate predictions and evaluate performance metrics.
#' The model to be used is the result of the `plr_estimation` function.
#'
#' @param df A data frame containing the new dataset for prediction.
#' @param y A character string specifying the column name of the target
#'   variable in `df`. It must be a binomial variable with values 0 or 1.
#'   Default is "y".
#' @param model_to_use A list containing the parameters from a model trained
#'   by the `plr_estimation` function, specifically the `betas_hat`,
#'   `c_hat`, and `model` components.
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
#' @param c_function_of_covariates A logical flag used internally to suppress
#'   printing of results if the output is handled by an external function
#'   (e.g., `covYI`). Default is `FALSE`.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#'
#' @return A list containing the predictions and a collection of performance
#'   metrics, including:
#'   \item{model}{The trained model object (e.g., from `glmnet` or `ncvreg`).}
#'   \item{model_type}{The model type ("logLasso", "logMCP", etc.).}
#'   \item{betas_hat}{The estimated vector of optimal beta coefficients.}
#'   \item{youden_index}{The Youden Index value at the optimal cut-point.}
#'   \item{sensitivity}{The sensitivity at the optimal cut-point.}
#'   \item{specificity}{The specificity at the optimal cut-point.}
#'   \item{geometric_mean}{The geometric mean of sensitivity and specificity.}
#'   \item{fdr}{The False Discovery Rate at the optimal cut-point.}
#'   \item{mcc}{The Matthews Correlation Coefficient.}
#'   \item{corrclass}{The overall correct classification rate.}
#'   \item{auc}{The Area Under the ROC Curve (AUC).}
#'   \item{lambda}{The penalization parameter used in the estimation.}
#'   \item{c_hat}{The optimal cut-point for the `z_hat` score.}
#'   \item{z_hat}{A data frame with the estimated linear combination of
#'     regressors for each observation.}
#'   \item{y_hat}{A data frame with the predicted binary outcomes (0 or 1).}
#'   \item{n_total_var}{The total number of regressors considered.}
#'   \item{n_predicted_non_zeros}{The number of estimated betas that are
#'     non-zero.}
#'   \item{TP, TN, FP, FN}{The counts from the confusion matrix.}
#'   \item{estimation_time}{The time taken for the prediction.}
#'
#' @examples
#' library(pye)
#' # --- 1. Simulate data ---
#' # Simulate a small dataset with 50 observations, 200 total features,
#' # and 20 covariates (not used here, but good for context).
#' set.seed(42) # for reproducibility
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df_train <- sim_data$train_df_scaled # scaled training data
#' df_test <- sim_data$test_df_scaled # scaled training data
#' X <- sim_data$X                  # names of regressors
#' y <- sim_data$y                  # name of target variable
#' regressors_betas <- sim_data$nregressors # true non-zero betas count
#'
#' # --- 2. Estimate using logLasso ---
#' model <- plr_estimation(df = df_train, X = X, y = y,
#'   model_type = "logMCP", lambda = 0.05,
#'   regressors_betas = regressors_betas, trace = 1)
#'
#' predictions <- plr_predict(df = df_test, y = y, model_to_use = model,
#'    trace = 1, regressors_betas = regressors_betas)
#'
#' print(predictions)
#'
#' @importFrom Rmpfr mpfr
#' @importFrom stats predict
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @export
plr_predict <- function (df,
                          y = "y", model_to_use, fold = NULL,
                          regressors_betas = NULL, trace = 1,
                          c_function_of_covariates = FALSE,
                          max.print = 10) {

  if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
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

  # Checks
  if (!is.logical(c_function_of_covariates)) {stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!is.null(fold) && (!is.numeric(fold) || length(fold) != 1)) {stop("The parameter 'fold' must be a single numeric value or NULL.")}
  if (!is.null(regressors_betas) && !is.numeric(regressors_betas)) {stop("The parameter 'regressors_betas' must be a numeric vector or NULL.")}
  if (!is.list(model_to_use) || !all(c("betas_hat", "c_hat", "model_type") %in% names(model_to_use))) {
    stop("The parameter 'model_to_use' must be a valid list result from the 'plr_estimation' function and contain 'betas_hat', 'c_hat', and 'model_type'.")
  }

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  # Extract model components
  betas_hat <- model_to_use$betas_hat
  c_hat <- model_to_use$c_hat
  model_type <- model_to_use$model_type
  lambda <- model_to_use$lambda
  X <- model_to_use$X_model

  # Ensure all variables required by the model (X) are present in the new data (df)
  if (length(X) == 0) {stop("The model contains no regressors (X). Cannot perform prediction.")}
  if (!all(X %in% names(df))) {
    missing_vars <- X[!X %in% names(df)]
    stop("The following required regressors from the trained model are missing in the new data: ", paste(missing_vars, collapse = ", "))
  }

  ID <- rownames(df)
  # The regressors are explicitly ordered according to X
  df1 <- cbind(ID, df[, c(y, X)], drop = FALSE)

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Prepare data matrices
  df_x_as_matrix <- as.matrix(df1[, X, drop = FALSE])

  # --- Predict z_hat score ---

  #alpha change on the bases of the used model
  alpha <- switch(model_type, logLasso = 1, logElasticNet = 0.5, 0)

  #compute z_hat
  valid_models <- c("logLasso", "logElasticNet", "logSCAD", "logMCP")
  if (model_type %in% c("logLasso", "logElasticNet")) {
    z_hat <- as.matrix(stats::predict(model_to_use$model, df_x_as_matrix, type = "response"))
  } else if (model_type %in% c("logSCAD", "logMCP")) {
    eta <- Rmpfr::mpfr(sweep(df_x_as_matrix %*% model_to_use$model$beta[-1], 2, model_to_use$model$beta[1], "+"), precBits = 120)
    z_hat <- as.matrix(as.numeric(exp(eta) / (1 + exp(eta))))
  } else {stop("The parameter 'model_type' must be one of: ", paste(valid_models, collapse = ", "))}

  colnames(z_hat) <- "z_hat"
  df_hat <- cbind(df1, z_hat)

  # --- Apply the estimated cut-off (c_hat) from the training phase ---
  c <- c_hat

  # Predict binary outcomes
  y_hat <- ifelse(df_hat["z_hat"] < c, 0, 1)
  colnames(y_hat) <- "y_hat"
  df_hat <- cbind(df_hat, y_hat)
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  # --- Calculate performance measures ---
  opt <- OptimalCutpoints::optimal.cutpoints(data = df_hat, X = "z_hat", status = y, methods = "Youden", tag.healthy = 0)
  auc <-  mean(opt$Youden$Global$measures.acc$AUC)
  youden_index <- mean(opt$Youden$Global$optimal.criterion)

  #confusion matix
  TP <- sum(ifelse(df_hat[y] == 1 & df_hat["z_hat"] >= c, 1, 0))
  TN <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] < c, 1, 0))
  FP <- sum(ifelse(df_hat[y] == 0 & df_hat["z_hat"] >= c, 1, 0))
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

  # if all the betas are zero, the measures are zeros
  if (sum(betas_hat) == 0) {
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    gm <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  if (trace %in% c(1, 2)) {
    #print the results only if c_function_of_covariates = FALSE
    if (c_function_of_covariates == FALSE) {
      cat("Prediction executed using the Penalized Logistic model type:", model_type, "; \n")
      cat("-> ")
      if (!is.null(fold)) {cat("fold =", fold, "; ")}
      if (!is.na(lambda)) {cat("lambda:", lambda, "; ")}
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

  return(list(model = model_to_use$model,
              model_type = model_type,
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



# Create S4 class for glmnet cross-validation output
setClass(Class = "plr_cross_validation_output",
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

#' Cross-validation helper for glmnet-based fits (robust, CRAN-ready)
#'
#' Internal helper used by plr_compute_cv. Fits glmnet models for each
#' lambda on the training fold, optionally fits covYI for each (lambda, tau),
#' and collects train/test performance measures for the current fold `k`.
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom methods new
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext file_ext
#' @noRd
#' @keywords internal
plr.cv <- function (df, X, y, C, lambda, tau,
                    w_g, c_to_use, alpha,
                    alpha_g, penalty_g,
										folds_i, k,
                    regressors_betas,
										model_type,
                    kernel_g, a1_g, a2_g, trend_g,
                    gamma_start_input,
										gamma_start_default,
                    regressors_gammas,
                    max_iter_g, delta_g, max_alpha_g,
                    stepsizeShrink_g,
										min_alpha_g, convergence_error_g,
                    auc, aauc, aYI,
										youden_index, sensitivity,
                    specificity, geometric_mean,
										fdr, mcc, corrclass,
                    n_betas, n_gammas, trace, used_cores,
                    c_function_of_covariates,
                    simultaneous, run_aauc, log_file) {

  # Basic validations
  if (!is.data.frame(df) || nrow(df) == 0L) {stop("`df` must be a non-empty data.frame.")}
  if (!is.character(X) || length(X) == 0L) stop("`X` must be a character vector with at least one name.")
  if (!is.character(y) || length(y) != 1L) stop("`y` must be a single column name (character).")
  if (!all(c(y, X) %in% names(df))) stop("Some columns in `X`/`y` are not present in `df`.")
  if (!is.numeric(lambda) || length(lambda) < 1L) stop("`lambda` must be a numeric vector of length >= 1.")
  if ((model_type != "logSCAD" && model_type != "logMCP") && (!is.numeric(alpha) || length(alpha) != 1L)) stop("`alpha` must be a single numeric value.")
  if (!is.numeric(used_cores) || length(used_cores) != 1L || used_cores < 1L) stop("`used_cores` must be a positive integer.")
  if (!is.numeric(k) || length(k) != 1L) stop("`k` must be a single integer specifying the fold index.")

  test_i <- which(folds_i == k)
  if (length(test_i) == 0L) stop("No test indices found for the specified fold `k`.")
  train_df <- df[-test_i, , drop = FALSE]
  test_df <- df[test_i, , drop = FALSE]

  greek <- if (c_function_of_covariates) "tau" else "lambda"
  if (trace %in% c(1, 2)) {
    cat("----------------------------------------------------------------\n")
    cat("|       starting with the", k, "-th fold for the CV of ", greek, "       |\n")
    cat("----------------------------------------------------------------\n")
  }

  #if ((c_function_of_covariates == TRUE) & (length(lambda) == 1)) {
  #  #in this case we have to estimate the combination z_hat over all the dataset
  #  train_df1 <- df
  #} else {
  #  train_df1 <- train_df
  #}

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

    if (trace %in% c(1, 2)) cat("Parallel computing for cross-validation started on", length(cl), "cores.")

    # Load the package on all workers
    parallel::clusterCall(cl, function() library(pye))

    # --- Fit primary models ---
    parallel::clusterExport(cl, c("train_df", "X", "y", "alpha", "model_type", "regressors_betas", "k", "trace", "lambda", "c_to_use"), envir = environment())
    fitted_models <- parallel::parLapply(cl, lambda, function(x) plr_estimation(df = train_df, X = X, y = y, alpha = alpha,
                                                                                   lambda = x,
                                                                                   c_to_use = c_to_use,
                                                                                   regressors_betas = regressors_betas,
                                                                                   model_type = model_type, fold = k,
                                                                                   trace = trace))
  } else {

    fitted_models <- lapply(lambda, function(x) plr_estimation(df = train_df, X = X, y = y, alpha = alpha,
                                                                  lambda = x,
                                                                  c_to_use = c_to_use,
                                                                  regressors_betas = regressors_betas,
                                                                  model_type = model_type, fold = k,
                                                                  trace = trace))
  }

  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat"))

  if (length(gamma_start_input) == 0) {
  # If gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
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
    for (i in seq_along(lambda)) {
      # Measures on the train set - these are multiple tables based on the number of considered TAUs
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
    # --- Parallel / Sequential Execution ---
    if (!is.null(cl)) { # Parallel execution (reusing the cluster)
      if (trace > 0) cat("Start parallel computing for cross-validation of covYI, I am using:", length(cl), "cores.\n")
      if (used_cores > max.cores) {
        warning("The number of specified cores (", used_cores, ") is larger than the number of physical cores available (", max.cores, ")!")
      }

      parallel::clusterExport(cl, c("train_df", "tau", "y", "C", "trace", "penalty_g", "w_g",
                                    "alpha_g", "a1_g", "a2_g", "trend_g", "kernel_g", "k",
                                    "gamma_start_input1", "gamma_start_default",
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

      fitted_gammas <- lapply(seq_along(z_hat), function (x) lapply(tau, function(t) covYI_KS_estimation(
                                                                    df = cbind(train_df, z_hat = z_hat[[x]]$z_hat),
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
  all_measures_test <- mapply(function(z) plr_predict(df = test_df,
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

    cov_results <- lapply(seq_along(lambda), function(xx) lapply(listing[[xx]],
                                                                function(x) covYI_KS(df = cbind(test_df1, z_hat = z_hat[[xx]]),
                                                                                      z = "z_hat", y = y, C = C1,
                                                                                      gammas = x[[1]],
                                                                                      tau = x[[2]],
																																											w = w_g,
                                                                                      kernel = kernel_g, alpha = alpha_g,
                                                                                      a1 = a1_g, a2 = a2_g, penalty = penalty_g,
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
        cat("youden_index:", all_measures_test["youden_index", ][[i]], "; sensitivity:", all_measures_test["sensitivity", ][[i]], "; specificity:", all_measures_test["specificity", ][[i]], "; geometric_mean:", all_measures_test["geometric_mean", ][[i]], "; fdr:", all_measures_test["fdr", ][[i]], "; mcc:", all_measures_test["mcc", ][[i]], "; auc:", all_measures_test["auc", ][[i]], "; corrclass:", all_measures_test["corrclass", ][[i]], "; \n")
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
  return(methods::new("plr_cross_validation_output", auc = auc,
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










#' @title Cross-Validation for Optimal glmnet Regularization Parameter Selection
#'
#' @description Performs k-fold cross-validation to select optimal penalization
#'   parameters for penalized logistic regression models. This function supports
#'   models fitted via the `glmnet` and `ncvreg` packages (Lasso, Elastic-Net,
#'   SCAD, MCP).
#'
#'   The primary goal is to identify the optimal `lambda` for the main
#'   regressors (`X`). Optionally, if `c_function_of_covariates` is `TRUE`,
#'   the function can also perform a search for the optimal `tau` parameter for
#'   a set of covariates (`C`) used to model a subject-specific cut-off point
#'   via the `covYI` method.
#'
#'   The hyperparameter tuning can be done in two ways:
#'   1.  Simultaneous: A grid search over all combinations of `lambda` and
#'       `tau` is performed.
#'   2.  Sequential: A two-step process where the best `lambda` is first
#'       selected (based on `measure_to_select_lambda`), and then, using that
#'       optimal `lambda`, the best `tau` is determined.
#'
#' @param n_folds Integer. The number of folds for cross-validation. Must be
#'   at least 2.
#' @param df A `data.frame` containing the complete dataset, including
#'   predictors, the target variable, and any covariates.
#' @param X A character vector of column names from `df` to be used as
#'   regressors (predictors). Defaults to all columns in `df` that are not `y`
#'   or `C`.
#' @param y A character string specifying the column name of the binary target
#'   variable (must be coded as 0 and 1). Default is `"y"`.
#' @param C A character vector of column names from `df` to be used as
#'   covariates for `covYI`. Default is `NULL`.
#' @param lambda A numeric vector of one or more non-negative penalization
#'   parameters for the regressors `X`.
#' @param tau A numeric vector of one or more non-negative penalization
#'   parameters for the covariates `C`. Used only if `c_function_of_covariates`
#'   is `TRUE`. Default is 0.
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param c_to_use character or numeric value specifying the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If set
#'   to `"Youden"`, the function computes the optimal cut-point using Youden's Index.
#'   If a numeric value is provided, it is used directly as the threshold.
#'   Default is "Youden".
#' @param alpha A numeric value between 0 and 1 for the Elastic-Net mixing
#'   parameter in `glmnet`. `alpha = 1` is the lasso (default), and `alpha = 0`
#'   is ridge regression.
#' @param model_type A character string specifying the penalized logistic
#'   regression model to use. Must be one of: `"logLasso"`, `"logElasticNet"`,
#'   `"logSCAD"`, or `"logMCP"`.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param trace Integer (0, 1, or 2). Controls the verbosity of the output.
#'   `0` = no output, `1` = summary results (default), `2` = detailed step-by-step
#'   output.
#' @param seed Integer. A random seed for ensuring reproducibility of the
#'   fold assignments. Default is 1.
#' @param used_cores Integer. The number of CPU cores to use for parallel
#'   processing of folds. If `1` (default), execution is sequential.
#' @param scaling Logical. If `TRUE`, the predictor and covariate columns of
#'   `df` are standardized before model fitting. Default is `FALSE`.
#' @param c_function_of_covariates Logical. If `TRUE`, the cut-off point is
#'   estimated as a function of covariates `C` using the `covYI` method.
#'   Default is `FALSE`.
#' @param simultaneous Logical. If `c_function_of_covariates` is `TRUE`, this
#'   determines the tuning strategy. If `TRUE`, a grid search over `lambda` and
#'   `tau` is performed. If `FALSE` (default), a sequential two-step search is
#'   performed.
#' @param measure_to_select_lambda Character. The performance metric used to
#'   select the optimal `lambda` in the first step of a sequential search
#'   (when `simultaneous = FALSE`). Must be one of `"auc"`, `"aauc"`, `"aYI"`,
#'   `"yi"`, `"sen"`, `"spc"`, `"gm"`, `"fdr"`, `"mcc"`, or `"ccr"`. Default is
#'   `"ccr"` (correct classification rate).
#' @param alpha_g Numeric. The Elastic-Net mixing parameter for the `covYI`
#'   model. Default is 0.5.
#' @param penalty_g Character. The penalty type for `covYI`. Must be one of
#'   `"L12"`, `"L1"`, `"EN"`, `"SCAD"`, `"MCP"`. Default is `"L1"`.
#' @param kernel_g Character. The kernel for density estimation in `covYI`.
#'   Default is `"gaussian"`.
#' @param a1_g,a2_g Numeric. Regularization parameters for SCAD and MCP
#'   penalties in `covYI`. Defaults are 3.7 and 3.0, respectively.
#' @param trend_g Character. Optimization algorithm for `covYI`: `"monotone"`
#'   (mmAPG) or `"nonmonotone"` (mnmAPG). Default is `"monotone"`.
#' @param gamma_start_input An optional numeric vector for the starting
#'   coefficients (`gammas`) in `covYI`.
#' @param gamma_start_default Character. Default starting point for `gammas` if
#'   `gamma_start_input` is `NULL`. `"zeros"` (default) or `"corr"`.
#' @param regressors_gammas An optional numeric vector of true gamma
#'   coefficients for simulation/evaluation. Default is `NULL`.
#' @param max_iter_g Integer. Maximum iterations for `covYI` optimization.
#'   Default is 10000.
#' @param delta_g,max_alpha_g,stepsizeShrink_g,min_alpha_g,convergence_error_g
#'   Numeric. Advanced optimization parameters for `covYI`. See
#'   `covYI_KS_estimation` for details.
#' @param run_aauc Logical. If `FALSE` (default), the adjusted-AUC (aAUC) and
#'   adjusted-YI (aYI) are not computed to save time.
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If `NULL`, output goes to the console.
#'   Default is `"log_logistic_models.txt"`.
#'
#' @return A list containing the aggregated results of the cross-validation.
#'   \item{model_type}{The primary model type used (e.g., `"logLasso"`).}
#'   \item{penalty_g}{The penalty type used for `covYI` (if applicable).}
#'   \item{cv_time}{Total time for the cross-validation process.}
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
#'   \item{lambda_hat_`[measure]`}{The optimal `lambda` value selected based on
#'     maximizing each performance measure (e.g., `lambda_hat_auc`,
#'     `lambda_hat_ccr`).}
#'   \item{tau_hat_`[measure]`}{The optimal `tau` value for each measure, if
#'     `c_function_of_covariates` was `TRUE`.}
#'   \item{c_function_of_covariates, simultaneous, measure_to_select_lambda}{
#'     Input parameters defining the CV strategy.}
#'   \item{lambda_star}{The optimal `lambda` chosen in the first step of a
#'     sequential search. `NA` otherwise.}
#'   \item{n_betas}{A matrix (`n_folds` x `length(lambda)`) showing the number
#'     of non-zero beta coefficients.}
#'   \item{n_betas_star}{Number of non-zero betas for `lambda_star` in a
#'     sequential search.}
#'   \item{n_gammas}{A list of matrices showing the number of non-zero gamma
#'     coefficients for each `lambda` and `tau` combination.}
#'   \item{betas}{A list of the estimated beta coefficients for each fold and
#'     `lambda`.}
#'   \item{betas_star}{A list of estimated beta coefficients for `lambda_star`
#'     in a sequential search.}
#'   \item{gammas}{A list of the estimated gamma coefficients for each fold,
#'     `lambda`, and `tau` combination.}
#'
#' @examples
#' library(pye)
#'
#' # 1. Simulate a small dataset
#' set.seed(123)
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#'
#' # 2. Define a small grid for hyperparameters
#' # In a real scenario, these would be calibrated (e.g., using calibrate_lambda_max)
#' # and the sequence would be longer.
#' lambda_seq <- c(0.1, 0.05)
#' tau_seq <- c(0.2, 0.1)
#'
#' # 3. Run cross-validation
#' # This example uses covariates in a sequential search.
#' # For a quick run, we use a small number of folds and iterations.
#' cv_results <- plr_compute_cv(
#'   n_folds = 2,
#'   df = df,
#'   X = X,
#'   y = y,
#'   C = C,
#'   model_type = "logLasso",
#'   lambda = lambda_seq,
#'   tau = tau_seq,
#'   trace = 0, # Use trace = 1 for more details
#'   used_cores = 1,
#'   c_function_of_covariates = TRUE,
#'   simultaneous = FALSE,
#'   measure_to_select_lambda = "auc", #"auc","aauc","aYI","ccr","yi","gm","pye")
#'   penalty_g = "L1",
#'   max_iter_g = 9 # Low for example speed
#' )
#'
#' # 4. Inspect the results
#' cat("Cross-validation finished in:", round(cv_results$cv_time, 2), "minutes\n")
#' cat("Optimal lambda based on AUC:", cv_results$lambda_hat_auc, "\n")
#' cat("Optimal tau based on AUC:", cv_results$tau_hat_ccr, "\n")
#'
#' # Show mean test AUC for each lambda/tau combination
#' # The structure is nested: results$auc$`lambda_...`$test
#' mean_test_auc <- sapply(cv_results$auc, function(res) colMeans(res$test))
#' print("Mean Test AUC matrix (rows = lambda, cols = tau):")
#' print(mean_test_auc)
#'
#' @export
plr_compute_cv <- function (n_folds, df, X = NULL, y = "y",
                            C = NULL, lambda,
                            tau = 0,
														w_g = 0.5,
														c_to_use = "Youden",
														alpha = 0.5,
														regressors_betas = NULL, trace = 1,
                            model_type = c("logLasso", "logElasticNet", "logSCAD", "logMCP"),
                            seed = 1, used_cores = 1, scaling = FALSE,
                            c_function_of_covariates = FALSE,
                            simultaneous = FALSE,
                            measure_to_select_lambda = "ccr",
                            alpha_g = 0.5, penalty_g = "L1",
														kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                            trend_g = "monotone", gamma_start_input = NULL,
														gamma_start_default = "zeros",
                            regressors_gammas = NULL, max_iter_g = 10000,
                            delta_g = 1e-5, max_alpha_g = 10000,
                            stepsizeShrink_g = 0.8, min_alpha_g = 1e-12,
														convergence_error_g = 1e-7,
                            run_aauc = FALSE, log_file = "log_logistic_models.txt") {

  # Start calculation of estimation time
  start_time <- Sys.time()

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input 'df' must be a non-empty data frame.")
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

  # Validate main model parameters
  model_type <- match.arg(model_type)
  if (!is.numeric(lambda) || length(lambda) < 1 || any(lambda < 0)) stop("Parameter 'lambda' must be a numeric vector of non-negative values.")
  # --- Model Type Validation ---
  valid_models <- c("logLasso", "logElasticNet", "logSCAD", "logMCP")
  if (!(model_type %in% valid_models)) {
    stop("The parameter 'model_type' must be one of: ", paste(valid_models, collapse = ", "))
  }
  # --- Alpha Validation and Standardization based on Model Type ---
  # 1. Check basic alpha properties (must be a single number/NA)
  if (!is.numeric(alpha) || length(alpha) != 1) {
    if (!(model_type %in% c("logSCAD", "logMCP") && is.na(alpha))) {
      stop("The parameter 'alpha' must be a single numeric value (or NA for SCAD/MCP).")
    }
  }
  if (model_type == "logLasso") {
    # Lasso is strictly L1, so alpha must be 1.
    if (alpha != 1) {
      if (trace %in% c(2)) {
        warning("For 'logLasso', 'alpha' is forced to 1, overriding the input value.")
      }
      alpha <- 1
    }
  } else if (model_type == "logElasticNet") {
    # Elastic-Net requires 0 < alpha < 1. Check range.
    if (alpha <= 0 || alpha >= 1) stop("For 'logElasticNet', 'alpha' must be a single numeric value between 0 (exclusive) and 1 (exclusive).")
  } else if (model_type %in% c("logSCAD", "logMCP")) {
    # SCAD/MCP do not use alpha mixing parameter. Set to NA for clarity in output.
    if (trace %in% c(2) && !is.na(alpha) && alpha!=0) {
      warning(paste0("For '", model_type, "', 'alpha' is not used for penalty mixing and is set to 0."))
			alpha <- 0
    }
  }

  # Validate execution and CV parameters
  if (!(trace %in% c(0, 1, 2))) stop("'trace' must be 0, 1, or 2.")
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {
    stop("'used_cores' must be a single positive integer.")
  }
  if (!is.logical(scaling)) stop("'scaling' must be a logical (TRUE/FALSE).")
  if (!is.logical(run_aauc)) stop("'run_aauc' must be a logical (TRUE/FALSE).")
  if (!is.numeric(seed) || length(seed) != 1) stop("'seed' must be a single numeric value.")

  # Validate parameters for covariate-adjusted cut-point (covYI)
  if (!is.logical(c_function_of_covariates)) stop("'c_function_of_covariates' must be a logical (TRUE/FALSE).")
  if (!is.logical(simultaneous)) stop("'simultaneous' must be a logical (TRUE/FALSE).")

  tau_initial <- 0 # Initialize the storage of the original tau

  # Check if tau exists when c_function_of_covariates = TRUE
  if (c_function_of_covariates) {
    if (is.null(tau) || length(tau) == 0) {stop("Parameter 'tau' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && length(tau) == 1 && tau == 0) {stop("Parameter 'tau' cannot be a single value of 0 if 'c_function_of_covariates' is TRUE.")}
    if (is.numeric(tau) && sum(tau == 0) == length(tau)) {stop("Parameter 'tau' cannot be a vector of all zeros if 'c_function_of_covariates' is TRUE.")}
    if (is.null(C) || length(C) == 0) { stop("Parameter 'C' cannot be NULL or empty if 'c_function_of_covariates' is TRUE.")}
    if (simultaneous == FALSE) { # In sequential search, tau is set to 0 for the first step
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
  # df1 <- cbind(ID, df[, (names(df) %in% c(y,X,C))]) #OLD
  df1 <- cbind(ID, df[, c(y, X, C)], drop = FALSE)

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty_g %in% valid_penalties)) {stop("Parameter 'penalty_g' must be one of: ", paste(valid_penalties, collapse = ", "))}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(max_iter_g) || length(max_iter_g) != 1 || max_iter_g < 2 || floor(max_iter_g) != max_iter_g) {stop("Parameter 'max_iter_g' must be an integer greater than or equal to 2.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(model_type %in% c("logLasso", "logElasticNet", "logSCAD", "logMCP"))) { stop("The parameter model_type is wrongly specified")}
  valid_measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
  if (!(measure_to_select_lambda %in% valid_measures)) {stop("Parameter 'measure_to_select_lambda' must be one of: ", paste(valid_measures, collapse = ", "))}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}
  if (!(is.character(c_to_use) && c_to_use == "Youden") && !(is.numeric(c_to_use) && length(c_to_use) == 1 && !is.na(c_to_use))) stop("Parameter 'c_to_use' must be either the string 'Youden' or a single numeric value.")

  #standardize df1
  if (scaling == TRUE) {
    df1 <- scaling_df_for_pye (df = df1, X = colnames(df1[, names(df1) %in% c(X, C)]), y = "y")$df_scaled
  }

  set.seed(seed)

  # Check if df1 is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1[[y]][df1[[y]] == 1]) < 2 * n_folds)) {stop("df1 contains too few 1s for this number of folds")
  } else if ((length(df1[[y]][df1[[y]] == 0]) < 2 * n_folds)) {stop("df1 contains too few 0s for this number of folds")}

  # Divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1, y)), 1:2]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 0, ])), replace = FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 1, ])), replace = FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0, fold_i_1))
  folds_i <- merge(df1[, 1:2], df_sort[, 1:3], by = 'ID', all = FALSE, sort = FALSE)[, 4]

  # Names of the columns
  #lambdanames <- unlist(lapply(lambda, function (x) paste("lambda", x, sep = "=")))
  lambdanames <- paste("lambda", lambda, sep = "=")
  taunames <- paste("tau", tau, sep = "=")
  foldnames <- paste("fold", 1:n_folds, sep = "=")
  # Accuracy measures (train and test)
  list_of_measures <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                        "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

  auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures) {
    assign(mes, lapply(lambda, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)),
                                                  test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))))
    eval(substitute(names(x) <- lambdanames, list(x = as.symbol(mes))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(foldnames, lambdanames))
  n_gammas <- lapply(lambda, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
  names(n_gammas) <- lambdanames

  if (trace > 0) {
    cat("Starting CV with model_type:", model_type, "\n")
  }

  # fill the matrices
  results <- mapply(function(k) plr.cv(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, C = C,
                                          lambda = lambda, tau = tau, w_g = w_g,
                                          c_to_use = c_to_use, alpha = alpha,
                                          alpha_g = alpha_g, penalty_g = penalty_g,
                                          folds_i = folds_i,
                                          k = k, regressors_betas = regressors_betas,
                                          trace = trace,
                                          model_type, auc = auc,
                                          aauc = aauc, aYI = aYI,
                                          youden_index = youden_index,
                                          sensitivity = sensitivity,
                                          specificity = specificity,
                                          geometric_mean = geometric_mean, fdr = fdr,
                                          mcc = mcc,
                                          corrclass = corrclass,
                                          n_betas = n_betas, n_gammas = n_gammas,
                                          used_cores = used_cores,
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
    assign(mes[], lapply(seq_along(lambda), function(i) list(train = wrapper(results, mes, i, "train", tau),
                                                             test = wrapper(results, mes, i, "test", tau))))
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

    # now we execute plr.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
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
    results <- mapply(function(k) plr.cv(df = df1[, names(df1) != "ID", drop = FALSE],
		                                     X = X, y = y, C = C,
																				 lambda = lambda_star, tau = tau, w_g = w_g,
                                         c_to_use = c_to_use, alpha = alpha,
                                         alpha_g = alpha_g, penalty_g = penalty_g,
                                         folds_i = folds_i,
                                         k = k, regressors_betas = regressors_betas,
                                         trace = trace,
                                         model_type, auc = auc,
                                         aauc = aauc, aYI = aYI,
                                         youden_index = youden_index,
                                         sensitivity = sensitivity,
                                         specificity = specificity,
                                         geometric_mean = geometric_mean, fdr = fdr,
                                         mcc = mcc,
                                         corrclass = corrclass,
                                         n_betas = n_betas_star, n_gammas = n_gammas,
                                         used_cores = used_cores,
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
    cat("----------------------> END OF THE CROSS-VALIDATION OF THE PEN. LOG. METHOD <----------------- \n")
    cat("------------------> For the whoole Cross-validation it took:", cv_time, "minutes <------------- \n")
  }
  return(list(model_type = model_type,
              c_to_use = c_to_use,
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
