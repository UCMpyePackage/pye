#https://cran.r-project.org/web/packages/penalizedSVM/penalizedSVM.pdf
#Becker, N., Werft, W., Toedt, G., Lichter, P. and Benner, A.(2009)
#PenalizedSVM: a R-package for feature selection SVM classification,
#Bioinformatics, 25(13), p 1711-1712
#https://cran.r-project.org/web/packages/sparseSVM/sparseSVM.pdf

# Support Vector Machine (SVM) classification with simultaneous feature
# selection using penalty functions is implemented.
# The smoothly clipped absolute deviation (SCAD), 'L1-norm',
# 'Elastic Net' ('L1-norm' and 'L2-norm') and 'Elastic SCAD'
# (SCAD and 'L2-norm') penalties are available. The tuning parameters
# can be found using either a fixed grid or a interval search.

#' @title Penalized svm Estimation for Coefficient and Feature Selection
#'
#' @description
#' Fit a support vector machine (SVM) classifier with optional feature
#' selection / regularization using multiple penalty families and kernels.
#' This function implements a consistent interface to: (i) SCAD penalized SVM
#' estimators, (ii) SCAD+L2 (elastic-SCAD) estimators, and (iii) sparse/elastic
#' sparse SVM via sparseSVM. The routine returns estimated coefficients,
#' sample decision scores, an estimated optimal cutpoint (Youden criterion)
#' and a comprehensive set of classification performance measures.
#'
#' @details
#' The function accepts a user supplied design matrix specification `X` and a
#' binary response `y` coded as 0/1 (or -1/1). Internally `y` is converted to
#' the form required by the selected estimation routine. For standard linear
#' SVMs the function will attempt to reconstruct linear predictor weights
#' (when available) from the fitted model object. For penalized estimators the
#' function extracts penalty-specific coefficients when present. The optimal
#' cutpoint for mapping decision scores to class labels is obtained by applying
#' the Youden index via OptimalCutpoints; when that routine fails the function
#' falls back to a conservative default. Output is returned as a named list
#' designed for subsequent evaluation or cross-validation.
#'
#' @param df data.frame containing the observations (rows) and variables (columns).
#' @param X character vector with names of predictors in `df`. If omitted the
#'   function will attempt to use all columns in `df` except `y`.
#' @param y single character string with the name of the binary outcome column
#'   in `df` (allowed codings: 0/1 or -1/1).
#' @param model_type character scalar selecting the estimation method. Allowed
#'   values include: "SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM".
#' @param c_to_use character or numeric value specifying the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If set
#'   to `"Youden"`, the function computes the optimal cut-point using Youden's Index.
#'   If a numeric value is provided, it is used directly as the threshold. 
#'   If `NULL`, a default is assigned based on `model_type`: `"SCADSVM"` uses 
#'   `0`, `"ElasticSCADSVM"`, `"l1SVM"` or `"enSVM"` use `"Youden"`.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param lambda numeric. Primary regularization parameter used by penalized
#'   estimators (SCAD, ElasticSCAD, sparseSVM). If NULL, methods that require
#'   a penalty will return an error.
#' @param lambda2 numeric. Secondary L2 penalty used by ElasticSCADSVM.
#' @param alpha numeric in (0, 1]. Elastic-net mixing parameter for sparseSVM
#'   (alpha = 1 corresponds to L1). Default is 0.5.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param trace integer verbosity level: 0 = silent, 1 = brief summary, 2 = verbose.
#' @param seed integer used to set RNG where supported by the back-end routines.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#'
#' @return A named list with the principal elements:
#'   \item{model}{The fitted model object returned by the selected backend
#'     (e.g. penalizedSVM object, sparseSVM object).}
#'   \item{model_type}{Character string with the selected model type.}
#'   \item{betas_hat}{Numeric named vector with estimated predictor coefficients
#'     where recoverable (zeros if unavailable). Names follow `X`.}
#'   \item{c_hat}{Numeric cutpoint computed by Youden's criterion on the
#'     training decision scores.}
#'   \item{z_hat}{Data.frame with columns `ID` and `z_hat` (decision score).}
#'   \item{y_hat}{Data.frame with columns `ID` and `y_hat` (predicted labels at `c_hat`).}
#'   \item{youden_index, auc, sensitivity, specificity, geometric_mean, fdr, mcc, 
#'     corrclass}{Numeric scalar performance summaries computed on the training data.}
#'   \item{TP, TN, FP, FN}{Confusion matrix cell counts at `c_hat`.}
#'   \item{lambda, lambda2, alpha}{Model and tuning parameters used.}
#'   \item{n_betas, n_total_var, n_predicted_zeros, n_predicted_non_zeros}{Variable selection 
#'     diagnostics.}
#'   \item{estimation_time}{Elapsed wall clock time (minutes).}
#'
#' @references
#' See references for penalizedSVM, sparseSVM and OptimalCutpoints
#' packages for algorithmic details.
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' regressors_betas<-sim_data$nregressors
#' model_type <- "SCADSVM"
#' lambda <- 0.1
#'
#' res <- psvm_estimation(df = df, X = X, y = y, model_type = model_type, lambda = lambda,
#'         regressors_betas = regressors_betas, trace = 1)
#'
#' print(res)
#'
#' @importFrom sparseSVM sparseSVM
#' @importFrom penalizedSVM scadsvc scad_L2.svc
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom stats setNames
#'
#' @export
psvm_estimation <- function (df,
                                X = NULL,
                                y = "y",
                                lambda = NULL,
                                lambda2 = 0.05,
                                alpha = 0.5,
                                model_type,
                                c_to_use = NULL,
                                regressors_betas = NULL,
                                fold = NULL,
                                trace = 1,
                                seed = 1,
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

  if (!is.null(lambda) && (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)) stop("'lambda' must be NULL or a non-negative numeric scalar.")
  if (!is.numeric(lambda2) || length(lambda2) != 1L || lambda2 < 0) stop("'lambda2' must be a non-negative numeric scalar.")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha > 1) stop("'alpha' must lie in (0, 1].")
  if (!is.character(model_type) || length(model_type) != 1L) stop("'model_type' must be a single character string.")
  if (!is.numeric(trace) || !(trace %in% c(0, 1, 2))) stop("'trace' must be 0, 1,2.")
  if (!is.numeric(seed) || length(seed) != 1L) stop("'seed' must be a single integer-like numeric.")
  # Enforce model_type choice
  model_type <- match.arg(model_type, choices = c("SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM"))

  if (is.null(c_to_use)) {c_to_use <- switch(model_type, "SCADSVM" = 0, "ElasticSCADSVM" = "Youden", "l1SVM" = "Youden", "enSVM" = "Youden")} # In PYE and covYI paper they were: "l1SVM" = "MCT", "enSVM" = "MCT")}
  if (!(is.character(c_to_use) && c_to_use == "Youden") && !(is.numeric(c_to_use) && length(c_to_use) == 1 && !is.na(c_to_use))) stop("Parameter 'c_to_use' must be either the string 'Youden' or a single numeric value.")

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  # ---- Prepare data ----
  ID <- rownames(df)
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1, -1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0, 1 or -1 (excluding NA).")}
  if (length(unique(df1[[y]])) != 2) {stop("The target variable 'y' must contain two unique values (0 and 1 or -1 and 1).")}

  df_x_as_matrix <- as.matrix(df1[, X, drop = FALSE])
  df_y_as_matrix <- as.matrix(factor(df1[[y]]))
  # Convert y in -1/1 if it is not
  if (!all((levels(factor(df_y_as_matrix)) == c("-1", "1")))) {
    df_y_as_matrix <- replace(df_y_as_matrix, df_y_as_matrix == "0", "-1")
  }

  # initialize outputs
  betas_hat <- setNames(rep(0, length(X)), X)
  c_hat <- 0
  z_hat <- data.frame(ID = df1$ID, z_hat = rep(0, nrow(df1)), stringsAsFactors = FALSE)
  y_hat <- data.frame(ID = df1$ID, y_hat = rep(0, nrow(df1)), stringsAsFactors = FALSE)
  model <- NULL

  # model_type can be: "SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM"
 if (model_type == "SCADSVM") {
    #scadsvc: for SCAD SVM
    tmp <- try( SCADSVM <- penalizedSVM::scadsvc(lambda1 = lambda,
                                                 x = df_x_as_matrix,
                                                 y = df_y_as_matrix,
                                                 a = 3.7,
                                                 maxIter = 100000,
                                                 verbose = FALSE,
                                                 seed = seed, tol = 10^(-8)),
                                                 silent = TRUE)

    count <- 1L
    lambda_new <- lambda
    df_x_as_matrix1 <- df_x_as_matrix

    while(inherits(tmp, "try-error")) {
      # Informational message to help debugging in long runs
      cat("SCADSVM attempt error number", count, " - attempting remedies.\n")
      if (!is.null(tmp)) print(tmp)

      # Strategy 1: slightly increase lambda to avoid numerical problems
      lambda_new <- lambda_new + 0.000001 #with this little modification it works and do not return error
      cat(" -> trying with lambda = ", format(lambda_new, scientific = TRUE), "\n")
      tmp <- try( SCADSVM <- penalizedSVM::scadsvc(lambda1 = lambda_new,
                                                   x = df_x_as_matrix,
                                                   y = df_y_as_matrix,
                                                   a = 3.7, maxIter = 100000,
                                                   verbose = FALSE,
                                                   seed = seed, tol = 10^(-5)),
                                                   silent = TRUE)
      count <- count + 1L
      if (count > 100L) {
        if (!inherits(tmp, "try-error"))  break

        # Strategy 2: progressively drop the first columns and retry on reduced matrix.
        # This tries to handle cases where some columns break the backend (defensive).
        column_start <- 1L
        lambda_new <- lambda # reset lambda perturbation when switching to column-dropping
        while (inherits(tmp, "try-error") && column_start < ncol(df_x_as_matrix)) {
          cat("SCADSVM attempt error number", count, " - attempting remedies.\n")
          if (!is.null(tmp)) print(tmp)
          column_start <- column_start + 1L
          cat(" -> dropping first", column_start - 1L, "columns and retrying (start column = ", column_start, ").\n")
          df_x_as_matrix1 <- df_x_as_matrix[, column_start:ncol(df_x_as_matrix), drop = FALSE]
          tmp <- try( SCADSVM <- penalizedSVM::scadsvc(lambda1 = lambda,
                                                        x = df_x_as_matrix1,
                                                        y = df_y_as_matrix,
                                                        a = 3.7, maxIter = 100000,
                                                        verbose = FALSE,
                                                        seed = seed, tol = 10^(-5)),
                                                        silent = TRUE)
          count <- count + 1L
          if (count > 200) {break}
        }
        # break outer loop whether succeeded or not to avoid infinite loops
        break
      }
    }

    SCADSVM <- penalizedSVM::scadsvc(lambda1 = lambda_new,
                                     x = df_x_as_matrix1,
                                     y = df_y_as_matrix,
                                     a = 3.7, maxIter = 100000, verbose = FALSE,
                                     seed = seed, tol = 10^(-5))

    #for predictions
    #pred = ifelse(df_x_as_matrix[, names(fit_SCADSVM$w)] %*% fit_SCADSVM$w + fit_SCADSVM$b > 0, 1, -1)

    if (SCADSVM[1] == "No variable selected.") {
      z_hat <- as.data.frame(rep(0, nrow(df_x_as_matrix)))
    } else {
      z_hat <- as.data.frame(SCADSVM$fitted)
      # Populate betas_hat using available weights (SCADSVM$w) when present.
      betas_hat[which(names(betas_hat) %in% attr(SCADSVM$w[!is.na(names(SCADSVM$w))], "names"))] <- SCADSVM$w[!is.na(names(SCADSVM$w))]
    }

    # Prepare the output
    model <- SCADSVM
    lambda <- lambda_new
    lambda2 <- NA
    alpha <- NA

  } else if (model_type == "ElasticSCADSVM") {
    #scad_L2.svc: for Elastic SCAD SVM
    tmp <- try(ElasticSCADSVM <- penalizedSVM::scad_L2.svc(lambda1 = lambda,
                                                           lambda2 = lambda2,
                                                           x = df_x_as_matrix,
                                                           y = df_y_as_matrix,
                                                           a = 3.7, maxIter = 100000,
                                                           verbose = FALSE, tol = 10^(-5)))

    count <- 1L
    lambda_new <- lambda
    df_x_as_matrix1 <- df_x_as_matrix

    while(inherits(tmp, "try-error")) {
      # Informational message to help debugging in long runs
      cat("ElasticSCADSVM attempt error number", count, " - attempting remedies.\n")
      if (!is.null(tmp)) print(tmp)

      # Strategy 1: slightly increase lambda to avoid numerical problems
      lambda_new <- lambda_new + 0.000001 #with this little modification it works and do not return error
      tmp <- try( ElasticSCADSVM <- penalizedSVM::scad_L2.svc(lambda = lambda_new,
                                                              lambda2 = lambda2,
                                                              x = df_x_as_matrix,
                                                              y = df_y_as_matrix,
                                                              a = 3.7, maxIter = 100000,
                                                              verbose = FALSE, tol = 10^(-5)), silent = TRUE)
      count <- count + 1L
      if (count > 100L) {
        if (!inherits(tmp, "try-error"))  break

        # Strategy 2: progressively drop the first columns and retry on reduced matrix.
        # This tries to handle cases where some columns break the backend (defensive).
        column_start <- 1L
        lambda_new <- lambda # reset lambda perturbation when switching to column-dropping
        while (inherits(tmp, "try-error") && column_start < ncol(df_x_as_matrix)) {
          cat("SCADSVM attempt error number", count, " - attempting remedies.\n")
          if (!is.null(tmp)) print(tmp)
          column_start <- column_start + 1L
          cat(" -> dropping first", column_start - 1L, "columns and retrying (start column = ", column_start, ").\n")
          df_x_as_matrix1 <- df_x_as_matrix[, column_start:ncol(df_x_as_matrix), drop = FALSE]
          tmp <- try( ElasticSCADSVM <- penalizedSVM::scad_L2.svc(lambda1 = lambda_new,
                                                                  lambda2 = lambda2,
                                                                  x = df_x_as_matrix1,
                                                                  y = df_y_as_matrix,
                                                                  a = 3.7, maxIter = 100000,
                                                                  verbose = FALSE, tol = 10^(-5)), silent = TRUE)
          count <- count + 1L
          if (count > 200) {break}
        }
        # break outer loop whether succeeded or not to avoid infinite loops
        break
      }
    }

    ElasticSCADSVM <- penalizedSVM::scad_L2.svc(lambda1 = lambda_new,
                                                 lambda2 = lambda2,
                                                 x = df_x_as_matrix1,
                                                 y = df_y_as_matrix,
                                                 a = 3.7, maxIter = 100000,
                                                 verbose = FALSE, tol = 10^(-5))

    if (ElasticSCADSVM[1] == "No variable selected.") {
      z_hat <- as.data.frame(rep(0, nrow(df_x_as_matrix)))
    } else {
      z_hat <- as.data.frame(ElasticSCADSVM$fitted)
      # Populate betas_hat using available weights (ElasticSCADSVM$w) when present.
      betas_hat[which(names(betas_hat) %in% attr(ElasticSCADSVM$w[!is.na(names(ElasticSCADSVM$w))], "names"))] <- ElasticSCADSVM$w[!is.na(names(ElasticSCADSVM$w))]
    }

    # Prepare the output
    model <- ElasticSCADSVM
    lambda <- lambda_new
    lambda2 <- lambda2
    alpha <- NA

  } else if (model_type == "l1SVM" || model_type == "enSVM") {

    # sparseSVM backends (l1 or elastic net)
    alpha_use <- if (model_type == "l1SVM") 1 else alpha

    l1SVM_or_enSVM <- sparseSVM::sparseSVM(X = df_x_as_matrix,
                                           y = df_y_as_matrix,
                                           alpha = alpha_use,
                                           gamma = 0.1,
                                           lambda = if (lambda == 1) {c(2, lambda)} else {c(1, lambda)}, #c(1, lambda)
                                           preprocess = "none", # c("standardize", "rescale", "none"),
                                           screen = "SR", # c("ASR", "SR", "none"),
                                           max.iter = 100000, eps = 1e-7,
                                           message = F)

    #ATTENTION: I invert the signs since in the original package (sparseSVM) if the model predicts -1 and not 1!
    if (l1SVM_or_enSVM$levels[1] == "-1") {
      # if the predicted value is -1
      z_hat <- as.data.frame(-stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[1] + (df_x_as_matrix %*% -stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[-1])) # coef (l1Sl1SVM_or_enSVMVM, lambda, exact = TRUE)[1])
    } else {
      # if the predicted value is 1 we apply the normal calculation
      z_hat <- as.data.frame(stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[1] + (df_x_as_matrix %*% stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[-1])) # coef (l1SVM_or_enSVM, lambda, exact = TRUE)[1])
    }
    # NB: after many controls, I understood that the estimation is biased in centering the estimation in 0.
    # This produces bad performances in CC and other measures that depend on the cut-off even if the AUC is good!

    # Extract coefficients
    if (l1SVM_or_enSVM$levels[1] == "-1") {
      betas_hat <- -stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[-1] #I dont save the intercept
    } else {
      betas_hat <- stats::coef (l1SVM_or_enSVM, lambda, exact = TRUE)[-1] #I dont save the intercept
    }

    # Prepare the output
    model <- l1SVM_or_enSVM
    lambda <- lambda
    lambda2 <- NA
    alpha <- alpha_use

  } else {stop("The parameter model_type is not valid!")}

  colnames(z_hat) <- "z_hat"
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
  c_hat <- c
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  if ((trace  %in% c(1, 2))) {
    cat("Estimation done using the Penalized SVM method type:", model_type, "; \n")
		cat("-> ")
    if (!is.null(fold)) {cat("fold =", fold, "; ")}
    if (length(lambda) != 0) {cat("lambda:", lambda, "; ")}
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

  return(list(model = model,
              model_type = model_type,
              betas_hat = betas_hat,
              X_model = X,
              youden_index = youden_index,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr, mcc = mcc,
              corrclass = corrclass,
              auc = auc,
              lambda = lambda, lambda2 = lambda2,
              alpha = alpha,
              c_hat = c_hat, z_hat = z_hat, y_hat = y_hat,
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









#' @title Prediction and Performance Evaluation for Penalized svm Models
#'
#' @description
#' Apply a previously fitted SVM model (from `psvm_estimation`) to new data
#' and generate predictions along with performance metrics.
#'
#' @details
#' This function takes a data frame and a fitted SVM model (produced by
#' `psvm_estimation`) to generate predictions on new data. It supports
#' various SVM model types, including linear, polynomial, radial, sigmoid,
#' SCADSVM, ElasticSCADSVM, l1SVM, and enSVM. The function calculates
#' decision scores (z_hat) and predicted binary outcomes (y_hat) based on a
#' cut-off value. It also computes and returns a comprehensive set of
#' classification performance measures.
#'
#' @param df data.frame. The input dataset containing observations (rows) and
#'   variables (columns) for prediction.
#' @param y character. The name of the target variable column in `df`. This
#'   column should contain binary outcomes (0 and 1). Default is "y".
#' @param model_to_use list. A list containing the fitted model and related
#'   information, as returned by the `psvm_estimation` function. This list
#'   must include elements named `betas_hat`, `c_hat`, `X_model`, and
#'   `model_type`.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param trace integer. Controls the level of verbosity during the function's
#'   execution.
#'   - 0: No output is printed.
#'   - 1: A summary of the results is printed.
#'   - 2: Detailed information about the calculations is printed.
#'   Default is 1.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#' @param c_function_of_covariates logical. If `TRUE`, indicates that the
#'   cut-off point (`c_hat`) should be estimated as a function of covariates
#'   (using a method like covYI). If `FALSE`, a single cut-off value is used
#'   for all predictions. This parameter primarily affects whether results are
#'   printed within this function or handled externally. Default is `FALSE`.
#'
#' @return A named list with the following elements:
#'   \item{model}{The fitted model object used for prediction (same as in
#'     `model_to_use`).}
#'   \item{model_type}{Character string indicating the type of SVM model used
#'     (same as in `model_to_use`).}
#'   \item{betas_hat}{Numeric vector of estimated beta coefficients from the
#'     fitted model (same as in `model_to_use`).}
#'   \item{youden_index}{Numeric. Youden's J statistic, calculated on the
#'     predicted values.}
#'   \item{sensitivity}{Numeric. Sensitivity (true positive rate) of the
#'     predictions.}
#'   \item{specificity}{Numeric. Specificity (true negative rate) of the
#'     predictions.}
#'   \item{geometric_mean}{Numeric. Geometric mean of sensitivity and
#'     specificity.}
#'   \item{fdr}{Numeric. False discovery rate.}
#'   \item{mcc}{Numeric. Matthews correlation coefficient.}
#'   \item{corrclass}{Numeric. Correct classification rate (accuracy).}
#'   \item{auc}{Numeric. Area under the ROC curve.}
#'   \item{lambda}{Numeric. Value of lambda used in the model (same as in
#'     `model_to_use`).}
#'   \item{lambda2}{Numeric. Value of lambda2 used in the model (same as in
#'     `model_to_use`).}
#'   \item{alpha}{Numeric. Value of alpha used in the model (same as in
#'     `model_to_use`).}
#'   \item{c_hat}{Numeric. The cut-off value used to convert decision scores
#'     into binary predictions (same as in `model_to_use`).}
#'   \item{z_hat}{Data.frame with columns "ID" and "z_hat" (decision scores).}
#'   \item{y_hat}{Data.frame with columns "ID" and "y_hat" (predicted binary
#'     outcomes).}
#'   \item{n_total_var}{Integer. Total number of variables (regressors) in the
#'     model.}
#'   \item{n_predicted_zeros}{Integer. Number of regressors with estimated
#'     coefficients of zero.}
#'   \item{n_predicted_non_zeros}{Integer. Number of regressors with estimated
#'     coefficients different from zero.}
#'   \item{n_caught_betas}{Integer. If `regressors_betas` is provided, the
#'     number of non-zero "true" betas that were also identified as non-zero in
#'     the model.}
#'   \item{n_non_caught_betas}{Integer. If `regressors_betas` is provided, the
#'     number of non-zero "true" betas that were incorrectly estimated as zero
#'     in the model.}
#'   \item{n_caught_zero}{Integer. If `regressors_betas` is provided, the
#'     number of zero "true" betas that were also estimated as zero in the
#'     model.}
#'   \item{n_zero_not_caught}{Integer. If `regressors_betas` is provided, the
#'     number of zero "true" betas that were incorrectly estimated as non-zero
#'     in the model.}
#'   \item{TP}{Integer. Number of true positives.}
#'   \item{TN}{Integer. Number of true negatives.}
#'   \item{FP}{Integer. Number of false positives.}
#'   \item{FN}{Integer. Number of false negatives.}
#'   \item{estimation_time}{difftime. The time elapsed during the prediction
#'     process.}
#'
#' @examples
#' library(pye)
#' # Generate synthetic data
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' train_df <- sim_data$train_df_scaled
#' test_df <- sim_data$test_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' regressors_betas <- sim_data$nregressors
#'
#' # Estimate an SVM model
#' model <- psvm_estimation(
#'   df = train_df, X = X, y = y, lambda = 0.1, model_type = "SCADSVM",
#'   regressors_betas = regressors_betas, trace = 0
#' )
#'
#' # Apply the model to the test data
#' predictions <- psvm_predict(df = test_df, y = y, model_to_use = model,
#'   trace = 1, regressors_betas = regressors_betas, c_function_of_covariates = FALSE
#' )
#'
#' # Print the predictions
#' print(predictions)
#'
#' @export
psvm_predict <- function (df,
                          y = "y",
                          model_to_use,
                          fold = NULL,
                          regressors_betas = NULL,
                          trace = 1,
                          max.print = 10,
                          c_function_of_covariates = FALSE) {

  # Start calculation of estimation time
  start_time <- Sys.time()
	
	if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
  # Set max.print temporarily
  old_options <- options(max.print = max.print)
  on.exit(options(old_options))

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
  X <- model_to_use$X_model
  model_type <- model_to_use$model_type

  ID <- rownames(df)
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Ensure all variables required by the model (X) are present in the new data (df)
  if (length(X) == 0) {stop("The model contains no regressors (X). Cannot perform prediction.")}
  if (!all(X %in% names(df))) {
    missing_vars <- X[!X %in% names(df)]
    stop("The following required regressors from the trained model are missing in the new data: ", paste(missing_vars, collapse = ", "))
  }

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1, -1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0, 1 or -1 (excluding NA).")}
  if (length(unique(df1[[y]])) != 2) {stop("The target variable 'y' must contain two unique values (0 and 1 or -1 and 1).")}

  # Convert the data frame to a matrix for prediction
  df_x_as_matrix <- as.matrix(df1[, X, drop = FALSE])
  df_y_as_matrix <- as.matrix(factor(df1[[y]]))

  #convert y in -1,1 if they are not:
  if (!all((levels(factor(df_y_as_matrix)) == c("-1", "1")))) {
    df_y_as_matrix <- replace(df_y_as_matrix, df_y_as_matrix == "0", "-1")
  }

  #model_type can be: "SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM"
  if (model_type %in% c("SCADSVM", "ElasticSCADSVM")) {
    # SCADSVM and ElasticSCADSVM prediction
    if (model_to_use$model[1] == "No variable selected.") {
      # All predictions are zero if no variable is selected
      z_hat <- as.data.frame(rep(0, nrow(df_x_as_matrix)))
    } else {
      # Calculate decision scores
      z_hat <- as.data.frame(as.matrix(df_x_as_matrix[, names(model_to_use$model$w[!is.na(names(model_to_use$model$w))])]) %*% model_to_use$model$w[!is.na(names(model_to_use$model$w))] + model_to_use$model$b)
    }

  } else if (model_type %in% c("l1SVM", "enSVM")) {
    # l1SVM and enSVM prediction
      # If the predicted value is -1
    if (model_to_use$model$levels[1] == "-1") {
      # Invert the signs since in the original package (sparseSVM) the model predicts -1 and not 1!
      z_hat <- as.data.frame(-stats::coef (model_to_use$model, model_to_use$lambda, exact = TRUE)[1] + df_x_as_matrix %*% -stats::coef (model_to_use$model, model_to_use$lambda, exact = TRUE)[-1])
    } else {
      # If the predicted value is 1 we apply the normal calculation
      z_hat <- as.data.frame(stats::coef (model_to_use$model, model_to_use$lambda, exact = TRUE)[1] + df_x_as_matrix %*% stats::coef (model_to_use$model, model_to_use$lambda, exact = TRUE)[-1])
    }
  } else {stop("The parameter model_type is not valid!")}

  colnames(z_hat) <- "z_hat"
  df_hat <- cbind(df1, z_hat)

  # Predict binary outcomes based on the cut-off
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
  z_hat <- df_hat[, c("ID", "z_hat")]
  y_hat <- df_hat[, c("ID", "y_hat")]

  if (trace %in% c(1, 2)) {
    #print the results only if c_function_of_covariates = FALSE
    if (c_function_of_covariates == FALSE) {
      cat("Predictions executed using the Penalizez SVM model type:", model_type, "; \n")
      cat("-> ")
      if (!is.null(fold)) {cat("fold =", fold, "; ")}
      if (!is.na(model_to_use$lambda)) {cat("lambda:", model_to_use$lambda, "; ")}
      cat("youden_index:", youden_index, "; sensitivity:", sensitivity, "; specificity:", spec, "; geometric_mean:", gm, "; fdr:", fdr, "; mcc:", mcc, "; auc:", auc, "; corrclass:", corrclass, "; \n")
      cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, "; betas_hat: \n")
      visualize_betas <- betas_hat[which(betas_hat != 0)]
      print( visualize_betas)
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
              geometric_mean = gm, fdr = fdr,
              mcc = mcc,
              corrclass = corrclass,
              auc = auc,
              lambda = model_to_use$lambda,
              lambda2 = model_to_use$lambda2,
              alpha = model_to_use$alpha,
              c_hat = model_to_use$c_hat,
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









#k-fold CV function to be used inide a loop with "k" the number of folds
#create the output class of the svm_cross_validation function
setClass(Class = "svm_cross_validation_output",
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
#' @noRd
#' @keywords internal
psvm.cv <- function (df, X, y, C, lambda, tau, 
                     w_g, c_to_use, alpha, alpha_g, 
										 penalty_g, folds_i, k, 
										 regressors_betas,
                     trace, model_type, kernel_g, 
										 a1_g, a2_g, trend_g, 
										 gamma_start_input, 
										 gamma_start_default,
                     regressors_gammas, 
										 max_iter_g, delta_g, 
										 max_alpha_g, 
										 stepsizeShrink_g, 
										 min_alpha_g,
                     convergence_error_g,
                     auc, aauc, aYI, youden_index, 
										 sensitivity, specificity,
                     geometric_mean, fdr, mcc, 
										 corrclass, n_betas, n_gammas, 
										 used_cores, c_function_of_covariates,
                     simultaneous, run_aauc, log_file) {


  # This function performs the cross-validation for a single fold.
  # It is designed to be used internally by the psvm_compute_cv function.
  # test_i: indices of the test set for the current fold
  test_i <- which(folds_i == k)
  # train_df: data frame containing the training data
  train_df <- df[-test_i, ]
  # test_df: data frame containing the test data
  test_df <- df[test_i, ]

  greek <- if (c_function_of_covariates) "tau" else "lambda"
  if (trace %in% c(1, 2)) {
    cat("----------------------------------------------------------------\n")
    cat("|       starting with the", k, "-th fold for the CV of ", greek, "      |\n")
    cat("----------------------------------------------------------------\n")
  }

  if (model_type %in% c("SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM")) {

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
      parallel::clusterExport(cl, c("train_df", "X", "y", "model_type", "regressors_betas", "k", "trace", "lambda", "c_to_use"), envir = environment())
      fitted_models <- parallel::parLapply(cl, lambda, function(x) psvm_estimation(df = train_df, X = X, y = y,
                                                                                      lambda = x, lambda2 = 0.05,
                                                                                      alpha = alpha,
                                                                                      c_to_use = c_to_use,
                                                                                      regressors_betas = regressors_betas,
                                                                                      model_type = model_type,
                                                                                      fold = k, trace = trace))
    } else {

      fitted_models <- lapply(lambda, function(x) psvm_estimation(df = train_df, X = X, y = y,
                                                                     lambda = x, lambda2 = 0.05,
                                                                     alpha = alpha,
                                                                     c_to_use = c_to_use,
                                                                     regressors_betas = regressors_betas,
                                                                     model_type = model_type,
                                                                     fold = k, trace = trace))
    }
  }

  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat"))

  if (length(gamma_start_input) == 0) {
    # If gamma_start_input is not present, use the optimal c of the betas estimation as the starting point of the constant
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
                                                gamma_start_default = gamma_start_default,
                                                trace = trace,
                                                alpha = alpha_g, a1 = a1_g, a2 = a2_g,
                                                penalty = penalty_g,
                                                max_iter = max_iter_g,
                                                min_alpha = min_alpha_g,
                                                convergence_error = convergence_error_g,
                                                regressors_gammas = regressors_gammas, fold = k,
                                                trend = trend_g,
                                                stepsizeShrink = stepsizeShrink_g,
                                                delta = delta_g, max_alpha = max_alpha_g,
                                                kernel = kernel_g,
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
  all_measures_test <- mapply(function(z) psvm_predict(df = test_df,
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

  return(methods::new("svm_cross_validation_output", auc = auc,
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





#' @title Cross-Validation for Optimal Penalized svm Regularization Parameter
#' Selection
#'
#' @description Performs k-fold cross-validation for penalized SVM models,
#'   optionally incorporating covariate-adjusted cut-point optimization via
#'   the covYI method. This function estimates the optimal values of the
#'   regularization parameters lambda (for feature selection) and tau (for
#'   covariate adjustment of the classification threshold).
#'
#' @param n_folds Integer. The number of folds to use for cross-validation
#'   (must be >= 2).
#' @param df Data frame. The input dataset containing the target variable,
#'   regressors, and covariates.
#' @param X Character vector. Names of the regressor variables (features).
#'   Defaults to all columns in "df" that are not "y" or "C". Can also be a
#'   data frame containing only the regressors.
#' @param y Character. Name of the binary target variable (outcome), coded as
#'   0 and 1. Defaults to "y". Can also be a data frame containing only the
#'   target variable.
#' @param C Character vector. Names of the covariate variables for
#'   covariate-adjusted cut-point estimation. Defaults to "NULL" (no
#'   covariates). Can also be a data frame containing only the covariates.
#' @param lambda Numeric vector. The regularization parameter(s) for the
#'   regressors "X".
#' @param tau Numeric vector. The regularization parameter(s) for the
#'   covariates "C" in "covYI". If "c_function_of_covariates = TRUE", this
#'   parameter cannot be "NULL" or contain only zeros. Default is 0 (no
#'   penalization).
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the 
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by 
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param alpha numeric in (0, 1]. Elastic-net mixing parameter for sparseSVM
#'   (alpha = 1 corresponds to L1). Default is 0.5.
#' @param c_to_use Character or numeric. Specifies the classification
#'   cut-point used to convert continuous predictions into binary outcomes. If
#'   set to "Youden", the function computes the optimal cut-point using
#'   Youden's Index. If a numeric value is provided, it is used directly as 
#'   the threshold. If "NULL", a default is assigned based on 
#'   "model_type": "SCADSVM" uses "0", "ElasticSCADSVM", "l1SVM" or "enSVM" 
#'   use "Youden".
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   "NULL".
#' @param trace Integer. Level of verbosity: 0 = silent, 1 = brief summary, 2
#'   = verbose. Default is 1.
#' @param model_type Character. The SVM model to use. Options: "SCADSVM",
#'   "ElasticSCADSVM", "l1SVM", "enSVM".
#' @param seed Integer. Random seed for reproducibility. Default is 1.
#' @param used_cores Integer. Number of cores to use for parallel processing.
#'   If 1, no parallelization is used. Default is 1.
#' @param scaling Logical. If "TRUE", the dataset is scaled. Default is
#'   "FALSE".
#' @param c_function_of_covariates Logical. If "TRUE", "covYI" is used to
#'   estimate the cut-off point as a function of the covariates. Default is
#'   "FALSE".
#' @param simultaneous Logical. If "TRUE" (and
#'   "c_function_of_covariates = TRUE"), gammas are estimated simultaneously
#'   with betas. If "FALSE", gammas are estimated as a second step. Default
#'   is "FALSE".
#' @param measure_to_select_lambda Character. The performance measure used to
#'   select the best lambda when "simultaneous = FALSE". Options: "auc",
#'   "aauc", "aYI", "ccr", "yi", "gm", "fdr", "mcc", "sen",
#'   "spc". Default is "ccr" (correct classification rate).
#' @param alpha_g Numeric. Elastic-Net mixing parameter for "covYI". Default
#'   is 0.5.
#' @param penalty_g Character. The penalty of "covYI". Options: "L12",
#'   "L1", "EN", "SCAD", "MCP". Default is "L1".
#' @param kernel_g Character. Kernel type for density estimation in "covYI".
#'   Default is "gaussian".
#' @param a1_g Numeric. Parameter for the SCAD and MCP penalties in "covYI".
#'   Default is 3.7.
#' @param a2_g Numeric. Parameter for the MCP penalty in "covYI". Default is
#'   3.0.
#' @param trend_g Character. For "covYI", if "monotone", mmAPG is used; if
#'   "nonmonotone", mnmAPG is used. Default is "monotone".
#' @param gamma_start_input Numeric vector. A specific starting point for
#'   gammas. Default is "NULL".
#' @param gamma_start_default Character. Sets the default starting point of
#'   gamma. If "zeros", it starts with all zero values; if "corr", it
#'   starts with the value of the correlation of every regressor with the
#'   target variable.
#' @param regressors_gammas Numeric vector. Optional vector containing the true
#'   gammas (if known).
#' @param max_iter_g Integer. Maximum number of iterations in the algorithms
#'   mmAPG and mnmAPG in "covYI". Default is 10000.
#' @param delta_g Numeric. Parameter for the convergence condition of the
#'   optimization algorithm of "covYI". Default is 1e-5.
#' @param max_alpha_g Numeric. Maximum value of the step-parameter alpha in
#'   "covYI". Default is 100.
#' @param stepsizeShrink_g Numeric. Parameter to adjust the step-size in the
#'   backtracking line-search, in the optimization of "covYI". Taking values
#'   between 0 and 1, the closer to 1, the more accurate the estimation will
#'   be, the longer it will take and vice versa. Default is 0.8.
#' @param min_alpha_g Numeric. Minimum value of the step-parameter alpha in
#'   "covYI". Default is 1e-12.
#' @param convergence_error_g Numeric. Error to accept for considering the
#'   algorithm converged in "covYI". Default is 1e-7.
#' @param run_aauc Logical. If "FALSE", the aAUC and aYI are not computed, to
#'   save estimation time if not requested. Default is "FALSE".
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If "NULL", output goes to the console.
#'   Defaults to "log_SVM_models.txt".
#'
#' @return A list containing the optimal values of lambda (and possibly tau) to
#'   estimate betas (and possibly gammas) and the value of the main accuracy
#'   measure for all the folds. The list includes:
#'   \item{model_type}{The SVM model type used.}
#'   \item{c_to_use}{The classification cut-point used.}
#'   \item{penalty_g}{The penalty type used for covYI.}
#'   \item{cv_time}{The time taken for cross-validation.}
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
#'   \item{lambda_hat_yi, lambda_hat_auc, lambda_hat_aauc, lambda_hat_aYI,
#'     lambda_hat_ccr, lambda_hat_sen, lambda_hat_spc, lambda_hat_gm}{The
#'     optimal lambda values selected based on different performance
#'     measures.}
#'   \item{tau_hat_yi, tau_hat_auc, tau_hat_aauc, tau_hat_aYI, tau_hat_ccr,
#'     tau_hat_sen, tau_hat_spc, tau_hat_gm}{The optimal tau values
#'     selected based on different performance measures (if
#'     `c_function_of_covariates = TRUE`).}
#'   \item{c_function_of_covariates}{Indicates whether the cut-off point was
#'     estimated as a function of covariates.}
#'   \item{simultaneous}{Indicates whether betas and gammas were estimated
#'     simultaneously.}
#'   \item{measure_to_select_lambda}{The measure used to select lambda.}
#'   \item{lambda_star}{The optimal lambda value when using a sequential
#'     search.}
#'   \item{n_betas}{The number of non-zero beta coefficients.}
#'   \item{n_betas_star}{The number of non-zero beta coefficients when using
#'     lambda_star.}
#'   \item{n_gammas}{The number of non-zero gamma coefficients.}
#'   \item{betas_star}{The estimated beta coefficients when using
#'     lambda_star.}
#'   \item{betas}{The estimated beta coefficients.}
#'   \item{gammas}{The estimated gamma coefficients.}
#'
#' @examples
#' library(pye)
#' # Generate a sample dataset
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_betas <- sim_data$nregressors
#' regressors_gammas <- sim_data$ncovariates
#'
#' # Define parameters
#' lambda <- c(0.1, 0.05)
#' c_function_of_covariates <- FALSE
#' model_type <- "SCADSVM"
#'
#' # Run cross-validation
#' result <- psvm_compute_cv(n_folds = 2,
#'                           df = df,
#'                           X = X,
#'                           y = y,
#'                           C = C,
#'                           lambda = lambda,
#'                           regressors_betas = regressors_betas,
#'                           regressors_gammas = regressors_gammas,
#'                           c_function_of_covariates = c_function_of_covariates,
#'                           model_type = model_type
#' )
#'
#' # Print the results
#' print(result$model_type)
#' print(result$lambda_hat_ccr)
#'
#' @export
psvm_compute_cv <- function (n_folds, df, X = NULL, y = "y", C = NULL, lambda,
                             tau = 0,
														 w_g = 0.5, 
														 c_to_use = NULL, 
														 regressors_betas = NULL, trace = 1,
														 alpha = 0.5,
                             model_type = c("SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM"),
                             seed = 1, used_cores = 1, 
														 scaling = FALSE,
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
                             delta_g = 1e-5, 
														 max_alpha_g = 10000,
                             stepsizeShrink_g = 0.8, 
														 min_alpha_g = 1e-12,
                             convergence_error_g = 1e-7,
                             run_aauc = FALSE, log_file = "log_SVM_models.txt") {

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

  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty_g %in% valid_penalties)) {stop("Parameter 'penalty_g' must be one of: ", paste(valid_penalties, collapse = ", "))}
	if (!is.numeric(w_g) || length(w_g) != 1 || w_g < 0 || w_g > 1) {stop("Parameter 'w_g' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(lambda) || length(lambda) < 1 || any(lambda < 0)) stop("Parameter 'lambda' must be a numeric vector of non-negative values.")
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter trace has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(model_type %in% c("SCADSVM", "ElasticSCADSVM", "l1SVM", "enSVM"))) { stop("The parameter model_type is wrongly specified")}
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
  if (is.null(c_to_use)) {c_to_use <- switch(model_type, "SCADSVM" = 0, "ElasticSCADSVM" = "Youden", "l1SVM" = "Youden", "enSVM" = "Youden")}
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
    assign(mes, lapply(lambda, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)),
                                                 test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))))
    eval(substitute(names(x) <- lambdanames, list(x = as.symbol(mes))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(foldnames, unlist(lapply(lambda, function (x) paste("lambda", x, sep = "=")))))
  n_gammas <- lapply(lambda, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
  names(n_gammas) <- lambdanames

  if (trace > 0) {
    cat("Starting CV with model_type:", model_type, "\n")
  }

  # fill the matrices
  results <- mapply(function(k) psvm.cv(df = df1[, names(df1) != "ID", drop = FALSE], 
	                                      X = X, y = y, C = C, 
																				lambda = lambda, tau = tau, w_g = w_g,
                                        alpha = alpha,
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

    # now we execute psvm.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
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
    results <- mapply(function(k) psvm.cv(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, C = C,
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
    list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity",
                           "geometric_mean", "fdr", "mcc", "corrclass")

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
    cat("------------------> END OF THE CROSS-VALIDATION OF THE PEN. SVM METHOD <-------------------- \n")
    cat("-------------> For the whoole Cross-validation it took:", cv_time, "minutes <---------------- \n")
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