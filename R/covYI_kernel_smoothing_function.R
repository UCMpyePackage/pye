#' @title Covariate-Adjusted Youden Index (covYI) method with Kernel Smoothing
#' Estimation of the Density Functions
#'
#' @description This function implements the Covariate-Adjusted Youden Index
#' with Kernel Smoothing (covYI_KS). It calculates the Youden Index,
#' incorporating covariate information through a kernel smooth density estimator,
#' and can apply various penalization types (L1 / 2, L1, Elastic-Net (EN), SCAD,
#' and MCP) to the covariate coefficients.
#' The function returns the penalized Youden Index values, along with key
#' classification accuracy measures, and the gradient of the Youden Index with
#' respect to the covariate coefficients.
#'
##' @param df A data frame containing the input data, including the latent
#'   variable/score (`z`), the binary target variable (`y`), and covariates (`C`).
#' @param z Character string. The column name in `df` representing the latent
#'   variable or score (e.g., from a primary prediction model). Default is "z_hat".
#' @param y Character string. The column name in `df` representing the binary
#'   target variable (0 or 1). Default is "y".
#' @param C Character vector. Column names from `df` to be used as covariate
#'   variables for the cut-point `c`. This vector must include "const" if an
#'   intercept is desired for the covariate combination.
#' @param gammas Numeric vector. The coefficients of the covariates (`C`) used
#'   to evaluate the cut-point `c`. The first element of this vector is assumed
#'   to correspond to the constant (intercept) term if "const" is in `C`.
#' @param tau Numeric. The penalization parameter for the covariates `C` in
#'   the `covYI` calculation. Default is 0, indicating no penalization. Must
#'   be a single non-negative value.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the 
#'   Weighted Youden Index. Sensitivity is weighted by `w` and specificity by 
#'   `1 - w`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param kernel Character string. The kernel type to use for the estimation
#'   of the density function in the Kernel Smooth method. Supported options
#'   include "gaussian", "normal", "uniform", "rectangular", "triangular",
#'   "epanechnikov", "biweight", "triweight", "tricube", "parzen", "cosine",
#'   and "optcosine". Default is "gaussian".
#' @param alpha Numeric. A value between 0 and 1, representing the elastic-net
#'   mixing parameter for the EN penalization term. Default is 0.5.
#' @param a1 Numeric. The `a` parameter for the SCAD penalization term.
#'   Default is 3.7.
#' @param a2 Numeric. The `a` parameter for the MCP penalization term.
#'   Default is 3.0.
#' @param penalty Character string. The type of penalization to apply. Must be
#'   one of "L12", "L1", "EN", "SCAD", or "MCP". Default is "L1".
#' @param h_exponent Numeric. The exponent used in the bandwidth calculation
#'   for the Kernel Smooth density estimation. Default is 0.2.
#' @param prediction Logical. If `TRUE`, the Youden Index (YI) returned is the
#'   empirical maximum Youden Index derived from `OptimalCutpoints` for the
#'   given `z` and `y`. If `FALSE`, it's the `yi` calculated from `f0 - f1`.
#'   Default is `FALSE`.
#' @param run_aauc Logical. If `FALSE`, the covariate-adjusted AUC (aAUC) and
#'   adjusted Youden Index (aYI) are not computed, which can save computation
#'   time if these measures are not required. Default is `FALSE`.
#'
#' @return A list containing the results of the `covYI_KS` calculation and
#'   related performance measures.
#'   \item{covYI_KS_L12}{Numeric. The Penalized Youden Index (pye) using the
#'     L1 / 2 penalty.}
#'   \item{covYI_KS_L1}{Numeric. The Penalized Youden Index (pye) using the
#'     L1 penalty.}
#'   \item{covYI_KS_EN}{Numeric. The Penalized Youden Index (pye) using the
#'     Elastic-Net penalty.}
#'   \item{covYI_KS_SCAD}{Numeric. The Penalized Youden Index (pye) using the
#'     SCAD penalty.}
#'   \item{covYI_KS_MCP}{Numeric. The Penalized Youden Index (pye) using the
#'     MCP penalty.}
#'   \item{gr_yi}{Numeric vector. The gradient of the Youden Index (`yi`) with
#'     respect to the `gammas` coefficients (excluding the constant term).}
#'   \item{youden_index}{Numeric. The unpenalized Youden Index calculated
#'     using the Kernel Smooth density estimator (or the empirical Youden Index
#'     if \code{prediction = TRUE}). It is weighted if w is not 0.5.}
#'   \item{sensitivity}{Numeric. The sensitivity (True Positive Rate) of the
#'     classification.}
#'   \item{specificity}{Numeric. The specificity (True Negative Rate) of the
#'     classification.}
#'   \item{geometric_mean}{Numeric. The geometric mean of sensitivity and specificity.}
#'   \item{fdr}{Numeric. The False Discovery Rate.}
#'   \item{mcc}{Numeric. The Matthews Correlation Coefficient.}
#'   \item{auc}{Numeric. The Area Under the ROC Curve (AUC) from `OptimalCutpoints`.}
#'   \item{aauc}{Numeric. The covariate-adjusted Area Under the ROC Curve (aAUC).
#'     This is 0 if `run_aauc` is `FALSE` or if computation fails.}
#'   \item{aYI}{Numeric. The adjusted Youden Index (aYI). This is 0 if
#'     `run_aauc` is `FALSE` or if computation fails.}
#'   \item{corrclass}{Numeric. The overall correct classification rate.}
#'   \item{z_hat}{Data frame. A data frame with 'ID' and the `z` values used
#'     in the calculation.}
#'   \item{c_hat}{Data frame. A data frame with 'ID' and the estimated cut-point
#'     values (`c_hat`) for each observation.}
#'   \item{y_hat}{Data frame. A data frame with 'ID' and the predicted binary
#'     outcomes (`y_hat`).}
#'   \item{TP}{Numeric. Number of True Positives.}
#'   \item{TN}{Numeric. Number of True Negatives.}
#'   \item{FP}{Numeric. Number of False Positives.}
#'   \item{FN}{Numeric. Number of False Negatives.}
#'   \item{input_data}{List. A list containing the input parameters used for
#'     this function call: `gammas`, `tau`, `w`, `alpha`, `a1`, `a2`, and `kernel`.}
#'
#' @examples
#' library(pye)
#'
#' # Simulate a sample dataset with covariates
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#'
#' # Extract necessary components from the simulated data
#' df <- sim_data$train_df_scaled
#' y <- sim_data$y
#' C <- sim_data$C
#'
#' # The combination of biomarkers (regressors) is known when we apply covYI
#' z <- sim_data$z
#'
#' #Merge 'df' and 'z' by row names, adding 'z' as a column and preserving row names.
#' df$ID <- rownames(df)
#' df1 <- merge(x = df, y = z, by = "ID", all.x = TRUE)
#' rownames(df1) <- df1$ID
#' df1$ID <- NULL
#'
#' # Input some variables
#' penalty <- "L12"
#' tau <- 0.1
#' gammas <- rep(1, length(C))
#' c <- 0
#'
#' # Run covYI_KS to evaluate pye for the given gammas and penalty
#' covYI_result <- covYI_KS(df = df1, z = "z", y = y, C = C, gammas = gammas, tau = tau,
#'   alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty)
#' print(covYI_result)
#'
#' # Example with a different penalty
#' penalty <- "SCAD"
#'
#' # Run covYI_KS to evaluate pye for the given gammas and penalty
#' covYI_result <- covYI_KS(df = df1, z = "z", y = y, C = C, gammas = gammas, tau = tau,
#'   alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty)
#' print(covYI_result)
#'
#' @importFrom evmix kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom ROCnReg AROC.sp compute.threshold.AROC
#' @importFrom plyr join
#' @importFrom stats sd IQR dnorm
#' @export
covYI_KS <- function(df, 
                     z = "z_hat", 
										 y = "y", 
										 C, 
										 gammas, 
										 tau,
										 w = 0.5,
										 kernel = "gaussian", 
										 alpha = 0.5, 
										 a1 = 3.7, 
										 a2 = 3,
                     penalty = "L1", 
										 h_exponent = 0.2, 
										 prediction = FALSE, 
										 run_aauc = FALSE) {

  # IMPORTANT: The first element of 'gammas' is assumed to be the constant term if 'C' includes "const".

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

  # Handle z parameter
  if (inherits(z, "data.frame")) {
    if (length(names(z)) > 1) {
      warning("Parameter 'z' is a data.frame with more than one column. Only the first column ('", names(z)[1], "') will be used.")
    }
    z <- names(z)[1]
  } else if (!inherits(z, "character") || length(z) != 1) {
    stop("Parameter 'z' must be a single column name (character).")
  }
  if (length(z) == 0) {
    stop("No 'z' specified or found in the data frame.")
  }
  if (!(z %in% names(df))) {
    stop("Variable 'z' ('", z, "') not found in the input data frame 'df'.")
  }

  # Handle y parameter
  if (inherits(y, "data.frame")) {
    if (length(names(y)) > 1) {
      warning("Parameter 'y' is a data.frame with more than one column. Only the first column ('", names(y)[1], "') will be used.")
    }
    y <- names(y)[1]
  } else if (!inherits(y, "character") || length(y) != 1) {
      stop("Parameter 'y' must be a single column name (character).")
  }
  if (!(y %in% names(df))) {
    stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")
  }

  # Handle C parameter
  if (inherits(C, "data.frame")) {
    C <- names(C)
  } else if (!inherits(C, "character")) {
    stop("Parameter 'C' must be a character vector of column names, a data.frame, or NULL.")
  }
  if (is.null(C) || length(C) == 0) {
    stop("C cannot be NULL or empty.")
  }

  # Validate target variable y (0, 1)
  if (is.factor(df[[y]]) || is.character(df[[y]])) {df[[y]] <- as.numeric(as.character(df[[y]]))}
  if (!all(sort(unique(df[[y]])) %in% c(0, 1)) || anyNA(df[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Validate tau
  if (length(tau) != 1 || !is.numeric(tau) || tau < 0) {stop("Parameter 'tau' must be a single non-negative numeric value.")}
	
	if (!is.numeric(w) || w < 0 || w > 1) {stop("The weight 'w' must be a numeric value between 0 and 1.")}

  if (!(penalty %in% c("L12", "L1", "EN", "SCAD", "MCP"))) {stop("A wrong value has been assigned to the parameter penalty.")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))) {
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(gammas) != length(C)) {
    stop("The number of element of gammas is different then the number of columns in C")
  }

  # Check for "ID" column conflict: 'ID' is used internally
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}
  ID <- rownames(df)

  # Check if all C column names (excluding 'const' if present) are in df
  C_without_const <- C[C != "const"]
  if (!all(C_without_const  %in% names(df))) {
    stop("Not all specified covariates ('C') are found in the input data frame 'df'.")
  }

  # Create "const" if needed
  if (("const" %in% C) && (!("const" %in% names(df)))) {
    const <- rep(1, length(ID)) #the constant for the coefficients gammas
    df1 <- cbind(ID, df[, c(y, z), drop = FALSE], const, df[, (names(df) %in% C), drop = FALSE])
  } else {
    df1 <- cbind(ID, df[, c(y, z, C), drop = FALSE])
  }

  # Divide control and diseased groups based on 'y'
  z_y0 <- as.matrix(df1[df1[[y]] == 0, (names(df1) %in% z), drop = FALSE])#drop = FALSE permits it to remain a data frame even if it is composed by a single column
  z_y1 <- as.matrix(df1[df1[[y]] == 1, (names(df1) %in% z), drop = FALSE])

  # Calculate c_hat for each observation
  c_y0 <- as.matrix(df1[df1[[y]] == 0, (names(df1) %in% C), drop = FALSE]) %*% gammas
  c_y1 <- as.matrix(df1[df1[[y]] == 1, (names(df1) %in% C), drop = FALSE]) %*% gammas

  # Kernel density estimation for f0 (healthy) and f1 (diseased)
  # For healthy group
  # if (kernel %in% c("gaussian", "epanechnikov", "triweight", "tricube", "biweight", "cosine") & sum(z_y0) != 0) {
  #   h0 = kedd::h.bcv(as.numeric(c-z_y0), kernel = kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h0 <- 0.9 * min(stats::sd(z_y0), stats::IQR(z_y0) / 1.34) * length(z_y0)^(-h_exponent)
  # }
  if (h0 < 0.1) {h0 <- 0.1} # Safeguard against NaN or too small bandwidth
  t0 <- (c_y0 - z_y0) / h0
  f0 <- sum(evmix::kpz(z = as.numeric(t0), kernel = kernel)) / length(z_y0)
  #f0 <- sum(pnorm(t0, 0, 1)) / length(z_y0)

  # For diseased group
  # if (kernel %in% c("gaussian", "epanechnikov", "triweight", "tricube", "biweight", "cosine") & sum(z_y1) != 0) {
  #   h1 <- kedd::h.bcv(as.numeric(z_y1), kernel = kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h1 <- 0.9 * min(stats::sd(z_y1), stats::IQR(z_y1) / 1.34) * length(z_y1)^(-h_exponent)
  # }
  if (h1 < 0.1) {h1 <- 0.1} # Safeguard against NaN or too small bandwidth
  t1 <- (c_y1 - z_y1) / h1
  f1 <- sum(evmix::kpz(z = as.numeric(t1), kernel = kernel)) / length(z_y1)
  #f1 <- sum(pnorm(t1, 0, 1)) / length(z_y1)

  # Youden Index based on kernel densities
  #yi <- f0 - f1
	yi <- 2 * (1 - w) * f0 - 2 * w * f1 + 2 * w - 1 #formula from Wang et. al. 2025

  # Evaluate the gradient of Youden Index with respect to gammas
  # If 'const' is the first element of C, its gradient is handled implicitly.
  # The 'gr_yi' will have length equal to 'length(gammas)'.
  #h0 <- 0.9 * min(stats::sd(z_y0), stats::IQR(z_y0) / 1.34) * (length(z_y0)^(-h_exponent))
  #t0 <- (c_y0 - z_y0) / h0
  # Gradient for f0 (healthy group)
  # Since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")) {
    gr0_gammas <- apply (df1[df1[[y]] == 0, (names(df1) %in% C), drop = FALSE], 2,
                        function(g) (t(stats::dnorm(t0, 0, 1)) %*% (g / h0)) / length(z_y0))
  } else {
    gr0_gammas <- apply (df1[df1[[y]] == 0, (names(df1) %in% C), drop = FALSE], 2,
            function(g) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% (g / h0)) / length(z_y0))
  }

  #h1=0.9 * min(stats::sd(z_y1), stats::IQR(z_y1) / 1.34) * (length(z_y1)^(-h_exponent))
  #t1 <- (c_y1 - z_y1) / h1
  # Gradient for f1 (diseased group)
  # Since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")) {
    gr1_gammas <- apply (df1[df1[[y]] == 1, (names(df1) %in% C), drop = FALSE], 2,
                        function(g) (t(stats::dnorm(t1, 0, 1)) %*% (g / h1)) / length(z_y1))
  } else {
    gr1_gammas <- apply (df1[df1[[y]] == 1, (names(df1) %in% C), drop = FALSE], 2,
           function(g) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% (g / h1)) / length(z_y1))
  }

  # Total gradient of Youden Index
  #gr_yi <- gr0_gammas - gr1_gammas
	gr_yi <- 2 * (1 - w) * gr0_gammas - 2 * w * gr1_gammas
	

  # Ensure gr_yi is a named numeric vector matching 'gammas' structure
  gr_yi <- as.numeric(gr_yi)
  names(gr_yi) <- C

  #if (abs(sum(gr_yi)) < 0.0001) {
  #  gr_yi <- gr_yi*(10^(round(abs(log10(abs(gr_yi)) + 1))))
  #}

  # Append calculated c to df1
  c <- as.data.frame(rbind(c_y0, c_y1))
  c['ID'] <- as.numeric(rownames(c))
  #c <- c[order(c$ID), ]
  names(c)[names(c) == "V1"] <- "c_hat"
  df1 <- plyr::join(df1, c, by = "ID", type = "left")
  rownames(df1) <- df1$ID # Restore original row names

  # Find the optimal cut-point using OptimalCutpoints for empirical measures
  #opt <- OptimalCutpoints::optimal.cutpoints(data = df1, X = z, status = y, methods = "Youden", tag.healthy = 0)
	opt <- OptimalCutpoints::optimal.cutpoints(data = df1, X = z, status = y, methods = "Youden", tag.healthy = 0, control = OptimalCutpoints::control.cutpoints(generalized.Youden = TRUE,  CFN = max(w, 1e-6), CFP = max(1 - w, 1e-6)))

  # Compute y_hat based on z and c_hat
  df1["y_hat"] <- ifelse(df1$z > df1$c_hat, 1, 0)

  # Confusion matrix components based on c
  TP <- sum(ifelse(z_y1 >= c_y1, 1, 0))
  TN <- sum(ifelse(z_y0 < c_y0, 1, 0))
  FP <- sum(ifelse(z_y0 >= c_y0, 1, 0))
  FN <- sum(ifelse(z_y1 < c_y1, 1, 0))

  # Compute standard classification measures
  spec <- TN / (TN + FP)
  fnr <- FN / (FN + TP)
  #yi1 <- spec - fnr

  # Sensitivity
  sensitivity <- 1 - fnr

  # Geometric mean
  gm <- sqrt(spec * sensitivity)

  # FDR
  fdr <- FP / (FP + TP)
  # Handle division by zero for FDR if FP + TP is zero
  if (is.nan(fdr) && (FP + TP) == 0) fdr <- 0

  # MCC
  mcc <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  # Correct classification:
  corrclass <- (TP + TN) / nrow(df1)

  # AUC and Youden Index from OptimalCutpoints for the overall dataset
  #try({test <- AUC::auc(AUC::roc(df1$z_hat,factor(df1[[y]])))}, silent = TRUE)
  auc <- mean(opt$Youden$Global$measures.acc$AUC)
  yi1 <- mean(opt$Youden$Global$optimal.criterion)

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.sp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                      p = seq(0, 1, l = 101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y", tag.healthy = 0,
  #                                     data = df1, p = seq(0, 1, l = 101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.bnp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1, p = seq(0, 1, l = 101),
  #                           compute.lpml = FALSE, compute.WAIC = FALSE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y", tag.h = 0, data = df1,
  #                              p = seq(0, 1, l = 101), B = 500)

  #start_time2 <- Sys.time()
  #this is longher then AROC.bsp
  #tmp2 <- try({aroc <- AROC::AROC.bnp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                   scale = TRUE, p = seq(0, 1, l = 101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                   a = 2, b = 0.5, L = 10, nsim = 5000, nburn = 1000)}, silent = TRUE)
  #cat("temp2:")
  #print(tmp2)
  #print(difftime(Sys.time(), start_time2 , units = "mins"))
  #start_time1 <- Sys.time()

  #this does not have the YI
  #aroc <- AROC::AROC.bsp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                scale = TRUE, p = seq(0, 1, l = 101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                a = 2, b = 0.5, nsim = 5000, nburn = 1500)

  # Adjusted AUC (aAUC) and adjusted Youden Index (aYI) using ROCnReg
  aauc <- 0
  aYI <- 0
  if (run_aauc) {
    tmp <- try({

      #AROC.sp function estimates the covariate-adjusted ROC curve (AROC) using the semiparametric approach
      #proposed by Janes and Pepe (2009).

      #if we have all zeros, the function returns a warning message. We suppress it just for this function:
      options(warn = -1)

      aroc <- ROCnReg::AROC.sp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1,
                               p = seq(0, 1, l = 101), B = 500)

      options(warn = 0)

      #NOt werking properly:
      #aroc <- ROCnReg::AROC.bnp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1,
      #                         p = seq(0, 1, l = 101))

      #aroc <- ROCnReg::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y",
      #                             tag.h = 0, data = df1, bw = "LS", regtype = "LC",
      #                             pauc = ROCnReg::pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
      #                             B = 500)

      ### Threshold values based on the YI
      th_AROC.sp <- ROCnReg::compute.threshold.AROC(aroc, criterion = "YI")
      aauc <- aroc$AUC[1]
      aYI <- th_AROC.sp$YI
      #summary.AROC(aroc)

    }, silent = TRUE)
    #print(difftime(Sys.time(), start_time1 , units = "mins"))
    if (inherits(tmp, "try-error")) {
      aauc <- 0
      aYI <- 0
    }
  }

  #if all the betas are zero, the measures are zeros
  #if (sum(gammas) == 0) {
  #  spec <- 0
  #  fnr <- 0
  #  sensitivity <- 0
  #  gm <- 0
  #  yi1 <- 0
  #  fdr <- 0
  #  mcc <- 0
  #  corrclass <- 0
  #}

  # Penalization terms
  # L_(1 / 2) penalization
  phi_L12_g <- tau * sum(abs(gammas)^(1 / 2))
  if (is.na(phi_L12_g)) {stop('Problem with the L1 / 2 penalization (phi_L12_g) in covYI_KS function. Check gammas or tau.')}

  # L1 penalization
  phi_L1_g <- tau * sum(abs(gammas))
  if (is.na(phi_L1_g)) {stop('Problem with the L1 penalization (phi_L1_g) in covYI_KS function. Check gammas or tau.')}

  # Elastic-Net penalization
  phi_EN_g <- tau * ((alpha * sum(abs(gammas))) + ((1 - alpha) / 2) * (sum(gammas^2)))
  if (is.na(phi_EN_g)) {stop('Problem with the Elastic-Net penalization (phi_EN_g) in covYI_KS function. Check gammas, tau, or alpha.')}

  # SCAD penalization
  phi_SCAD_g <- SCAD_function(gammas, tau, a = a1)
  if (is.na(phi_SCAD_g)) {stop('Problem with the SCAD penalization (phi_SCAD_g) in covYI_KS function. Check gammas, tau, or a1.')}

  # MCP penalization
  phi_MCP_g <- MCP_function(gammas, tau, a = a2)
  if (is.na(phi_MCP_g)) {stop('Problem with the MCP penalization (phi_MCP_g) in covYI_KS function. Check gammas, tau, or a2.')}

  # Join yi and the penalty functions to get covYI values
  covYI_KS_L12 <- yi - phi_L12_g
  covYI_KS_L1 <- yi - phi_L1_g
  covYI_KS_EN <- yi - phi_EN_g
  covYI_KS_SCAD <- yi - phi_SCAD_g
  covYI_KS_MCP <- yi - phi_MCP_g

  # If in prediction mode, use the Youden Index from OptimalCutpoints
  if (prediction == TRUE) {
    yi <- yi1
  }

  return(list(covYI_KS_L12 = covYI_KS_L12,
              covYI_KS_L1 = covYI_KS_L1,
              covYI_KS_EN = covYI_KS_EN,
              covYI_KS_SCAD = covYI_KS_SCAD,
              covYI_KS_MCP = covYI_KS_MCP,
              gr_yi = gr_yi,
              youden_index = yi,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr,
              mcc = mcc,
              auc = auc,
              aauc = aauc,
              aYI = aYI,
              corrclass = corrclass,
              z_hat = df1[, c("ID", z), drop = FALSE],
              c_hat = df1[, c("ID", "c_hat"), drop = FALSE] ,
              y_hat = df1[, c("ID", "y_hat"), drop = FALSE],
              TP = TP, TN = TN, FP = FP, FN = FN,
              input_data = list(gammas = gammas, tau = tau, w = w, alpha = alpha, a1 = a1, a2 = a2, prediction = prediction, kernel = kernel)))
}









#' @title Estimation of Optimal Covariate Coefficients through Penalized
#' Covariate-Adjusted Youden Index
#'
#' @description This function estimates the optimal values of the `gammas`
#' coefficients by maximizing the Penalized Youden Index (pye) function. To
#' achieve this, it uses either the monotonic (mmAPG) or non-monotonic (mnmAPG)
#' Accelerated Proximal Gradient algorithms for nonconvex programming.
#'
#' @param df The input dataset as a data frame.
#' @param z A single regressor or a known combination of regressors. It should
#'   be a character string specifying the column name in `df`. Default is "z_hat".
#' @param y The target variable, which must be a binomial (0 or 1) variable. It
#'   should be a character string specifying the column name in `df`. Default is "y".
#' @param C A character vector of column names from `df` to be used as covariate
#'   variables. It is mandatory and cannot be `NULL` or empty.
#' @param tau The penalization parameter applied to the covariates. Default is 0,
#'   which corresponds to no penalization.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the 
#'   Weighted Youden Index. Sensitivity is weighted by `w` and specificity by 
#'   `1 - w`. Default is 0.5.
#' @param penalty The type of penalization to apply. Must be one of "L12", "L1",
#'   "EN", "SCAD", or "MCP". Default is "L1".
#' @param gamma_start_input A numeric vector of specific starting points for the
#'   `gammas` coefficients. Default is `NULL`.
#' @param gamma_start_default A character string to set the default starting
#'   point for `gammas`. If "zeros", it starts with all zero values. If "corr",
#'   it starts with the correlation of each regressor with the target variable.
#'   Default is "zeros".
#' @param alpha The elastic-net mixing parameter, for use with the "EN" penalty.
#'   A value between 0 and 1. Default is 0.5.
#' @param a1 The `a` parameter for the SCAD penalization term. Default is 3.7.
#' @param a2 The `a` parameter for the MCP penalization term. Default is 3.0.
#' @param regressors_gammas A vector containing the true gamma coefficients (if known).
#'   Used for diagnostic purposes. Default is `NULL`.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param max_iter The maximum number of iterations for the optimization algorithms.
#'   Default is 10000.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#' @param trend A character string specifying the optimization algorithm to use.
#'   If "monotone", the mmAPG algorithm is used. If "nonmonotone", mnmAPG is used.
#'   Default is "monotone".
#' @param delta The parameter for the convergence condition of the optimization
#'   algorithm. Default is 1e-5.
#' @param max_alpha The maximum value of the step-size parameter alpha.
#'   Default is 10000.
#' @param stepsizeShrink The parameter to adjust the step-size in the
#'   backtracking line search. Taking values between 0 and 1, the closer to 1,
#'   the more accurate the estimation will be. Default is 0.8.
#' @param min_alpha The minimum value of the step-size parameter alpha.
#'   Default is 1e-10.
#' @param convergence_error The error tolerance for considering the algorithm to
#'   have converged. Default is 1e-7.
#' @param trace An integer controlling the level of output. 2: visualize all
#'   steps, 1: visualize just the final result, 0: visualize nothing. Default is 1.
#' @param seed An integer to fix the random seed. Default is 1.
#' @param kernel The kernel type to use for the estimation of the density
#'   function. Currently only "gaussian" is fully supported. Default is "gaussian".
#' @param run_aauc A logical value. If `FALSE`, the aAUC and aYI metrics are not
#'   computed to save estimation time. Default is `FALSE`.
#'
#' @return A list containing the optimal value of gammas and other key metrics.
#'   \item{gammas_hat_penalty}{The vector of optimal gamma coefficients found by 
#'     the algorithm.}
#'   \item{covYI_KS_penalty}{The final maximum penalized Youden Index (pye) value.}
#'   \item{gr_yi}{The gradient of the Youden Index with respect to the gammas 
#'     at the optimal solution.}
#'   \item{tau}{The input penalization parameter.}
#'   \item{penalty}{The type of penalty used.}
#'   \item{gammas_start}{The starting point for the gamma coefficients.}
#'   \item{kernel}{The kernel type used.}
#'   \item{c_hat}{A data frame with the estimated cut-point for each observation.}
#'   \item{z_hat}{A data frame with the regressor values used in the calculation.}
#'   \item{y_hat}{A data frame with the predicted binary outcomes.}
#'   \item{youden_index}{The Youden Index value at the optimal solution.}
#'   \item{sensitivity}{The sensitivity (True Positive Rate) at the optimal 
#'     cut-point.}
#'   \item{specificity}{The specificity (True Negative Rate) at the optimal 
#'     cut-point.}
#'   \item{geometric_mean}{The geometric mean of sensitivity and specificity.}
#'   \item{fdr}{The False Discovery Rate.}
#'   \item{mcc}{The Matthews Correlation Coefficient.}
#'   \item{auc}{The Area Under the ROC Curve (AUC).}
#'   \item{aauc}{The covariate-adjusted AUC (aAUC). Only computed if `run_aauc` 
#'     is `TRUE`.}
#'   \item{aYI}{The adjusted Youden Index (aYI). Only computed if `run_aauc` 
#'     is `TRUE`.}
#'   \item{corrclass}{The overall correct classification rate.}
#'   \item{TP}{Number of True Positives.}
#'   \item{TN}{Number of True Negatives.}
#'   \item{FP}{Number of False Positives.}
#'   \item{FN}{Number of False Negatives.}
#'   \item{n_gammas}{The number of non-zero gamma coefficients in the final 
#'     estimation.}
#'   \item{n_total_var_gammas}{The total number of gamma coefficients estimated.}
#'   \item{n_predicted_zeros_gammas}{The number of estimated gammas that are zero.}
#'   \item{n_predicted_non_zeros_gammas}{The number of estimated gammas that are 
#'     non-zero.}
#'   \item{n_caught_gammas}{Number of true non-zero gammas correctly identified as 
#'     non-zero (if `regressors_gammas` is provided).}
#'   \item{n_non_caught_gammas}{Number of true non-zero gammas incorrectly identified 
#'     as zero.}
#'   \item{n_caught_zero_gammas}{Number of true zero gammas correctly identified as 
#'     zero.}
#'   \item{n_zero_not_caught_gammas}{Number of true zero gammas incorrectly identified 
#'     as non-zero.}
#'   \item{input_parameters}{\code{character vector}. A list containing the 
#'     input parameters.}
#'   \item{estimation_time}{The time taken for the estimation, in minutes.}
#'   \item{niter}{The number of iterations performed by the optimization algorithm.}
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_gammas <- sim_data$ncovariates
#'
#' # The combination of biomarkers (regressors) is known when we apply covYI
#' z <- sim_data$z
#'
#' #Merge 'df' and 'z' by row names, adding 'z' as a column and preserving row names.
#' df$ID <- rownames(df)
#' df1 <- merge(x = df, y = z, by = "ID", all.x = TRUE)
#' rownames(df1) <- df1$ID
#' df1$ID <- NULL
#'
#' tau <- 0.04
#'
#' covYI_estimation_result <- covYI_KS_estimation(df = df1, z = "z", y = y,
#'   C = C, tau = tau, penalty = "SCAD", trace = 2,
#'   gamma_start_default = "zeros", regressors_gammas = regressors_gammas,
#'   max_iter = 5, run_aauc = TRUE)
#'
#' print(covYI_estimation_result)
#'
#' @importFrom stats setNames cor
#' @export
#Estimation of the parameter using covYI_KS
covYI_KS_estimation <- function(df, z = "z_hat", y = "y", C, tau, w = 0.5, penalty = "L1", gamma_start_input = NULL,
                                gamma_start_default = "zeros", alpha = 0.5, a1 = 3.7, a2 = 3,
                                regressors_gammas = NULL, fold = NULL, max_iter = 10000, max.print = 10,
                                trend = "monotone", delta = 1e-5, max_alpha = 10000, stepsizeShrink = 0.8,
                                min_alpha = 1e-10, convergence_error = 1e-7, trace = 1, seed = 1, kernel = "gaussian",
                                run_aauc = FALSE) {


  start_time <- Sys.time()

	if (!is.numeric(max.print) || length(max.print) != 1 || max.print <= 0) {
    stop("The parameter 'max.print' must be a single positive integer.")
  }
	old_options <- options(max.print = max.print)
  on.exit(options(old_options))

  # --- Input Parameter Validation and Standardization ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")

  # Handle z parameter
  if (inherits(z, "data.frame")) {
    if (length(names(z)) > 1) {
      warning("Parameter 'z' is a data.frame with more than one column. Only the first column ('", names(z)[1], "') will be used.")
    }
    z <- names(z)[1]
  } else if (!inherits(z, "character") || length(z) != 1) {
    stop("Parameter 'z' must be a single column name (character).")
  }
  if (length(z) == 0) {
    stop("No 'z' specified or found in the data frame.")
  }
  if (!(z %in% names(df))) {
    stop("Variable 'z' ('", z, "') not found in the input data frame 'df'.")
  }

  # Handle y parameter
  if (inherits(y, "data.frame")) {
    if (length(names(y)) > 1) {
      warning("Parameter 'y' is a data.frame with more than one column. Only the first column ('", names(y)[1], "') will be used.")
    }
    y <- names(y)[1]
  } else if (!inherits(y, "character") || length(y) != 1) {
      stop("Parameter 'y' must be a single column name (character).")
  }
  if (!(y %in% names(df))) {
    stop("The target variable 'y' ('", y, "') is not found in the input data frame 'df'.")
  }

  # Handle C parameter
  if (inherits(C, "data.frame")) {
    C <- names(C)
  } else if (!inherits(C, "character")) {
    stop("Parameter 'C' must be a character vector of column names, a data.frame, or NULL.")
  }
  if (is.null(C) || length(C) == 0) {
    stop("C cannot be NULL or empty.")
  }

  #check if const already exists in the dataset
  if ("const" %in% colnames(df)) {stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal const")}
  if ("const" %in% colnames(C)) {stop("const already exists as column in C! Please delete or rename this column since I need this name to set the internal const")}

  # Check if all C column names are in df
  if (!all(C %in% names(df))) {
    stop("Not all specified covariates ('C') are found in the input data frame 'df'.")
  }

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  const <- rep(1, length(ID)) #the constant for the coefficients gammas
  df1 <- cbind(ID, df[, c(y, z), drop = FALSE], const, df[, C, drop = FALSE])

  #adding "const" in C
  C1 <- c("const", C)

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Validate other parameters
  if (length(tau) != 1 || !is.numeric(tau) || tau < 0) {stop("Parameter 'tau' must be a single non-negative numeric value.")}
	if (!is.numeric(w) || w < 0 || w > 1) {stop("The weight 'w' must be a numeric value between 0 and 1.")}
  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 2 || !is.integer(as.integer(max_iter))) {stop("Parameter 'max_iter' needs to be an integer and at least 2.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter 'trace' has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel %in% valid_kernels)) {stop(paste0("Parameter 'kernel' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}

  gammas1_initial_corr <- NULL

  names_gammas <- C1
  gammas_start <- NULL

  #initializing the parameters (gammas)
  if (is.null(gamma_start_input) || length(gamma_start_input) == 0) {
    if (gamma_start_default == "zeros") {
      gammas_start <- rep(0, length(C1))
    } else if (gamma_start_default == "corr") {

      #compute the corr between every x and y
      corr.xy <- data.frame(matrix(ncol = ncol(df1[, (names(df1) %in% C1)]), nrow = 0))
      names(corr.xy) <- names_gammas
      corr.xy[1:3, ] <- t(sapply(c("pearson", "kendall", "spearman"),
                                function (h) c(0, unlist(lapply(C1[-1], function (x) stats::cor(df1[[y]], df1[, x],  method = h))))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      gammas1_initial_corr <- c(as.numeric(mean.corr.xy))
      gammas_start <- gammas1_initial_corr

    } else {stop("The parameter gamma_start_default is wrongly defined. It can be only equal to 'zeros', 'corr'.")}
  } else {
    if (length(gamma_start_input) != length(C1)) {
      stop("The length of 'gamma_start_input' must be equal to the number of ",
           "covariates plus the constant (", length(C1), ").Remeber that the first parameter to include is the constant.")
    }
    gammas_start <- gamma_start_input
  }

  names(gammas_start) <- names_gammas
  gammas1 <- gammas_start

  a <- if (penalty == "SCAD") a1 else a2
  prox_penalty <- get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  delta_fx_gammas <- function(gammas) {
    result <- -covYI_KS(df = df1[, !(names(df1) %in% c("ID")), drop = FALSE], z = z, y = y, C = C1[C1 %in% names(gammas)], gammas = gammas,
                       tau = tau, w = w, kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty)$gr_yi

    result <- result[names(result) %in% C1]
    return(result)
  }

  proxx_gammas <- function(x, eta) {
    #old version:
    #result <- c(x[(names(x) == "const")], prox_penalty(betas = x[!(names(x) == "const")], lambda = eta * tau, alpha = alpha, a = a)) #-1 since I dont penalize the constant
    #new version (let the const enter in the proximal operator)
    result <- prox_penalty(betas = x, lambda = eta * tau, alpha = alpha, a = a)
    return(result)
  }

  Fx <- function(aa) {
    result <- -getElement(covYI_KS(df = df1[, !(names(df1) %in% c("ID")), drop = FALSE], z = z, y = y, C = C1[C1 %in% names(aa)], gammas = aa, tau = tau,
                            w = w, kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty), paste0("covYI_KS_", penalty))
    return(result)
  }

  # Run the optimization algorithm
  if (trend == "monotone") {
    estim <- mmAPG(x0 = gammas1, c_pos = NULL, delta_fx = delta_fx_gammas, proxx = proxx_gammas, Fx = Fx,
                    lambda = tau, penalty = penalty, fold = fold, stepsizeShrink = stepsizeShrink, delta = delta,
                    max_iter = max_iter, min_alpha = min_alpha, convergence_error = convergence_error,
                    max_alpha = max_alpha, trace = trace, seed = seed, max.print = max.print, zeros_stay_zeros_from_iteration = 20)
  } else if (trend == "nonmonotone") {
    estim <- mnmAPG(x0 = gammas1, c_pos = NULL, delta_fx = delta_fx_gammas, proxx = proxx_gammas, Fx = Fx,
                     lambda = tau, penalty = penalty, fold = fold, stepsizeShrink = stepsizeShrink, delta = delta,
                     max_iter = max_iter, min_alpha = min_alpha, convergence_error = convergence_error,
                     max_alpha = max_alpha, trace = trace, seed = seed, max.print = max.print, zeros_stay_zeros_from_iteration = 20)
  }

  gammas1 <- estim$x1

  covYI_KS_value <- covYI_KS(df = df1[, !(names(df1) %in% c("ID"))], z = z, y = y, C = C1,
                           gammas = gammas1, tau = tau, w = w, kernel = kernel,
                           alpha = alpha, a1 = a1, a2 = a2, penalty = penalty, run_aauc = run_aauc)

  if (trace %in% c(1, 2)) {
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend, "version ; \n")
    if (!is.null(fold)) {cat("fold =", fold, "; ")}
    visualize_gammas <- c(gammas1[which(gammas1 != 0)])
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters, "; tau:", tau, "; weight:", w, "; penalty:", penalty, "; covYI_KS:", getElement(covYI_KS_value,  paste0("covYI_KS_", penalty)), "; youden_index:", covYI_KS_value$youden_index, "; aYI:", covYI_KS_value$aYI, "; sensitivity:", covYI_KS_value$sensitivity, "; specificity:", covYI_KS_value$specificity, "; geometric_mean:", covYI_KS_value$geometric_mean, "; fdr:", covYI_KS_value$fdr, "; mcc:", covYI_KS_value$mcc, "; auc:", covYI_KS_value$auc, "; aauc:", covYI_KS_value$aauc, "; corrclass:", covYI_KS_value$corrclass, ";\n")
    cat("TP:", covYI_KS_value$TP, "; TN:", covYI_KS_value$TN, "; FP:", covYI_KS_value$FP, "; FN:", covYI_KS_value$FN, "; gammas: \n")
    print(visualize_gammas)
  }

  #total number of variables
  n_total_var_gammas <- length(gammas1)

  #n of predicted zeros
  n_predicted_zeros_gammas <- sum(gammas1 == 0)

  #n of predicted gammas different from 0
  n_predicted_non_zeros_gammas <- sum(gammas1 != 0)

  n_caught_gammas <- NA
  n_non_caught_gammas <- NA
  n_caught_zero_gammas <- NA
  n_zero_not_caught_gammas <- NA

  #compute other measures
  if (!is.null(regressors_gammas) && length(gammas1[-1]) == length(regressors_gammas)) {

    #n of gammas caught
    n_caught_gammas <- sum(gammas1[-1][regressors_gammas != 0] != 0)

    #n of gammas not caught
    n_non_caught_gammas <- sum(gammas1[-1][regressors_gammas != 0] == 0)

    #n of zeros caught
    n_caught_zero_gammas <- sum(gammas1[-1][regressors_gammas == 0] == 0)

    #n of zeros not caught
    n_zero_not_caught_gammas <- sum(gammas1[-1][regressors_gammas == 0] != 0)
  } else if (!is.null(regressors_gammas)) {
    stop("The length of 'regressors_gammas' does not match the number of covariates.")
  }

  results <- list("t2" = gammas1,
                 "t3" = getElement(covYI_KS_value, paste0("covYI_KS_", penalty)),
                 gr_yi = covYI_KS_value$gr_yi,
                 tau = tau,
                 penalty = penalty,
                 gammas_start = gammas_start,
                 kernel = kernel,
                 c_hat = covYI_KS_value$c_hat,
                 z_hat = covYI_KS_value$z_hat,
                 y_hat = covYI_KS_value$y_hat,
                 youden_index = covYI_KS_value$youden_index,
                 sensitivity = covYI_KS_value$sensitivity,
                 specificity = covYI_KS_value$specificity,
                 geometric_mean = covYI_KS_value$geometric_mean,
                 fdr = covYI_KS_value$fdr,
                 mcc = covYI_KS_value$mcc,
                 auc = covYI_KS_value$auc,
                 aauc = covYI_KS_value$aauc,
                 aYI = covYI_KS_value$aYI,
                 corrclass = covYI_KS_value$corrclass,
                 TP = covYI_KS_value$TP,
                 TN = covYI_KS_value$TN,
                 FP = covYI_KS_value$FP,
                 FN = covYI_KS_value$FN,
                 n_gammas = length(gammas1[which(gammas1 != 0)]),
                 n_total_var_gammas = n_total_var_gammas,
                 n_predicted_zeros_gammas = n_predicted_zeros_gammas,
                 n_predicted_non_zeros_gammas = n_predicted_non_zeros_gammas,
                 n_caught_gammas = n_caught_gammas,
                 n_non_caught_gammas = n_non_caught_gammas,
                 n_caught_zero_gammas = n_caught_zero_gammas,
                 n_zero_not_caught_gammas = n_zero_not_caught_gammas)

  names(results)[1:2] <- paste0(c("gammas_hat_", "covYI_KS_"), penalty)

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
	if (trace %in% c(1, 2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  results <- c(results ,
    input_parameters = list(
      df = df, 
      z = z, 
      y = y, 
      C = C, 
      tau = tau, 
      w = w, 
      penalty = penalty, 
      gamma_start_input = gamma_start_input, 
      gamma_start_default = gamma_start_default, 
      alpha = alpha, 
      a1 = a1, 
      a2 = a2, 
      regressors_gammas = regressors_gammas, 
      fold = fold, 
      max_iter = max_iter, 
      max.print = max.print, 
      trend = trend, 
      delta = delta, 
      max_alpha = max_alpha, 
      stepsizeShrink = stepsizeShrink, 
      min_alpha = min_alpha, 
      convergence_error = convergence_error, 
      trace = trace, 
      seed = seed, 
      kernel = kernel, 
      run_aauc = run_aauc
    ),
		estimation_time = estimation_time, 
		niter = estim$tot_iters)
  class(results) <- ("covYI_KS_estimation")
  return(results)
}