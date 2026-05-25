#' @title Penalized Youden Index (pye) function via Kernel Smoothing Density
#'
#' @description The Penalized Youden Index (pye) function based on the Kernel
#' Smooth density estimator. It not only performs value of pye but it returns
#' all the necessary for the estimation process, like measure of fit and
#' derivatives. It works for all the considered penalties (L12, L1,
#' EN, SCAD and MCP)
#'
#' @param df A `data.frame` containing the input dataset. Must include the
#'   columns specified in `X` and `y`.
#' @param X A `character` vector specifying the names of the regressor variables
#'   (biomarkers) to consider for the pye evaluation. Alternatively, a
#'   `data.frame` containing only these regressors can be provided, in which case
#'   their column names will be used. Default is all columns in `df` not
#'   specified as `y`.
#' @param y A `character` string specifying the name of the target binary
#'   variable (outcome). This variable must contain only values 0 and 1.
#'   Alternatively, a single-column `data.frame` containing only the target
#'   variable can be provided, in which case its column name will be used.
#'   Default is "y".
#' @param betas A `numeric` vector representing the coefficients of the
#'   biomarker combination (Z = X \%*\% betas) used for the evaluation of pye.
#' @param lambda A `numeric` value specifying the penalization parameter for
#'   the regressors `X`.
#' @param c A `numeric` value representing the cut-off point used to classify
#'   observations based on the combined biomarker score Z. Can be a single
#'   value (applied to all observations) or a numeric vector of length `nrow(df)`
#'   for observation-specific cut-offs. Default is 0.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index. Sensitivity is weighted by `w` and specificity by
#'   `1 - w`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param kernel A `character` string indicating the kernel type to use for the
#'   estimation of the density function. Currently, only "gaussian" is tested
#'   and supported. Default is "gaussian".
#' @param alpha A `numeric` value between 0 and 1, used as the mixing parameter
#'   for the Elastic-Net penalization term. Only relevant if `penalty` is "EN".
#'   Default is 0.5.
#' @param a1 A `numeric` value specifying the regularization parameter 'a' for
#'   the SCAD penalty (as defined in the original paper by Fan and Li). Only
#'   relevant if `penalty` is "SCAD". Default is 3.7.
#' @param a2 A `numeric` value specifying the regularization parameter 'gamma'
#'   for the MCP penalty (as defined in the original paper by Zhang). Only
#'   relevant if `penalty` is "MCP". Default is 3.0.
#' @param penalty A `character` string indicating the type of penalty considered
#'   for the pye calculation. Must be one of "L12", "L1", "EN" (Elastic-Net),
#'   "SCAD", or "MCP". Default is "L1".
#' @param h_exponent A `numeric` value for the exponent of the bandwidth `h` in
#'   the Kernel Smooth density estimation.
#' @param use_opt_c A `logical` value. If `TRUE`, the function internally
#'   estimates and uses the optimal cut-off point 'c' based on the given betas.
#'   If `FALSE` (default), the `c` parameter provided directly is used.
#' @param prediction A `logical` value. If `TRUE`, the empirical maximum Youden
#'   Index is returned as the Youden index performance measure. This is useful
#'   for evaluating predictive performance on new data.
#' @param print.CDF.plot A `logical` value. If `TRUE`, a plot of the Cumulative
#'   Distribution Functions (CDFs) for cases and controls, based on the
#'   combined regressor scores Z, will be printed. Default is `FALSE`.
#'
#' @return A list with the following components:
#'   \item{pye_L12}{\code{numeric}. The Penalized Youden Index value using the
#'     \code{L1 / 2} penalty.}
#'   \item{pye_L1}{\code{numeric}. The Penalized Youden Index value using the
#'     L1 penalty.}
#'   \item{pye_EN}{\code{numeric}. The Penalized Youden Index value using the
#'     Elastic-Net penalty.}
#'   \item{pye_SCAD}{\code{numeric}. The Penalized Youden Index value using the
#'     SCAD penalty.}
#'   \item{pye_MCP}{\code{numeric}. The Penalized Youden Index value using the
#'     MCP penalty.}
#'   \item{gr_yi}{\code{numeric vector}. The gradient of the Youden Index with
#'     respect to the betas and the cut-off point 'c'. The last element
#'     corresponds to the gradient with respect to 'c'.}
#'   \item{youden_index}{\code{numeric}. The unpenalized Youden Index calculated
#'     using the Kernel Smooth density estimator (or the empirical Youden Index
#'     if \code{prediction = TRUE}). It is weighted if w is not 0.5.}
#'   \item{sensitivity}{\code{numeric}. The sensitivity (True Positive Rate) for
#'     the given 'betas' and 'c'.}
#'   \item{specificity}{\code{numeric}. The specificity (True Negative Rate) for
#'     the given 'betas' and 'c'.}
#'   \item{geometric_mean}{\code{numeric}. The geometric mean of sensitivity and
#'     specificity.}
#'   \item{fdr}{\code{numeric}. The False Discovery Rate (FP / (FP + TP)).}
#'   \item{mcc}{\code{numeric}. The Matthews Correlation Coefficient.}
#'   \item{auc}{\code{numeric}. The Area Under the ROC Curve (AUC) for the
#'     given 'z_hat' scores.}
#'   \item{corrclass}{\code{numeric}. The overall correct classification rate
#'     ((TP + TN) / N).}
#'   \item{empir_suggest_yi_c}{\code{numeric}. The optimal cut-off point 'c'
#'     suggested by the 'Youden' method from \code{OptimalCutpoints} package
#'     based on the empirical Youden Index.}
#'   \item{corrclass_with_YI_c}{\code{numeric}. The correct classification rate
#'     using the empirically suggested optimal cut-off point
#'     (\code{empir_suggest_yi_c}).}
#'   \item{yi1_yi_c}{\code{numeric}. The Youden Index obtained using the
#'     empirically suggested optimal cut-off point (\code{empir_suggest_yi_c}).}
#'   \item{c_hat}{\code{data.frame} with columns \code{ID} and \code{c_hat}. The
#'     cut-off point(s) used for each observation (either a single value repeated
#'     or input 'c').}
#'   \item{z_hat}{\code{data.frame} with columns \code{ID} and \code{z_hat}. The
#'     combined biomarker scores (\code{Z = X * betas}) for each observation.}
#'   \item{y_hat}{\code{data.frame} with columns \code{ID} and \code{y_hat}. The
#'     predicted binary outcomes (0 or 1) based on 'z_hat' and 'c_hat'.}
#'   \item{TP}{\code{numeric}. Number of True Positives.}
#'   \item{TN}{\code{numeric}. Number of True Negatives.}
#'   \item{FP}{\code{numeric}. Number of False Positives.}
#'   \item{FN}{\code{numeric}. Number of False Negatives.}
#'   \item{input_data}{\code{list}. A list containing the input parameters used
#'     for the calculation:
#'     \itemize{
#'       \item \code{beta} (\code{numeric vector}): The input 'betas' coefficients.
#'       \item \code{w} (\code{numeric vector}): The input 'w' weighting parameter.
#'       \item \code{lambda} (\code{numeric}): The input 'lambda' penalization parameter.
#'       \item \code{alpha} (\code{numeric}): The input 'alpha' parameter for Elastic-Net.
#'       \item \code{a1} (\code{numeric}): The input 'a1' parameter for SCAD.
#'       \item \code{a2} (\code{numeric}): The input 'a2' parameter for MCP.
#'       \item \code{prediction} (\code{logical}): The input 'prediction' flag.
#'       \item \code{c_hat} (\code{numeric} or \code{vector}): The cut-off point(s)
#'         that were ultimately used (could be the input 'c' or
#'         'empir_suggest_yi_c' if \code{use_opt_c = TRUE}).
#'       \item \code{kernel} (\code{character}): The kernel type used for
#'         density estimation.
#'     }
#'   }
#'
#' @examples
#' library(pye)
#' sim_data <- create_sample_with_covariates(
#'     rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' penalty <- "L12"
#' lambda <- 0.1
#' betas <- rep(1, length(X))
#' c <- 0
#'
#' pye_result <- pye_KS(df = df[, names(df) %in% c(X,y)], X = X, y = y, betas = betas,
#'   lambda = lambda, c = c, alpha = 0.5, a1 = 3.7, a2 = 3, penalty = penalty)
#' print(pye_result)
#'
#' @importFrom evmix kpz kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom plyr join
#' @importFrom stats ecdf
#' @importFrom ggplot2 ggplot geom_line aes labs theme_minimal
#' @export
pye_KS <- function(df,
                   X = NULL,
                   y = "y",
                   betas,
                   lambda,
                   c = 0,
									 w = 0.5,
                   kernel = "gaussian",
                   alpha = 0.5,
                   a1 = 3.7,
                   a2 = 3,
                   penalty = "L1",
                   h_exponent = 0.2,
                   use_opt_c = FALSE,
                   prediction = FALSE,
                   print.CDF.plot = FALSE) {

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

  if (length(c) != nrow(df)) {
    if (length(c) == 1) {
      c <- rep(c, nrow(df))
    } else {stop("c can only be of class numeric on length 1 or equal to the number rows of df")}
  }

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  # df1 <- cbind(ID, df[, (names(df) %in% c(y, X)), drop = FALSE]) #OLD
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}
	if (!is.numeric(w) || w < 0 || w > 1) {stop("The weight 'w' must be a numeric value between 0 and 1.")}

  if (length(lambda) > 1) {stop("pye_KS stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12", "L1", "EN", "SCAD", "MCP"))) {stop("A wrong value has been assigned to the parameter penalty.")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))) {
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(betas) != length(X)) {
    stop("The number of element of betas is different then the number of columns in X")
  }

  #divide control and diseased and multiply by betas
  z_y0 <- as.matrix(df1[df1[[y]] == 0, X, drop = FALSE]) %*% betas #drop = FALSE permits it to remain a data frame even if it is composed by a single column
  z_y1 <- as.matrix(df1[df1[[y]] == 1, X, drop = FALSE]) %*% betas

  #divide c
  c_y0 <- as.matrix(c[df1[[y]] == 0])
  c_y1 <- as.matrix(c[df1[[y]] == 1])

  #append data
  z <- as.data.frame(rbind(z_y0, z_y1))
  z['ID'] <- as.numeric(rownames(z))
  z <- z[order(z$ID), ]
  names(z)[names(z) == 'V1'] <- 'z_hat'
  df2 <- plyr::join(df1, z, by = 'ID')
  #print(cbind(df2["y"], df2["z_hat"]))
  rownames(df2) <- df2$ID

  #find the optimal cut-point
  opt <- OptimalCutpoints::optimal.cutpoints(data = df2, X = "z_hat", status = y, methods = "Youden", tag.healthy = 0, control = OptimalCutpoints::control.cutpoints(generalized.Youden = TRUE,  CFN = max(w, 1e-6), CFP = max(1 - w, 1e-6)))
  empir_suggest_yi_c <- mean(opt$Youden$Global$optimal.cutoff$cutoff)

  #try the Corr Class with the suggested c
  TP_yi_c <- sum(ifelse(z_y1 >= empir_suggest_yi_c, 1, 0))
  TN_yi_c <- sum(ifelse(z_y0 < empir_suggest_yi_c, 1, 0))
  FP_yi_c <- sum(ifelse(z_y0 >= empir_suggest_yi_c, 1, 0))
  FN_yi_c <- sum(ifelse(z_y1 < empir_suggest_yi_c, 1, 0))

  #compute the Youden index
  spec_yi_c <- TN_yi_c / (TN_yi_c + FP_yi_c)
  fnr_yi_c <- FN_yi_c / (FN_yi_c + TP_yi_c)
  #yi1_yi_c <- spec_yi_c - fnr_yi_c
	yi1_yi_c <- 2 * (1 - w) * spec_yi_c + 2 * w * (1 - fnr_yi_c) - 1
  corrclass_with_YI_c <- (TP_yi_c + TN_yi_c) / nrow(df1)

  #if we want to use the empirical c
  if (use_opt_c == TRUE) {
    c <- empir_suggest_yi_c
  }

  #kernel from the YI paper
  #healthy
  # if (kernel %in% c("gaussian", "epanechnikov", "triweight", "tricube", "biweight", "cosine") & sum(z_y0) != 0) {
  #   h0 = kedd::h.bcv(as.numeric(c-z_y0), kernel = kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h0 <- 0.9 * min(stats::sd(z_y0), stats::IQR(z_y0) / 1.34) * length(z_y0)^(-h_exponent)
  # }
  if (h0 < 0.1) {h0 <- 0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t0 <- as.numeric((c_y0 - z_y0) / h0)
  cum0 <- evmix::kpz(z = t0, kernel = kernel)
  f0 <- sum(cum0) / length(z_y0)
  #f0 <- sum(pnorm(t0, 0, 1)) / length(z_y0)

  #diseased
  # if (kernel %in% c("gaussian", "epanechnikov", "triweight", "tricube", "biweight", "cosine") & sum(z_y1) != 0) {
  #   h1 = kedd::h.bcv(as.numeric(z_y1), kernel = kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h1 <- 0.9 * min(stats::sd(z_y1), stats::IQR(z_y1) / 1.34) * length(z_y1)^(-h_exponent)
  # }
  if (h1 < 0.1) {h1 <- 0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t1 <- as.numeric((c_y1 - z_y1) / h1)
  cum1 <- evmix::kpz(z = t1, kernel = kernel)
  f1 <- sum(cum1) / length(z_y1)
  #f1 <- sum(pnorm(t1, 0, 1)) / length(z_y1)

  #yi <- f0 - f1
	yi <- 2 * (1 - w) * f0 - 2 * w * f1 + 2 * w - 1 #formula from Wang et. al. 2025

  if (print.CDF.plot == TRUE) {

    # Compute the empirical distribution functions using the Gaussian kernel
    z_y0_ord <- z_y0[order(z_y0)]
    z_y1_ord <- z_y1[order(z_y1)]
    ecdf0 <- stats::ecdf (z_y0_ord)
    ecdf1 <- stats::ecdf (z_y1_ord)

    # Create a data frame for the plot
    plot_data <- data.frame(x_axis = c(z_y0_ord, z_y1_ord),
                            y_axis = c(ecdf0(z_y0_ord), ecdf1(z_y1_ord)),
                            groups_to_consider = c(rep("CDF_y0", length(z_y0_ord)), rep("CDF_y1", length(z_y1_ord))))

    # Crea il plot
    x_axis <- y_axis <- groups_to_consider <- NULL
    ggplot2::ggplot(plot_data, ggplot2::aes(x = x_axis, y = y_axis, color = groups_to_consider)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Z", y = "CDF of Z for case and conntrol patients", color = "Groups") +
      ggplot2::theme_minimal()
  }

  #evaluate the gradient
  #h0 <- 0.9 * min(stats::sd(z_y0), stats::IQR(z_y0) / 1.34) * (length(z_y0)^(-h_exponent))
  #t0 <- (c_y0 - z_y0) / h0
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")) {
    gr0_betas <- apply (df1[df1[[y]] == 0, X, drop = FALSE], 2,
                       function(x) (t(stats::dnorm(t0, 0, 1)) %*% ((-x) / h0)) / length(z_y0))
    gr0_c <- sum(stats::dnorm(t0, 0, 1)) / (length(z_y0) * h0)
  } else {
    gr0_betas <- apply (df1[df1[[y]] == 0, X, drop = FALSE], 2,
           function(x) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% ((-x) / h0)) / length(z_y0))
    gr0_c <- sum(evmix::kdz(z = as.numeric(t0), kernel = kernel)) / (length(z_y0) * h0)
  }
  names(gr0_c) <- "c"

  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")) {
    gr1_betas <- apply (df1[df1[[y]] == 1, X, drop = FALSE], 2,
                        function(x) (t(stats::dnorm(t1, 0, 1)) %*% ((-x) / h1)) / length(z_y1))
    gr1_c <- sum(stats::dnorm(t1, 0, 1)) / (length(z_y1) * h1)
  } else {
    gr1_betas <- apply (df1[df1[[y]] == 1, X, drop = FALSE], 2,
           function(x) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% ((-x) / h1)) / length(z_y1))
    gr1_c <- sum(evmix::kdz(z = as.numeric(t1), kernel = kernel)) / (length(z_y1) * h1)
  }
  names(gr1_c) <- "c"

  #gradient
  #gr_yi <- c(gr0_betas, gr0_c) - c(gr1_betas, gr1_c)
	gr_yi <- 2 * (1 - w) * c(gr0_betas, gr0_c) - 2 * w * c(gr1_betas, gr1_c)

  #append data
  z <- as.data.frame(rbind(z_y0, z_y1))
  z['ID'] <- as.numeric(rownames(z))
  #z <- z[order(z$ID), ]
  names(z)[names(z) == 'V1'] <- 'z_hat'
  df1 <- plyr::join(df1, z, by = 'ID')
  # df1<-merge(df1, z, by = 'ID', all = TRUE)
  rownames(df1) <- df1$ID

  #assign the ID to c
  c2 <- cbind(ID = df1["ID"], c_hat = c)
  rownames(c2) <- c2$ID
  #append data (c)
  df1 <- plyr::join(df1, c2, by = 'ID')
  # df1<-merge(df1, c, by = 'ID', all = TRUE)
  rownames(df1) <- df1$ID

  #the following few lines would be the way to estimate c separately with respect of
  #the betas, but we try to estimate all at the same time
  # if (mode == "estimation") {
  #   c_hat <- if (length(c$Youden$Global$optimal.cutoff$cutoff) == 1) {c$Youden$Global$optimal.cutoff$cutoff} else {c$Youden$Global$optimal.cutoff$cutoff[round(length(c$Youden$Global$optimal.cutoff$cutoff) / 2)]}
  # } else if (mode != "prediction") {stop("The parameter mode is not valid!")}
  #

  df1["y_hat"] <- ifelse(df1$z_hat >= df1$c_hat, 1, 0)

  # Calculate performance measures
  # AUC e YI
  auc <-  mean(opt$Youden$Global$measures.acc$AUC)

	# For the YI - Extract sensitivity and specificity from the opt object
  se_emp <- mean(opt$Youden$Global$measures.acc$Se)
  sp_emp <- mean(opt$Youden$Global$measures.acc$Sp)
	#yi1 <- mean(opt$Youden$Global$optimal.criterion)
  yi1 <- 2 * w * se_emp + 2 * (1 - w) * sp_emp - 1


  # Confusion matrix components
  TP <- sum(ifelse(z_y1 >= c_y1, 1, 0))
  TN <- sum(ifelse(z_y0 < c_y0, 1, 0))
  FP <- sum(ifelse(z_y0 >= c_y0, 1, 0))
  FN <- sum(ifelse(z_y1 < c_y1, 1, 0))

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
  corrclass <- (TP + TN) / nrow(df1)

  #if all the betas are zero, the measures are zeros
  if (sum(betas) == 0) {
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    spec <- 0
    yi <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  #L_(1 / 2) penalization
  phi_L12 <- lambda * sum(abs(betas)^(1 / 2))#NB:we DO NOT add the intercept!!!
  #gr_phi_L12 <- lambda*betas/(2 * abs(betas)^(3/2))
  #gr_phi_L12[is.na(gr_phi_L12)] <- 0
  if (is.na(phi_L12)) {
    stop('problem with the penalization phi_L12 in "pye" function')}

  #L1 penalization
  phi_L1 <- lambda * sum(abs(betas))#NB:we DO NOT add the intercept!!!
  #gr_phi_L1 <- lambda * sign(betas)
  #gr_phi_L1[is.na(gr_phi_L1)] <- 0
  if (is.na(phi_L1)) {
    stop('problem with the penalization phi_L1 in "pye" function')}

  #Elastic-Net penalization

  phi_EN <- lambda * ((alpha) * (sum(abs(betas))) + ((1 - alpha) / 2) * (sum(betas^2)))
  #gr_phi_EN <- lambda * ((alpha)*sign(betas) + (1 - alpha)*betas)
  #gr_phi_EN[is.na(gr_phi_EN)] <- 0
  if (is.na(phi_EN)) {
    stop('problem with the penalization phi_EN in "pye" function')}

  #SCAD penalization
  phi_SCAD <- SCAD_function(betas, lambda, a = a1)
  #gr_phi_SCAD <- SCAD_derivative(betas, lambda, a = a1)
  #gr_phi_SCAD[is.na(gr_phi_SCAD)] <- 0
  if (is.na(phi_SCAD)) {
    stop('problem with the penalization phi_SCAD in "pye" function')}

  #MCP penalization
  phi_MCP <- MCP_function(betas, lambda, a = a2)
  #gr_phi_MCP <- MCP_derivative(betas, lambda,a = a2)
  #gr_phi_MCP[is.na(gr_phi_MCP)] <- 0
  if (is.na(phi_MCP)) {
    stop('problem with the penalization phi_MCP in "pye" function')}

  #join yi and the penalty functions
  pye_L12 <- yi - phi_L12
  pye_L1 <- yi - phi_L1
  pye_EN <- yi - phi_EN
  pye_SCAD <- yi - phi_SCAD
  pye_MCP <- yi - phi_MCP

  #gr_pye_L12 <- c(gr_yi[-length(gr_yi)] - gr_phi_L12, gr_yi[length(gr_yi)])
  #gr_pye_L1 <- c(gr_yi[-length(gr_yi)] - gr_phi_L1, gr_yi[length(gr_yi)])
  #gr_pye_EN <- c(gr_yi[-length(gr_yi)] - gr_phi_EN, gr_yi[length(gr_yi)])
  #gr_pye_SCAD <- c(gr_yi[-length(gr_yi)] - gr_phi_SCAD, gr_yi[length(gr_yi)])
  #gr_pye_MCP <- c(gr_yi[-length(gr_yi)] - gr_phi_MCP, gr_yi[length(gr_yi)])

  #if we are in prediction this is the right yi
  if (prediction == TRUE) {
    yi <- yi1
  }

  return(list(pye_L12 = pye_L12,
              pye_L1 = pye_L1,
              pye_EN = pye_EN,
              pye_SCAD = pye_SCAD,
              pye_MCP = pye_MCP,
              # gr_pye_L12 = gr_pye_L12, gr_pye_L1 = gr_pye_L1, gr_pye_EN = gr_pye_EN, gr_pye_SCAD = gr_pye_SCAD, gr_pye_MCP = gr_pye_MCP,
              # gr_phi_L12 = gr_phi_L12, gr_phi_L1 = gr_phi_L1, gr_phi_EN = gr_phi_EN, gr_phi_SCAD = gr_phi_SCAD, gr_phi_MCP = gr_phi_MCP,
              gr_yi = gr_yi,
              youden_index = yi,
              sensitivity = sensitivity,
              specificity = spec,
              geometric_mean = gm,
              fdr = fdr,
              mcc = mcc,
              auc = auc,
              corrclass = corrclass,
              empir_suggest_yi_c = empir_suggest_yi_c,
              corrclass_with_YI_c = corrclass_with_YI_c,
              yi1_yi_c = yi1_yi_c,
              c_hat = c,
              z_hat = df1[, c("ID", "z_hat")],
              y_hat = df1[, c("ID", "y_hat")],
              TP = TP,
              TN = TN,
              FP = FP,
              FN = FN,
              input_data = list(beta = betas, lambda = lambda, w = w, alpha = alpha,
                                a1 = a1, a2 = a2, prediction = prediction, c_hat = c, kernel = kernel)
              ))
}





#' @title pye KS Estimation for Coefficient and Feature Selection
#'
#' @description function to estimate the optimal value of betas and c maximizing
#' the pye function. To find the optimum the mmAPG and mnmAPG algorithms are
#' used.
#'
#' @param df A `data.frame` containing the input dataset. Must include the
#'   columns specified in `X` and `y`.
#' @param X A `character` vector specifying the names of the regressor variables
#'   (biomarkers) to consider in the estimation. Alternatively, a `data.frame`
#'   containing only these regressors can be provided, in which case their
#'   column names will be used. Default is all columns in `df` not specified as
#'   `y`.
#' @param y A `character` string specifying the name of the target binary
#'   variable (outcome). This variable must contain only values 0 and 1.
#'   Alternatively, a single-column `data.frame` containing only the target
#'   variable can be provided, in which case its column name will be used.
#'   Default is "y".
#' @param lambda A `numeric` value specifying the penalization parameter for
#'   the regressors `X`. Must be a single non-negative value.
#' @param penalty A `character` string indicating the type of penalty to apply.
#'   Must be one of "L12", "L1", "EN" (Elastic-Net), "SCAD", or "MCP". Default
#'   is "L1".
#' @param w A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index. Sensitivity is weighted by `w` and specificity by
#'   `1 - w`. Default is 0.5.
#' @param beta_start_input A `numeric` vector of specific starting points for
#'   the beta coefficients. If `NULL` (default), the starting points are
#'   determined by `beta_start_default`. If provided, its length must match the
#'   number of regressors in `X`.
#' @param beta_start_default A `character` string defining the default starting
#'   points for betas if `beta_start_input` is `NULL`. If "zeros" (default),
#'   betas start with a vector of all zeros. If "corr", betas start with the
#'   absolute value of the correlation of each regressor with the target
#'   variable `y`.
#' @param max.print The number of elements to show when printing results.
#'   Default is 10.
#' @param alpha A `numeric` value between 0 and 1, used as the mixing parameter
#'   for the Elastic-Net penalization term. Only relevant if `penalty` is "EN".
#'   Default is 0.5.
#' @param a1 A `numeric` value specifying the regularization parameter 'a' for
#'   the SCAD penalty (as defined in the original paper by Fan and Li). Only
#'   relevant if `penalty` is "SCAD". Default is 3.7.
#' @param a2 A `numeric` value specifying the regularization parameter 'gamma'
#'   for the MCP penalty (as defined in the original paper by Zhang). Only
#'   relevant if `penalty` is "MCP". Default is 3.0.
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param fold numeric. An optional fold number, used when the function is
#'   called within a cross-validation loop. This is primarily for tracking
#'   and reporting purposes. Default is `NULL`.
#' @param max_iter A `numeric` value specifying the maximum number of iterations
#'   for the optimization algorithms (mmAPG and mnmAPG). Default is 10000.
#' @param trend A `character` string indicating the optimization algorithm
#'   variant to use. If "monotone" (default), the Monotone Accelerated Proximal
#'   Gradient (mmAPG) algorithm is used. If "nonmonotone", the Non-Monotone
#'   Accelerated Proximal Gradient (mnmAPG) algorithm is used.
#' @param delta A `numeric` value specifying the parameter for the convergence
#'   condition of the optimization algorithm. Default is 1e-5.
#' @param max_alpha A `numeric` value representing the maximum allowed value for
#'   the step-size parameter `alpha` in the backtracking line-search. Default
#'   is 10000.
#' @param stepsizeShrink A `numeric` value between 0 and 1, used to adjust the
#'   step-size in the backtracking line-search within the optimization of pye.
#'   Values closer to 1 lead to higher accuracy but longer estimation times.
#'   Default is 0.8.
#' @param min_alpha A `numeric` value representing the minimum allowed value for
#'   the step-size parameter `alpha`. Default is 1e-10.
#' @param convergence_error A `numeric` value specifying the error tolerance for
#'   considering the optimization algorithm converged. Default is 1e-7.
#' @param trace A `numeric` value controlling the verbosity of the output. 0: no
#'   visualization, 1: visualize only the final result (default), 2: visualize
#'   all steps during optimization.
#' @param seed A `numeric` value to fix the random seed for reproducibility.
#'   Default is 1.
#' @param kernel A `character` string indicating the kernel type to use for the
#'   estimation of the density function. Currently, only "gaussian" is tested
#'   and supported. Default is "gaussian".
#' @param c_zero_fixed A `logical` value. If `TRUE`, the estimation process
#'   considers the cut-off point `c` as fixed and equal to zero, which can
#'   reduce estimation complexity. If `FALSE` (default), `c` can vary and is
#'   estimated by the pye optimization.
#' @param zeros_stay_zeros_from_iteration A `numeric` value specifying the
#'   iteration number from which parameters that have reached zero in the
#'   estimation cannot change (i.e., remain zero) anymore. This helps preserve
#'   the sparsity of the solution. Default is 20.
# #' @param long_suffix A `character` string. If provided, it's used as a prefix
# #'   to identify longitudinal variables (e.g., "visit_"). For proper management
# #'   of longitudinal variables or coefficients referring to an original
# #'   longitudinal dimension, `long_suffix` is expected to be followed by a
# #'   number (e.g., "visit_1", "visit_2"). If `NULL` (default), no special
# #'   handling for longitudinal variables is applied, presuming the input `df`
# #'   does not contain such variables.
#'
#' @return A list with the following components:
#'   \item{betas_hat_penalty}{\code{numeric vector}. The estimated optimal
#'     coefficients (betas) for the regressors X, with \code{penalty} replaced
#'     by the actual penalty used (e.g., \code{betas_hat_L1}). Non-zero values
#'     indicate selected features.}
#'   \item{pye_KS_penalty}{\code{numeric}. The optimized Penalized Youden Index
#'     (pye) value achieved for the given \code{penalty} (e.g.,
#'     \code{pye_KS_L1}).}
#'   \item{gr_yi}{\code{numeric vector}. The gradient of the Youden Index at the
#'     estimated optimal point.}
#'   \item{lambda}{\code{numeric}. The penalization parameter used in the
#'     estimation.}
#'   \item{penalty}{\code{character}. The type of penalty used for
#'     regularization (e.g., "L1", "SCAD").}
#'   \item{betas_start}{\code{numeric vector}. The initial starting point for
#'     the betas used in the optimization algorithm.}
#'   \item{kernel}{\code{character}. The kernel type used for density estimation
#'     (e.g., "gaussian").}
#'   \item{c_hat}{\code{numeric}. The estimated optimal cut-off point 'c'.}
#'   \item{z_hat}{\code{data.frame} with columns \code{ID} and \code{z_hat}. The
#'     combined biomarker scores (\code{Z = X * betas}) for each observation
#'     based on the estimated \code{betas_hat}.}
#'   \item{y_hat}{\code{data.frame} with columns \code{ID} and \code{y_hat}. The
#'     predicted binary outcomes (0 or 1) based on \code{z_hat} and
#'     \code{c_hat}.}
#'   \item{youden_index}{\code{numeric}. The unpenalized Youden Index calculated
#'     using the Kernel Smooth density estimator at the estimated optimal point.}
#'   \item{sensitivity}{\code{numeric}. The sensitivity (True Positive Rate) at
#'     the estimated optimal point.}
#'   \item{specificity}{\code{numeric}. The specificity (True Negative Rate) at
#'     the estimated optimal point.}
#'   \item{geometric_mean}{\code{numeric}. The geometric mean of sensitivity and
#'     specificity at the estimated optimal point.}
#'   \item{fdr}{\code{numeric}. The False Discovery Rate (FP / (FP + TP)) at the
#'     estimated optimal point.}
#'   \item{mcc}{\code{numeric}. The Matthews Correlation Coefficient at the
#'     estimated optimal point.}
#'   \item{auc}{\code{numeric}. The Area Under the ROC Curve (AUC) for the
#'     \code{z_hat} scores at the estimated optimal point.}
#'   \item{corrclass}{\code{numeric}. The overall correct classification rate
#'     ((TP + TN) / N) at the estimated optimal point.}
#'   \item{n_betas}{\code{numeric}. The number of non-zero (selected) betas in
#'     the final estimation.}
#'   \item{n_total_var}{\code{numeric}. The total number of regressors
#'     considered (length of X).}
#'   \item{n_predicted_zeros}{\code{numeric}. The number of estimated betas that
#'     are exactly zero.}
#'   \item{n_predicted_non_zeros}{\code{numeric}. The number of estimated betas
#'     that are non-zero.}
#'   \item{n_caught_betas}{\code{numeric} or \code{NA}. If \code{regressors_betas}
#'     is provided, the number of true non-zero betas correctly estimated as
#'     non-zero. Otherwise \code{NA}.}
#'   \item{n_non_caught_betas}{\code{numeric} or \code{NA}. If
#'     \code{regressors_betas} is provided, the number of true non-zero betas
#'     incorrectly estimated as zero. Otherwise \code{NA}.}
#'   \item{n_caught_zero}{\code{numeric} or \code{NA}. If \code{regressors_betas}
#'     is provided, the number of true zero betas correctly estimated as zero.
#'     Otherwise \code{NA}.}
#'   \item{n_zero_not_caught}{\code{numeric} or \code{NA}. If
#'     \code{regressors_betas} is provided, the number of true zero betas
#'     incorrectly estimated as non-zero. Otherwise \code{NA}.}
#'   \item{input_parameters}{\code{character vector}. A list containing the
#'     input parameters.}
#'   \item{estimation_time}{\code{difftime}. The total time taken for the
#'     estimation process, in minutes.}
#'   \item{niter}{\code{numeric}. The total number of iterations performed by
#'     the optimization algorithm (mmAPG or mnmAPG).}
#'
#' @examples
#' # Load the package
#' library(pye)
#'
#' # 1. Simulate data for the example
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' regressors_betas <- sim_data$nregressors # True betas for evaluation
#'
#' # 2. Set estimation parameters
#' penalty_type <- "SCAD"
#' lambda_val <- 0.1
#' trend_algo <- "monotone" # or "nonmonotone"
#' start_betas <- "zeros"
#' alpha_val <- 0.5
#' c_fixed <- FALSE
#'
#' # 3. Run the pye estimation
#' pye_estimation_result <- pye_KS_estimation(
#'   df = df,
#'   X = X,
#'   y = y,
#'   penalty = penalty_type,
#'   trend = trend_algo,
#'   trace = 2, # Show full trace for example
#'   beta_start_default = start_betas,
#'   beta_start_input = NULL, # No specific starting point
#'   lambda = lambda_val,
#'   alpha = alpha_val,
#'   a1 = 3.7, # Default for SCAD
#'   a2 = 3.0, # Default for MCP (not used for SCAD here, but good to include if relevant)
#'   regressors_betas = regressors_betas,
#'   c_zero_fixed = c_fixed,
#'   max_iter = 5 # Keep iterations low for quick example run
#' )
#'
#' # 4. Print results
#' print(pye_estimation_result)
#'
#'
#' @importFrom stats setNames
#' @export
#Estimation of the parameter using pye_KS
pye_KS_estimation <- function(df, X = NULL, y = "y",
                               lambda, penalty = "L1", w = 0.5, beta_start_input = NULL,
                               beta_start_default = "zeros", max.print = 10,
                               alpha = 0.5, a1 = 3.7, a2 = 3, regressors_betas = NULL, fold = NULL, max_iter = 10000,
                               trend = "monotone", delta = 1e-5, max_alpha = 10000, stepsizeShrink = 0.8,
                               min_alpha = 1e-10, convergence_error = 1e-7,
                               trace = 1, seed = 1, c_zero_fixed = FALSE, kernel = "gaussian",
                               zeros_stay_zeros_from_iteration = 20) {

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

  # Handle X and y parameters
  if (is.null(X)) X <- setdiff(names(df), y)
  if (inherits(X, "data.frame")) X <- names(X)
  if (!is.character(X) || length(X) == 0) stop("'X' must be a character vector of column names or a data.frame.")
  if (!all(X %in% names(df))) stop("Not all regressors in 'X' are found in 'df'.")

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  ID <- rownames(df)
  #df1 <- cbind(ID, df[, (names(df) %in% c(y, X)), drop = FALSE]) #OLD
  df1 <- cbind(ID, df[, c(y, X), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Further input validation
  if (length(lambda) != 1 || !is.numeric(lambda) || lambda < 0) {stop("Parameter 'lambda' must be a single non-negative numeric value.")}
  valid_penalties <- c("L12", "L1", "EN", "SCAD", "MCP")
  if (!(penalty %in% valid_penalties)) {stop(paste0("A wrong value has been assigned to the parameter 'penalty'. Must be one of: ", paste(valid_penalties, collapse = ", "), "."))}
	if (!is.numeric(w) || length(w) != 1 || w < 0 || w > 1) {stop("Parameter 'w' must be a single numeric value between 0 and 1.")}
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 2 || !is.integer(as.integer(max_iter))) {stop("Parameter 'max_iter' needs to be an integer and at least 2.")}
  if (!(trace %in% c(0, 1, 2))) {stop("The parameter 'trace' has been wrongly assigned. It can be 0 (no print), 1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))) {stop("The parameter 'trend' has been wrongly assigned. It can be 'monotone' or 'nonmonotone'.")}
  if (!is.logical(c_zero_fixed)) {stop("Parameter 'c_zero_fixed' must be a logical (TRUE/FALSE).")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel %in% valid_kernels)) {stop(paste0("Parameter 'kernel' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}

  betas1_initial_zeros <- rep(0, length(X))
  betas1_initial_corr <- NULL

  #define c
  c <- 0
  names(c) <- "c"

  #initializing the parameters (c always starts from zero)
  if (length(beta_start_input) == 0) {
    if (beta_start_default == "zeros") {

      betas_start <- betas1_initial_zeros
      names(betas_start) <- X
      betas1 <- betas_start

    } else if (beta_start_default == "corr") {

      #compute the corr between every x and y
      corr.xy <- data.frame(matrix(ncol = length(X), nrow = 0))
      names(corr.xy) <- X
      corr.xy[1:3, ] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(lapply(X, function (x) stats::cor(df1[[y]], df1[, x],  method = h)))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      betas1_initial_corr <- as.numeric(mean.corr.xy)
      betas_start <- betas1_initial_corr
      names(betas_start) <- X
      betas1 <- betas_start

    } else {stop("The parameter beta_start_default can only be equal to 'zeros', 'corr'.")}
  } else {
    if (length(beta_start_input) == length(X)) {
      betas_start <- beta_start_input
      names(betas_start) <- X
      betas1 <- betas_start
    } else if (length(beta_start_input) != 0) {
      stop("The length of the parameter beta_start_input is not equal to ", length(X) , "\n")
    }
  }
  #add c to the betas
  betas1 <- c(betas1[!(names(betas1) %in% names(c))], c) #merge betas1 and c

  a <- if (penalty == "SCAD") a1 else a2
  prox_penalty <- get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  #wrappers
  delta_fx <- function(x) {
    if (c_zero_fixed == TRUE) {
      x[names(x) == "c"] <- 0
    }		
    result <- -pye_KS(df = df1[, names(df1) != "ID", drop = FALSE], X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")], lambda = lambda,
                      c = x[(names(x) == "c")], kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty, w = w)$gr_yi
    if (c_zero_fixed == TRUE) {
      result[names(result) == "c"] <- 0
    }
    return(result)
  }

  proxx <- function(x, eta) {
    if (c_zero_fixed == TRUE) {
      x[names(x) == "c"] <- 0
    }		
    result <- c(prox_penalty(betas = x[!(names(x) == "c")], lambda = eta * lambda, alpha = alpha, a = a), x[(names(x) == "c")])
    #if (length(long_suffix) == 0) {
    #  result <- c(prox_penalty(betas = x[!(names(x) == "c")], lambda = eta * lambda, alpha = alpha, a = a), x[(names(x) == "c")])
    #} else {
    #  #we have to apply the penalty to all the betas if the variables referring to the same longitudinal biomarkers
    #  #to allow the "all or nothing" criteria to improve interpretability of the result
    #  #identify which variables are from longitudinal data
    #  longitudinal_vars <- unique(sub(paste0(long_suffix, ".*"), "", grep(paste0(long_suffix, ".*"), names(x), value = TRUE)))
    #  vars_dedup <- unique(sub(paste0(long_suffix, ".*"), "", names(x)))
    #  dgrfree <- as.numeric(sub(paste0(".*", long_suffix, "(\\d+)"), "\\1", grep(paste0(long_suffix, ".*"), names(x), value = TRUE)))
    #  dgrfree <- max(dgrfree[!is.na(dgrfree)])
    #
    #  result <- matrix(NA, nrow=length(vars_dedup), ncol = dgrfree+1)
    #  rownames(result) <- c(vars_dedup)
    #  colnames(result) <- paste0("t",0:dgrfree)
    #
    #  ordered_x <- x[order(names(x))]
    #  # Fill the result matrix
    #  for(row in 1:length(rownames(result))) {
    #    if (length(grep(paste0("^",rownames(result)[row], long_suffix, ".*"), names(ordered_x)))>0) {
    #      vec_index <- grep(paste0("^",rownames(result)[row], long_suffix, ".*"), names(ordered_x))
    #    } else {
    #      vec_index <- grep(paste0("^",rownames(result)[row], "$"), names(ordered_x))
    #    }
    #    result[row, 1:length(vec_index)] <- ordered_x[vec_index]
    #  }
    #
    #  x2 <- rowMeans(abs(result), na.rm = TRUE)
    #  temp <- c(prox_penalty(betas = x2[-length(x2)], lambda = eta * lambda, alpha = alpha, a = a), x2[length(x2)])
    #  #result <- c(prox_penalty(betas = x[-c_pos], lambda = eta * lambda, alpha = alpha, a = a), x[c_pos])
    #  #result[rownames(result) %in% names(temp[temp==0]), ] <- result[rownames(result) %in% names(temp[temp==0]), ]*0
    #  result2 <- (temp/x2)*result
    #  result2[is.nan(result2)] <- 0
    #  vec <- as.vector(result2)
    #  names(vec) <- rep(rownames(result2), times = ncol(result2))
    #  long_names <- expand.grid(unique(names(vec[names(vec) %in% longitudinal_vars])), paste0(long_suffix, 0:dgrfree))
    #  names(vec)[names(vec) %in% longitudinal_vars] <- paste(long_names$Var1, long_names$Var2, sep = "")
    #  vec <- vec[!is.na(vec)]
    #  vec <- vec[match(names(x), names(vec))] #reorder as of x
    #  result <- vec
    #}
    return(result)
  }

  Fx <- function(x) {
    if (c_zero_fixed == TRUE) {
      x[(names(x) == "c")] <- 0
    }		
    result <- -getElement(pye_KS(df = df1[, names(df1) != "ID", drop = FALSE], X = X[X %in% names(x)], y = y, betas = x[!(names(x) == "c")],
                                 lambda = lambda, c = x[(names(x) == "c")], kernel = kernel,
                                 alpha = alpha, a1 = a1, a2 = a2, penalty = penalty, w = w), paste0("pye_", penalty))
    return(result)
  }


  # Run the optimization algorithm
  if (trend == "monotone") {
    estim <- mmAPG(x0 = betas1, c_pos = length(betas1), delta_fx = delta_fx, proxx = proxx, Fx = Fx, lambda = lambda, penalty = penalty,
                   fold = fold, stepsizeShrink = stepsizeShrink, delta = delta, max_alpha = max_alpha,
                   max_iter = max_iter, min_alpha = min_alpha, convergence_error = convergence_error,
                   trace = trace, seed = seed, max.print = max.print, zeros_stay_zeros_from_iteration = zeros_stay_zeros_from_iteration)
  } else if (trend == "nonmonotone") {
    estim <- mnmAPG(x0 = betas1, c_pos = length(betas1), delta_fx = delta_fx, proxx = proxx, Fx = Fx, lambda = lambda, penalty = penalty,
                    fold = fold, stepsizeShrink = stepsizeShrink, delta = delta, max_alpha = max_alpha,
                    max_iter = max_iter, min_alpha = min_alpha, convergence_error = convergence_error,
                    trace = trace, seed = seed, max.print = max.print, zeros_stay_zeros_from_iteration = zeros_stay_zeros_from_iteration)
  }

  #divide betas1 and c
  betas_hat <- estim$x1[!(names(estim$x1) == "c")]
  c_hat <- estim$x1[(names(estim$x1) == "c")]


  #lambda == 1234 is a secret code to test the performance of the correlation alone without optimization
  #if (lambda == 1234) {
  #  #compute the corr between every x and y
  #  corr.xy <- data.frame(matrix(ncol = length(X), nrow = 0))
  #  names(corr.xy) <- X
  #  corr.xy[1:3, ] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(lapply(X, function (x) stats::cor(df1[[y]], df1[, x],  method = h)))))
  #  row.names(corr.xy) <- c("pearson", "kendall", "spearman")
  #  #mean of the corrs and sort names
  #  mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
  #  betas1_initial_corr <- as.numeric(mean.corr.xy)
  #  betas_start <- betas1_initial_corr
  #  names(betas_start) <- X
  #
  #  betas_hat <- betas_start
  #  c_hat <- c
  #}

  #compute z_hat
  pye_KS_value <- pye_KS(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, betas = betas_hat, lambda = lambda,
                         c = c_hat, kernel = kernel, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty, w = w)
  pye_value_for_current_penalty <- getElement(pye_KS_value, paste0("pye_", penalty))

  #z_hat <- pye_KS_value$z_hat$z_hat
  #df_hat <- cbind(df1, z_hat)

  youden_index <- pye_KS_value$youden_index
  sensitivity <- pye_KS_value$sensitivity
  specificity <- pye_KS_value$specificity
  geometric_mean <- pye_KS_value$geometric_mean
  fdr <- pye_KS_value$fdr
  mcc <- pye_KS_value$mcc
  auc <- pye_KS_value$auc
  corrclass <- pye_KS_value$corrclass

  if (trace %in% c(1, 2)) {
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend, "version ; \n")
    if (!is.null(fold)) {cat("fold =", fold, "; \n")}
    visualize_betas <- c(betas_hat[which(betas_hat != 0)], c_hat)
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters , "; lambda:", lambda, "; weight:", w, "; penalty:", penalty, "; pye_KS:", pye_value_for_current_penalty, "; youden_index:", pye_KS_value$youden_index, "; sensitivity:", pye_KS_value$sensitivity, "; specificity:", pye_KS_value$specificity, "; geometric_mean:", pye_KS_value$geometric_mean, "; fdr:", pye_KS_value$fdr, "; mcc:", pye_KS_value$mcc, "; auc:", pye_KS_value$auc, "; corrclass:", pye_KS_value$corrclass, "; \n")
    cat("TP:", pye_KS_value$TP, "; TN:", pye_KS_value$TN, "; FP:", pye_KS_value$FP, "; FN:", pye_KS_value$FN, "; betas: \n")
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


  results <- list("t1" = betas_hat,
                 "t2" = pye_value_for_current_penalty,
                 gr_yi = pye_KS_value$gr_yi,
                 lambda = lambda,
                 penalty = penalty,
                 betas_start = betas_start,
                 kernel = kernel,
                 c_hat = c_hat,
                 z_hat = pye_KS_value$z_hat,
                 y_hat = pye_KS_value$y_hat,
                 youden_index = youden_index,
                 sensitivity = sensitivity,
                 specificity = specificity,
                 geometric_mean = geometric_mean,
                 fdr = fdr,
                 mcc = mcc,
                 auc = auc,
                 corrclass = corrclass,
                 n_betas = n_betas,
                 n_total_var = n_total_var,
                 n_predicted_zeros = n_predicted_zeros,
                 n_predicted_non_zeros = n_predicted_non_zeros,
                 n_caught_betas = n_caught_betas,
                 n_non_caught_betas = n_non_caught_betas,
                 n_caught_zero = n_caught_zero,
                 n_zero_not_caught = n_zero_not_caught,
								 input_parameters = list(
                   df = df,
									 X = X,
									 y = y,
									 lambda = lambda,
									 penalty = penalty,
									 w = w,
                   beta_start_input = beta_start_input,
									 beta_start_default = beta_start_default,
                   max.print = max.print,
									 alpha = alpha,
									 a1 = a1,
									 a2 = a2,
                   regressors_betas = regressors_betas,
									 fold = fold,
									 max_iter = max_iter,
                   trend = trend,
									 delta = delta,
									 max_alpha = max_alpha,
                   stepsizeShrink = stepsizeShrink,
									 min_alpha = min_alpha,
                   convergence_error = convergence_error,
									 trace = trace,
									 seed = seed,
                   c_zero_fixed = c_zero_fixed,
									 kernel = kernel,
                   zeros_stay_zeros_from_iteration = zeros_stay_zeros_from_iteration
                 )
               )

  names(results)[1:2] <- paste0(c("betas_hat_", "pye_KS_"), penalty)

  estimation_time <- difftime(Sys.time(), start_time, units = "mins")
  if (trace %in% c(1, 2)) {
    cat("Estimation time:", format(estimation_time, units = "mins"), "\n\n\n")
  }

  results <- c(results , estimation_time = estimation_time, niter = estim$tot_iters)

  return(results)
}










#---------------> cross-validation of pye KS <----------------

#create the output class of the pye.cv function
setClass(Class = "pye_KS_CV_output",
         representation(
           penalty = "character",
           penalty_covariates = "character",
           kernel = "character",
           pye_KS_L12 = "ANY",
           pye_KS_L1 = "ANY",
           pye_KS_EN = "ANY",
           pye_KS_SCAD = "ANY",
           pye_KS_MCP = "ANY",
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
           betas = "list",
					 c_pye = "list",
           gammas = "list"
         )
)

# to extract just part of the estim of pye
subset_pye_KS <- function(df, X, y, betas, lambda, c, w, fold, alpha, trace, a1, a2, penalty,
                           cv_time, niter, kernel, c_function_of_covariates = FALSE) {

  #ID <- rownames(df)
  pye_KS_result <- pye_KS(df = df[, names(df) != "ID", drop = FALSE], X = X, y = y, betas = betas, lambda = lambda, c = c,
                          w = w, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty, prediction = TRUE, kernel = kernel)

  #z_hat <- pye_KS_result$z_hat$z_hat

  #optimal cutpoint (the one estimated in the training)
  TP <- pye_KS_result$TP
  TN <- pye_KS_result$TN
  FP <- pye_KS_result$FP
  FN <- pye_KS_result$FN
  c_hat <- c

  if (trace %in% c(1, 2)) {
    #print the results only if c_function_of_covariates = FALSE
    if (c_function_of_covariates == FALSE) {
      cat("-> Results on the TEST SET of pye\n")
      cat("-> algorithm: pye_KS_proximal_gradient_method ; ")
      if (!is.null(fold)) {cat("fold =", fold, "; \n")}
      visualize_betas <- c(betas[which(betas != 0)], c_hat)
      cat("lambda:", lambda, "; weight:", w, "; penalty:", penalty, "; pye_KS:", getElement(pye_KS_result,  paste0("pye_", penalty)), "; youden_index:", pye_KS_result$youden_index, "; sensitivity:", pye_KS_result$sensitivity, "; specificity:", pye_KS_result$specificity, "; geometric_mean:", pye_KS_result$geometric_mean, "; fdr:", pye_KS_result$fdr, "; mcc:", pye_KS_result$mcc, "; auc:", pye_KS_result$auc, "; corrclass:", pye_KS_result$corrclass, " \n")
      cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, ";  betas: \n")
      print(visualize_betas)
      cat("Cross-validation time:", cv_time, "; Number of iterations:", niter, "\n\n\n")
    }
  }

  return(list(pye_KS_L12 = pye_KS_result$pye_L12,
              pye_KS_L1 = pye_KS_result$pye_L1,
              pye_KS_EN = pye_KS_result$pye_EN,
              pye_KS_SCAD = pye_KS_result$pye_SCAD,
              pye_KS_MCP = pye_KS_result$pye_MCP,
              youden_index = pye_KS_result$youden_index,
              sensitivity = pye_KS_result$sensitivity,
              specificity = pye_KS_result$specificity,
              geometric_mean = pye_KS_result$geometric_mean,
              fdr = pye_KS_result$fdr,
              mcc = pye_KS_result$mcc,
              auc = pye_KS_result$auc,
              corrclass = pye_KS_result$corrclass,
              TP = TP,
              TN = TN,
              FP = FP,
              FN = FN,
              z_hat = pye_KS_result$z_hat))
}

#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster
#' @importFrom methods new
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext file_ext
#' @noRd
#' @keywords internal
pye_KS.cv <- function (df, X, y, C, lambda, tau, w, w_g,
                       trace = 1, alpha, alpha_g, a1, a2,
											 penalty, penalty_g, folds_i, k,
                       regressors_betas = NULL,
											 kernel_g, a1_g, a2_g,
                       trend_g, gamma_start_input,
											 gamma_start_default,
                       regressors_gammas = NULL, max_iter_g,
											 delta_g, max_alpha_g, stepsizeShrink_g,
                       min_alpha_g, convergence_error_g,
                       pye_KS_L12, pye_KS_L1, pye_KS_EN, pye_KS_SCAD, pye_KS_MCP,
                       auc, aauc, aYI, youden_index,
											 sensitivity, specificity, geometric_mean, fdr,
                       mcc, corrclass, n_betas, n_gammas, used_cores,
											 trend, delta, max_alpha, kernel,
                       beta_start_input, beta_start_default,
											 c_zero_fixed, max_iter, min_alpha, convergence_error,
                       stepsizeShrink, c_function_of_covariates,
											 simultaneous, run_aauc, log_file) { #, long_suffix) {

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  greek <- if (c_function_of_covariates) "tau" else "lambda"
  if (trace %in% c(1, 2)) {
    cat("----------------------------------------------------------------\n")
    cat("|       starting with the", k, "-th fold for the CV of ", greek,"      |\n")
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

    if (trace %in% c(1, 2)) cat("Parallel computing for cross-validation started on ", length(cl), "cores.\n")

    # Load the package on all workers
    parallel::clusterCall(cl, function() library(pye))

    # --- Fit primary models ---
    parallel::clusterExport(cl, c("train_df", "X", "y", "C", "trace", "alpha", "alpha_g", "penalty", "penalty_g",
                                  "w", "a1_g", "a2_g", "trend_g",
                                  "regressors_betas", "k", "trend", "delta", "c_function_of_covariates",
                                  "max_alpha", "kernel", "beta_start_input", "beta_start_default",
                                  "c_zero_fixed", "a2", "a1", #"long_suffix",
                                  "max_iter", "min_alpha", "convergence_error", "stepsizeShrink"), envir = environment())

    fitted_models <- parallel::parLapply(cl, lambda, function(x) pye_KS_estimation(df = train_df, X = X, y = y, lambda = x,
		                                                                                w = w,
                                                                                    beta_start_input = beta_start_input,
                                                                                    beta_start_default = beta_start_default,
                                                                                    trace = trace,
                                                                                    alpha = alpha,
                                                                                    a1 = a1, a2 = a2,
                                                                                    penalty = penalty, max_iter = max_iter,
                                                                                    convergence_error = convergence_error,
                                                                                    regressors_betas = regressors_betas, fold = k,
                                                                                    trend = trend,
                                                                                    stepsizeShrink = stepsizeShrink,
                                                                                    delta = delta, max_alpha = max_alpha,
                                                                                    min_alpha = min_alpha,
                                                                                    kernel = kernel,
                                                                                    c_zero_fixed = c_zero_fixed)) #, long_suffix = long_suffix))
  } else {

    fitted_models <- lapply(lambda, function(x) pye_KS_estimation(df = train_df, X = X, y = y,
                                                                   lambda = x, w = w,
                                                                   beta_start_input = beta_start_input,
                                                                   beta_start_default = beta_start_default, trace = trace,
                                                                   alpha = alpha,
                                                                   a1 = a1, a2 = a2, penalty = penalty, max_iter = max_iter,
                                                                   min_alpha = min_alpha,
                                                                   convergence_error = convergence_error,
                                                                   regressors_betas = regressors_betas, fold = k,
                                                                   trend = trend,
                                                                   stepsizeShrink = stepsizeShrink,
                                                                   delta = delta, max_alpha = max_alpha,
                                                                   kernel = kernel,
                                                                   c_zero_fixed = c_zero_fixed)) #, long_suffix = long_suffix))
  }


  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat"))
  #if ((c_function_of_covariates == TRUE) & (length(lambda) == 1)) {
  #  #in this case we have to estimate the combination z_hat over all the dataset
  #  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat")[-test_i, ])
  #} else {
  #  z_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "z_hat"))
  #}

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

  #here I have to put the covYI estimation
  #computing c
  if (c_function_of_covariates == FALSE) {

    #measures on the train set - if we don't compute c with covariates, we are only dependent of lambda:
    temp_pye <- get(paste0("pye_KS_", penalty))

    for (i in seq_along(lambda)) {
      #measures on the train set
      #these are multiple tables based on the number of considered TAUs
      temp_pye[[i]]$train[k, ] <- getElement(fitted_models[[i]], paste0("pye_KS_", penalty))
      auc[[i]]$train[k, ]  <- getElement(fitted_models[[i]], "auc")
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

    gammas_hat <- lapply(seq_along(lambda), function(y) NA)
    names(gammas_hat) <- lambda

  } else {

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

    temp_pye <- get(paste0("pye_KS_", penalty))

    for (i in seq_along(lambda)) {
      # Measures on the train set
      temp_pye[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("covYI_KS_", penalty_g))))
      auc[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "auc")))
      aauc[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "aauc")))
      aYI[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "aYI")))
      youden_index[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "sensitivity")))
      specificity[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "fdr")))
      mcc[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "mcc")))
      corrclass[[i]]$train[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "corrclass")))

      n_gammas[[i]][k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], "n_gammas")))
      gammas_hat[[i]] <- lapply(seq_along(tau), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("gammas_hat_", penalty_g)))
      names(gammas_hat[[i]]) <- paste("tau", tau, sep = "=")
    }
  }

  # n_betas is only based on lambda, not tau!
  n_betas[k, ] <- unlist(lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "n_betas")))

  # create a list of the results - it is only based on lambda, not tau!
  betas_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], paste0("betas_hat_", penalty)))
  names(betas_hat) <- paste("lambda", lambda, sep = "=")

  # c_hat is the fixed value of c in case we don't use the covariates to estimate a patient's specific cut-off value
  c_hat <- lapply(seq_along(lambda), function(x) getElement(fitted_models[[x]], "c_hat"))
  names(c_hat) <- paste("lambda", lambda, sep = "=")

  betas <- betas_hat
  gammas <- gammas_hat
	c_pye <- c_hat

  # --- Test Set Evaluation ---
  all_measures_test <- mapply(function(x, xx, z, zz) subset_pye_KS(df = test_df, X = X, y = y, betas = x, lambda = z,
                                                           c = xx, w = w, fold = k, alpha = alpha, a1 = a1, a2 = a2, penalty = penalty,
                                                           cv_time = fitted_models[[zz]]$estimation_time,
                                                           niter = fitted_models[[zz]]$niter, kernel = kernel, trace = trace,
                                                           c_function_of_covariates = c_function_of_covariates),
                              betas_hat, c_hat, lambda, seq_along(lambda))

  # --- Organize test set results ---
  if (c_function_of_covariates == FALSE) {
    # Measures on the test set - if we don't compute c with covariates, we are only dependent of lambda:
    for (i in seq_along(lambda)) {
      temp_pye[[i]]$test[k, ] <- all_measures_test[paste0("pye_KS_", penalty), ][[i]]
      auc[[i]]$test[k, ]  <- all_measures_test["auc", ][[i]]
      aauc[[i]]$test[k, ]  <- NA
      aYI[[i]]$test[k, ]  <- NA
      youden_index[[i]]$test[k, ] <- all_measures_test["youden_index", ][[i]]
      sensitivity[[i]]$test[k, ] <- all_measures_test["sensitivity", ][[i]]
      specificity[[i]]$test[k, ]  <- all_measures_test["specificity", ][[i]]
      geometric_mean[[i]]$test[k, ]  <- all_measures_test["geometric_mean", ][[i]]
      fdr[[i]]$test[k, ]  <- all_measures_test["fdr", ][[i]]
      mcc[[i]]$test[k, ]  <- all_measures_test["mcc", ][[i]]
      corrclass[[i]]$test[k, ] <- all_measures_test["corrclass", ][[i]]
    }
    assign(paste0("pye_KS_", penalty), temp_pye)
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
                                                        z = "z_hat", y = y, C = C1, w = w_g,
                                                        gammas = x[[1]], tau = x[[2]], kernel = kernel_g, alpha = alpha_g,
                                                        a1 = a1_g, a2 = a2_g, penalty = penalty_g, prediction = TRUE,
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
        cat("With lambda:", lambda[i], ",  penalty:", penalty, "and weight:", w,  ". Accuracy measures only pye: \n")
        visualize_betas <- c(betas_hat[[i]][which(betas_hat[[i]] != 0)], stats::setNames(c_hat[[i]], "c"))
        cat("pye_KS:", all_measures_test[paste0("pye_KS_", penalty), ][[i]], "; youden_index:", all_measures_test["youden_index", ][[i]], "; sensitivity:", all_measures_test["sensitivity", ][[i]], "; specificity:", all_measures_test["specificity", ][[i]], "; geometric_mean:", all_measures_test["geometric_mean", ][[i]], "; fdr:", all_measures_test["fdr", ][[i]], "; mcc:", all_measures_test["mcc", ][[i]], "; auc:", all_measures_test["auc", ][[i]], "; corrclass:", all_measures_test["corrclass", ][[i]], "; \n")
        cat("TP:", all_measures_test["TP", ][[i]], "; TN:", all_measures_test["TN", ][[i]], "; FP:", all_measures_test["FP", ][[i]], "; FN:", all_measures_test["FN", ][[i]], ";  betas: \n")
        print(visualize_betas)
        cat("Cross-validation time:", fitted_models[[i]]$estimation_time, "; Number of iterations:", fitted_models[[i]]$niter, "\n\n")

        for (ii in seq_along(cov_results[[i]])) {
          cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend, "version ; \n")
          visualize_gammas <- c(gammas_hat[[i]][[ii]][which(gammas_hat[[i]][[ii]] != 0)])
          cat("tau:", tau[ii], "; weight:", w_g, "; penalty:", penalty_g, "; ", paste0("covYI_KS_", penalty_g), ":" , getElement(cov_results[[i]][[ii]], paste0("covYI_KS_", penalty_g)), "; youden_index:", cov_results[[i]][[ii]]$youden_index, "; aYI:", cov_results[[i]][[ii]]$aYI, "; sensitivity:", cov_results[[i]][[ii]]$sensitivity, "; specificity:", cov_results[[i]][[ii]]$specificity, "; geometric_mean:", cov_results[[i]][[ii]]$geometric_mean, "; fdr:", cov_results[[i]][[ii]]$fdr, "; mcc:", cov_results[[i]][[ii]]$mcc, "; auc:", cov_results[[i]][[ii]]$auc, "; aauc:", cov_results[[i]][[ii]]$aauc, "; corrclass:", cov_results[[i]][[ii]]$corrclass, "; \n")
          cat("TP:", cov_results[[i]][[ii]]$TP, "; TN:", cov_results[[i]][[ii]]$TN, "; FP:", cov_results[[i]][[ii]]$FP, "; FN:", cov_results[[i]][[ii]]$FN, "; gammas: \n")
          print(visualize_gammas)
          cat("Cross-validation time:", fitted_gammas[[i]][[ii]]$estimation_time, "; Number of iterations:", fitted_gammas[[i]][[ii]]$niter, "\n\n\n")
        }
				cat("\n")
      }

      # Store test measures
      #NB: when c_function_of_covariates = TRUE, the measure pye_KS coincides with the covYI_KS measure
      temp_pye[[i]]$test[k, ] <- unlist(lapply(seq_along(tau), function(yy) getElement(cov_results[[i]][[yy]], paste0("covYI_KS_", penalty_g))))
      assign(paste0("pye_KS_", penalty) , temp_pye)
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

  return(methods::new("pye_KS_CV_output", penalty = penalty,
                                          penalty_covariates = penalty_g,
                                          kernel = kernel,
                                          pye_KS_L12 = pye_KS_L12,
                                          pye_KS_L1 = pye_KS_L1,
                                          pye_KS_EN = pye_KS_EN,
                                          pye_KS_SCAD = pye_KS_SCAD,
                                          pye_KS_MCP = pye_KS_MCP,
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
                                          n_betas = n_betas,
                                          n_gammas = n_gammas,
                                          betas = betas,
																					c_pye = c_pye,
                                          gammas = gammas))
}








#' @title Cross-Validation for Optimal pye KS Regularization Parameter
#' Selection
#'
#' @description This function performs cross-validation to select the optimal
#' lambda (and potentially tau) values for estimating biomarker coefficients
#' (betas) and the cut-off point (c), using the Penalized Youden Index (pye)
#' based on Kernel Smooth density estimation. It can also integrate covariate
#' information via covYI.
#'
#' @param n_folds A `numeric` value specifying the number of folds for the
#'   cross-validation.
#' @param df A `data.frame` containing the input dataset. It must include the
#'   columns specified in `X` and `y`, and optionally `C`.
#' @param X A `character` vector specifying the names of the regressor
#'   variables (biomarkers) to consider in the estimation. Alternatively, a
#'   `data.frame` containing only these regressors can be provided, in which
#'   case their column names will be used. Defaults to all columns in `df`
#'   not specified as `y` or `C`.
#' @param y A `character` string specifying the name of the target binary
#'   variable (outcome). This variable must contain only values 0 and 1.
#'   Alternatively, a single-column `data.frame` containing only the target
#'   variable can be provided, in which case its column name will be used.
#'   Default is "y".
#' @param C A `character` vector specifying the names of covariate variables.
#'   Alternatively, a `data.frame` containing these covariates can be provided,
#'   in which case their column names will be used. Default is `NULL`,
#'   indicating no covariates are used.
#' @param lambda A `numeric` vector of penalization parameters for the
#'   regressors `X`. The cross-validation will evaluate models across these
#'   lambda values.
#' @param w A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in pye. Sensitivity is weighted by `w` and specificity by
#'   `1 - w`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param w_g A `numeric` value between 0 and 1 specifying the weight for the
#'   Weighted Youden Index in covYI. Sensitivity is weighted by `w_g` and specificity by
#'   `1 - w_g`. Default is 0.5 (which corresponds to the standard Youden Index).
#' @param trace A `numeric` value controlling the verbosity of the output.
#'   `0`: no visualization; `1`: visualize only the final result (default);
#'   `2`: visualize all steps during optimization.
#' @param alpha A `numeric` value between 0 and 1, used as the mixing parameter
#'   for the Elastic-Net penalization term applied to `X`. Only relevant if
#'   `penalty` is "EN". Default is 0.5.
#' @param alpha_g A `numeric` value between 0 and 1, used as the mixing
#'   parameter for the Elastic-Net penalization term applied to covariates `C`
#'   in covYI. Only relevant if `penalty_g` is "EN". Default is 0.5.
#' @param tau A `numeric` vector of penalization parameters for the covariates
#'   `C` in covYI. If 0 (default), no penalization term is applied to covariates.
#'   If `c_function_of_covariates` is `TRUE`, `tau` cannot be `NULL` or contain
#'   only zeros.
#' @param a1 A `numeric` value specifying the regularization parameter 'a' for
#'   the SCAD penalty (as defined by Fan and Li). Only relevant if `penalty` is
#'   "SCAD". Default is 3.7.
#' @param a2 A `numeric` value specifying the regularization parameter 'gamma'
#'   for the MCP penalty (as defined by Zhang). Only relevant if `penalty` is
#'   "MCP". Default is 3.0.
#' @param penalty A `character` string indicating the type of penalty considered
#'   for the pye calculation. Must be one of "L12", "L1", "EN" (Elastic-Net),
#'   "SCAD", or "MCP". Default is "L1".
#' @param regressors_betas numeric vector. An optional vector containing the
#'   "true" beta coefficients, if known. This is used to calculate additional
#'   performance measures related to variable selection accuracy. Default is
#'   `NULL`.
#' @param seed A `numeric` value to fix the random seed for reproducibility.
#'   Default is 1.
#' @param used_cores A `numeric` value indicating the number of cores to use
#'   for parallelization. If equal to 1 (default), no parallelization is
#'   adopted.
#' @param trend A `character` string indicating the optimization algorithm
#'   variant to use for pye estimation. If "monotone" (default), the Monotone
#'   Accelerated Proximal Gradient (mmAPG) algorithm is used. If "nonmonotone",
#'   the Non-Monotone Accelerated Proximal Gradient (mnmAPG) algorithm is used.
#' @param delta A `numeric` value specifying the parameter for the convergence
#'   condition of the optimization algorithm for pye. Default is 1e-5.
#' @param max_alpha A `numeric` value representing the maximum allowed value for
#'   the step-size parameter `alpha` in the backtracking line-search for pye
#'   optimization. Default is 10000.
#' @param kernel A `character` string indicating the kernel type to use for the
#'   estimation of the density function for pye. Currently, only "gaussian" is
#'   tested and supported. Default is "gaussian".
#' @param beta_start_input A `numeric` vector of specific starting points for
#'   the beta coefficients. If `NULL` (default), the starting points are
#'   determined by `beta_start_default`. If provided, its length must match the
#'   number of regressors in `X`.
#' @param max_iter A `numeric` value specifying the maximum number of iterations
#'   for the optimization algorithms (mmAPG and mnmAPG) for pye. Default is 10000.
#' @param min_alpha A `numeric` value representing the minimum allowed value for
#'   the step-size parameter `alpha` for pye optimization. Default is 1e-10.
#' @param convergence_error A `numeric` value specifying the error tolerance for
#'   considering the pye optimization algorithm converged. Default is 1e-7.
#' @param stepsizeShrink A `numeric` value between 0 and 1, used to adjust the
#'   step-size in the backtracking line-search within the pye optimization.
#'   Values closer to 1 lead to higher accuracy but longer estimation times.
#'   Default is 0.8.
#' @param beta_start_default A `character` string defining the default starting
#'   points for betas if `beta_start_input` is `NULL`. If "zeros" (default),
#'   betas start with a vector of all zeros. If "corr", betas start with the
#'   absolute value of the correlation of each regressor with the target
#'   variable `y`.
#' @param scaling A `logical` value. If `TRUE`, the input dataset `df` is
#'   scaled internally. If `FALSE` (default), no scaling is performed.
#' @param c_zero_fixed A `logical` value. If `TRUE`, the estimation process
#'   considers the cut-off point `c` as fixed and equal to zero, which can
#'   reduce estimation complexity. If `FALSE` (default), `c` can vary and is
#'   estimated by the pye optimization.
#' @param c_function_of_covariates A `logical` value. If `TRUE`, covYI is used
#'   to estimate the cut-off point `c` as a function of the covariate
#'   information (`C`). If `FALSE` (default), covariate information is ignored
#'   for the cut-off estimation.
#' @param simultaneous A `logical` value. If `c_function_of_covariates` is
#'   `TRUE`, this parameter defines if `gammas` (covariate coefficients) need
#'   to be estimated simultaneously with `betas` or as a second step. Default
#'   is `FALSE`, meaning a sequential estimation.
#' @param measure_to_select_lambda A `character` string specifying the measure
#'   used to select `lambda` if `simultaneous` is `FALSE` (i.e., when the
#'   cross-validation process selects `lambda` first and then `tau`). Must be
#'   one of "auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", or
#'   "ccr". Default is "ccr" (correct classification rate).
#' @param penalty_g A `character` string indicating the type of penalty for
#'   covYI (when `c_function_of_covariates` is `TRUE`). Must be one of "L12",
#'   "L1", "EN" (Elastic-Net), "SCAD", or "MCP". Default is "L1".
#' @param kernel_g A `character` string indicating the kernel type to use for
#'   the density estimation in covYI. Currently, only "gaussian" is tested and
#'   supported. Default is "gaussian".
#' @param a1_g A `numeric` value specifying the regularization parameter 'a' for
#'   the SCAD penalty in covYI. Only relevant if `penalty_g` is "SCAD".
#'   Default is 3.7.
#' @param a2_g A `numeric` value specifying the regularization parameter 'gamma'
#'   for the MCP penalty in covYI. Only relevant if `penalty_g` is "MCP".
#'   Default is 3.0.
#' @param trend_g A `character` string indicating the optimization algorithm
#'   variant to use for covYI. If "monotone" (default), mmAPG is used. If
#'   "nonmonotone", mnmAPG is used.
#' @param gamma_start_input A `numeric` vector of specific starting points for
#'   the gamma coefficients in covYI. If `NULL` (default), starting points are
#'   determined by `gamma_start_default`.
#' @param gamma_start_default A `character` string defining the default starting
#'   point of gamma coefficients. If "zeros" (default), they start with all
#'   zero values. If "corr", they start with the absolute correlation of each
#'   covariate with the target variable.
#' @param regressors_gammas A `numeric` vector containing the true gamma
#'   coefficients (if known, for simulation/testing purposes). Default is `NULL`.
#' @param max_iter_g A `numeric` value specifying the maximum number of
#'   iterations for the optimization algorithms (mmAPG and mnmAPG) in covYI.
#'   Default is 10000.
#' @param delta_g A `numeric` value specifying the parameter for the convergence
#'   condition of the optimization algorithm of covYI. Default is 1e-5.
#' @param max_alpha_g A `numeric` value representing the maximum allowed value
#'   for the step-size parameter `alpha` in covYI. Default is 10000.
#' @param stepsizeShrink_g A `numeric` value between 0 and 1, used to adjust the
#'   step-size in the backtracking line-search within the covYI optimization.
#'   Values closer to 1 lead to higher accuracy but longer estimation times.
#'   Default is 0.8.
#' @param min_alpha_g A `numeric` value representing the minimum allowed value
#'   for the step-size parameter `alpha` in covYI. Default is 1e-12.
#' @param convergence_error_g A `numeric` value specifying the error tolerance
#'   for considering the covYI algorithm converged. Default is 1e-7.
#' @param run_aauc A `logical` value. If `FALSE` (default), the aAUC and aYI
#'   measures are not computed, saving estimation time if not requested.
# #' @param long_suffix A `character` string. If provided, it's used as a prefix
# #'   to identify longitudinal variables (e.g., "visit_"). For proper management
# #'   of longitudinal variables or coefficients referring to an original
# #'   longitudinal dimension, `long_suffix` is expected to be followed by a
# #'   number (e.g., "visit_1", "visit_2"). If `NULL` (default), no special
# #'   handling for longitudinal variables is applied, presuming the input `df`
# #'   does not contain such variables.
#' @param log_file Character. Path to a file for logging output from parallel
#'   workers. If `NULL`, output goes to the console.
#'   Default is "log_pye_ks_models.txt".
#'
#' @return A `list` with the following components:
#'   \item{penalty}{\code{character}. The type of penalty used for regressors X.}
#'   \item{penalty_g}{\code{character}. The type of penalty used for covariates C
#'   (if applicable).}
#'   \item{kernel}{\code{character}. The kernel type used for density estimation of pye.}
#'   \item{cv_time}{\code{difftime}. The total time taken for the cross-validation
#'     process, in minutes.}
#'   \item{pye_KS_L12, pye_KS_L1, pye_KS_EN, pye_KS_SCAD, pye_KS_MCP}{\code{list}.
#'     For each penalty type, a list containing `train` and `test` matrices.
#'     These matrices store the Penalized Youden Index values for each fold
#'     and each combination of `lambda` and `tau` (if applicable). Their
#'     dimensions are `n_folds` rows by `length(tau)` columns (or 1 if `length(tau)`
#'     is 1).}
#'   \item{auc_first_step, aauc_first_step, aYI_first_step, youden_index_first_step,
#'     sensitivity_first_step, specificity_first_step, geometric_mean_first_step,
#'     fdr_first_step, mcc_first_step, corrclass_first_step, auc, aauc, aYI, youden_index,
#'     sensitivity, specificity, geometric_mean,fdr, mcc, corrclass}{\code{list}.
#'     Similar to `pye_KS_penalty` components, these lists contain `train` and
#'     `test` matrices for each respective performance measure (AUC, adaptive AUC,
#'     adaptive Youden Index, unpenalized Youden Index, sensitivity, specificity,
#'     geometric mean, false discovery rate, Matthews Correlation Coefficient, and
#'     correct classification rate) of the first-step method and of the final method
#'     (which coincide if c_function_of_covariates is FALSE), across all folds and
#'     `lambda`/`tau` combinations.}
#'   \item{lambda_hat_yi, lambda_hat_auc, lambda_hat_aauc,
#'     lambda_hat_aYI, lambda_hat_ccr, lambda_hat_sen,
#'     lambda_hat_spc, lambda_hat_gm, lambda_hat_pye}{\code{numeric}.
#'     The optimal `lambda` value selected based on maximizing the corresponding
#'     performance measure (e.g., `_yi` for Youden Index, `_auc` for AUC). If no
#'     optimal lambda is found, it will be `NA`.}
#'   \item{tau_hat_yi, tau_hat_auc, tau_hat_aauc,
#'     tau_hat_aYI, tau_hat_ccr, tau_hat_sen,
#'     tau_hat_spc, tau_hat_gm, tau_hat_pye}{\code{numeric}.
#'     The optimal `tau` value selected based on maximizing the corresponding
#'     performance measure. This is only relevant and will be a value other
#'     than `NA` if `c_function_of_covariates` is `TRUE`. Otherwise, it will be `NA`.}
#'   \item{c_function_of_covariates}{\code{logical}. Indicates whether the cut-off
#'     point was estimated as a function of covariates.}
#'   \item{simultaneous}{\code{logical}. Indicates whether betas and gammas were
#'     estimated simultaneously (if `c_function_of_covariates` was `TRUE`).}
#'   \item{measure_to_select_lambda}{\code{character}. The measure used to select
#'     the optimal `lambda` (if `simultaneous` is `FALSE`).}
#'   \item{lambda_star}{\code{numeric}. The single optimal lambda value chosen
#'     (when `c_function_of_covariates` is `TRUE` and `simultaneous` is `FALSE`),
#'     based on `measure_to_select_lambda`. `NA` otherwise.}
#'   \item{n_betas}{\code{matrix}. A matrix showing the number of non-zero (selected)
#'     betas for each fold and each `lambda` value. Its dimensions are `n_folds`
#'     rows by `length(lambda)` columns.}
#'   \item{n_betas_star}{\code{matrix}. If a two-step optimization (sequential
#'     lambda then tau) is performed (`c_function_of_covariates = TRUE` and
#'     `simultaneous = FALSE`), this matrix shows the number of non-zero betas for
#'     each fold using `lambda_star`. Dimensions are `n_folds` rows by 1 column.}
#'   \item{n_gammas}{\code{list}. A list, where each element corresponds to a `lambda`
#'     value and contains a matrix showing the number of non-zero (selected)
#'     gammas for each fold and each `tau` value. Dimensions are `n_folds` rows
#'     by `length(tau)` columns.}
#'   \item{betas}{\code{matrix} or \code{list}. The estimated beta coefficients.
#'     If `length(lambda)` is 1, a matrix with `n_folds` rows. Otherwise, a list
#'     of beta vectors for each fold and each lambda.}
#'   \item{c_pye}{\code{list}. A list where each element corresponds to a
#'     \code{lambda} value and contains a matrix showing the estimated optimal
#'     cut-off point \eqn{c} for each fold and each \code{lambda} value. If
#'     \code{c_zero_fixed} is \code{TRUE}, all entries will be 0.}
#'   \item{betas_star}{\code{matrix} or \code{list}. If a two-step optimization is
#'     performed, the estimated beta coefficients using `lambda_star`. A matrix
#'     with `n_folds` rows if `length(lambda_star)` is 1, otherwise a list of
#'     beta vectors.}
#'   \item{gammas}{\code{list}. A list, where each element corresponds to a `lambda`
#'     value (or `lambda_star` if sequential) and contains a matrix or vector
#'     of estimated gamma coefficients for each fold and `tau` value.}
#'
#' @examples
#' # Load the package
#' library(pye)
#'
#' # 1. Simulate data for the example
#' \donttest{
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$train_df_scaled
#' X <- sim_data$X
#' y <- sim_data$y
#' C <- sim_data$C
#' regressors_betas <- sim_data$nregressors # True betas for evaluation
#' regressors_gammas <- sim_data$ncovariates # True gammas for evaluation
#'
#' # 2. Set cross-validation parameters
#' penalty <- "SCAD" # Penalty for betas in pye estimation
#' penalty_g <- "L12" # Penalty for gammas in covYI estimation
#' trend <- "monotone" # Trend for the KS estimation
#' alpha <- 0.5
#' c_function_of_covariates <- TRUE # Use covariates for 'c' estimation
#' used_cores <- 1 # For this example, no parallelization
#' max_iter <- 10 # Keep iterations low for a quick example run
#'
#' # 3. Calibrate lambda_max and lambda_min for betas (pye estimation)
#' lambda_seq <- create_lambda(n = 3, lmax = 1.5, lmin = 0.05)
#' lambda_seq <- as.numeric(formatC(lambda_seq, format = "e", digits = 9))
#'
#' # 4. Calibrate tau_max and tau_min for gammas (covYI estimation), if
#' tau_seq <- create_lambda(n = 3, lmax = 0.15, lmin = 0.005)
#' tau_seq <- as.numeric(formatC(tau_seq, format = "e", digits = 9))
#'
#' # 5. Run the cross-validation
#' pye_cv_result <- pye_KS_compute_cv(
#'   n_folds = 3,
#'   df = df,
#'   X = X,
#'   y = y,
#'   C = C,
#'   lambda = lambda_seq,
#'   tau = tau_seq,
#'   trace = 1, # Show final results
#'   alpha = alpha,
#'   penalty = penalty,
#'   regressors_betas = regressors_betas,
#'   regressors_gammas = regressors_gammas,
#'   seed = 1,
#'   used_cores = used_cores,
#'   trend = trend,
#'   max_iter = max_iter,
#'   c_function_of_covariates = c_function_of_covariates,
#'   measure_to_select_lambda = "ccr",
#'   penalty_g = penalty_g,
#'   trend_g = trend,
#'   max_iter_g = max_iter
#' )
#'
#' # 6. Print results and access optimal lambda/tau
#' cat("\nOptimal Lambda (based on CCR):", pye_cv_result$lambda_hat_ccr, "\n")
#' if (c_function_of_covariates == TRUE) {
#'   cat("Optimal Tau (based on CCR):", pye_cv_result$tau_hat_ccr, "\n")
#' }
#'
#' # You can access other results like:
#' pye_cv_result$auc # AUC values
#' pye_cv_result$n_betas # Number of non-zero betas for each lambda
#' pye_cv_result$n_gammas # Number of non-zero gammas for each tau
#' }
#'
#' @export
pye_KS_compute_cv <- function (n_folds, df, X = NULL, y = "y", C = NULL, lambda,
                                w = 0.5, w_g = 0.5, trace = 1, alpha = 0.5, alpha_g = 0.5, tau = 0,
                                a1 = 3.7, a2 = 3, penalty = "L1", regressors_betas = NULL,
                                seed = 1, used_cores = 1, trend = "monotone", delta = 1e-5, max_alpha = 10000,
                                kernel = "gaussian", beta_start_input = NULL,
                                max_iter = 10000, #reduced for the CV
                                min_alpha = 1e-10,
                                convergence_error = 1e-7,
                                stepsizeShrink = 0.8,
                                beta_start_default = "zeros", scaling = FALSE, c_zero_fixed = FALSE,
                                c_function_of_covariates = FALSE,
                                simultaneous = FALSE,
                                measure_to_select_lambda = "ccr",
                                penalty_g = "L1", kernel_g = "gaussian", a1_g = 3.7, a2_g = 3,
                                trend_g = "monotone", gamma_start_input = NULL, gamma_start_default = "zeros",
                                regressors_gammas = NULL, max_iter_g = 10000,
                                delta_g = 1e-5, max_alpha_g = 10000,
                                stepsizeShrink_g = 0.8, min_alpha_g = 1e-12, convergence_error_g = 1e-7,
                                run_aauc = FALSE, log_file = "log_pye_ks_models.txt") { #, long_suffix=NULL) {

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
    cat("Setting 'tau' equal to zero since 'c_function_of_covariates' is FALSE \n")
  }

  # Check for "ID" column conflict
  if ("ID" %in% colnames(df)) {stop("The column name 'ID' already exists in 'df'. Please rename or remove it, as it's used internally.")}

  # Create a working data frame with ID and selected columns
  ID <- rownames(df)
  #df1 <- cbind(ID=ID, y = df[, (names(df) %in% c(y))], df[, (names(df) %in% c(X, C))]) #OLD
  df1 <- cbind(ID = ID, y = df[, y, drop = FALSE], df[, c(X, C), drop = FALSE])

  # Validate target variable y (0, 1)
  if (is.factor(df1[[y]]) || is.character(df1[[y]])) {df1[[y]] <- as.numeric(as.character(df1[[y]]))}
  if (!all(sort(unique(df1[[y]])) %in% c(0, 1)) || anyNA(df1[[y]])) {stop("The target variable 'y' must contain only values 0 and 1 (excluding NA).")}
  if (length(unique(df1[[y]])) < 2) {stop("The target variable 'y' must contain at least two unique values (0 and 1).")}

  # Further input validation
  if (!is.numeric(lambda) || length(lambda) < 1 || any(lambda < 0)) stop("Parameter 'lambda' must be a numeric vector of non-negative values.")
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
  if (!is.logical(scaling)) {stop("Parameter 'scaling' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_zero_fixed)) {stop("Parameter 'c_zero_fixed' must be a logical (TRUE/FALSE).")}
  if (!is.logical(c_function_of_covariates)) {stop("Parameter 'c_function_of_covariates' must be a logical (TRUE/FALSE).")}
  if (!is.logical(run_aauc)) {stop("Parameter 'run_aauc' must be a logical (TRUE/FALSE).")}
  valid_kernels <- c("gaussian", "normal", "uniform", "rectangular", "triangular", "epanechnikov",
                     "biweight", "triweight", "tricube", "parzen", "cosine", "optcosine")
  # NB: kernels: "normal", "uniform", "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "tricube", "parzen",
  # "cosine", "optcosine", have not been deeply tested. Most of the work has been done with "gaussian" kernel
  if (!(kernel %in% valid_kernels)) {stop(paste0("Parameter 'kernel' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  if (!(kernel_g %in% valid_kernels)) {stop(paste0("Parameter 'kernel_g' is not supported. Must be one of: ", paste(valid_kernels, collapse = ", "), "."))}
  valid_measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
  if (!(measure_to_select_lambda %in% valid_measures)) {stop("Parameter 'measure_to_select_lambda' must be one of: ", paste(valid_measures, collapse = ", "))}
  if (!is.numeric(used_cores) || length(used_cores) != 1 || used_cores <= 0 || floor(used_cores) != used_cores) {stop("The parameter 'used_cores' must be a single positive integer.")}

  set.seed(seed)

  #standardize df1
  if (scaling == TRUE) {
    df1 <- scaling_df_for_pye (df = df1, X = colnames(df1[, names(df1) %in% c(X, C)]), y = "y")$df_scaled
  }

  #check if df is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1[[y]][df1[[y]] == 1]) < 2 * n_folds)) {stop("df contains too few 1s for this number of folds")
  } else if ((length(df1[[y]][df1[[y]] == 0]) < 2 * n_folds)) {stop("df contains too few 0s for this number of folds")}

  # Divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1, y)), c("ID", y)]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 0, ])), replace = FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort[[y]] == 1, ])), replace = FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0, fold_i_1))
  folds_i <- merge(df1[, 1:2], df_sort[, 1:3], by = 'ID', all = FALSE, sort = FALSE)[, 4]

  # Names of the columns
  lambdanames <- paste("lambda", lambda, sep = "=")
  taunames <- paste("tau", tau, sep = "=")
  foldnames <- paste("fold", 1:n_folds, sep = "=")
  # Accuracy measures (train and test)
  list_of_measures <- c("pye_KS_L12", "pye_KS_L1", "pye_KS_EN", "pye_KS_SCAD", "pye_KS_MCP",
                        "auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                        "fdr", "mcc", "corrclass")
  pye_KS_L12 <- pye_KS_L1 <- pye_KS_EN <- pye_KS_SCAD <- pye_KS_MCP <- NULL
  auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures) {
    #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(foldnames, taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(c(1:n_folds), taunames)))
    assign(mes, lapply(lambda, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)), test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))))
    eval(substitute(names(x) <- lambdanames, list(x = as.symbol(mes))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames = list(foldnames, lambdanames))
  n_gammas <- lapply(lambda, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
  names(n_gammas) <- lambdanames

  if (trace > 0) {
    cat("Starting CV with the following estimation/convergence parameters: \n")
    cat("max_iter:", max_iter, "\n")
    cat("min_alpha:", min_alpha, "\n")
    cat("convergence_error:", convergence_error, "\n")
    cat("stepsizeShrink:", stepsizeShrink, "\n")
    cat("c_function_of_covariates:", c_function_of_covariates, "\n")
  }

  # fill the matrices
  results <- mapply(function(k) pye_KS.cv(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, C = C,
                                          lambda = lambda, tau = tau,
																					w = w, w_g = w_g,
																					trace = trace,
                                          alpha = alpha, a1 = a1, a2 = a2,
                                          penalty = penalty, alpha_g = alpha_g,
                                          penalty_g = penalty_g,
                                          folds_i, k, regressors_betas,
                                          pye_KS_L12 = pye_KS_L12, pye_KS_L1 = pye_KS_L1,
                                          pye_KS_EN = pye_KS_EN, pye_KS_SCAD = pye_KS_SCAD,
                                          pye_KS_MCP = pye_KS_MCP,
                                          auc = auc, aauc = aauc, aYI = aYI,
                                          youden_index = youden_index,
                                          sensitivity = sensitivity,
                                          specificity = specificity,
                                          geometric_mean = geometric_mean, fdr = fdr,
                                          mcc = mcc, corrclass = corrclass,
                                          n_betas = n_betas, n_gammas = n_gammas,
                                          max_iter = max_iter, min_alpha = min_alpha,
                                          convergence_error = convergence_error,
                                          stepsizeShrink = stepsizeShrink,
                                          used_cores = used_cores, kernel = kernel,
                                          trend = trend, delta = delta,
                                          max_alpha = max_alpha,
                                          beta_start_input = beta_start_input,
                                          beta_start_default = beta_start_default,
                                          c_zero_fixed = c_zero_fixed,
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
                                          run_aauc = run_aauc, log_file = log_file), #, long_suffix = long_suffix),
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

  # Prepare the result
  temp_pye <- get(paste0("pye_KS_", penalty))
  temp_pye[] <- lapply(seq_along(lambda), function(i) list(train = wrapper(results, paste0("pye_KS_", penalty), i, "train", tau) ,
                                                           test = wrapper(results, paste0("pye_KS_", penalty), i, "test", tau)))
  assign(paste0("pye_KS_", penalty) , temp_pye)

  for (mes in list_of_measures) {
    assign(mes[], lapply(seq_along(lambda), function(i) list(train = wrapper(results, mes, i, "train", tau),
                                                             test = wrapper(results, mes, i, "test", tau))))
    eval(substitute(names(x) <- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep = "="))), list(x = as.symbol(mes))))
  }

  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ])))
  n_gammas[] <- lapply(seq_along(lambda), function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ]))))

  betas <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas)))
	c_pye <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@c_pye)))
  gammas <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

  if (length(lambda) == 1) {
    n_betas[] <- t(t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ]))))
    n_gammas[] <- lapply(seq_along(lambda), function(i) t(t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ])))))

    betas <- t(t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas))))
		c_pye <- t(t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@c_pye))))
    gammas <- t(t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas))))
  }

  rownames(betas) <- rownames(gammas) <- rownames(c_pye) <- foldnames

  measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye")
  list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity",
                        "geometric_mean", "fdr", "mcc", "corrclass", paste0("pye_KS_", penalty))

  #t(sapply(c(1:n_folds), function(k) sapply(seq_along(lambda), function(i) results[[k]]@auc[[i]]$train[k])))

  lambda_hat_yi <- lambda_hat_auc <- lambda_hat_aauc <- lambda_hat_aYI <- lambda_hat_ccr <- NULL
  lambda_hat_sen <- lambda_hat_spc <- lambda_hat_gm <- lambda_hat_pye <- NULL
  tau_hat_yi <- tau_hat_auc <- tau_hat_aauc <- tau_hat_aYI <- tau_hat_ccr <- NULL
  tau_hat_sen <- tau_hat_spc <- tau_hat_gm <- tau_hat_pye <- NULL
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

    # now we execute pye_KS.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
    lambda_star <- get(paste0("lambda_hat_", measure_to_select_lambda))
    if (is.na(lambda_star)) {stop("Optimal lambda for sequential search was NA. Skipping tau search.")}

    # Re-initialize result structures for the second run
    taunames <- paste("tau", tau, sep = "=")
    lambdanames_star <- paste("lambda", lambda_star, sep = "=")
    # Accuracy measures (train and test)
    list_of_measures <- c("pye_KS_L12", "pye_KS_L1", "pye_KS_EN", "pye_KS_SCAD", "pye_KS_MCP",
                          "auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                          "fdr", "mcc", "corrclass")
    pye_KS_L12 <- pye_KS_L1 <- pye_KS_EN <- pye_KS_SCAD <- pye_KS_MCP <- NULL
    auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL

    for (mes in list_of_measures) {
      #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(foldnames, taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(c(1:n_folds), taunames)))
      assign(mes, lapply(lambda_star, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)), test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(c(1:n_folds), taunames)))))
      eval(substitute(names(x) <- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep = "="))), list(x = as.symbol(mes))))
    }

    n_betas_star <- matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames = list(foldnames, lambdanames_star))
    n_gammas <- lapply(lambda_star, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames = list(foldnames, taunames)))
    names(n_gammas) <- lambdanames_star

    if (trace > 0) {
      cat(" \n Starting the CV of tau using as lambda:", lambda_star, ", that is the best value of lambda as per:", measure_to_select_lambda, "\n")
    }

    # Re-run CV with optimal lambda and full tau grid
    results <- mapply(function(k) pye_KS.cv(df = df1[, names(df1) != "ID", drop = FALSE], X = X, y = y, C = C,
		                                        lambda = lambda_star, tau = tau,
		                                        w = w, w_g = w_g,
																						trace = trace,
                                            alpha = alpha, a1 = a1, a2 = a2,
                                            penalty = penalty, alpha_g = alpha_g,
                                            penalty_g = penalty_g,
                                            folds_i = folds_i, k = k,
                                            regressors_betas = regressors_betas,
                                            pye_KS_L12 = pye_KS_L12, pye_KS_L1 = pye_KS_L1,
                                            pye_KS_EN = pye_KS_EN, pye_KS_SCAD = pye_KS_SCAD,
                                            pye_KS_MCP = pye_KS_MCP,
                                            auc = auc, aauc = aauc, aYI = aYI,
                                            youden_index = youden_index,
                                            sensitivity = sensitivity,
                                            specificity = specificity,
                                            geometric_mean = geometric_mean, fdr = fdr,
                                            mcc = mcc, corrclass = corrclass,
                                            n_betas = n_betas_star, n_gammas = n_gammas,
                                            max_iter = max_iter, min_alpha = min_alpha,
                                            convergence_error = convergence_error,
                                            stepsizeShrink = stepsizeShrink,
                                            used_cores = used_cores, kernel = kernel, trend = trend,
                                            delta = delta, max_alpha = max_alpha,
                                            beta_start_input = beta_start_input,
                                            beta_start_default = beta_start_default,
                                            c_zero_fixed = c_zero_fixed,
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
                                            run_aauc = run_aauc, log_file = log_file), #, long_suffix = long_suffix),
                                            seq(1:n_folds))

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

    # Prepare the result
    temp_pye <- get(paste0("pye_KS_", penalty))
    temp_pye[] <- lapply(1, function(i) list(train = wrapper(results, paste0("pye_KS_", penalty), i, "train", tau) , test = wrapper(results, paste0("pye_KS_", penalty), i, "test", tau)))
    #temp_pye[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]], paste0("pye_KS_", penalty))$train[k, ]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]],paste0("pye_KS_", penalty))$test[k, ]))))
    assign(paste0("pye_KS_", penalty) , temp_pye)

    for (mes in list_of_measures) {
      assign(mes, lapply(1, function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
      eval(substitute(names(x) <- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep = "="))), list(x = as.symbol(mes))))
    }

    n_betas_star[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k, ])))
    n_gammas[] <- lapply(1, function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k, ]))))

    betas_star <- matrix(sapply(c(1:n_folds), function(k) results[[k]]@betas), byrow = TRUE, dimnames = list(foldnames, lambdanames_star))
    rownames(betas_star) <- foldnames
    gammas <- lapply(1, function(i) { x <- t(sapply(c(1:n_folds), function(k) results[[k]]@gammas[[i]])); rownames(x) <- foldnames; x})
    names(gammas) <- lambdanames_star

    measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye")
    list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity",
                           "geometric_mean", "fdr", "mcc", "corrclass", paste0("pye_KS_", penalty))

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
            assign(paste0("tau_hat", measures[i]), tau)
          }
        }
      }
    }
  }

  cv_time <- difftime(Sys.time(), start_time, units = "mins")

  if (trace %in% c(1, 2)) {
    cat("----------------------> END OF THE CROSS-VALIDATION OF THE PYE METHOD <----------------- \n")
    cat("-------------> For the whoole Cross-validation it took:", cv_time, "minutes <------------ \n")
  }

  return(list(penalty = penalty,
              penalty_g = penalty_g,
							w = w, w_g = w_g,
              kernel = kernel,
              cv_time = cv_time,
              pye_KS_L12 = pye_KS_L12,
              pye_KS_L1 = pye_KS_L1,
              pye_KS_EN = pye_KS_EN,
              pye_KS_SCAD = pye_KS_SCAD,
              pye_KS_MCP = pye_KS_MCP,
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
              lambda_hat_pye = lambda_hat_pye,
              tau_hat_yi = tau_hat_yi,
              tau_hat_auc = tau_hat_auc,
              tau_hat_aauc = tau_hat_aauc,
              tau_hat_aYI = tau_hat_aYI,
              tau_hat_ccr = tau_hat_ccr,
              tau_hat_sen = tau_hat_sen,
              tau_hat_spc = tau_hat_spc,
              tau_hat_gm = tau_hat_gm,
              tau_hat_pye = tau_hat_pye,
              c_function_of_covariates = c_function_of_covariates,
              simultaneous = simultaneous,
              measure_to_select_lambda = measure_to_select_lambda,
              lambda_star = lambda_star,
              n_betas = n_betas,
              n_betas_star = n_betas_star,
              n_gammas = n_gammas,
              betas = betas,
							c_pye = c_pye,
              betas_star = betas_star,
              gammas = gammas))
}
