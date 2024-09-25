#' @title pye_KS
#'
#' @description The Penalizad Youden Index (PYE) function based on the Kernel
#' Smooth density estimator. It not only performs value of PYE but it returns
#' all the necessary for the estimation process, like measure of fit and
#' derivatives. It works for all the considered penalties (L12, L1,
#' EN, SCAD and MCP)
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param betas the coefficients of the biomarker combination used for the
#' evaluation of PYE
#' @param lambda the penalization parameter of the regressors X
#' @param c the cut-off point
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param h_exponent parameter of the bandwidth of the Kernel Smooth
#' @param use_opt_c if TRUE the function uses the estimate of the best value of
#' c with the given beta, FALSE otherwise
#' @param prediction if TRUE the empirical maximum Youden index is returned as
#' Youden index performance measure
#' @param print.CDF.plot if TRUE it prints also the plot of the CDFs of cases
#' and controls for the compination of regressors Z. Default is FALSE
#'
#' @return a list containing the value of PYE for the given penalty, the value
#' of the main accuracy measure and the gradient of PYE computer the the given
#' point.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' X <- simMicroarrayData_cov02_dim50_covariates$X
#' y <- simMicroarrayData_cov02_dim50_covariates$y
#' C <- simMicroarrayData_cov02_dim50_covariates$C
#' penalty <- "L12"
#' lambda <- 0.1
#' betas <- rep(1, length(X))
#' c <- 0
#'
#' PYE_result <- pye_KS(df=df[,names(df) %in% c(X,y)], X=X, y=y, betas=betas,
#'   lambda=lambda, c=c, alpha=0.5, a1=3.7, a2=3, penalty=penalty)
#' print(PYE_result)
#'
#' @importFrom evmix kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom plyr join
#' @importFrom stats ecdf
#' @importFrom ggplot2 geom_line aes labs theme_minimal
#' @export
pye_KS <- function(df, X=names(df[,!(names(df) == y)]), y="y", betas, lambda, c=0,
                   kernel="gaussian", alpha=0.5, a1=3.7, a2=3,
                   penalty="L1", h_exponent=0.2, use_opt_c=FALSE,
                   prediction=FALSE, print.CDF.plot=FALSE){

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (!inherits(class(X), "character")){stop("X can only be of class character or data.frame.")}
  if (!inherits(class(y), "character")){stop("y can only be of class character or data.frame.")}
  if (length(c) != nrow(df)){
    if (length(c) == 1){
      c <- rep(c, nrow(df))
    } else {stop("c can only be of class numeric on length 1 or equal to the number rows of df")}
  }

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X)), drop=FALSE])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("pye_KS stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))){
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(betas) != length(X)){
    stop("The number of element of betas is different then the number of columns in df")
  }

  #divide control and diseased and multiply by betas
  z_y0<- as.matrix(df1[df1$y == 0, names(df1) %in% X, drop=FALSE])%*%betas #drop=FALSE permits it to remain a data frame even if it is composed by a single column
  z_y1<- as.matrix(df1[df1$y == 1, names(df1) %in% X, drop=FALSE])%*%betas

  #divide c
  c_y0 <- as.matrix(c[df1$y==0])
  c_y1 <- as.matrix(c[df1$y==1])

  #append data
  z<-as.data.frame(rbind(z_y0, z_y1))
  z['ID']=as.numeric(rownames(z))
  z <- z[order(z$ID), ]
  names(z)[names(z) == 'V1'] <- 'z_hat'
  df2<-plyr::join(df1, z, by='ID')
  #print(cbind(df2["y"],df2["z_hat"]))
  rownames(df2) <- df2$ID

  #find the optimal cut-point
  opt <- OptimalCutpoints::optimal.cutpoints(data= df2, X = "z_hat", status = "y", methods = "Youden", tag.healthy =0)
  empir_suggest_yi_c <- if(length(opt$Youden$Global$optimal.cutoff$cutoff)==1){opt$Youden$Global$optimal.cutoff$cutoff} else {opt$Youden$Global$optimal.cutoff$cutoff[round(length(opt$Youden$Global$optimal.cutoff$cutoff)/2)]}

  #try the Corr Class with the suggested c
  TP_yi_c <- sum(ifelse(z_y1 >= empir_suggest_yi_c ,1,0))
  TN_yi_c <- sum(ifelse(z_y0 < empir_suggest_yi_c ,1,0))
  FP_yi_c <- sum(ifelse(z_y0 >= empir_suggest_yi_c ,1,0))
  FN_yi_c <- sum(ifelse(z_y1 < empir_suggest_yi_c ,1,0))

  #compute the Youden index
  spec_yi_c <- TN_yi_c/(TN_yi_c+FP_yi_c)
  fnr_yi_c <- FN_yi_c/(FN_yi_c+TP_yi_c)
  yi1_yi_c <- spec_yi_c - fnr_yi_c
  corrclass_with_YI_c <- (TP_yi_c + TN_yi_c)/nrow(df1)

  #if we want to use the empirical c
  if (use_opt_c==TRUE){
    c = empir_suggest_yi_c
  }

  #kernel from the YI paper
  #healthy
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y0)!=0){
  #   h0 = kedd::h.bcv(as.numeric(c-z_y0), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
    h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*length(z_y0)^(-h_exponent)
  # }
  if(h0<0.1) {h0=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t0 <- as.numeric((c_y0-z_y0)/h0)
  cum0 <- evmix::kpz(z = t0, kernel = kernel)
  f0 <- sum(cum0)/length(z_y0)
  #f0 <- sum(pnorm(t0, 0, 1))/length(z_y0)

  #diseased
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y1)!=0){
  #   h1= kedd::h.bcv(as.numeric(z_y1), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
    h1=0.9*min(stats::sd(z_y1), stats::IQR(z_y1)/1.34)*length(z_y1)^(-h_exponent)
  # }
  if(h1<0.1) {h1=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t1 <- as.numeric((c_y1-z_y1)/h1)
  cum1 <- evmix::kpz(z = t1, kernel = kernel)
  f1 <- sum(cum1)/length(z_y1)
  #f1 <- sum(pnorm(t1, 0, 1))/length(z_y1)

  yi <- f0 - f1

  if(print.CDF.plot == TRUE){

    # Calcola le funzioni di distribuzione empirica con il kernel gaussiano
    z_y0_ord <- z_y0[order(z_y0)]
    z_y1_ord <- z_y1[order(z_y1)]
    ecdf0 <- stats::ecdf(z_y0_ord)
    ecdf1 <- stats::ecdf(z_y1_ord)

    # Crea un frame di dati per il plot
    plot_data <- data.frame(x = c(z_y0_ord, z_y1_ord),
                            y = c(ecdf0(z_y0_ord), ecdf1(z_y1_ord)),
                            group = c(rep("CDF_y0",length(z_y0_ord)), rep("CDF_y1", length(z_y1_ord))))

    # Crea il plot
    ggplot2::ggplot(plot_data, ggplot2::aes(x="x", y="y", color="group")) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Z", y = "CDF of Z for case and conntrol patients", color= "Groups") +
      ggplot2::theme_minimal()
  }

  #evaluate the gradient
  #h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*(length(z_y0)^(-h_exponent))
  #t0<-(c_y0-z_y0)/h0
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr0_betas <- apply (df1[df1$y == 0, names(df1) %in% X, drop=FALSE], 2,
                       function(x) (t(stats::dnorm(t0, 0, 1)) %*% ((-x)/h0))/length(z_y0))
    gr0_c <- sum(stats::dnorm(t0,0,1))/(length(z_y0)*h0)
  } else {
    gr0_betas <- apply (df1[df1$y == 0, names(df1) %in% X, drop=FALSE], 2,
           function(x) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% ((-x)/h0))/length(z_y0))
    gr0_c <- sum(evmix::kdz(z = as.numeric(t0), kernel = kernel))/(length(z_y0)*h0)
  }
  names(gr0_c)="c"

  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr1_betas <- apply (df1[df1$y == 1, names(df1) %in% X, drop=FALSE], 2,
                        function(x) (t(stats::dnorm(t1, 0, 1)) %*% ((-x)/h1))/length(z_y1))
    gr1_c <- sum(stats::dnorm(t1,0,1))/(length(z_y1)*h1)
  } else {
    gr1_betas <- apply (df1[df1$y == 1, names(df1) %in% X, drop=FALSE], 2,
           function(x) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% ((-x)/h1))/length(z_y1))
    gr1_c <- sum(evmix::kdz(z = as.numeric(t1), kernel = kernel))/(length(z_y1)*h1)
  }
  names(gr1_c)="c"

  #gradient
  gr_yi <- c(gr0_betas, gr0_c) - c(gr1_betas, gr1_c)

  #append data
  z<-as.data.frame(rbind(z_y0, z_y1))
  z['ID']=as.numeric(rownames(z))
  #z <- z[order(z$ID), ]
  names(z)[names(z) == 'V1'] <- 'z_hat'
  df1<-plyr::join(df1, z, by='ID')
  #df1<-merge(df1, z, by='ID', all = TRUE)
  rownames(df1) <- df1$ID

  #assign the ID to c
  c2<-cbind(ID=df1["ID"], c_hat=c)
  rownames(c2) <- c2$ID
  #append data (c)
  df1<-plyr::join(df1, c2, by='ID')
  #df1<-merge(df1, c, by='ID', all = TRUE)
  rownames(df1) <- df1$ID

  #the following few lines would be the way to estimate c separately with respect of
  #the betas, but we try to estimate all at the same time
  # if (mode == "estimation"){
  #   c_hat <- if(length(c$Youden$Global$optimal.cutoff$cutoff)==1){c$Youden$Global$optimal.cutoff$cutoff} else {c$Youden$Global$optimal.cutoff$cutoff[round(length(c$Youden$Global$optimal.cutoff$cutoff)/2)]}
  # } else if (mode != "prediction"){stop("The parameter mode is not valid!")}
  #

  df1["y_hat"] = ifelse(df1$z_hat >= df1$c_hat ,1,0)

  #confusion matrix
  TP <- sum(ifelse(z_y1 >= c_y1 ,1,0))
  TN <- sum(ifelse(z_y0 < c_y0 ,1,0))
  FP <- sum(ifelse(z_y0 >= c_y0 ,1,0))
  FN <- sum(ifelse(z_y1 < c_y1 ,1,0))

  #compute the Youden index
  spec <- TN/(TN+FP)
  fnr <- FN/(FN+TP)

  #sensitivity
  sensitivity = 1-fnr

  #geometric mean
  gm <- sqrt(spec*sensitivity)

  #compute other measures: FDR:
  fdr <- FP/(FP + TP)

  #compute other measures: FDR:
  mcc <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  #compute the correct classification:
  corrclass <- (TP + TN)/nrow(df1)

  #AUC e YI
  auc <-  opt$Youden$Global$measures.acc$AUC[1]
  yi1 <- opt$Youden$Global$optimal.criterion

  #if all the betas are zero, the measures are zeros
  if(sum(betas)==0){
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    gm <- 0
    yi1 <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
    auc <- 0
  }

  #L_(1/2) penalization
  phi_L12 <- lambda*sum(abs(betas)^(1/2))#NB:we DO NOT add the intercept!!!
  #gr_phi_L12 <- lambda*betas/(2*abs(betas)^(3/2))
  #gr_phi_L12[is.na(gr_phi_L12)] <- 0
  if(is.na(phi_L12)){
    stop('problem with the penalization phi_L12 in "pye" function')}

  #L1 penalization
  phi_L1 <- lambda*sum(abs(betas))#NB:we DO NOT add the intercept!!!
  #gr_phi_L1 <- lambda*sign(betas)
  #gr_phi_L1[is.na(gr_phi_L1)] <- 0
  if(is.na(phi_L1)){
    stop('problem with the penalization phi_L1 in "pye" function')}

  #Elastic-Net penalization

  phi_EN <- lambda*((alpha)*(sum(abs(betas))) + ((1-alpha)/2)*(sum(betas^2)))
  #gr_phi_EN <- lambda*((alpha)*sign(betas) + (1-alpha)*betas)
  #gr_phi_EN[is.na(gr_phi_EN)] <- 0
  if(is.na(phi_EN)){
    stop('problem with the penalization phi_EN in "pye" function')}

  #SCAD penalization
  phi_SCAD <- SCAD_function(betas, lambda, a=a1)
  #gr_phi_SCAD <- SCAD_derivative(betas, lambda, a=a1)
  #gr_phi_SCAD[is.na(gr_phi_SCAD)] <- 0
  if(is.na(phi_SCAD)){
    stop('problem with the penalization phi_SCAD in "pye" function')}

  #MCP penalization
  phi_MCP <- MCP_function(betas, lambda, a=a2)
  #gr_phi_MCP <- MCP_derivative(betas, lambda,a=a2)
  #gr_phi_MCP[is.na(gr_phi_MCP)] <- 0
  if(is.na(phi_MCP)){
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
  if (prediction==TRUE){
    yi=yi1
  }

  return(list("pye_L12"=pye_L12, "pye_L1"=pye_L1, "pye_EN" = pye_EN, "pye_SCAD" = pye_SCAD, "pye_MCP" = pye_MCP,
              #"gr_pye_L12"=gr_pye_L12, "gr_pye_L1"=gr_pye_L1, "gr_pye_EN"=gr_pye_EN, "gr_pye_SCAD"=gr_pye_SCAD, "gr_pye_MCP"=gr_pye_MCP,
              #"gr_phi_L12"=gr_phi_L12, "gr_phi_L1"=gr_phi_L1, "gr_phi_EN"=gr_phi_EN, "gr_phi_SCAD"=gr_phi_SCAD, "gr_phi_MCP"=gr_phi_MCP,
              "gr_yi"=gr_yi,
              "youden_index"=yi, "sensitivity"=sensitivity, "specificity"=spec,
              "geometric_mean"=gm, "fdr"=fdr, "mcc"=mcc, "auc"=auc, "corrclass"=corrclass,
              "empir_suggest_yi_c"=empir_suggest_yi_c, "corrclass_with_YI_c"=corrclass_with_YI_c,
              "yi1_yi_c"=yi1_yi_c, "c_hat"=c, "z_hat"=df1[,c("ID","z_hat")],
              "y_hat"=df1[,c("ID","y_hat")], "TP"=TP, "TN"=TN, "FP"=FP, "FN"=FN,
              "input_data" = list("beta"=betas, "lambda"=lambda, "alpha"=alpha,
                                  "a1"=a1, "a2"=a2, "prediction"=prediction, "c_hat"=c, "kernel"=kernel)
              ))
}





#' @title pye_KS_estimation
#'
#' @description function to estimate the optimal value of betas and c maximizing
#' the PYE function. To find the optimum the mmAPG and mnmAPG algorithms are
#' used.
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param lambda the penalization parameter of the regressors X
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
#' @param max.print number of elements to show if printing the results. Default
#' is 10
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param fold number of the fold if the function is called in cross-validation.
#' Just for visualization purposes. Default is NULL
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 10000
#' @param trend if "monotone", mmAPG is used, if "nonmonotone", mnmAPG is used.
#' Default is "monotone"
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param seed fix the seed. Default is 1
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param c_zero_fixed if TRUE the estimation process considers c, the cut-off
#' point, as fixed and equal to zero, to reduce the complexity of the estimation.
#' If FALSE c can vary and be different from zero, estimated by pye. Default
#' is FALSE
#' @param long_suffix is the standard prefix to identify the longitudinal
#' variables. For a proper management of the longitudinal variables
#' (or coefficient referring to an original longitudinal dimension),
#' long_suffix needs to be followed by a number. Default is NULL since we
#' presume that the starting df does not contain any longitudinal variable.
#'
#' @return a list containing the optimal value of betas and c, the value
#' of the main accuracy measure.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' X <- simMicroarrayData_cov02_dim50_covariates$X
#' y <- simMicroarrayData_cov02_dim50_covariates$y
#' C <- simMicroarrayData_cov02_dim50_covariates$C
#' regressors_betas<-simMicroarrayData_cov02_dim50_covariates$nregressors
#' penalty <- "SCAD"
#' lambda <- 0.1
#' trend <- "monotone" #or "nonmonotone"
#' beta_start_default <- "zeros"
#' alpha <- 0.5
#' c_zero_fixed <- FALSE
#'
#' PYE_estimation_result <- pye_KS_estimation(df=df, X=X, y=y, penalty=penalty,
#'   trend = trend, trace=2, beta_start_default=beta_start_default,
#'   beta_start_input=NULL, lambda=lambda, alpha=alpha, a1=3.7, a2=3,
#'   regressors_betas=regressors_betas, c_zero_fixed=c_zero_fixed, max_iter=5)
#' print(PYE_estimation_result)
#'
#' @importFrom stats setNames
#' @export
#Estimation of the parameter using pye_KS
pye_KS_estimation <- function(df, X=names(df[,!(names(df) == y)]), y="y",
                               lambda, penalty="L1", beta_start_input=NULL,
                               beta_start_default="zeros", max.print=10,
                               alpha=0.5, a1=3.7, a2=3, regressors_betas=NULL, fold=NULL, max_iter=10000,
                               trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                               min_alpha=1e-10, convergence_error=1e-7,
                               trace=1, seed=1, c_zero_fixed = FALSE, kernel = "gaussian",
                               zeros_stay_zeros_from_iteration=20, long_suffix=NULL){

  start_time <- Sys.time()

  options(max.print = max.print)

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors_betas to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (!inherits(class(X), "character")){stop("X can only be of class character or data.frame.")}
  if (!inherits(class(y), "character")){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X)), drop=FALSE])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("pye_KS_estimation stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}

  betas1_initial_zeros <- rep(0, length(X))
  betas1_initial_corr <- NULL

  #define c
  c<-0
  names(c)<-"c"

  #initializing the parameters (c always starts from zero)
  if (length(beta_start_input)==0){
    if (beta_start_default=="zeros"){

      betas_start <- betas1_initial_zeros
      names(betas_start) <- X
      betas1 <- betas_start

    } else if (beta_start_default=="corr") {

      #compute the corr between every x and y
      corr.xy <- data.frame(matrix(ncol = length(X), nrow = 0))
      names(corr.xy) <- X
      corr.xy[1:3,] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(lapply(X, function (x) stats::cor(df1$y, df1[,x],  method = h)))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      betas1_initial_corr <- as.numeric(mean.corr.xy)
      betas_start <- betas1_initial_corr
      names(betas_start) <- X
      betas1 <- betas_start

    } else {stop("The parameter beta_start_default can only be equal to 'zeros', 'corr'.")}
  } else {
    if (length(beta_start_input) == length(X)){
      betas_start <- beta_start_input
      names(betas_start) <- X
      betas1 <- betas_start
    } else if (length(beta_start_input) != 0){
      stop("The length of the parameter beta_start_input is not equal to ", length(X) , "\n")
    }
  }
  #add c to the betas
  betas1 <- c(betas1[!(names(betas1) %in% names(c))], c) #merge betas1 and c

  if (penalty == "SCAD") {a=a1} else {a=a2}
  prox_penalty = get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  #wrappers
  delta_fx <- function(x){
    if (c_zero_fixed==TRUE){
      x[names(x) == "c"]<-0
    }
    result <- -pye_KS(df=df1[,names(df1)!="ID", drop=FALSE], X=X[X %in% names(x)], y=y, betas=x[!(names(x) == "c")], lambda=lambda,
                      c=x[(names(x) == "c")], kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty)$gr_yi
    if (c_zero_fixed==TRUE){
      result[names(result) == "c"] <-0
    }
    return(result)
  }

  proxx <- function(x, eta){
    if (c_zero_fixed==TRUE){
      x[names(x) == "c"]<-0
    }
    if (length(long_suffix)==0){
      result <- c(prox_penalty(betas=x[!(names(x) == "c")], lambda=eta*lambda, alpha=alpha, a=a), x[(names(x) == "c")])
    } else {
      #we have to apply the penalty to all the betas if the variables referring to the same longitudinal biomarkers
      #to allow the "all or nothing" criteria to improve interpretability of the result
      #identify which variables are from longitudinal data
      longitudinal_vars <- unique(sub(paste0(long_suffix, ".*"), "", grep(paste0(long_suffix, ".*"), names(x), value = TRUE)))
      vars_dedup <- unique(sub(paste0(long_suffix, ".*"), "", names(x)))
      dgrfree <- as.numeric(sub(paste0(".*", long_suffix, "(\\d+)"), "\\1", grep(paste0(long_suffix, ".*"), names(x), value = TRUE)))
      dgrfree <- max(dgrfree[!is.na(dgrfree)])

      result <- matrix(NA, nrow=length(vars_dedup), ncol=dgrfree+1)
      rownames(result) <- c(vars_dedup)
      colnames(result) <- paste0("t",0:dgrfree)

      ordered_x <- x[order(names(x))]
      # Fill the result matrix
      for(row in 1:length(rownames(result))){
        if(length(grep(paste0("^",rownames(result)[row], long_suffix, ".*"), names(ordered_x)))>0){
          vec_index <- grep(paste0("^",rownames(result)[row], long_suffix, ".*"), names(ordered_x))
        } else {
          vec_index <- grep(paste0("^",rownames(result)[row],"$"), names(ordered_x))
        }
        result[row, 1:length(vec_index)] <- ordered_x[vec_index]
      }

      x2 <- rowMeans(abs(result), na.rm = TRUE)
      temp <- c(prox_penalty(betas=x2[-length(x2)], lambda=eta*lambda, alpha=alpha, a=a), x2[length(x2)])
      #result <- c(prox_penalty(betas=x[-c_pos], lambda=eta*lambda, alpha=alpha, a=a), x[c_pos])
      #result[rownames(result) %in% names(temp[temp==0]),] <- result[rownames(result) %in% names(temp[temp==0]),]*0
      result2 <- (temp/x2)*result
      result2[is.nan(result2)] <- 0
      vec <- as.vector(result2)
      names(vec) <- rep(rownames(result2), times=ncol(result2))
      long_names <- expand.grid(unique(names(vec[names(vec) %in% longitudinal_vars])), paste0(long_suffix, 0:dgrfree))
      names(vec)[names(vec) %in% longitudinal_vars] <- paste(long_names$Var1, long_names$Var2, sep = "")
      vec <- vec[!is.na(vec)]
      vec <- vec[match(names(x), names(vec))] #reorder as of x
      result <- vec
    }
    return(result)
  }

  Fx <- function(x){
    if (c_zero_fixed==TRUE){
      x[(names(x) == "c")]<-0
    }
    result <- -getElement(pye_KS(df=df1[,names(df1)!="ID", drop=FALSE], X=X[X %in% names(x)], y=y, betas=x[!(names(x) == "c")],
                                 lambda=lambda, c=x[(names(x) == "c")], kernel=kernel,
                                 alpha=alpha, a1=a1, a2=a2, penalty=penalty), paste0("pye_", penalty))
    return(result)
  }


  if (trend == "monotone"){
    estim <- mmAPG(x0=betas1, c_pos=length(betas1), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda, penalty=penalty,
                                  fold=fold, stepsizeShrink=stepsizeShrink, delta=delta, max_alpha= max_alpha,
                                  max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                  trace=trace, seed=seed, max.print=max.print,
                                  zeros_stay_zeros_from_iteration=zeros_stay_zeros_from_iteration)
  } else if (trend == "nonmonotone"){
    estim <- mnmAPG(x0=betas1, c_pos=length(betas1), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda, penalty=penalty,
                                      fold=fold, stepsizeShrink=stepsizeShrink, delta=delta, max_alpha=max_alpha,
                                      max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                      trace=trace, seed=seed, max.print=max.print,
                                      zeros_stay_zeros_from_iteration=zeros_stay_zeros_from_iteration)
  }

  #divide betas1 and c
  betas_hat <- estim$x1[!(names(estim$x1) == "c")]
  c_hat <- estim$x1[(names(estim$x1) == "c")]


  #lambda == 1234 is a secret code to test the performance of the correlation alone without optimization
  if(lambda == 1234){
    #compute the corr between every x and y
    corr.xy <- data.frame(matrix(ncol = length(X), nrow = 0))
    names(corr.xy) <- X
    corr.xy[1:3,] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(lapply(X, function (x) stats::cor(df1$y, df1[,x],  method = h)))))
    row.names(corr.xy) <- c("pearson", "kendall", "spearman")
    #mean of the corrs and sort names
    mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
    betas1_initial_corr <- as.numeric(mean.corr.xy)
    betas_start <- betas1_initial_corr
    names(betas_start) <- X

    betas_hat <- betas_start
    c_hat <- c
  }

  #compute z_hat
  pye_KS_value <- pye_KS(df=df1[,names(df1)!="ID"], X=X, y=y, betas=betas_hat, lambda=lambda,
                         c=c_hat, kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty)

  z_hat <- pye_KS_value$z_hat$z_hat
  df_hat <- cbind(df1, z_hat)

  youden_index <- pye_KS_value$youden_index
  sensitivity <- pye_KS_value$sensitivity
  specificity <- pye_KS_value$specificity
  geometric_mean <- pye_KS_value$geometric_mean
  fdr <- pye_KS_value$fdr
  mcc <- pye_KS_value$mcc
  auc <- pye_KS_value$auc
  corrclass <- pye_KS_value$corrclass

  if (trace %in% c(1,2)){
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend, "version. ; \n")
    if (length(fold)!=0) {cat("fold =", fold, "; \n")}
    visualize_betas <- c(betas_hat[which(betas_hat!=0)], c_hat)
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters , "; lambda:", lambda, "; penalty:", penalty, "; pye_KS:", getElement(pye_KS_value,  paste0("pye_",penalty)), "; youden_index:", pye_KS_value$youden_index, "; sensitivity:", pye_KS_value$sensitivity, "; specificity:", pye_KS_value$specificity, "; geometric_mean:", pye_KS_value$geometric_mean, "; fdr:", pye_KS_value$fdr, "; mcc:", pye_KS_value$mcc, "; auc:", pye_KS_value$auc, "; corrclass:", pye_KS_value$corrclass, "; \n")
    cat("TP:", pye_KS_value$TP, "; TN:", pye_KS_value$TN, "; FP:", pye_KS_value$FP, "; FN:", pye_KS_value$FN, ";  betas: \n")
    print(visualize_betas)
    cat("\n")
  }

  #total number of variables
  n_total_var= length(betas_hat)

  #n of predicted zeros
  n_predicted_zeros= sum(betas_hat==0)

  #n of predicted betas different from 0
  n_predicted_non_zeros = sum(betas_hat!=0)

  #compute other measures
  if (length(regressors_betas)!=0){

    #n of beta catch
    n_caught_betas = sum(betas_hat[regressors_betas!=0]!=0)

    #n of beta not catch
    n_non_caught_betas = sum(betas_hat[regressors_betas!=0]==0)

    #n of zeros catch
    n_caught_zero = sum(betas_hat[regressors_betas==0]==0)

    #n of zeros not catch
    n_zero_not_caught = sum(betas_hat[regressors_betas==0]!=0)
  } else {
    n_caught_betas=NA
    n_non_caught_betas=NA
    n_caught_zero=NA
    n_zero_not_caught=NA
  }

  name1 = paste0("betas_hat_", penalty)
  name2 = paste0("pye_KS_", penalty)

  results <-list("t1" = betas_hat,
                 "t2" = getElement(pye_KS_value, paste0("pye_", penalty)),
                 gr_yi = pye_KS_value$gr_yi,
                 lambda=lambda,
                 penalty = penalty,
                 betas_start = betas_start,
                 kernel=kernel,
                 c_hat=c_hat,
                 z_hat=pye_KS_value$z_hat,
                 y_hat=pye_KS_value$y_hat,
                 youden_index= youden_index,
                 sensitivity=sensitivity,
                 specificity=specificity,
                 geometric_mean=geometric_mean,
                 fdr=fdr,
                 mcc=mcc,
                 auc=auc,
                 corrclass=corrclass,
                 n_betas=length(betas_hat[which(betas_hat!=0)]),
                 n_total_var=n_total_var,
                 n_predicted_zeros=n_predicted_zeros,
                 n_predicted_non_zeros=n_predicted_non_zeros,
                 n_caught_betas=n_caught_betas,
                 n_non_caught_betas=n_non_caught_betas,
                 n_caught_zero=n_caught_zero,
                 n_zero_not_caught=n_zero_not_caught,
                 input_parameters=c(X=X, y=y, alpha=alpha)
                 )

  name_all = names(results)
  name_all[c(1,2)] = c(name1, name2)
  results <- stats::setNames(results, name_all)

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
  }

  results = c(results , estimation_time = estimation_time, niter=estim$tot_iters)

  return(results)
}










#---------------> cross-validation of pye KS <----------------

#create the output class of the pye.cv function
setClass(Class="pye_KS_CV_output",
         representation(
           penalty="character",
           penalty_covariates="character",
           kernel="character",
           pye_KS_L12= "ANY",
           pye_KS_L1= "ANY",
           pye_KS_EN= "ANY",
           pye_KS_SCAD= "ANY",
           pye_KS_MCP= "ANY",
           auc= "ANY",
           aauc= "ANY",
           aYI= "ANY",
           youden_index="ANY",
           sensitivity="ANY",
           specificity="ANY",
           geometric_mean="ANY",
           fdr="ANY",
           mcc="ANY",
           corrclass="ANY",
           n_betas="ANY",
           n_gammas="ANY",
           betas="list",
           gammas="list"
         )
)

#to extract just part of the estim of pye
subset_pye_KS <- function(df, X, y, betas, lambda, c, fold, alpha, trace, a1, a2, penalty,
                           cv_time, niter, kernel, c_function_of_covariates=FALSE){

  ID <- rownames(df)
  pye_KS_result <- pye_KS(df=df[,names(df)!="ID"], X=X, y=y, betas=betas, lambda=lambda, c=c,
                          alpha=alpha, a1=a1, a2=a2, penalty=penalty, prediction=TRUE, kernel=kernel)

  #z_hat <- pye_KS_result$z_hat$z_hat

  #optimal cutpoint (the one estimated in the training)
  youden_index <- pye_KS_result$youden_index
  sensitivity <- pye_KS_result$sensitivity
  specificity <- pye_KS_result$specificity
  geometric_mean <- pye_KS_result$geometric_mean
  fdr <- pye_KS_result$fdr
  mcc <- pye_KS_result$mcc
  auc <- pye_KS_result$auc
  aauc <- NULL
  aYI <- NULL
  corrclass <- pye_KS_result$corrclass
  TP <- pye_KS_result$TP
  TN <- pye_KS_result$TN
  FP <- pye_KS_result$FP
  FN <- pye_KS_result$FN
  c_hat <- c

  if (trace %in% c(1,2)) {
    #print the results only if c_function_of_covariates = FALSE
    if(c_function_of_covariates==FALSE){
      cat("-> Results on the TEST SET of PYE\n")
      cat("-> algorithm: pye_KS_proximal_gradient_method ; ")
      if (length(fold)!=0) {cat("fold =", fold, "; \n")}
      visualize_betas <- c(betas[which(betas!=0)], c_hat)
      cat("lambda:", lambda, "; penalty:", penalty, "; pye_KS:", getElement(pye_KS_result,  paste0("pye_",penalty)), "; youden_index:", pye_KS_result$youden_index, "; sensitivity:", pye_KS_result$sensitivity, "; specificity:", pye_KS_result$specificity, "; geometric_mean:", pye_KS_result$geometric_mean, "; fdr:", pye_KS_result$fdr, "; mcc:", pye_KS_result$mcc, "; auc:", pye_KS_result$auc, "; corrclass:", pye_KS_result$corrclass, " \n")
      cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, ";  betas: \n")
      print(visualize_betas)
      cat("Cross-validation time:", cv_time, "; Number of iterations:", niter, "\n")
      cat("\n")
    }
  }

  return(list(pye_KS_L12=pye_KS_result$pye_L12, pye_KS_L1=pye_KS_result$pye_L1, pye_KS_EN=pye_KS_result$pye_EN,
              pye_KS_SCAD= pye_KS_result$pye_SCAD, pye_KS_MCP=pye_KS_result$pye_MCP,
              youden_index=pye_KS_result$youden_index, sensitivity=pye_KS_result$sensitivity,
              specificity=pye_KS_result$specificity, geometric_mean=pye_KS_result$geometric_mean,
              fdr=pye_KS_result$fdr, mcc=pye_KS_result$mcc, auc=pye_KS_result$auc,
              corrclass=pye_KS_result$corrclass, TP=TP, TN=TN, FP=FP, FN=FN,
              z_hat=pye_KS_result$z_hat))
}

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom stats setNames
pye_KS.cv <- function (df, X, y, C, lambda, tau, trace=1, alpha, alpha_g, a1, a2, penalty, penalty_g, folds_i, k,
                        regressors_betas=NULL, kernel_g, a1_g, a2_g,
                        trend_g, gamma_start_input, gamma_start_default,
                        regressors_gammas=NULL, max_iter_g, delta_g, max_alpha_g, stepsizeShrink_g,
                        min_alpha_g, convergence_error_g,
                        pye_KS_L12, pye_KS_L1, pye_KS_EN, pye_KS_SCAD, pye_KS_MCP,
                        auc, aauc, aYI, youden_index, sensitivity, specificity, geometric_mean, fdr,
                        mcc, corrclass, n_betas, n_gammas, used_cores, trend, delta, max_alpha, kernel,
                        beta_start_input, beta_start_default, c_zero_fixed, max_iter, min_alpha, convergence_error,
                        stepsizeShrink, c_function_of_covariates, simultaneous, run_aauc, long_suffix){

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  if(c_function_of_covariates==TRUE){greek<-"tau  "} else {greek<-"lambda"}
  if(trace %in% c(1,2)){
    cat(" ---------------------------------------------------------\n |  starting with the", k ,"-th fold for the CV of ", greek,"  |\n ---------------------------------------------------------\n")
  }

  if((c_function_of_covariates == TRUE) & (length(lambda) == 1)) {
    #in this case we have to estimate the combination z_hat over all the dataset
    train_df1 <- df
  } else {
    train_df1 <- train_df
  }

  #start cores
  #if no cores are defined lets use the 70% of the available cores
  if (length(used_cores)==0){
    if (parallel::detectCores()>1){used_cores=round(parallel::detectCores()*0.7,0)} else {used_cores=1}
  }

  if (used_cores > 1) {
    max.cores <- parallel::detectCores()
    if (used_cores > max.cores) {
      #stop("The number of cores specified (", used_cores,
      warning("The number of cores specified (", used_cores,
           ") is larger than the number of available cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_pye_ks_models.txt", setup_strategy = "sequential")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for cross-validation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("train_df", "X", "y", "C", "trace", "alpha", "alpha_g", "penalty", "penalty_g",
                                  "a1_g", "a2_g", "trend_g",
                                  "regressors_betas", "k", "trend", "delta", "c_function_of_covariates",
                                  "max_alpha", "kernel", "beta_start_input", "beta_start_default",
                                  "c_zero_fixed", "a2", "a1", "long_suffix",
                                  "max_iter", "min_alpha", "convergence_error", "stepsizeShrink"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))
    fitted_models <- parallel::parLapply(cl, lambda, function(x) pye_KS_estimation(df=train_df1, X=X, y=y, lambda=x,
                                                                                    beta_start_input=beta_start_input,
                                                                                    beta_start_default=beta_start_default,
                                                                                    trace=trace,
                                                                                    alpha=alpha,
                                                                                    a1=a1, a2=a2,
                                                                                    penalty=penalty, max_iter=max_iter,
                                                                                    convergence_error=convergence_error,
                                                                                    regressors_betas=regressors_betas, fold=k,
                                                                                    trend=trend,
                                                                                    stepsizeShrink=stepsizeShrink,
                                                                                    delta=delta, max_alpha=max_alpha,
                                                                                    min_alpha=min_alpha,
                                                                                    kernel=kernel,
                                                                                    c_zero_fixed=c_zero_fixed,
                                                                                    long_suffix=long_suffix))


    on.exit(parallel::stopCluster(cl))

  } else {

    fitted_models <- lapply(lambda, function(x) pye_KS_estimation(df=train_df1, X=X, y=y,
                                                                   lambda=x,
                                                                   beta_start_input=beta_start_input,
                                                                   beta_start_default=beta_start_default, trace=trace,
                                                                   alpha=alpha,
                                                                   a1=a1, a2=a2, penalty=penalty, max_iter=max_iter,
                                                                   min_alpha=min_alpha,
                                                                   convergence_error=convergence_error,
                                                                   regressors_betas=regressors_betas, fold=k,
                                                                   trend=trend,
                                                                   stepsizeShrink=stepsizeShrink,
                                                                   delta=delta, max_alpha=max_alpha,
                                                                   kernel=kernel,
                                                                   c_zero_fixed=c_zero_fixed,
                                                                   long_suffix=long_suffix
                                                                   ))
  }


  #z_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "z_hat"))
  if((c_function_of_covariates == TRUE) & (length(lambda) == 1)) {
    #in this case we have to estimate the combination z_hat over all the dataset
    z_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "z_hat")[-test_i, ])
  } else {
    z_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "z_hat"))
  }

  if(length(gamma_start_input) == 0){
	#if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
    gamma_start_input1 <- lapply(c(1:length(lambda)), function(x) {
														gamma_start_input1 <- c(getElement(fitted_models[[x]], "c_hat"), rep(0, length(C)))
														names(gamma_start_input1) <- c("const", C)
														gamma_start_input1
														})
  } else {
	  gamma_start_input1 <- lapply(c(1:length(lambda)), function(x) {
														gamma_start_input1 <- gamma_start_input
														names(gamma_start_input1) <- c("const", C)
														gamma_start_input1
														})
  }

  #here I have to put the covYI estimation
  #computing c
  if(c_function_of_covariates == FALSE) {

    #measures on the train set - if we don't compute c with covariates, we are only dependent of lambda:
    temp_pye = get(paste0("pye_KS_", penalty))

    for (i in 1:length(lambda)){
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

    gammas_hat <- lapply(c(1:length(lambda)), function(y) NA)
    names(gammas_hat) <- lambda

  } else {

  	#start cores
    #if no cores are defined lets use the 70% of the available cores
    if (length(used_cores)==0){
      if (parallel::detectCores()>1){used_cores=round(parallel::detectCores()*0.7,0)} else {used_cores=1}
    }

    #use covYI to compute the value of c
    if (used_cores > 1) {
      max.cores <- parallel::detectCores()
      if (used_cores > max.cores) {
        #stop("The number of cores specified (", used_cores,
        warning("The number of cores specified (", used_cores,
                ") is larger than the number of available cores (",
                max.cores, ")!")
      }
      #cl <- parallel::makeCluster(used_cores, outfile="log_pye_ks_models.txt", setup_strategy = "sequential")
      #if (!("cluster" %in% class(cl)))
      #  stop("cl is not of class 'cl'; see ?makeCluster")

      cat("Start parallel computing for cross-validation, I am using:", length(cl),"cores \n")
      parallel::clusterExport(cl, c("train_df", "tau", "y", "C", "trace", "penalty_g",
                                    "alpha_g", "a1_g", "a2_g", "trend_g", "kernel_g", "k",
                                    "gamma_start_input1", "gamma_start_default", "z_hat",
                                    "regressors_gammas", "max_iter_g", "delta_g", "max_alpha_g", "stepsizeShrink_g",
                                    "min_alpha_g", "convergence_error_g", "run_aauc"), envir = environment())
      parallel::clusterCall(cl, function() require(pye))
      if(length(lambda)>=length(tau)){
        print("The parallel loop is on lambda")
        fitted_gammas <- parallel::parLapply(cl, c(1:length(z_hat)), function (x) lapply(tau, function(t) covYI_KS_estimation(df=cbind(train_df, z_hat=z_hat[[x]]$z_hat),
                                              																														z="z_hat", y=y, C=C, tau=t,
                                              																														gamma_start_input=gamma_start_input1[[x]],
                                              																														gamma_start_default=gamma_start_default, trace=trace,
                                              																														alpha=alpha_g, a1=a1_g, a2=a2_g, penalty=penalty_g,
                                              																														max_iter=max_iter_g,
                                              																														min_alpha=min_alpha_g,
                                              																														convergence_error=convergence_error_g,
                                              																														regressors_gammas=regressors_gammas, fold=k,
                                              																														trend=trend_g,
                                              																														stepsizeShrink=stepsizeShrink_g,
                                              																														delta=delta_g, max_alpha=max_alpha_g, kernel=kernel_g,
                                              																														run_aauc=run_aauc)))
      } else {
        print("The parallel loop is on tau")
        fitted_gammas <- lapply(c(1:length(z_hat)), function (x) parallel::parLapply(cl, tau, function(t) covYI_KS_estimation(df=cbind(train_df, z_hat=z_hat[[x]]$z_hat),
                                                                                                                 z="z_hat", y=y, C=C, tau=t,
                                                                                                                 gamma_start_input=gamma_start_input1[[x]],
                                                                                                                 gamma_start_default=gamma_start_default, trace=trace,
                                                                                                                 alpha=alpha_g, a1=a1_g, a2=a2_g, penalty=penalty_g,
                                                                                                                 max_iter=max_iter_g,
                                                                                                                 min_alpha=min_alpha_g,
                                                                                                                 convergence_error=convergence_error_g,
                                                                                                                 regressors_gammas=regressors_gammas, fold=k,
                                                                                                                 trend=trend_g,
                                                                                                                 stepsizeShrink=stepsizeShrink_g,
                                                                                                                 delta=delta_g, max_alpha=max_alpha_g, kernel=kernel_g,
                                                                                                                 run_aauc=run_aauc)))
      }

      on.exit(parallel::stopCluster(cl))

    } else {

      fitted_gammas <- lapply(c(1:length(z_hat)), function (x) lapply(tau, function(t) covYI_KS_estimation(df=cbind(train_df, z_hat=z_hat[[x]]$z_hat),
																								  z="z_hat", y=y, C=C, tau=t,
																								  gamma_start_input=gamma_start_input1[[x]],
																								  gamma_start_default=gamma_start_default, trace=trace,
																								  alpha=alpha_g, a1=a1_g, a2=a2_g, penalty=penalty_g,
																								  max_iter=max_iter_g,
																								  min_alpha=min_alpha_g,
																								  convergence_error=convergence_error_g,
																								  regressors_gammas=regressors_gammas, fold=k,
																								  trend=trend_g,
																								  stepsizeShrink=stepsizeShrink_g,
																								  delta=delta_g, max_alpha=max_alpha_g, kernel=kernel_g,
																								  run_aauc=run_aauc)))
    }

    names(fitted_gammas) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
    for (zy in 1:length(lambda)){
      names(fitted_gammas[[zy]]) <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
    }

    #organize the result
    names <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
    gammas_hat <- lapply(lambda, function(x) sapply(names, function(xx) NULL))
    names(gammas_hat)<- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep="=")))

    temp_pye = get(paste0("pye_KS_", penalty))

    for (i in 1:length(lambda)){
      #measures on the train set
      #these are multiple tables based on the number of considered TAUs
      temp_pye[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("covYI_KS_", penalty_g))))
      auc[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "auc")))
      aauc[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "aauc")))
      aYI[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "aYI")))
      youden_index[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "sensitivity")))
      specificity[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "fdr")))
      mcc[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "mcc")))
      corrclass[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "corrclass")))

      n_gammas[[i]][k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "n_gammas")))
      gammas_hat[[i]] <- lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("gammas_hat_", penalty_g)))
      names(gammas_hat[[i]])<- lapply(tau, function (xx) paste("tau", xx, sep="="))
    }
  }

  #n_betas is only based on lambda, not tau!
  n_betas[k, ] <- unlist(lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "n_betas")))

  #create a list of the results (NB: betas_hat includes c_hat) - it is only based on lambda, not tau!
  betas_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], paste0("betas_hat_", penalty)))
  names(betas_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  #c_hat is the fixed value of c in case we don't use the covariates to estimate a patient's specific cut-off value
  c_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "c_hat"))
  names(c_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  betas <- betas_hat
  names(betas) <- lambda
  gammas <- gammas_hat

  #then I have to evaluate the different taus with subset_pye_KS and choose the best tau
  #based on the performance on the test set.

  #measures on the test set
  all_measures_test <- mapply(function(x, xx, z, zz) subset_pye_KS(df=test_df, X=X, y=y, betas=x,lambda=z,
                                                           c=xx, fold = k, alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                                           cv_time=fitted_models[[zz]]$estimation_time,
                                                           niter=fitted_models[[zz]]$niter, kernel=kernel, trace=trace,
                                                           c_function_of_covariates=c_function_of_covariates),
                              betas_hat, c_hat, lambda, 1:length(lambda))

  #organize the results
  if(c_function_of_covariates == FALSE) {

    #measures on the test set - if we don't compute c with covariates, we are only dependent of lambda:
    for (i in 1:length(lambda)){
      #measures on the train set
      #these are multiple tables based on the number of considered TAUs
      temp_pye[[i]]$test[k, ] <- all_measures_test[paste0("pye_KS_", penalty),][[i]]
      auc[[i]]$test[k, ]  <- all_measures_test["auc",][[i]]
      aauc[[i]]$test[k, ]  <- NA
      aYI[[i]]$test[k, ]  <- NA
      youden_index[[i]]$test[k, ] <- all_measures_test["youden_index",][[i]]
      sensitivity[[i]]$test[k, ] <- all_measures_test["sensitivity",][[i]]
      specificity[[i]]$test[k, ]  <- all_measures_test["specificity",][[i]]
      geometric_mean[[i]]$test[k, ]  <- all_measures_test["geometric_mean",][[i]]
      fdr[[i]]$test[k, ]  <- all_measures_test["fdr",][[i]]
      mcc[[i]]$test[k, ]  <- all_measures_test["mcc",][[i]]
      corrclass[[i]]$test[k, ] <- all_measures_test["corrclass",][[i]]
    }

    assign(paste0("pye_KS_", penalty), temp_pye)

  } else {

    #put "const" in C if not already there
    if(!("const" %in% C)){
      C1 <- c("const", C)
    }

    if("const" %in% names(test_df)){
      test_df1 <- cbind(test_df[y], test_df[,"const"], test_df[,(names(test_df) %in% C)])
    } else {
      const <- rep(1, nrow(test_df)) #the constant for the coefficients gammas
      test_df1 <- cbind(test_df[y], const, test_df[,(names(test_df) %in% C)])
    }

    #use covYI to compute the value of c
    listing <- lapply(lambda, function(xx) lapply(tau, function(x) NULL))
    for (l in 1:length(lambda)){for (t in 1:length(tau)){ listing[[l]][[t]] <- list(gammas_hat[[l]][[t]],tau[[t]])}}

    z_hat <- lapply(c(1:length(lambda)), function(x) getElement(all_measures_test["z_hat",][[x]], "z_hat"))

    cov_results <- lapply(1:length(lambda), function(xx) lapply(listing[[xx]], function(x) covYI_KS(df=cbind(test_df1,
                                                        z_hat=z_hat[[xx]]), z="z_hat", y=y, C=C1,
                                                        gammas=x[[1]], tau=x[[2]], kernel=kernel_g, alpha=alpha_g,
                                                        a1=a1_g, a2=a2_g, penalty=penalty_g, prediction=TRUE,
                                                        run_aauc=run_aauc)))

    #name the lambdas
    names(cov_results) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
    for (i in 1:length(lambda)){
      #name the taus
      names(cov_results[[i]]) <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))

      if (trace %in% c(1,2)){
        #print the estimation
        cat("\n")
        cat("Final results on TEST SET of covYI for the FOLD:", k, ": \n")
        cat("With lambda:", lambda[i], " and  penalty:", penalty, ". Accuracy measures only PYE: \n")
        visualize_betas <- c(betas_hat[[i]][which(betas_hat[[i]]!=0)], stats::setNames(c_hat[[i]], "c"))
        cat("pye_KS:", all_measures_test[paste0("pye_KS_", penalty),][[i]], "; youden_index:", all_measures_test["youden_index",][[i]], "; sensitivity:", all_measures_test["sensitivity",][[i]], "; specificity:", all_measures_test["specificity",][[i]], "; geometric_mean:", all_measures_test["geometric_mean",][[i]], "; fdr:", all_measures_test["fdr",][[i]], "; mcc:", all_measures_test["mcc",][[i]], "; auc:", all_measures_test["auc",][[i]], "; corrclass:", all_measures_test["corrclass",][[i]], "; \n")
        cat("TP:", all_measures_test["TP",][[i]], "; TN:", all_measures_test["TN",][[i]], "; FP:", all_measures_test["FP",][[i]], "; FN:", all_measures_test["FN",][[i]], ";  betas: \n")
        print(visualize_betas)
		    cat("Cross-validation time:", fitted_models[[i]]$estimation_time, "; Number of iterations:", fitted_models[[i]]$niter, "\n")
        cat("\n")


        for (ii in 1:length(cov_results[[i]])){
          cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend,"version. ; \n")
          visualize_gammas = c(gammas_hat[[i]][[ii]][which(gammas_hat[[i]][[ii]]!=0)])
          cat("tau:", tau[ii], "; penalty:", penalty_g, "; covYI_KS_", penalty_g, ":" , getElement(cov_results[[i]][[ii]], paste0("covYI_KS_", penalty_g)), "; youden_index:", cov_results[[i]][[ii]]$youden_index, "; aYI:", cov_results[[i]][[ii]]$aYI, "; sensitivity:", cov_results[[i]][[ii]]$sensitivity, "; specificity:", cov_results[[i]][[ii]]$specificity, "; geometric_mean:", cov_results[[i]][[ii]]$geometric_mean, "; fdr:", cov_results[[i]][[ii]]$fdr, "; mcc:", cov_results[[i]][[ii]]$mcc, "; auc:", cov_results[[i]][[ii]]$auc, "; aauc:", cov_results[[i]][[ii]]$aauc, "; corrclass:", cov_results[[i]][[ii]]$corrclass, "; \n")
          cat("TP:", cov_results[[i]][[ii]]$TP, "; TN:", cov_results[[i]][[ii]]$TN, "; FP:", cov_results[[i]][[ii]]$FP, "; FN:", cov_results[[i]][[ii]]$FN, "; gammas: \n")
          cat("\n")
          print(visualize_gammas)
		      cat("Cross-validation time:", fitted_gammas[[i]][[ii]]$estimation_time, "; Number of iterations:", fitted_gammas[[i]][[ii]]$niter, "\n")
          cat("\n")
        }
      }

      #measures on the test set
      #these are multiple tables based on the number of considered TAUs
      #NB: when c_function_of_covariates = TRUE, the measure pye_KS coincides with the covYI_KS measure
      temp_pye[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], paste0("covYI_KS_", penalty_g))))
      assign(paste0("pye_KS_", penalty) , temp_pye)
      auc[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "auc")))
      aauc[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "aauc")))
      aYI[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "aYI")))
      youden_index[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "sensitivity")))
      specificity[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "fdr")))
      mcc[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "mcc")))
      corrclass[[i]]$test[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(cov_results[[i]][[yy]], "corrclass")))
    }
  }

  return(new("pye_KS_CV_output", penalty=penalty, penalty_covariates=penalty_g, kernel=kernel,
             pye_KS_L12=pye_KS_L12, pye_KS_L1=pye_KS_L1, pye_KS_EN=pye_KS_EN, pye_KS_SCAD=pye_KS_SCAD, pye_KS_MCP=pye_KS_MCP,
             auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index, sensitivity=sensitivity,
             specificity=specificity, geometric_mean=geometric_mean, fdr=fdr, mcc=mcc, corrclass=corrclass,
             n_betas=n_betas, n_gammas=n_gammas, betas=betas, gammas=gammas))
}


#' @title pye_KS_compute_cv
#'
#' @description function to perform the cross-validation to select the best
#' value of lambda (and possibly tau) for the estimation of betas and c (and
#' possibly gammas) using pye (and possibly covYI).
#'
#' @param n_folds number of fold of the cross validation
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y and C
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param C covariate variables. It can be of type dataframe, containing the
#' same covariates included in df, or just a vector of character. Default is NULL
#' @param lambda the penalization parameter of the regressors X
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param alpha_g parameter for the Elastic-Net penalization term in covYI.
#' Default is 0.5
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param seed fix the seed. Default is 1
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param trend if "monotone", mmAPG is used, if "nonmonotone", mnmAPG is used.
#' Default is "monotone"
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 10000
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
#' @param scaling if TRUE, the dataset is scaled. FALSE otherwise, Default is
#' FALSE.
#' @param c_zero_fixed if TRUE the estimation process considers c, the cut-off
#' point, as fixed and equal to zero, to reduce the complexity of the estimation.
#' If FALSE c can vary and be different from zero, estimated by pye. Default
#' is FALSE
#' @param c_function_of_covariates if TRUE, covYI is used to estimate the
#' cut-off point as function of the convariate information. If FALSE, the
#' covariate information is ignored. Default is FALSE
#' @param simultaneous in case of c_function_of_covariates=TRUE, it defines if
#' gammas needs to be estimated simultaneously or as a second step with respect
#' of betas. Default is FALSE, i.e. not simultaneously
#' @param measure_to_select_lambda the measure used to select lambda if
#' simultaneous=FALSE, i.e. when the cross-validation proecess select fist
#' lambda and then tau. Default is "ccr", i.e. the correct classification rate
#' @param penalty_g the considered penalty of covYI. To be chosen among L12, L1,
#' EN, SCAD and MCP. Default if "L1"
#' @param kernel_g the kernel type to use for the estimation of the density
#' function (tested only for "gaussian") in covYI.  Default is "gaussian"
#' @param a1_g parameter for the SCAD and MCP penalization term in covYI. Default is 3.7
#' @param a2_g parameter for the MCP penalization term in covYI. Default is 3.0
#' @param trend_g for covYI. If "monotone", mmAPG is used, if "nonmonotone",
#' mnmAPG is used. Default is "monotone"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point of gamma.
#' If "zeros", it starts with all zero values, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable. Default
#' is "zeros"
#' @param regressors_gammas a vector containing the real gammas (if known).
#' Default is NULL
#' @param max_iter_g maximum number of iterations in the algorithms mmAPG and
#' mnmAPG in covYI. Default is 10000
#' @param delta_g parameter for the convergence condition of the optimization
#' algorithm of covYI. Default is 1e-5
#' @param max_alpha_g maximum value of the step-parameter alpha in covYI.
#' Default is 100
#' @param stepsizeShrink_g parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of covYI. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param min_alpha_g minimum value of the step-parameter alpha in covYI.
#' Default is 1e-12
#' @param convergence_error_g in covYI, it is error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
#'
#' @return a list containing the optimal value of lambda to estimate betas and
#' c, the value of the main accuracy measure for all the folds.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' X <- simMicroarrayData_cov02_dim50_covariates$X
#' y <- simMicroarrayData_cov02_dim50_covariates$y
#' C <- simMicroarrayData_cov02_dim50_covariates$C
#' regressors_betas<-simMicroarrayData_cov02_dim50_covariates$nregressors
#' regressors_gammas<-simMicroarrayData_cov02_dim50_covariates$ncovariates
#' penalty <- "SCAD"
#' trend = "monotone" #or "nonmonotone"
#' pye_starting_point = "zeros"
#' alpha = 0.5
#' c_zero_fixed <- FALSE
#' c_function_of_covariates <- TRUE
#' used_cores <- 1
#' used_penalty_pye <- c("L1", "MCP") #c("L12", "L1", "EN", "SCAD", "MCP")
#' max_iter <- 10
#' n_folds <- 3
#'
#' #pye Gaussian (and others) Kernel Smooth
#' for (p in used_penalty_pye){
#'
#'   name <- paste0("param_estimate_PYE_KS_covYI_", p)
#'
#'   #wrapper of the function
#'   wrapper_lambda <- function(lambda){
#'     return(pye_KS_estimation(df=df, X=X, penalty=p, trend=trend, trace=1,
#'       beta_start_default=pye_starting_point, beta_start_input=NULL,
#'       lambda=lambda, alpha=alpha, a1=3.7, a2=3,
#'       regressors_betas=regressors_betas, c_zero_fixed=c_zero_fixed,
#'       max_iter=max_iter))
#'   }
#'
#'   #lambda to max and min
#'   lambda_max <- calibrate_lambda_max (function_to_run=wrapper_lambda,
#'     var_to_check=paste0("betas_hat_",p), lambda_start = 5)
#'   lambda_min <- calibrate_lambda_min (function_to_run=wrapper_lambda,
#'     var_to_check=paste0("betas_hat_",p), lambda_start = 0.0001, max_var=100)
#'
#'	 cat("\n Pye penalty:", p, "; lambda_max=", lambda_max, "; lambda_min=",
#'	  lambda_min, "\n")
#'
#'   #create a suited lambda
#'   lambda <- create_lambda(n=3, lmax=lambda_max, lmin=lambda_min)
#'   lambda <- as.numeric(formatC(lambda, format = "e", digits = 9))
#'
#'   if(c_function_of_covariates == TRUE){
#'     #wrapper of the function for tau
#'     lambda_star <- (lambda_max+lambda_min)/2
#'     penalty_g <- p
#'     alpha_g <- alpha
#'     trend_g <- trend
#'     wrapper_tau <- function(tau){
#'       z_hat <- pye_KS_estimation(df=df, X=X, penalty=p, trend=trend, trace=1,
#'         beta_start_default=pye_starting_point, beta_start_input=NULL,
#'         lambda=lambda_star, alpha=alpha, a1=3.7, a2=3,
#'         regressors_betas=regressors_betas, c_zero_fixed=c_zero_fixed,
#'         max_iter=max_iter)$z_hat$z_hat
#'       return(covYI_KS_estimation(df=cbind(df, z_hat=z_hat), z="z_hat", y=y,
#'         C=C, tau=tau, gamma_start_input=NULL,
#'         gamma_start_default=pye_starting_point, trace=1, alpha=alpha,
#'         penalty=p, regressors_gammas=regressors_gammas, trend=trend,
#'         max_iter=max_iter, run_aauc=FALSE))
#'     }
#'
#'     #tau to max and min
#'     tau_max <- calibrate_lambda_max (function_to_run=wrapper_tau,
#'       var_to_check=paste0("gammas_hat_",p), lambda_start = 5, n_min_var=1)
#'     tau_min <- calibrate_lambda_min (function_to_run=wrapper_tau,
#'       var_to_check=paste0("gammas_hat_",p), lambda_start = 0.0001,
#'       max_var=(length(C)-1))
#'
#'	   cat("\n tau_max=", tau_max, "; tau_min=", tau_min, "\n")
#'
#'     #create a suited tau
#'     tau <- create_lambda(n=3, lmax=tau_max, lmin=tau_min)
#'     tau <- as.numeric(formatC(tau, format = "e", digits = 9))
#'   }
#'
#'   #start cross-validation
#'   assign(name, pye_KS_compute_cv(penalty=p, df=df,
#'     X=names(df[,!(names(df) %in% c(y,C))]), y=y, C=C, trace=1,
#'     beta_start_default=pye_starting_point, trend = trend,
#'     beta_start_input=NULL, n_folds=n_folds, lambda=lambda, tau=tau,
#'     alpha_g=alpha, alpha=alpha, a1=3.7, a2=3,
#'     regressors_betas=regressors_betas, regressors_gammas=regressors_gammas,
#'     used_cores=used_cores, kernel="gaussian", c_zero_fixed=c_zero_fixed,
#'     c_function_of_covariates=c_function_of_covariates, penalty_g=p,
#'     max_iter=max_iter, max_iter_g=max_iter))
#'
#'   #take the best lambda per measures (in case of same measure for diff lambdas,
#'   #we take the one associated to less betas)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_yi"), get(name)$lambda_hat_pye_KS_yi)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_auc"), get(name)$lambda_hat_pye_KS_auc)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_aauc"), get(name)$lambda_hat_pye_KS_aauc)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_aYI"), get(name)$lambda_hat_pye_KS_aYI)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_ccr"), get(name)$lambda_hat_pye_KS_ccr)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_gm"), get(name)$lambda_hat_pye_KS_gm)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_pye"), get(name)$lambda_hat_pye_KS_pye)
#'
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_yi"), get(name)$tau_hat_pye_KS_yi)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_auc"), get(name)$tau_hat_pye_KS_auc)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_aauc"), get(name)$tau_hat_pye_KS_aauc)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_aYI"), get(name)$tau_hat_pye_KS_aYI)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_ccr"), get(name)$tau_hat_pye_KS_ccr)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_gm"), get(name)$tau_hat_pye_KS_gm)
#'   assign(paste0("tau_hat_PYE_KS_covYI_", p ,"_pye"), get(name)$tau_hat_pye_KS_pye)
#'
#' }
#'
#'
#' @export
pye_KS_compute_cv <- function (n_folds, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C=NULL, lambda, trace=1,
                                alpha=0.5, alpha_g=0.5, tau=0,
                                a1=3.7, a2=3, penalty="L1", regressors_betas=NULL,
                                seed=1, used_cores=1, trend = "monotone", delta = 1e-5, max_alpha=10000,
                                kernel="gaussian", beta_start_input=NULL,
                                max_iter=10000, #reduced for the CV
                                min_alpha=1e-10,
                                convergence_error=1e-7,
                                stepsizeShrink=0.8,
                                beta_start_default="zeros", scaling=FALSE, c_zero_fixed=FALSE,
                                c_function_of_covariates=FALSE,
                                simultaneous=FALSE,
                                measure_to_select_lambda="ccr",
                                penalty_g="L1", kernel_g="gaussian", a1_g=3.7, a2_g=3,
                                trend_g="monotone", gamma_start_input=NULL, gamma_start_default="zeros",
                                regressors_gammas=NULL, max_iter_g=10000,
                                delta_g=1e-5, max_alpha_g=10000,
                                stepsizeShrink_g=0.8, min_alpha_g=1e-12, convergence_error_g=1e-7,
								                run_aauc=FALSE, long_suffix=NULL) {

  start_time <- Sys.time()

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors_betas to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (!inherits(class(X), "character")){stop("X can only be of class character or data.frame.")}
  if (!inherits(class(y), "character")){stop("y can only be of class character or data.frame.")}
  if(c_function_of_covariates==TRUE){
    if (!inherits(class(C), "character")){stop("C can only be of class character or data.frame.")}
  }

  tau2 <- 0

  #check if tau exists when c_function_of_covariates=TRUE
  if(c_function_of_covariates==TRUE){
    if(is.null(tau)){stop("tau cannot be NULL if c_function_of_covariates=TRUE.")}
    if(length(tau)==1){if(tau==0){stop("tau cannot be 0 if c_function_of_covariates=TRUE.")}}
    if(simultaneous == FALSE){
      c_function_of_covariates <- FALSE
      tau2 <- tau #if tau2 != 0 then c_function_of_covariates = TRUE, and we have to re-change it later
      tau <- 0
    }
  } else if(c_function_of_covariates ==FALSE) {
    tau <- 0
  } else {stop("c_function_of_covariates shoould be equal to TRUE or FALSE.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, y=df[,(names(df) %in% c(y))], df[,(names(df) %in% c(X,C))])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}
  if (!(measure_to_select_lambda %in% c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye"))){stop("A wrong value has been assigned to the parameter measure_to_select_lambda")}

  set.seed(seed)

  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye (df=df1, X=colnames(df1[, names(df1) %in% c(X,C)]), y="y")$df_scaled
  }

  #check if df is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1$y[df1$y==1])<2*n_folds)){stop("df contains too few 1s for this number of folds")
  } else if ((length(df1$y[df1$y==0])<2*n_folds)){stop("df contains too few 0s for this number of folds")}

  #divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1,y)), c("ID", y)]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==0,])), replace=FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==1,])), replace=FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0,fold_i_1))
  folds_i <- merge(df1[,1:2], df_sort[,1:3], by='ID', all = FALSE, sort = FALSE)[,4]

  #names of the columns
  lambdanames <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
  taunames <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
  foldnames <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
  #measures (train and test)
  list_of_measures <- c("pye_KS_L12", "pye_KS_L1", "pye_KS_EN", "pye_KS_SCAD", "pye_KS_MCP",
                        "auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                        "fdr", "mcc", "corrclass")
  pye_KS_L12 <- pye_KS_L1 <- pye_KS_EN <- pye_KS_SCAD <- pye_KS_MCP <- NULL
  auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures){
    #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(c(1:n_folds),taunames)))
    assign(mes, lapply(lambda, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(foldnames,taunames)))))
    eval(substitute(names(x)<- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames, unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))))
  n_gammas <- lapply(lambda, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(foldnames, taunames)))
  names(n_gammas) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))

  betas <- vector(mode="list", length=n_folds)
  names(betas) <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
  gammas <- vector(mode="list", length=n_folds)
  names(gammas) <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))

  results = list(pye_KS_L12=pye_KS_L12, pye_KS_L1=pye_KS_L1, pye_KS_EN= pye_KS_EN, pye_KS_SCAD=pye_KS_SCAD,
                 pye_KS_MCP= pye_KS_MCP, auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index,
                 sensitivity=sensitivity, specificity=specificity,
                 geometric_mean=geometric_mean, fdr=fdr, mcc=mcc,
                 corrclass=corrclass, n_betas=n_betas, n_gammas=n_gammas,
                 betas= betas, gammas=gammas)

  cat("Starting CV with the following estimation/convergence parameters: \n")
  cat("max_iter:", max_iter, "\n")
  cat("min_alpha:", min_alpha, "\n")
  cat("convergence_error:", convergence_error, "\n")
  cat("stepsizeShrink:", stepsizeShrink, "\n")

  #fill the matrices
  results <- mapply(function(k) pye_KS.cv(df=df1[,names(df1)!="ID"], X=X, y=y, C=C, lambda=lambda, tau=tau, trace=trace,
                                          alpha=alpha, a1=a1, a2=a2, penalty=penalty, alpha_g=alpha_g, penalty_g=penalty_g,
                                          folds_i, k, regressors_betas, pye_KS_L12=pye_KS_L12, pye_KS_L1=pye_KS_L1,
                                          pye_KS_EN=pye_KS_EN, pye_KS_SCAD=pye_KS_SCAD, pye_KS_MCP=pye_KS_MCP,
                                          auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index,
                                          sensitivity=sensitivity, specificity=specificity,
                                          geometric_mean=geometric_mean, fdr=fdr,
                                          mcc=mcc, corrclass=corrclass,
                                          n_betas=n_betas, n_gammas=n_gammas,
                                          max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                          stepsizeShrink=stepsizeShrink,
                                          used_cores=used_cores, kernel=kernel, trend=trend,
                                          delta=delta, max_alpha=max_alpha, beta_start_input=beta_start_input,
                                          beta_start_default=beta_start_default, c_zero_fixed=c_zero_fixed,
                                          c_function_of_covariates=c_function_of_covariates,
                                          kernel_g=kernel_g, a1_g=a1_g, a2_g=a2_g, trend_g=trend_g,
                                          gamma_start_input=gamma_start_input,
                                          gamma_start_default=gamma_start_default,
                                          regressors_gammas=regressors_gammas, max_iter_g=max_iter_g, delta_g=delta_g,
                                          max_alpha_g=max_alpha_g, stepsizeShrink_g=stepsizeShrink_g,
                                          min_alpha_g=min_alpha_g, convergence_error_g=convergence_error_g,
                                          run_aauc=run_aauc, long_suffix=long_suffix),
                                          seq(1:n_folds))

  wrapper <- function(results, mes, i, dataset, tau){
    m <- as.matrix(sapply(c(1:n_folds), function(k) getElement(getElement(results[[k]],mes)[[i]],dataset)[k,,drop=FALSE]))
    if(length(tau) > 1){
      m2 <- t(m)
    } else{
      m2 <- m
    }
    colnames(m2) <- taunames
    rownames(m2) <- foldnames
    return(m2)
  }

  #prepare the results
  temp_pye <- get(paste0("pye_KS_", penalty))
  temp_pye[] <- lapply(1:length(lambda), function(i) list(train = wrapper(results, paste0("pye_KS_", penalty), i, "train", tau) , test = wrapper(results, paste0("pye_KS_", penalty), i, "test", tau)))
  #temp_pye[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]], paste0("pye_KS_", penalty))$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]],paste0("pye_KS_", penalty))$test[k,]))))
  assign(paste0("pye_KS_", penalty) , temp_pye)

  for (mes in list_of_measures){
    #aauc[] <- lapply(lambda, function(i) list(train = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$test[k,])))))
    #names(aauc) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
    assign(mes, lapply(1:length(lambda), function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
    eval(substitute(names(x)<- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
  }

  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))
  n_gammas[] <- lapply(1:length(lambda), function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k,]))))

  measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye")
  list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity",
                        "geometric_mean", "fdr", "mcc", "corrclass", paste0("pye_KS_", penalty))

  #t(sapply(c(1:n_folds), function(k) sapply(1:length(lambda), function(i) results[[k]]@auc[[i]]$train[k])))

  lambda_hat_pye_KS_yi <- lambda_hat_pye_KS_auc <- lambda_hat_pye_KS_aauc <- lambda_hat_pye_KS_aYI <- lambda_hat_pye_KS_ccr <- NULL
  lambda_hat_pye_KS_sen <- lambda_hat_pye_KS_spc <- lambda_hat_pye_KS_gm <- lambda_hat_pye_KS_pye <- NULL
  tau_hat_pye_KS_yi <- tau_hat_pye_KS_auc <- tau_hat_pye_KS_aauc <- tau_hat_pye_KS_aYI <- tau_hat_pye_KS_ccr <- NULL
  tau_hat_pye_KS_sen <- tau_hat_pye_KS_spc <- tau_hat_pye_KS_gm <- tau_hat_pye_KS_pye <- NULL
  for (i in 1:length(measures)){
    #in general, create tau_hat equal to NA
    assign(paste0("tau_hat_pye_KS", measures[i]), NA)

    if(length(tau) != 1){
      measures_matrix <- t(sapply(1:length(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    } else {
      measures_matrix <- as.matrix(sapply(1:length(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    }

    rownames(measures_matrix) <- lambdanames
    colnames(measures_matrix) <- taunames

    if (!is.na(max(measures_matrix))){

      max_measures <- which(measures_matrix == max(measures_matrix), arr.ind = TRUE)[1,]
      assign(paste0("lambda_hat_pye_KS_", measures[i]), lambda[max_measures[1]])

      if (c_function_of_covariates==TRUE){
        if (length(tau)>1){
          assign(paste0("tau_hat_pye_KS_", measures[i]), tau[max_measures[2]])
        } else {
          assign(paste0("tau_hat_pye_KS_", measures[i]), tau)
        }
      }
    } else {
      assign(paste0("lambda_hat_pye_KS_", measures[i]), NA)
    }
  }

  #if tau2 is not zero, but the condition is valued for vectors as well
  if(length(tau2) > 1 || sum(tau2) > 0) {
    c_function_of_covariates <- TRUE
    tau <- tau2

    #re-create the output tables
    #names of the columns
    taunames <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
	  foldnames <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
    #measures (train and test)
    list_of_measures <- c("pye_KS_L12", "pye_KS_L1", "pye_KS_EN", "pye_KS_SCAD", "pye_KS_MCP",
                          "auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                          "fdr", "mcc", "corrclass")
    pye_KS_L12 <- pye_KS_L1 <- pye_KS_EN <- pye_KS_SCAD <- pye_KS_MCP <- NULL
    auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
	
	#now we execute pye_KS.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
    lambda_star <- get(paste0("lambda_hat_pye_KS_", measure_to_select_lambda))
	
    for (mes in list_of_measures){
      #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames= list(c(1:n_folds),taunames)))
      assign(mes, lapply(lambda_star, function(x) list (train = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(c(1:n_folds),taunames)))))
      eval(substitute(names(x)<- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
    }

    n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda_star), dimnames= list(foldnames, unlist(lapply(lambda_star, function (x) paste("lambda", x, sep="=")))))
    n_gammas <- lapply(lambda_star, function(x) matrix(NA, nrow = n_folds, ncol = length(tau), dimnames= list(foldnames, taunames)))
    names(n_gammas) <- unlist(lapply(lambda_star, function (x) paste("lambda", x, sep="=")))

    betas <- vector(mode="list", length=n_folds)
    names(betas) <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
    gammas <- vector(mode="list", length=n_folds)
    names(gammas) <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))

	cat(" \n Starting the CV of tau using as lambda:", lambda_star, ", that is the best value of lambda as per:", measure_to_select_lambda, "\n")
	
    #re-fill the matrices
    results <- mapply(function(k) pye_KS.cv(df=df1[,names(df1)!="ID"], X=X, y=y, C=C, lambda=lambda_star, tau=tau, trace=trace,
                                            alpha=alpha, a1=a1, a2=a2, penalty=penalty, alpha_g=alpha_g, penalty_g=penalty_g,
                                            folds_i, k, regressors_betas, pye_KS_L12=pye_KS_L12, pye_KS_L1=pye_KS_L1,
                                            pye_KS_EN=pye_KS_EN, pye_KS_SCAD=pye_KS_SCAD, pye_KS_MCP=pye_KS_MCP,
                                            auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index,
                                            sensitivity=sensitivity, specificity=specificity,
                                            geometric_mean=geometric_mean, fdr=fdr,
                                            mcc=mcc, corrclass=corrclass,
                                            n_betas=n_betas, n_gammas=n_gammas,
                                            max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                            stepsizeShrink=stepsizeShrink,
                                            used_cores=used_cores, kernel=kernel, trend=trend,
                                            delta=delta, max_alpha=max_alpha, beta_start_input=beta_start_input,
                                            beta_start_default=beta_start_default, c_zero_fixed=c_zero_fixed,
                                            c_function_of_covariates=c_function_of_covariates,
                                            kernel_g=kernel_g, a1_g=a1_g, a2_g=a2_g, trend_g=trend_g,
                                            gamma_start_input=gamma_start_input,
                                            gamma_start_default=gamma_start_default,
                                            regressors_gammas=regressors_gammas, max_iter_g=max_iter_g, delta_g=delta_g,
                                            max_alpha_g=max_alpha_g, stepsizeShrink_g=stepsizeShrink_g,
                                            min_alpha_g=min_alpha_g, convergence_error_g=convergence_error_g,
											run_aauc=run_aauc, long_suffix=long_suffix),
                                            seq(1:n_folds))

    wrapper <- function(results, mes, i, dataset, tau){
      m <- as.matrix(sapply(c(1:n_folds), function(k) getElement(getElement(results[[k]],mes)[[i]],dataset)[k,,drop=FALSE]))
      if(length(tau) > 1){
        m2 <- t(m)
      } else{
        m2 <- m
      }
      colnames(m2) <- taunames
	  rownames(m2) <- foldnames
      return(m2)
    }

    #prepare the results
    temp_pye <- get(paste0("pye_KS_", penalty))
    temp_pye[] <- lapply(1, function(i) list(train = wrapper(results, paste0("pye_KS_", penalty), i, "train", tau) , test = wrapper(results, paste0("pye_KS_", penalty), i, "test", tau)))
    #temp_pye[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]], paste0("pye_KS_", penalty))$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]],paste0("pye_KS_", penalty))$test[k,]))))
    assign(paste0("pye_KS_", penalty) , temp_pye)

    for (mes in list_of_measures){
      #aauc[] <- lapply(lambda, function(i) list(train = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$test[k,])))))
      #names(aauc) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
      assign(mes, lapply(1, function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
      eval(substitute(names(x)<- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
    }

    n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))
    n_gammas[] <- lapply(1, function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k,]))))

    betas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas)))
    gammas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

    measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye")
    list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity",
                           "geometric_mean", "fdr", "mcc", "corrclass", paste0("pye_KS_", penalty))

    for (i in 1:length(measures)){
      #in general, create tau_hat equal to NA
      assign(paste0("tau_hat_pye_KS", measures[i]), NA)
      if(length(tau) != 1){
        measures_matrix <- t(sapply(1, function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
      } else {
        measures_matrix <- as.matrix(sapply(1, function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
      }

      rownames(measures_matrix) <- paste0("lambda=", lambda_star)
      colnames(measures_matrix) <- taunames

      if (!is.na(max(measures_matrix))){

        #if the number of gammas is zero, we set the measure to zero
        mean_n_gammas <- t(sapply(n_gammas, function(i) colMeans(i)))
        measures_matrix[mean_n_gammas==0] <- 0

        max_measures <- which(measures_matrix == max(measures_matrix), arr.ind = TRUE)[1,]

        if (c_function_of_covariates==TRUE){
          if (length(tau)>1){
            assign(paste0("tau_hat_pye_KS_", measures[i]), tau[max_measures[2]])
          } else {
            assign(paste0("tau_hat_pye_KS", measures[i]), tau)
          }
        }
      }
    }
  } else {
    lambda_star <- NA
  }

  cv_time = difftime(Sys.time(), start_time , units = "mins")

  if (trace %in% c(1,2)) {
	cat("----------------> END OF THE CROSS-VALIDATION OF THE PYE METHODS <----------------- \n")
    cat("----------------> For the whoole Cross-validation it took:", cv_time, "minutes \n")
  }

  return(list(penalty=penalty, penalty_g=penalty_g, kernel=kernel, cv_time=cv_time, pye_KS_L12=pye_KS_L12,
              pye_KS_L1=pye_KS_L1, pye_KS_EN=pye_KS_EN, pye_KS_SCAD=pye_KS_SCAD, pye_KS_MCP=pye_KS_MCP,
              auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index, sensitivity=sensitivity,
              specificity=specificity, geometric_mean=geometric_mean,
              fdr=fdr, mcc=mcc, corrclass=corrclass,
              lambda_hat_pye_KS_yi=lambda_hat_pye_KS_yi, lambda_hat_pye_KS_auc=lambda_hat_pye_KS_auc,
              lambda_hat_pye_KS_aauc=lambda_hat_pye_KS_aauc, lambda_hat_pye_KS_aYI=lambda_hat_pye_KS_aYI,
              lambda_hat_pye_KS_ccr=lambda_hat_pye_KS_ccr, lambda_hat_pye_KS_sen=lambda_hat_pye_KS_sen,
              lambda_hat_pye_KS_spc=lambda_hat_pye_KS_spc,lambda_hat_pye_KS_gm=lambda_hat_pye_KS_gm,
              lambda_hat_pye_KS_pye=lambda_hat_pye_KS_pye,
              tau_hat_pye_KS_yi=tau_hat_pye_KS_yi, tau_hat_pye_KS_auc=tau_hat_pye_KS_auc,
              tau_hat_pye_KS_aauc=tau_hat_pye_KS_aauc, tau_hat_pye_KS_aYI=tau_hat_pye_KS_aYI,
              tau_hat_pye_KS_ccr=tau_hat_pye_KS_ccr, tau_hat_pye_KS_sen=tau_hat_pye_KS_sen,
              tau_hat_pye_KS_spc=tau_hat_pye_KS_spc, tau_hat_pye_KS_gm=tau_hat_pye_KS_gm,
              tau_hat_pye_KS_pye=tau_hat_pye_KS_pye, c_function_of_covariates=c_function_of_covariates,
              simultaneous=simultaneous, measure_to_select_lambda=measure_to_select_lambda,
              lambda_star=lambda_star, n_betas=n_betas, n_gammas=n_gammas, betas=betas, gammas=gammas))
}

