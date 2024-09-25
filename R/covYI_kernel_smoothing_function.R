#install.packages("KernSmooth")
#library(KernSmooth)

#' create the output class of covYI
# covYI_output_KS <- setClass(Class="covYI_output_KS",
#                        representation(
#                          covYI_L12 ="numeric",
#                          covYI_L1 ="numeric",
#                          covYI_EN ="numeric",
#                          covYI_SCAD ="numeric",
#                          covYI_MCP="numeric",
#                          gr_covYI_L12 ="numeric",
#                          gr_covYI_L1 ="numeric",
#                          gr_covYI_EN ="numeric",
#                          gr_covYI_SCAD ="numeric",
#                          gr_covYI_MCP ="numeric",
#                          gr_yi ="numeric",
#                          gr_phi_L12 ="numeric",
#                          gr_phi_L1 ="numeric",
#                          gr_phi_EN ="numeric",
#                          gr_phi_SCAD ="numeric",
#                          gr_phi_MCP ="numeric",
#                          youden_index ="numeric",
#                          sensitivity ="numeric",
#                          specificity ="numeric",
#                          fdr ="numeric",
#                          mcc ="numeric",
#                          auc ="numeric",
#                          aauc ="numeric",
#                          aYI ="numeric",
#                          corrclass="numeric",
#                          empir_suggest_yi_c="numeric",
#                          c_hat="numeric",
#                          z_hat="data.frame",
#                          y_hat="data.frame",
#                          TP = "numeric",
#                          TN = "numeric",
#                          FP = "numeric",
#                          FN = "numeric",
#                          input_data="list"
#                        )
# )

#' @title covYI_KS
#'
#' @description The Penalizad Youden Index (PYE) function based on the Kernel
#' Smooth density estimator. It not only performs value of PYE but it returns
#' all the necessary for the estimation process, like measure of fit and
#' derivatives. It works for all the considered penalties (L12, L1,
#' EN, SCAD and MCP)
#'
#' @param df the input dataset
#' @param z a single or a known combination of regressors. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "z_hat"
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param gammas the coefficients of the covariates combination used for the
#' evaluation of the cut-point c (the first element is the constant!!)
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param C covariate variables. It can be of type dataframe, containing the
#' same covariates included in df, or just a vector of character. Default is NULL
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param h_exponent parameter of the bandwidth of the Kernel Smooth
#' @param prediction if TRUE the empirical maximum Youden index is returned as
#' Youden index performance measure
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
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
#' #The combination of biomarkers (regressors) is known when we apply covYI
#' z <- simMicroarrayData_cov02_dim50_covariates$z
#' colnames(z) <- "z"
#' df1 <- transform(merge(x = df, y = z, by = 0, all.x = TRUE),
#'   row.names=Row.names, Row.names=NULL)
#' penalty <- "L12"
#' tau <- 0.1
#' gammas <- rep(1, length(C))
#' c <- 0
#' prox_penalty = get(paste0("proximal_operator_", penalty))
#'
#' PYE_result <- covYI_KS(df=df1, z="z", y=y, C=C, gammas=gammas, tau=tau,
#'   alpha=0.5, a1=3.7, a2=3, penalty=penalty)
#' print(PYE_result)
#'
#' @importFrom evmix kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom ROCnReg AROC.sp compute.threshold.AROC
#' @importFrom plyr join
#' @export
covYI_KS <- function(df, z="z_hat", y="y", C, gammas, tau, kernel="gaussian", alpha=0.5, a1=3.7, a2=3,
                     penalty="L1", h_exponent=0.2, prediction=FALSE, run_aauc=FALSE){

  #IMPORTANT: the first element of gammas is the constant term!!!!!!

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (inherits(class(C), "data.frame")){
    C <- names(C)
  } else if (!inherits(class(C), "character")){stop("C can only be of class character or data.frame.")}
  if (inherits(class(z), "data.frame")){
    z <- names(z)
  } else if (!inherits(class(z), "character")){stop("z can only be of class character or data.frame.")}
  if (inherits(class(y), "data.frame")){
    y <- names(y)
  } else if (!inherits(class(y), "character")){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  #check if const already exists in the dataset
  #if("const" %in% colnames(df)){stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  #Here I dont create the constant since it should already exist, but I check if it needs to be created
  if(("const" %in% C) & (!("const" %in% names(df)))){
    const <- rep(1, length(ID)) #the constant for the coefficients gammas
    df1 <- cbind(ID, df[, (names(df) %in% c(y,z)), drop=FALSE], const, df[, (names(df) %in% C), drop=FALSE])
  } else {
    df1 <- cbind(ID, df[, (names(df) %in% c(y,z)), drop=FALSE], df[, (names(df) %in% C), drop=FALSE])
  }

  #put "const" in C
  #C1 <- c("const", C) #C Already contains the constant

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(tau)>1){stop("covYI_KS stopped because tau needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))){
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(gammas) != length(C)){
    stop("The number of element of gammas is different then the number of columns in C")
  }

  #divide control and diseased
  z_y0 <- as.matrix(df1[df1$y==0, (names(df1) %in% z), drop=FALSE])#drop=FALSE permits it to remain a data frame even if it is composed by a single column
  z_y1 <- as.matrix(df1[df1$y==1, (names(df1) %in% z), drop=FALSE])

  #compute c
  #IMPORTANT: the first element of gammas is the constant term!

  c_y0 <- as.matrix(df1[df1$y==0, (names(df1) %in% C), drop=FALSE])%*%gammas
  c_y1 <- as.matrix(df1[df1$y==1, (names(df1) %in% C), drop=FALSE])%*%gammas

  #kernel from the YI paper
  #healthy
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y0)!=0){
  #   h0 = kedd::h.bcv(as.numeric(c-z_y0), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*length(z_y0)^(-h_exponent)
  # }
  if(h0<0.1) {h0=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t0 <- (c_y0-z_y0)/h0
  f0 <- sum(evmix::kpz(z = as.numeric(t0), kernel = kernel))/length(z_y0)
  #f0 <- sum(pnorm(t0, 0, 1))/length(z_y0)

  #diseased
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y1)!=0){
  #   h1= kedd::h.bcv(as.numeric(z_y1), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
  h1=0.9*min(stats::sd(z_y1), stats::IQR(z_y1)/1.34)*length(z_y1)^(-h_exponent)
  # }
  if(h1<0.1) {h1=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  t1 <- (c_y1-z_y1)/h1
  f1 <- sum(evmix::kpz(z = as.numeric(t1), kernel = kernel))/length(z_y1)
  #f1 <- sum(pnorm(t1, 0, 1))/length(z_y1)

  yi <- f0 - f1

  #evaluate the gradient. NOTE: the variables to which consider the gradient
  #are betas and gammas elements
  #h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*(length(z_y0)^(-h_exponent))
  #t0 <- (c_y0-z_y0)/h0
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr0_gammas <- apply (df1[df1$y==0,(names(df1) %in% C), drop=FALSE], 2,
                        function(g) (t(stats::dnorm(t0, 0, 1)) %*% (g/h0))/length(z_y0))
  } else {
    gr0_gammas <- apply (df1[df1$y==0,(names(df1) %in% C), drop=FALSE], 2,
            function(g) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% (g/h0))/length(z_y0))
  }

  #h1=0.9*min(stats::sd(z_y1), stats::IQR(z_y1)/1.34)*(length(z_y1)^(-h_exponent))
  #t1 <- (c_y1-z_y1)/h1
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr1_gammas <- apply (df1[df1$y == 1,(names(df1) %in% C), drop=FALSE], 2,
                        function(g) (t(stats::dnorm(t1, 0, 1)) %*% (g/h1))/length(z_y1))
  } else {
    gr1_gammas <- apply (df1[df1$y == 1,(names(df1) %in% C), drop=FALSE], 2,
           function(g) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% (g/h1))/length(z_y1))
  }

  #gradient
  gr_yi <- gr0_gammas - gr1_gammas #note that in gammas we are not including the constant here

  #if (abs(sum(gr_yi))<0.0001){
  #  gr_yi <- gr_yi*(10^(round(abs(log10(abs(gr_yi))+1))))
  #}

  #append data (c)
  c<-as.data.frame(rbind(c_y0, c_y1))
  c['ID'] <- as.numeric(rownames(c))
  #c <- c[order(c$ID), ]
  names(c)[names(c) == 'V1'] <- 'c_hat'
  df1<-plyr::join(df1, c, by='ID')
  #df1<-merge(df1, c, by='ID', all = TRUE)
  rownames(df1) <- df1$ID

  #find the optimal cut-point
  opt <- OptimalCutpoints::optimal.cutpoints(data= df1, X = z, status = "y", methods = "Youden", tag.healthy =0)

  #compute y_hat
  df1["y_hat"] <- ifelse(df1$z_hat > df1$c_hat ,1,0)

  #confusion matrix
  TP <- sum(ifelse(z_y1 >= c_y1 ,1,0))
  TN <- sum(ifelse(z_y0 < c_y0 ,1,0))
  FP <- sum(ifelse(z_y0 >= c_y0 ,1,0))
  FN <- sum(ifelse(z_y1 < c_y1 ,1,0))

  #compute the Youden index
  spec <- TN/(TN+FP)
  fnr <- FN/(FN+TP)
  #yi1 <- spec - fnr

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
  #try({test <- AUC::auc(AUC::roc(df1$z_hat,factor(df1$y)))}, silent = TRUE)
  auc <- opt$Youden$Global$measures.acc$AUC[1]
  yi1 <- opt$Youden$Global$optimal.criterion[1]

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.sp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                      p = seq(0,1,l=101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y", tag.healthy = 0,
  #                                     data = df1, p = seq(0,1,l=101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.bnp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1, p = seq(0,1,l=101),
  #                           compute.lpml = FALSE, compute.WAIC = FALSE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y", tag.h = 0, data = df1,
  #                              p = seq(0,1,l=101), B = 500)

  #start_time2 <- Sys.time()
  #this is longher then AROC.bsp
  #tmp2 <- try({aroc <- AROC::AROC.bnp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                   scale = TRUE, p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                   a = 2, b = 0.5, L = 10, nsim = 5000, nburn = 1000)}, silent = TRUE)
  #cat("temp2:")
  #print(tmp2)
  #print(difftime(Sys.time(), start_time2 , units = "mins"))
  #start_time1 <- Sys.time()

  #this does not have the YI
  #aroc <- AROC::AROC.bsp(formula.healthy = z_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                scale = TRUE, p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                a = 2, b = 0.5, nsim = 5000, nburn = 1500)

  if (run_aauc){
    tmp <- try({

      #AROC.sp function estimates the covariate-adjusted ROC curve (AROC) using the semiparametric approach
      #proposed by Janes and Pepe (2009).

      #if we have all zeros, the function returns a warning message. We suppress it just for this function:
      options(warn=-1)

      aroc <- ROCnReg::AROC.sp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1,
                               p = seq(0,1,l=101), B = 500)

      options(warn=0)

      #NOt werking properly:
      #aroc <- ROCnReg::AROC.bnp(formula.h = z_hat ~ c_hat, group = "y", tag.h = 0, data = df1,
      #                         p = seq(0,1,l=101))

      #aroc <- ROCnReg::AROC.kernel(marker = "z_hat", covariate = "c_hat", group = "y",
      #                             tag.h = 0, data = df1, bw = "LS", regtype = "LC",
      #                             pauc = ROCnReg::pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
      #                             B = 500)

      ### Threshold values based on the YI
      th_AROC.sp <- ROCnReg::compute.threshold.AROC(aroc, criterion = "YI")

      #summary.AROC(aroc)

    }, silent = TRUE)


    #print(difftime(Sys.time(), start_time1 , units = "mins"))

    if(inherits(tmp, "try-error")){
      aauc <- 0
      aYI <- 0
    }else{
      #summary(m3)
      aauc <- aroc$AUC[1]
      aYI <- th_AROC.sp$YI
      #print(aauc)
      #print(aYI)
    }
  } else {
    aauc <- 0
    aYI <- 0
  }

  #if all the betas are zero, the measures are zeros
  #if(sum(gammas)==0){
  #  spec <- 0
  #  fnr <- 0
  #  sensitivity <- 0
  #  gm <- 0
  #  yi1 <- 0
  #  fdr <- 0
  #  mcc <- 0
  #  corrclass <- 0
  #}

  #L_(1/2) penalization
  phi_L12_g <- tau*sum(abs(gammas[1:length(gammas)])^(1/2))#NB:we DO HAVE the intercept!!!
  if(is.na(phi_L12_g)){
    stop('problem with the penalization phi_L12_g in covYI function')}


  #L1 penalization
  phi_L1_g <- tau*sum(abs(gammas[1:length(gammas)]))#NB:we DO HAVE the intercept!!!
  if(is.na(phi_L1_g)){
    stop('problem with the penalization phi_L1_g in covYI function')}

  #Elastic-Net penalization
  phi_EN_g <- tau*((alpha)*(sum(abs(gammas[1:length(gammas)]))) + ((1-alpha)/2)*(sum(gammas[1:length(gammas)]^2)))
  if(is.na(phi_EN_g)){
    stop('problem with the penalization phi_EN_g in covYI function')}

  #SCAD penalization
  phi_SCAD_g <- SCAD_function(gammas[1:length(gammas)], tau, a=a1)
  if(is.na(phi_SCAD_g)){
    stop('problem with the penalization phi_SCAD_g in covYI function')}

  #MCP penalization
  phi_MCP_g <- MCP_function(gammas[1:length(gammas)], tau, a=a2)
  if(is.na(phi_MCP_g)){
    stop('problem with the penalization phi_MCP_g in covYI function')}

  #join yi and the penalty functions
  covYI_KS_L12 <- yi - phi_L12_g
  covYI_KS_L1 <- yi - phi_L1_g
  covYI_KS_EN <- yi - phi_EN_g
  covYI_KS_SCAD <- yi - phi_SCAD_g
  covYI_KS_MCP <- yi - phi_MCP_g

  #if we are in prediction this is the right yi
  if (prediction==TRUE){
    yi=yi1
  }

  return(list(covYI_KS_L12=covYI_KS_L12, covYI_KS_L1=covYI_KS_L1, covYI_KS_EN = covYI_KS_EN, covYI_KS_SCAD = covYI_KS_SCAD,
              covYI_KS_MCP = covYI_KS_MCP, gr_yi=gr_yi, youden_index=yi, sensitivity=sensitivity,
              specificity=spec, geometric_mean=gm, fdr=fdr, mcc=mcc, auc=auc, aauc=aauc, aYI=aYI,
              corrclass=corrclass, z_hat=df1[,c("ID", z)], c_hat=df1[,c("ID","c_hat")] ,
              y_hat=df1[,c("ID","y_hat")], TP=TP, TN=TN, FP=FP, FN=FN,
              input_data = list(gammas=gammas, tau=tau, alpha=alpha, a1=a1, a2=a2, mode=mode, kernel=kernel)))
}



#' @title covYI_KS_estimation
#'
#' @description function to estimate the optimal value of gammas and c maximizing
#' the PYE function. To find the optimum the mmAPG and mnmAPG algorithms are
#' used.
#'
#' @param df the input dataset
#' @param z a single or a known combination of regressors. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "z_hat"
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param C covariate variables. It can be of type dataframe, containing the
#' same covariates included in df, or just a vector of character. Default is NULL
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point of gamma.
#' If "zeros", it starts with all zero values, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable. Default
#' is "zeros"
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param regressors_gammas a vector containing the real gammas (if known).
#' Default is NULL
#' @param fold number of the fold if the function is called in cross-validation.
#' Just for visualization purposes. Default is NULL
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 10000
#' @param max.print number of elements to show if printing the results. Default
#' is 10
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
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
#'
#' @return a list containing the optimal value of gammas, the value
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
#' regressors_gammas<-simMicroarrayData_cov02_dim50_covariates$ncovariates
#' #The combination of biomarkers (regressors) is known when we apply covYI
#' z <- simMicroarrayData_cov02_dim50_covariates$z
#' colnames(z) <- "z"
#' df1 <- transform(merge(x = df, y = z, by = 0, all.x = TRUE),
#'   row.names=Row.names, Row.names=NULL)
#' penalty <- "SCAD"
#' tau <- 0.1
#' betas <- rep(1, length(X))
#' c <- 0
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#' gamma_start_default <- "zeros"
#' alpha <- 0.5
#'
#' covYI_estimation_result <- covYI_KS_estimation(df=df1, z="z", y=y, C=C, tau=tau,
#'   penalty=penalty, trace=2, gamma_start_default=gamma_start_default,
#'   gamma_start_input=NULL, alpha=alpha, regressors_gammas=regressors_gammas,
#'   max_iter=5, run_aauc=TRUE)
#' print(covYI_estimation_result)
#'
#' @importFrom stats setNames
#' @export
#Estimation of the parameter using covYI_KS
covYI_KS_estimation <- function(df, z="z_hat", y="y", C, tau, penalty="L1", gamma_start_input=NULL,
                                gamma_start_default="zeros", alpha=0.5, a1=3.7, a2=3,
                                regressors_gammas=NULL, fold=NULL, max_iter=10000, max.print=10,
                                trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                                min_alpha=1e-10, convergence_error=1e-7, trace=1, seed=1, kernel="gaussian",
								                run_aauc=FALSE){


  start_time <- Sys.time()

  options(max.print = max.print)

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (!inherits(class(C), "character")){stop("C can only be of class character or data.frame.")}
  if (!inherits(class(z), "character")){stop("X can only be of class character or data.frame.")}
  if (!inherits(class(y), "character")){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  #check if const already exists in the dataset
  if("const" %in% colnames(df)){stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal const")}
  if("const" %in% colnames(C)){stop("const already exists as column in C! Please delete or rename this column since I need this name to set the internal const")}

  if(sum(names(df) %in% C)==0){stop("No covariates (c) present in the dataframe (df).")}

  ID <- rownames(df)
  const <- rep(1, length(ID)) #the constant for the coefficients gammas
  df1 <- cbind(ID, df[,(names(df) %in% c(y,z)), drop=FALSE], const, df[,(names(df) %in% C)])

  #adding "const" in C
  C1 <- c("const", C)

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(tau)>1){stop("covYI_KS_estimation stopped because tau needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}


  gammas1_initial_zeros <- rep(0, length(C1))
  gammas1_initial_corr <- NULL

  names_gammas <- C1

  #initializing the parameters (gammas)
  if (length(gamma_start_input)==0){
    if (gamma_start_default=="zeros"){

      gammas_start <- gammas1_initial_zeros
      names(gammas_start) <- names_gammas
      gammas1 <- gammas_start

    } else if (gamma_start_default=="corr") {

      #compute the corr between every x and y
      corr.xy <- data.frame(matrix(ncol = ncol(df1[,(names(df1) %in% C1)]), nrow = 0))
      names(corr.xy) <- names_gammas
      corr.xy[1:3,] <- t(sapply(c("pearson", "kendall", "spearman"),
                                function (h) c(0,unlist(lapply(C1[-1], function (x) stats::cor(df1$y, df1[,x],  method = h))))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      gammas1_initial_corr <- c(as.numeric(mean.corr.xy))
      gammas_start <- gammas1_initial_corr
      names(gammas_start) <- names_gammas
      gammas1 <- gammas_start

    } else {stop("The parameter gamma_start_default is wrongly defined. It can be only equal to 'zeros', 'corr'.")}
  } else {
    if (length(gamma_start_input) == ncol(df1[,(names(df1) %in% C1)])){
      gammas_start <- gamma_start_input
      names(gammas_start) <- names_gammas
      gammas1 <- gammas_start
    } else if (length(gamma_start_input) != 0){
      stop("The length of the parameter gamma_start_input is not equal to ", ncol(df1[,(names(df1) %in% C1)]), ". Remeber that the first parameter to include is the constant.")
    }
  }

  if (penalty == "SCAD") {a=a1} else {a=a2}
  prox_penalty = get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  delta_fx_gammas <- function(gammas){
    result <- -covYI_KS(df=df1[,!(names(df1) %in% c("ID")), drop=FALSE], z=z, y=y, C=C1[C1 %in% names(gammas)], gammas=gammas,
                       tau=tau, kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty)$gr_yi
    result <- result[names(result) %in% C1]
    return(result)
  }

  proxx_gammas <- function(x, eta){
    #old version:
    #result <- c(x[(names(x) == "const")], prox_penalty(betas=x[!(names(x) == "const")], lambda=eta*tau, alpha=alpha, a=a)) #-1 since I dont penalize the constant
    #new version (let the const enter in the proximal operator)
    result <- prox_penalty(betas=x, lambda=eta*tau, alpha=alpha, a=a)
    return(result)
  }

  Fx <- function(aa){
    result <- -getElement(covYI_KS(df=df1[,!(names(df1) %in% c("ID")), drop=FALSE], z=z, y=y, C=C1[C1 %in% names(aa)], gammas=aa, tau=tau,
                            kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty), paste0("covYI_KS_", penalty))
    return(result)
  }

  if (trend == "monotone"){
    estim <- mmAPG(x0=gammas1, c_pos=NULL, delta_fx=delta_fx_gammas, proxx=proxx_gammas, Fx=Fx,
                    lambda=tau, penalty=penalty, fold=fold, stepsizeShrink=stepsizeShrink, delta=delta,
                    max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                    max_alpha=max_alpha, trace=trace, seed=seed, zeros_stay_zeros_from_iteration=20)
  } else if (trend == "nonmonotone"){
    estim <- mnmAPG(x0=gammas1, c_pos=NULL, delta_fx=delta_fx_gammas, proxx=proxx_gammas, Fx=Fx,
                     lambda=tau, penalty=penalty, fold=fold, stepsizeShrink=stepsizeShrink, delta=delta,
                     max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                     max_alpha=max_alpha, trace=trace, seed=seed, zeros_stay_zeros_from_iteration=20)
  }

  gammas1 <- estim$x1

  covYI_KS_value <- covYI_KS(df=df1[,!(names(df1) %in% c("ID"))], z=z, y=y, C=C1,
                           gammas=gammas1, tau=tau, kernel=kernel,
                           alpha=alpha, a1=a1, a2=a2, penalty=penalty, run_aauc=run_aauc)

  if (trace %in% c(1,2)){
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend,"version. ; \n")
    if (length(fold)!=0) {cat("fold =", fold, "; ")}
    visualize_gammas = c(gammas1[which(gammas1!=0)])
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters , "; tau:", tau, "; penalty:", penalty, "; covYI_KS:", getElement(covYI_KS_value,  paste0("covYI_KS_",penalty)), "; youden_index:", covYI_KS_value$youden_index, "; aYI:", covYI_KS_value$aYI, "; sensitivity:", covYI_KS_value$sensitivity, "; specificity:", covYI_KS_value$specificity, "; geometric_mean:", covYI_KS_value$geometric_mean, "; fdr:", covYI_KS_value$fdr, "; mcc:", covYI_KS_value$mcc, "; auc:", covYI_KS_value$auc, "; aauc:", covYI_KS_value$aauc, "; corrclass:", covYI_KS_value$corrclass, "; \n")
    cat("TP:", covYI_KS_value$TP, "; TN:", covYI_KS_value$TN, "; FP:", covYI_KS_value$FP, "; FN:", covYI_KS_value$FN, "; gammas: \n")
    cat("\n")
    print(visualize_gammas)
    cat("\n")
  }

  #total number of variables
  n_total_var_gammas= length(gammas1)

  #n of predicted zeros
  n_predicted_zeros_gammas= sum(gammas1==0)

  #n of predicted gammas different from 0
  n_predicted_non_zeros_gammas = sum(gammas1!=0)

  #compute other measures
  if (length(regressors_gammas)!=0){

    #n of gammas catch
    n_caught_gammas = sum(gammas1[regressors_gammas!=0]!=0)

    #n of gammas not catch
    n_non_caught_gammas = sum(gammas1[regressors_gammas!=0]==0)

    #n of zeros catch
    n_caught_zero_gammas = sum(gammas1[regressors_gammas==0]==0)

    #n of zeros not catch
    n_zero_not_caught_gammas = sum(gammas1[regressors_gammas==0]!=0)
  } else {
    n_caught_gammas=NA
    n_non_caught_gammas=NA
    n_caught_zero_gammas=NA
    n_zero_not_caught_gammas=NA
  }

  name1 = paste0("gammas_hat_", penalty)
  name2 = paste0("covYI_KS_", penalty)

  results <-list("t2" = gammas1,
                 "t3" = getElement(covYI_KS_value, paste0("covYI_KS_", penalty)),
                 gr_yi = covYI_KS_value$gr_yi,
                 tau=tau,
                 penalty = penalty,
                 gammas_start = gammas_start,
                 kernel=kernel,
                 c_hat=covYI_KS_value$c_hat,
                 z_hat=covYI_KS_value$z_hat,
                 y_hat=covYI_KS_value$y_hat,
                 youden_index= covYI_KS_value$youden_index,
                 sensitivity=covYI_KS_value$sensitivity,
                 specificity=covYI_KS_value$specificity,
                 geometric_mean=covYI_KS_value$geometric_mean,
                 fdr=covYI_KS_value$fdr,
                 mcc=covYI_KS_value$mcc,
                 auc=covYI_KS_value$auc,
                 aauc=covYI_KS_value$aauc,
                 aYI=covYI_KS_value$aYI,
                 corrclass=covYI_KS_value$corrclass,
                 TP=covYI_KS_value$TP,
                 TN=covYI_KS_value$TN,
                 FP=covYI_KS_value$FP,
                 FN=covYI_KS_value$FN,
                 n_gammas=length(gammas1[which(gammas1!=0)]),
                 n_total_var_gammas=n_total_var_gammas, n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                 n_predicted_non_zeros_gammas = n_predicted_non_zeros_gammas, n_caught_gammas=n_caught_gammas,
                 n_non_caught_gammas=n_non_caught_gammas, n_caught_zero_gammas=n_caught_zero_gammas,
                 n_zero_not_caught_gammas=n_zero_not_caught_gammas)
  name_all = names(results)
  name_all[c(1,2)] = c(name1, name2)
  results <- stats::setNames(results, name_all)

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  #if (trace %in% c(1,2)){
  #  print(estimation_time)
  #}

  results = c(results , estimation_time = estimation_time, niter=estim$tot_iters)
  class(results) <- ("covYI_KS_estimation")
  return(results)
}

