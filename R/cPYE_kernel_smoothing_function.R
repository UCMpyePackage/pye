#install.packages("KernSmooth")
#library(KernSmooth)

#' create the output class of cPYE
# cPYE_output_KS <- setClass(Class="cPYE_output_KS",
#                        representation(
#                          cPYE_L12 ="numeric",
#                          cPYE_L1 ="numeric",
#                          cPYE_EN ="numeric",
#                          cPYE_SCAD ="numeric",
#                          cPYE_MCP="numeric",
#                          gr_cPYE_L12 ="numeric",
#                          gr_cPYE_L1 ="numeric",
#                          gr_cPYE_EN ="numeric",
#                          gr_cPYE_SCAD ="numeric",
#                          gr_cPYE_MCP ="numeric",
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

#' @title cPYE_KS
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
#' of just a vector of character. Default is all not present in y and C
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param betas the coefficients of the biomarker combination used for the
#' evaluation of PYE
#' @param gammas the coefficients of the covariates combination used for the
#' evaluation of the cut-point c (the first element is the constant!!)
#' @param lambda the penalization parameter of the regressors X
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param C covariate variables. It can be of type dataframe, containing the
#' same covariates included in df, or just a vector of character.
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
#' penalty <- "L12"
#' lambda <- 0.1
#' tau <- 0.1
#' betas <- rep(1, length(X))
#' gammas <- rep(1, length(C))
#' c <- 0
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' PYE_result <- cPYE_KS(df=df, X=X, C=C, betas=betas, gammas=gammas,
#' lambda=lambda, tau=tau, alpha=0.5, a1=3.7, a2=3, penalty=penalty)
#' print(PYE_result)
#'
#' @importFrom evmix kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom ROCnReg AROC.sp compute.threshold.AROC
#' @importFrom plyr join
#' @export
cPYE_KS <- function(df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C, betas, gammas, lambda, tau,
                    kernel="gaussian", alpha=0.5, a1=3.7, a2=3, penalty="L1", h_exponent=0.2, prediction=FALSE,
                    run_aauc=FALSE){

  #IMPORTANT: the first element of gammas is the constant term!!!!!!

  #check df:
  if (nrow(df)==0){stop("df has no rows")}
  
  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(C) != "character"){stop("C can only be of class character or data.frame.")}
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  #check if const already exists in the dataset
  #if("const" %in% colnames(df)){stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  #Here I dont create the constant since it should already exist, but I check if it needs to be created
  if(("const" %in% C) & (!("const" %in% names(df)))){
    const <- rep(1, length(ID)) #the constant for the coefficients gammas
    df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], const, df[,(names(df) %in% C)])
  } else {
    df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], df[,(names(df) %in% C)])
  }
  #put "const" in C
  #C1 <- c("const", C) #C Already contains the constant

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("cPYE_KS stopped because lambda needs to have length 1")}
  if (!is.numeric(lambda)){stop("cPYE_KS stopped because lambda needs be numeric")}
  if (length(tau)>1){stop("cPYE_KS stopped because tau needs to have length 1")}
  if (!is.numeric(tau)){stop("cPYE_KS stopped because tau needs be numeric")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(run_aauc %in% c("TRUE","FALSE"))){stop("run_aauc can be only TRUE or FALSE. Default is FALSE")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))){
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(betas) != length(X)){
    stop("The number of element of betas is different then the number of columns in X")
  }
  if (length(gammas) != length(C)){
    stop("The number of element of gammas is different then the number of columns in C")
  }

  #divide control and diseased and multiply by betas
  z_y0<- as.matrix(df1[df1$y==0, (names(df1) %in% X), drop=FALSE])%*%betas#drop=FALSE permits it to remain a data frame even if it is composed by a single column
  z_y1<- as.matrix(df1[df1$y==1, (names(df1) %in% X), drop=FALSE])%*%betas

  #compute c
  #IMPORTANT: the first element of gammas is the constant term!
  c_y0 <- as.matrix(df1[df1$y==0,(names(df1) %in% C), drop=FALSE])%*%gammas
  c_y1 <- as.matrix(df1[df1$y==1,(names(df1) %in% C), drop=FALSE])%*%gammas

  #kernel from the YI paper
  #healthy
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y0)!=0){
  #   h0 = kedd::h.bcv(as.numeric(c-z_y0), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
    h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*length(z_y0)^(-h_exponent)
  # }
  if(h0==0) {h0=0.1}
  t0 <- (c_y0-z_y0)/h0
  f0 <- sum(evmix::kpz(z = as.numeric(t0), kernel = kernel))/length(z_y0)
  #f0 <- sum(pnorm(t0, 0, 1))/length(z_y0)

  #diseased
  # if (kernel %in% c("gaussian","epanechnikov", "triweight","tricube","biweight","cosine") & sum(z_y1)!=0){
  #   h1= kedd::h.bcv(as.numeric(z_y1), kernel=kernel, deriv.order = 0)$h #using "kedd" package only those kernel are avail.
  # } else {
    h1=0.9*min(stats::sd(z_y1), stats::IQR(z_y1)/1.34)*length(z_y1)^(-h_exponent)
  # }
  if(h1==0) {h1=0.1} #cannot divide for 0!
  t1 <- (c_y1-z_y1)/h1
  f1 <- sum(evmix::kpz(z = as.numeric(t1), kernel = kernel))/length(z_y1)
  #f1 <- sum(pnorm(t1, 0, 1))/length(z_y1)

  yi <- f0 - f1

  #evaluate the gradient. NOTE: the variables to which consider the gradient
  #are betas and gammas elements
  #h0=0.9*min(stats::sd(z_y0), stats::IQR(z_y0)/1.34)*(length(z_y0)^(-h_exponent))
  t0 <- (c_y0-z_y0)/h0
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr0_betas <- apply (df1[df1$y==0,(names(df1) %in% X)], 2,
                       function(b) (t(stats::dnorm(t0, 0, 1)) %*% ((-b)/h0))/length(z_y0))
    gr0_gammas <- apply (df1[df1$y==0,(names(df1) %in% C)], 2,
                        function(g) (t(stats::dnorm(t0, 0, 1)) %*% (g/h0))/length(z_y0))
  } else {
    gr0_betas <- apply (df1[df1$y == 0,(names(df1) %in% X)], 2,
           function(b) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% ((-b)/h0))/length(z_y0))
    gr0_gammas <- apply (df1[df1$y == 0,(names(df1) %in% C)], 2,
            function(g) (t(evmix::kdz(z = as.numeric(t0), kernel = kernel)) %*% (g/h0))/length(z_y0))
  }

  #h1=0.9*min(stats::sd(z_y1), stats::IQR(z_y1)/1.34)*(length(z_y1)^(-h_exponent))
  t1= (c_y1-z_y1)/h1
  #since dnorm is way faster I leave it implemented for the gaussian kernel:
  if (kernel %in% c("gaussian", "normal")){
    gr1_betas <- apply (df1[df1$y == 1,(names(df1) %in% X)], 2,
                        function(b) (t(stats::dnorm(t1, 0, 1)) %*% ((-b)/h1))/length(z_y1))
    gr1_gammas <- apply (df1[df1$y == 1,(names(df1) %in% C)], 2,
                        function(g) (t(stats::dnorm(t1, 0, 1)) %*% (g/h1))/length(z_y1))
  } else {
    gr1_betas <- apply (df1[df1$y == 1,(names(df1) %in% X)], 2,
           function(b) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% ((-b)/h1))/length(z_y1))
    gr1_gammas <- apply (df1[df1$y == 1,(names(df1) %in% C)], 2,
           function(g) (t(evmix::kdz(z = as.numeric(t1), kernel = kernel)) %*% (g/h1))/length(z_y1))
  }

  #gradient
  gr_yi <- c(gr0_betas, gr0_gammas) - c(gr1_betas, gr1_gammas) #note that in gammas we are not inluding the constant here

  #append data (z)
  zx<-as.data.frame(rbind(z_y0, z_y1))
  zx['ID']=as.numeric(rownames(zx))
  #zx <- zx[order(zx$ID), ]
  names(zx)[names(zx) == 'V1'] <- 'zx_hat'
  df1<-plyr::join(df1, zx, by='ID')
  #df1<-merge(df1, zx, by='ID', all = TRUE)
  rownames(df1) <- df1$ID

  #append data (c)
  c<-as.data.frame(rbind(c_y0, c_y1))
  c['ID']=as.numeric(rownames(c))
  #c <- c[order(c$ID), ]
  names(c)[names(c) == 'V1'] <- 'c_hat'
  df1<-plyr::join(df1, c, by='ID')
  #df1<-merge(df1, c, by='ID', all = TRUE)
  rownames(df1) <- df1$ID

  #test: sum zx and c
  #df1["z_hat"] = df1$zx_hat - df1$c_hat
  #or, normally:
  df1["z_hat"] <- df1$zx_hat

  #find the optimal cut-point
  opt <- OptimalCutpoints::optimal.cutpoints(data= df1, X = "z_hat", status = "y", methods = "Youden", tag.healthy =0)

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
  #try({test <- AUC::auc(AUC::roc(df1$zx_hat,factor(df1$y)))}, silent = TRUE)
  auc <- opt$Youden$Global$measures.acc$AUC[1]
  yi1 <- opt$Youden$Global$optimal.criterion

  #df_global <<- df1 # "<<-" is to save in the global environment

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.sp(formula.healthy = zx_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                      p = seq(0,1,l=101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- try({aroc <- AROC::AROC.kernel(marker = "zx_hat", covariate = "c_hat", group = "y", tag.healthy = 0,
  #                                     data = df1, p = seq(0,1,l=101), B = 500)}, silent = TRUE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.bnp(formula.h = zx_hat ~ c_hat, group = "y", tag.h = 0, data = df1, p = seq(0,1,l=101),
  #                           compute.lpml = FALSE, compute.WAIC = FALSE)

  #the following is too long:
  #tmp <- ROCnReg::AROC.kernel(marker = "zx_hat", covariate = "c_hat", group = "y", tag.h = 0, data = df1,
  #                              p = seq(0,1,l=101), B = 500)

  #start_time2 <- Sys.time()
  #this is longher then AROC.bsp
  #tmp2 <- try({aroc <- AROC::AROC.bnp(formula.healthy = zx_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                   scale = TRUE, p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                   a = 2, b = 0.5, L = 10, nsim = 5000, nburn = 1000)}, silent = TRUE)
  #cat("temp2:")
  #print(tmp2)
  #print(difftime(Sys.time(), start_time2 , units = "mins"))
  #start_time1 <- Sys.time()

  #this does not have the YI
  #aroc <- AROC::AROC.bsp(formula.healthy = zx_hat ~ c_hat, group = "y", tag.healthy = 0, data = df1,
  #                                scale = TRUE, p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
  #                                a = 2, b = 0.5, nsim = 5000, nburn = 1500)

  if (run_aauc==TRUE){
    tmp <- try({

      #AROC.sp function estimates the covariate-adjusted ROC curve (AROC) using the semiparametric approach
      #proposed by Janes and Pepe (2009).

      #if we have all zeros, the function returns a warning message. We suppress it just for this function:
      options(warn=-1)

      aroc <- ROCnReg::AROC.sp(formula.h = zx_hat ~ c_hat, group = "y", tag.h = 0, data = df1,
                              p = seq(0,1,l=101), B = 500)

      options(warn=0)

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
    aauc <- 0
    aYI <- 0
  }

  #L_(1/2) penalization
  phi_L12_b <- lambda*sum(abs(betas)^(1/2))#NB:we DO NOT add the intercept!!!
  phi_L12_g <- tau*sum(abs(gammas[2:length(gammas)])^(1/2))#NB:we DO HAVE the intercept!!!
  if(is.na(phi_L12_b)){
    stop('problem with the penalization phi_L12_b in cPYE function')}
  if(is.na(phi_L12_g)){
    stop('problem with the penalization phi_L12_g in cPYE function')}


  #L1 penalization
  phi_L1_b <- lambda*sum(abs(betas))#NB:we DO NOT add the intercept!!!
  phi_L1_g <- tau*sum(abs(gammas[2:length(gammas)]))#NB:we DO HAVE the intercept!!!
  if(is.na(phi_L1_b)){
    stop('problem with the penalization phi_L1_b in cPYE function')}
  if(is.na(phi_L1_g)){
    stop('problem with the penalization phi_L1_g in cPYE function')}

  #Elastic-Net penalization
  phi_EN_b <- lambda*((alpha)*(sum(abs(betas))) + ((1-alpha)/2)*(sum(betas^2)))
  phi_EN_g <- tau*((alpha)*(sum(abs(gammas[2:length(gammas)]))) + ((1-alpha)/2)*(sum(gammas[2:length(gammas)]^2)))
  if(is.na(phi_EN_b)){
    stop('problem with the penalization phi_EN_b in cPYE function')}
  if(is.na(phi_EN_g)){
    stop('problem with the penalization phi_EN_g in cPYE function')}

  #SCAD penalization
  phi_SCAD_b <- SCAD_function(betas, lambda, a=a1)
  phi_SCAD_g <- SCAD_function(gammas[2:length(gammas)], tau, a=a1)
  if(is.na(phi_SCAD_b)){
    stop('problem with the penalization phi_SCAD_b in cPYE function')}
  if(is.na(phi_SCAD_g)){
    stop('problem with the penalization phi_SCAD_g in cPYE function')}

  #MCP penalization
  phi_MCP_b <- MCP_function(betas, lambda, a=a2)
  phi_MCP_g <- MCP_function(gammas[2:length(gammas)], tau, a=a2)
  if(is.na(phi_MCP_b)){
    stop('problem with the penalization phi_MCP_b in cPYE function')}
  if(is.na(phi_MCP_g)){
    stop('problem with the penalization phi_MCP_g in cPYE function')}

  #join yi and the penalty functions
  cPYE_L12 <- yi - phi_L12_b - phi_L12_g
  cPYE_L1 <- yi - phi_L1_b - phi_L1_g
  cPYE_EN <- yi - phi_EN_b - phi_EN_g
  cPYE_SCAD <- yi - phi_SCAD_b - phi_SCAD_g
  cPYE_MCP <- yi - phi_MCP_b - phi_MCP_g

  #if we are in prediction this is the right yi
  if (prediction==TRUE){
    yi=yi1
  }

  return(list("cPYE_L12"=cPYE_L12, "cPYE_L1"=cPYE_L1, "cPYE_EN" = cPYE_EN, "cPYE_SCAD" = cPYE_SCAD,
              "cPYE_MCP" = cPYE_MCP, "gr_yi"=gr_yi, "youden_index"=yi, "sensitivity"=sensitivity,
              "specificity"=spec, "geometric_mean"=gm, "fdr"=fdr, "mcc"=mcc, "auc"=auc, "aauc"=aauc,
              "aYI"=aYI, "corrclass"=corrclass, "zx_hat"=df1[,c("ID","zx_hat")],
              "c_hat"=df1[,c("ID","c_hat")], "z_hat"=df1[,c("ID","z_hat")],
              "y_hat"=df1[,c("ID","y_hat")], "TP"=TP, "TN"=TN, "FP"=FP, "FN"=FN,
              "input_data" = list("beta"=betas, "gammas"=gammas, "lambda"=lambda,
                                  "tau"=tau, "alpha"=alpha, "a1"=a1, "a2"=a2, "mode"=mode, "kernel"=kernel)
              ))
}





#' @title cPYE_KS_estimation
#'
#' @description function to estimate the optimal value of betas and c maximizing
#' the PYE function. To find the optimum the mmAPG and mnmAPG algorithms are
#' used.
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y and C
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param C covariate variables. It can be of type dataframe, containing the
#' same covariates included in df, or just a vector of character.
#' @param lambda the penalization parameter of the regressors X
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If "zeros", it
#' starts with all zero values, if "corr" it starts with the value of the
#' correlation of every regressor with the target variable. Default is "zeros"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point of gammas. If "zeros", it
#' starts with all zero values, if "corr" it starts with the value of the
#' correlation of every regressor with the target variable. Default is "zeros"
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
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
#' regressors_gammas<-simMicroarrayData_cov02_dim50_covariates$ncovariates
#' penalty <- "SCAD"
#' lambda <- 0.1
#' tau <- 0.1
#' c <- 0
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#' param_start_default <- "zeros"
#' alpha <- 0.5
#'
#' cPYE_estimation_result <- cPYE_KS_estimation(df=df, X=X, y=y, C=C,
#'   penalty=penalty, trace=2,
#'   beta_start_default=param_start_default,
#'	 gamma_start_default=param_start_default,
#'   lambda=lambda, tau=tau, alpha=alpha, regressors_betas=regressors_betas,
#'   regressors_gammas=regressors_gammas)
#'
#' print(cPYE_estimation_result)
#'
#' @importFrom stats setNames
#' @export
#Estimation of the parameter using cPYE_KS
cPYE_KS_estimation <- function(df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C, lambda, tau, penalty="L1",
                              beta_start_input=NULL, beta_start_default="zeros",
							  gamma_start_input=NULL, gamma_start_default="zeros",
                              alpha=0.5, a1=3.7, a2=3, regressors_betas=NULL, regressors_gammas=NULL,
                              fold=NULL, max_iter=1000, max.print=10,
                              trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                              min_alpha=1e-10, convergence_error=1e-7,
                              trace=1, seed=1, kernel="gaussian", run_aauc=FALSE){


  start_time <- Sys.time()

  options(max.print = max.print)
  
  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(C) != "character"){stop("C can only be of class character or data.frame.")}
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  #check if const already exists in the dataset
  if("const" %in% colnames(df)){stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal const")}
  if("const" %in% colnames(C)){stop("const already exists as column in C! Please delete or rename this column since I need this name to set the internal const")}

  ID <- rownames(df)
  const <- rep(1, length(ID)) #the constant for the coefficients gammas
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], const, df[,(names(df) %in% C)])
  #put "const" in C
  C1 <- c("const", C)

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("cPYE_KS_estimation stopped because lambda needs to have length 1")}
  if (length(tau)>1){stop("cPYE_KS_estimation stopped because tau needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}


  betas1_initial_zeros <- rep(0, length(X))
  betas1_initial_corr <- NULL

  gammas1_initial_zeros <- rep(0, length(C1))
  gammas1_initial_corr <- NULL

  names_betas <- X
  names_gammas <- C1

  #initializing the parameters (betas)
  if (length(beta_start_input)==0){
    if (beta_start_default=="zeros"){

      betas_start <- betas1_initial_zeros
      names(betas_start) <- names_betas
      betas1 <- betas_start

    } else if (beta_start_default=="corr") {

      #compute the corr between every x and y
      corr.xy <- data.frame(matrix(ncol = ncol(df1[,(names(df1) %in% X)]), nrow = 0))
      names(corr.xy) <- names_betas
      corr.xy[1:3,] <- t(sapply(c("pearson", "kendall", "spearman"),
                                function (h) unlist(lapply(X, function (x) stats::cor(df1$y, df1[,x],  method = h)))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      betas1_initial_corr <- as.numeric(mean.corr.xy)
      betas_start <- betas1_initial_corr
      names(betas_start) <- names_betas
      betas1 <- betas_start

    } else {stop("The parameter beta_start_default can only be equal to 'zeros', 'corr'.")}
  } else {
    if (length(beta_start_input) == ncol(df1[,(names(df1) %in% X)])){
      betas_start <- beta_start_input
      names(betas_start) <- names_betas
      betas1 <- betas_start
    } else if (length(beta_start_input) != 0){
      stop("The length of the parameter beta_start_input is not equal to ", ncol(df1[,(names(df1) %in% X)]), ". Remeber that the last parameter it has to include is c.")
    }
  }

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
                                function (h) c(0,unlist(lapply(C1[-1], function (x) stats::cor(df1$y, df1[,(names(df1) %in% C1)][,x],  method = h))))))
      row.names(corr.xy) <- c("pearson", "kendall", "spearman")
      #mean of the corrs and sort names
      mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
      gammas1_initial_corr <- c(0,as.numeric(mean.corr.xy)) #the initial 0 is the constant
      gammas_start <- gammas1_initial_corr
      names(gammas_start) <- names_gammas
      gammas1 <- gammas_start

    } else {stop("The parameter gamma_start_default is defines. It can be only equal to 'zeros', 'corr'.")}
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
  prox_penalty <- get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  #wrappers
  delta_fx_betas <- function(betas, gammas){
    result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X[X %in% names(betas)], y=y, C=C1[C1 %in% names(gammas)], betas=betas, gammas=gammas, lambda=lambda,
                       tau=tau, kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty, run_aauc=run_aauc)$gr_yi
    result <- result[names(result) %in% X]
    return(result)
  }

  delta_fx_gammas <- function(betas, gammas){
    result <- -cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X[X %in% names(betas)], y=y, C=C1[C1 %in% names(gammas)], betas=betas, gammas=gammas, lambda=lambda,
                       tau=tau, kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty, run_aauc=run_aauc)$gr_yi
    result <- result[names(result) %in% C1]
    return(result)
  }

  proxx_betas <- function(x, eta){
    result <- prox_penalty(betas=x, lambda=eta*lambda, alpha=alpha, a=a)
    return(result)
  }

  proxx_gammas <- function(x, eta){
    #old version:
    #result <- c(x[(names(x) == "const")], prox_penalty(betas=x[!(names(x) == "const")], lambda=eta*tau, alpha=alpha, a=a)) #-1 since I dont penalize the constant
    #new version (let c enter in the proximal operator)
    result <- prox_penalty(betas=x, lambda=eta*tau, alpha=alpha, a=a) #deleted the -"const" since I dont penalize the constant
    return(result)
  }

  Fx <- function(x,z){
    result <- -getElement(cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X[X %in% names(x)], y=y, C=C1[C1 %in% names(z)], betas=x, gammas=z, lambda=lambda, tau=tau,
                            kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty, run_aauc=run_aauc), paste0("cPYE_", penalty))
    return(result)
  }

  if (trend == "monotone"){
    estim <- mmAPG2(x0=betas1, g0=gammas1, delta_fx_betas=delta_fx_betas, delta_fx_gammas=delta_fx_gammas,
                    proxx_betas=proxx_betas, proxx_gammas=proxx_gammas, Fx=Fx,
                    lambda=lambda, tau=tau, penalty=penalty, fold=fold, stepsizeShrink=stepsizeShrink, delta=delta,
                    max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                    max_alpha= max_alpha, trace=trace, seed=seed)
  } else if (trend == "nonmonotone"){
    estim <- mnmAPG2(x0=betas1, g0=gammas1, delta_fx_betas=delta_fx_betas, delta_fx_gammas=delta_fx_gammas,
                     proxx_betas=proxx_betas, proxx_gammas=proxx_gammas, Fx=Fx,
                     lambda=lambda, tau=tau, penalty=penalty, fold=fold, stepsizeShrink=stepsizeShrink, delta=delta,
                     max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                     max_alpha= max_alpha, trace=trace, seed=seed)
  }

  betas1 <- estim$x1
  gammas1 <- estim$g1

  cPYE_KS_value <- cPYE_KS(df=df1[,!(names(df1) %in% c("ID"))], X=X, y=y, C=C1, betas=betas1,
                           gammas=gammas1, lambda=lambda, tau=tau, kernel=kernel,
                           alpha=alpha, a1=a1, a2=a2, penalty=penalty, run_aauc=run_aauc)

  if (trace %in% c(1,2)){
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend,"version. ; \n")
    if (length(fold)!=0) {cat("fold =", fold, "; ")}
    visualize_betas <- c(betas1[which(betas1!=0)])
    visualize_gammas <- c(gammas1[which(gammas1!=0)])
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters , "; lambda:", lambda, "; tau:", tau, "; penalty:", penalty, "; cPYE_KS:", getElement(cPYE_KS_value,  paste0("cPYE_",penalty)), "; youden_index:", cPYE_KS_value$youden_index, "; aYI:", cPYE_KS_value$aYI, "; sensitivity:", cPYE_KS_value$sensitivity, "; specificity:", cPYE_KS_value$specificity, "; geometric_mean:", cPYE_KS_value$geometric_mean, "; fdr:", cPYE_KS_value$fdr, "; mcc:", cPYE_KS_value$mcc, "; auc:", cPYE_KS_value$auc, "; aauc:", cPYE_KS_value$aauc, "; corrclass:", cPYE_KS_value$corrclass, "; \n")
    cat("TP:", cPYE_KS_value$TP, "; TN:", cPYE_KS_value$TN, "; FP:", cPYE_KS_value$FP, "; FN:", cPYE_KS_value$FN, ";  betas: \n")
    print(visualize_betas)
    cat("\n gammas: \n")
    print(visualize_gammas)
    cat("\n")
  }

  #total number of variables
  n_total_var_betas= length(betas1)
  n_total_var_gammas= length(gammas1)

  #n of predicted zeros
  n_predicted_zeros_betas= sum(betas1==0)
  n_predicted_zeros_gammas= sum(gammas1==0)

  #n of predicted betas different from 0
  n_predicted_non_zeros_betas = sum(betas1!=0)
  n_predicted_non_zeros_gammas = sum(gammas1!=0)

  #compute other measures
  if (length(regressors_gammas)!=0){

    #n of beta catch
    n_caught_betas = sum(betas1[regressors_betas!=0]!=0)
    n_caught_gammas = sum(gammas1[regressors_gammas!=0]!=0)

    #n of beta not catch
    n_non_caught_betas = sum(betas1[regressors_betas!=0]==0)
    n_non_caught_gammas = sum(gammas1[regressors_gammas!=0]==0)

    #n? of zeros catch
    n_caught_zero_betas = sum(betas1[regressors_betas==0]==0)
    n_caught_zero_gammas = sum(gammas1[regressors_gammas==0]==0)

    #n of zeros not catch
    n_zero_not_caught_betas = sum(betas1[regressors_betas==0]!=0)
    n_zero_not_caught_gammas = sum(gammas1[regressors_gammas==0]!=0)
  } else {
    n_caught_betas=NA
    n_caught_gammas=NA
    n_non_caught_betas=NA
    n_non_caught_gammas=NA
    n_caught_zero_betas=NA
    n_caught_zero_gammas=NA
    n_zero_not_caught_betas=NA
    n_zero_not_caught_gammas=NA
  }

  name1 = paste0("betas_hat_", penalty)
  name2 = paste0("gammas_hat_", penalty)
  name3 = paste0("cPYE_KS_", penalty)
  name4 = paste0("gr_cPYE_KS_", penalty)

  results <-list("t1" = betas1,
                 "t2" = gammas1,
                 "t3" = getElement(cPYE_KS_value, paste0("cPYE_", penalty)),
                 "grad" = getElement(cPYE_KS_value, paste0("gr_cPYE_", penalty)),
                 lambda=lambda,
                 tau=tau,
                 penalty = penalty,
                 betas_start = betas_start,
                 gammas_start = gammas_start,
                 kernel=kernel,
                 c_hat=cPYE_KS_value$c_hat,
                 z_hat=cPYE_KS_value$z_hat,
                 y_hat=cPYE_KS_value$y_hat,
                 youden_index= cPYE_KS_value$youden_index,
                 sensitivity=cPYE_KS_value$sensitivity,
                 specificity=cPYE_KS_value$specificity,
                 geometric_mean=cPYE_KS_value$geometric_mean,
                 fdr=cPYE_KS_value$fdr,
                 mcc=cPYE_KS_value$mcc,
                 auc=cPYE_KS_value$auc,
                 aauc=cPYE_KS_value$aauc,
                 aYI=cPYE_KS_value$aYI,
                 corrclass=cPYE_KS_value$corrclass,
                 n_betas=length(betas1[which(betas1!=0)]),
                 n_gammas=length(gammas1[which(gammas1!=0)]),
                 n_total_var_betas=n_total_var_betas, n_predicted_zeros_betas=n_predicted_zeros_betas,
                 n_predicted_non_zeros_betas = n_predicted_non_zeros_betas, n_caught_betas=n_caught_betas,
                 n_non_caught_betas=n_non_caught_betas, n_caught_zero_betas=n_caught_zero_betas,
                 n_zero_not_caught_betas=n_zero_not_caught_betas,
                 n_total_var_gammas=n_total_var_gammas, n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                 n_predicted_non_zeros_gammas = n_predicted_non_zeros_gammas, n_caught_gammas=n_caught_gammas,
                 n_non_caught_gammas=n_non_caught_gammas, n_caught_zero_gammas=n_caught_zero_gammas,
                 n_zero_not_caught_gammas=n_zero_not_caught_gammas)
  name_all = names(results)
  name_all[c(1,2,3,4)] = c(name1, name2, name3, name4)
  results <- stats::setNames(results, name_all)

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
  }

  results = c(results , estimation_time = estimation_time, niter=estim$tot_iters)
  class(results) <- ("cPYE_KS_estimation")
  return(results)
}










#---------------> cross-validation of cPYE KS <----------------

#create the output class of the cPYE.cv function
setClass(Class="cPYE_KS_cross_validation_output",
         representation(
           penalty="character",
           kernel="character",
           cv_time="ANY",
           cPYE_KS_L12= "ANY",
           cPYE_KS_L1= "ANY",
           cPYE_KS_EN= "ANY",
           cPYE_KS_SCAD= "ANY",
           cPYE_KS_MCP= "ANY",
           auc = "ANY",
           aauc = "ANY",
           aYI = "ANY",
           youden_index ="ANY",
           sensitivity ="ANY",
           specificity ="ANY",
           geometric_mean ="ANY",
           fdr ="ANY",
           mcc="ANY",
           corrclass ="ANY",
           n_betas="ANY",
           n_gammas="ANY",
           betas="list",
           gammas="list"
         )
)

#to extract just part of the estim of cPYE
subset_cPYE_KS <- function(df, X, y, C1, betas, gammas, lambda, tau, fold, alpha, a1, a2, penalty, cv_time, niter, kernel, run_aauc){
  cPYE_KS_result <- cPYE_KS(df=df[,names(df)!=c("ID")], X=X, y=y, C=C1, betas=betas, gammas=gammas, lambda=lambda, tau=tau, alpha=alpha, a1=a1, a2=a2,
                           penalty=penalty, prediction=TRUE, kernel=kernel, run_aauc=run_aauc)

  cat("-> Results on the TEST SET \n")
  cat("-> algorithm: cPYE_KS_proximal_gradient_method ; ")
  if (length(fold)!=0) {cat("fold =", fold, "; ")}
  visualize_betas <- betas[which(betas!=0)]
  visualize_gammas <- gammas[which(gammas!=0)]
  cat("lambda:", lambda, "tau:", tau, "; penalty:", penalty, "; cPYE_KS:", getElement(cPYE_KS_result,  paste0("cPYE_",penalty)), "; youden_index:", cPYE_KS_result$youden_index, "; aYI:", cPYE_KS_result$aYI, "; sensitivity:", cPYE_KS_result$sensitivity, "; specificity:", cPYE_KS_result$specificity, "; geometric_mean:", cPYE_KS_result$geometric_mean, "; fdr:", cPYE_KS_result$fdr, "; mcc:", cPYE_KS_result$mcc, "; auc:", cPYE_KS_result$auc, "; aauc:", cPYE_KS_result$aauc, "; corrclass:", cPYE_KS_result$corrclass, " \n")
  cat("TP:", cPYE_KS_result$TP, "; TN:", cPYE_KS_result$TN, "; FP:", cPYE_KS_result$FP, "; FN:", cPYE_KS_result$FN, ";  betas: \n")
  print(visualize_betas)
  cat("\n gammas: \n")
  print(visualize_gammas)
  cat("Cross-validation time:", cv_time, "; Number of iterations:", niter)
  cat("\n")

  return(list("cPYE_KS_L12" =cPYE_KS_result$cPYE_L12, "cPYE_KS_L1"=cPYE_KS_result$cPYE_L1, "cPYE_KS_EN"=cPYE_KS_result$cPYE_EN, "cPYE_KS_SCAD" = cPYE_KS_result$cPYE_SCAD, "cPYE_KS_MCP" =cPYE_KS_result$cPYE_MCP, "youden_index"=cPYE_KS_result$youden_index, "sensitivity"=cPYE_KS_result$sensitivity, "specificity"=cPYE_KS_result$specificity, "geometric_mean"=cPYE_KS_result$geometric_mean, "fdr"=cPYE_KS_result$fdr, "mcc"=cPYE_KS_result$mcc, "auc"=cPYE_KS_result$auc, "aauc"=cPYE_KS_result$aauc, "aYI"=cPYE_KS_result$aYI, "corrclass"=cPYE_KS_result$corrclass))
}

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
cPYE_KS.cv <- function (df, X, y, C, lambda_tau, trace=1, alpha, a1, a2, penalty, folds_i, k, regressors_betas=NULL,
                        regressors_gammas=NULL, cPYE_KS_L12, cPYE_KS_L1, cPYE_KS_EN, cPYE_KS_SCAD, cPYE_KS_MCP, auc,
                        aauc, aYI, youden_index, sensitivity, specificity, geometric_mean,
                        fdr, mcc, corrclass, n_betas, n_gammas, used_cores, trend, delta, max_alpha, kernel,
                        max_iter, beta_start_input, beta_start_default, gamma_start_input, gamma_start_default, min_alpha,
						convergence_error, stepsizeShrink, run_aauc){

  #put "const" in C
  C1 <- c("const", C)

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  cat(" ---------------------------------------------------------------\n |  starting with the", k ,"-th fold for the CV of lambda and tau  |\n ---------------------------------------------------------------\n")

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
           ") is larger than the number of avaiable cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_CV_cPYE_KS.txt")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for cross-validation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("train_df", "X", "y", "C", "trace", "alpha", "penalty", "regressors_betas",
                                  "regressors_gammas", "k", "trend", "delta", "max_alpha", "kernel", "max_iter",
                                  "beta_start_input", "beta_start_default", "gamma_start_input", "gamma_start_default",
								  "a2", "a1", "min_alpha", "convergence_error", "stepsizeShrink"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))
    fitted_models <- parallel::parLapply(cl, lambda_tau, function(x) cPYE_KS_estimation(df=train_df[,names(train_df) != "const"],
                                                            X=X, y=y, C=C, lambda=x[1], tau=x[2],
                                                            beta_start_input=beta_start_input,
															beta_start_default=beta_start_default,
                                                            gamma_start_input=gamma_start_input,
                                                            gamma_start_default=gamma_start_default, trace=trace,
                                                            alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                                            regressors_betas=regressors_betas, max_iter=max_iter,
                                                            min_alpha=min_alpha, convergence_error=convergence_error,
                                                            regressors_gammas=regressors_gammas, fold=k, trend=trend,
                                                            stepsizeShrink=stepsizeShrink,
                                                            delta=delta, max_alpha=max_alpha, kernel=kernel,
															run_aauc=run_aauc))
    #measures on the train set
    temp_cPYE = get(paste0("cPYE_KS_", penalty))
    parallel::clusterExport(cl, list("fitted_models", paste0("cPYE_KS_", penalty)), envir = environment())
    temp_cPYE$train[k, ] <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("cPYE_KS_", penalty))))
    auc$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "auc")))
    aauc$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "aauc")))
    aYI$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "aYI")))
    youden_index$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "youden_index")))
    sensitivity$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "sensitivity")))
    specificity$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "specificity")))
    geometric_mean$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:length(lambda_tau)), function(x) getElement(fitted_models[[x]], "geometric_mean")))
    fdr$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "fdr")))
    mcc$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "mcc")))
    corrclass$train[k, ]  <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "corrclass")))
    n_betas[k, ] <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "n_betas")))
    n_gammas[k, ] <- unlist(parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "n_gammas")))

    #create a list of the results
    betas_hat <- parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("betas_hat_", penalty)))
    gammas_hat <- parallel::parLapply(cl, c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("gammas_hat_", penalty)))

    on.exit(parallel::stopCluster(cl))

  } else {

    fitted_models <- lapply(lambda_tau, function(x) cPYE_KS_estimation(df=train_df[,names(train_df) != "const"], X=X, y=y, C=C,
                                                       lambda=x[1], tau=x[2], beta_start_input=beta_start_input,
													   beta_start_default=beta_start_default,
                                                       gamma_start_input=gamma_start_input,
                                                       gamma_start_default=gamma_start_default, trace=trace,
                                                       alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                                       min_alpha=min_alpha, convergence_error=convergence_error,
                                                       stepsizeShrink=stepsizeShrink,
                                                       regressors_betas=regressors_betas, max_iter=max_iter,
                                                       regressors_gammas=regressors_gammas, fold=k, trend=trend,
                                                       delta=delta, max_alpha=max_alpha, kernel=kernel,
													   run_aauc=run_aauc))

    #measures on the train set
    temp_cPYE = get(paste0("cPYE_KS_", penalty))
    temp_cPYE$train[k, ] <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("cPYE_KS_", penalty))))
    auc$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "auc")))
    aauc$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "aauc")))
    aYI$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "aYI")))
    youden_index$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "youden_index")))
    sensitivity$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "sensitivity")))
    specificity$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "specificity")))
    geometric_mean$train[k, ]  <- unlist(lapply(c(1:length(lambda_tau)), function(x) getElement(fitted_models[[x]], "geometric_mean")))
    fdr$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "fdr")))
    mcc$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "mcc")))
    corrclass$train[k, ]  <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "corrclass")))
    n_betas[k, ] <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "n_betas")))
    n_gammas[k, ] <- unlist(lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], "n_gammas")))

    #create a list of the results
    betas_hat <- lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("betas_hat_", penalty)))
    gammas_hat <- lapply(c(1:(length(lambda_tau))), function(x) getElement(fitted_models[[x]], paste0("gammas_hat_", penalty)))

  }

  betas <- list(mapply(function(x) list(x), betas_hat))
  names(betas[[1]]) <- c(unlist(lapply(lambda_tau[1], function (x) paste0(paste("lambda", x, sep="="),"_", paste("tau", mean(lambda_tau[2]), sep="=")))), unlist(lapply(lambda_tau[2], function (x) paste0(paste("lambda", mean(lambda_tau[1]), sep="="),"_", paste("tau", x, sep="=")))))
  gammas <- list(mapply(function(x) list(x), gammas_hat))
  names(gammas[[1]]) <- c(unlist(lapply(lambda_tau[1], function (x) paste0(paste("lambda", x, sep="="),"_", paste("tau", mean(lambda_tau[2]), sep="=")))), unlist(lapply(lambda_tau[2], function (x) paste0(paste("lambda", mean(lambda_tau[1]), sep="="),"_", paste("tau", x, sep="=")))))

  #measures on the test set
  all_measures_test <- mapply(function(xx,xy,z,zz) subset_cPYE_KS(df=test_df, X=X, y=y, C1=C1,
                                                                  betas=xx, gammas=xy, lambda=z[1], tau=z[2],
                                                                  fold = k, alpha=alpha, a1=a1, a2=a2,
                                                                  penalty=penalty,
                                                                  cv_time=fitted_models[[zz]]$estimation_time,
                                                                  niter=fitted_models[[zz]]$niter, kernel=kernel,
																  run_aauc=run_aauc),
                              betas_hat, gammas_hat, lambda_tau, 1:length(lambda_tau))

  temp_cPYE$test[k, ] <- unlist(all_measures_test[paste0("cPYE_KS_", penalty),])
  assign(paste0("cPYE_KS_", penalty), temp_cPYE)
  auc$test[k, ] <- unlist(all_measures_test["auc",])
  aauc$test[k, ] <- unlist(all_measures_test["aauc",])
  aYI$test[k, ] <- unlist(all_measures_test["aYI",])
  youden_index$test[k, ] <- unlist(all_measures_test["youden_index",])
  sensitivity$test[k, ] <- unlist(all_measures_test["sensitivity",])
  specificity$test[k, ] <- unlist(all_measures_test["specificity",])
  geometric_mean$test[k, ] <- unlist(all_measures_test["geometric_mean",])
  fdr$test[k, ] <- unlist(all_measures_test["fdr",])
  mcc$test[k, ] <- unlist(all_measures_test["mcc",])
  corrclass$test[k, ] <- unlist(all_measures_test["corrclass",])

  return(new("cPYE_KS_cross_validation_output", "penalty"=penalty, "kernel"=kernel, "cPYE_KS_L12"=cPYE_KS_L12,
             "cPYE_KS_L1"=cPYE_KS_L1, "cPYE_KS_EN"=cPYE_KS_EN, "cPYE_KS_SCAD"=cPYE_KS_SCAD,
             "cPYE_KS_MCP"=cPYE_KS_MCP, "auc"=auc, "aauc"=aauc,
             "aYI"=aYI, "youden_index"=youden_index,
             "sensitivity"=sensitivity, "specificity"=specificity,
             "geometric_mean"=geometric_mean, "fdr"=fdr,
             "mcc"=mcc, "corrclass"=corrclass, "n_betas"=n_betas,
             "betas"=betas, "n_gammas"=n_gammas, "gammas"=gammas))
}


#' @title cPYE_KS_compute_cv
#'
#' @description function to perform the cross-validation to select the best
#' value of lambda and tau for the estimation of betas and c using PYE.
#'
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
#' same covariates included in df, or just a vector of character.
#' @param lambda the penalization parameter of the regressors X
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param regressors_gammas a vector containing the real gammas (if known).
#' Default is NULL
#' @param seed fix the seed. Default is 1
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param trend if "monotone", mmAPG is used, if "nonmonotone", mnmAPG is used.
#' Default is "monotone"
#' @param delta parameter for the convergence condition of the optimization
#' algorithm. Default is 1e-5
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 10000
#' @param max_alpha maximum value of the step-parameter alpha. Default is 1000
#' @param min_alpha minimum value of the step-parameter alpha. Default is 1e-12
#' @param convergence_error error to accept for considering the algorithm
#' converged. Default is 1e-5
#' @param stepsizeShrink parameter to adjust the step-size in the backtracking
#' line-search, in the optimization of pye. Taking values between 0 and 1,
#' the closer to 1, the more accurate the estimation will be, the longer it
#' will take and viceversa. Default is 0.8
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If "zeros", it
#' starts with all zero values, if "corr" it starts with the value of the
#' correlation of every regressor with the target variable. Default is "zeros"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point for gammas. If "zeros", it
#' starts with all zero values, if "corr" it starts with the value of the
#' correlation of every regressor with the target variable. Default is "zeros"
#' @param scaling if TRUE, the dataset is scaled. FALSE otherwise, Default is
#' FALSE.
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
#'
#' @return a list containing the optimal value of lambda and tau to estimate betas and
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
#' betas <- rep(1, length(X))
#' c <- 0
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#' param_start_default <- "zeros"
#' alpha <- 0.5
#' tau_star <- 0.5 #a starting value of tau, to calibrate lambda
#' n_folds <- 3
#' used_cores <- 1
#' max_iter <- 5
#' used_penalty_cPYE <- c("L12", "SCAD") #c("L12", "L1", "EN", "SCAD", "MCP")
#'
#' #cPYE Gaussian (and others) Kernel Smooth
#' for (p in used_penalty_cPYE){
#'
#'   name <- paste0("lambda_estimate_cPYE_KS_", p)
#'
#'   #wrapper of the function
#'   wrapper_lambda <- function(lambda){
#'     return(cPYE_KS_estimation(df=df, X=X, y=y, C=C, penalty=p, trace=2,
#'	   beta_start_default=param_start_default,
#'     gamma_start_default=param_start_default, lambda=lambda, tau=tau_star,
#'     alpha=alpha, regressors_betas=regressors_betas,
#'     regressors_gammas=regressors_gammas, max_iter=max_iter, run_aauc=FALSE))
#'   }
#'   #lambda to max and min
#'   lambda_max <- calibrate_lambda_max (function_to_run=wrapper_lambda,
#'     var_to_check=paste0("betas_hat_",p), lambda_start = 5)
#'   lambda_min <- calibrate_lambda_min (function_to_run=wrapper_lambda,
#'     var_to_check=paste0("betas_hat_",p), lambda_start = 0.0001, max_var = 100)
#'
#'   #create a suited lambda
#'   lambda <- create_lambda(n=3, lmax=lambda_max, lmin=lambda_min)
#'   lambda <- as.numeric(formatC(lambda, format = "e", digits = 9))
#'
#'   #create a suited tau
#'	 lambda_star <- (lambda_max+lambda_min)/2
#'   wrapper_tau <- function(tau){
#'     return(cPYE_KS_estimation(df=df, X=X, y=y, C=C, penalty=p, trace=2,
#'	   beta_start_default=param_start_default,
#'     gamma_start_default=param_start_default, lambda=lambda_star, tau=tau,
#'     alpha=alpha, regressors_betas=regressors_betas,
#'     regressors_gammas=regressors_gammas, max_iter=max_iter, run_aauc=FALSE))
#'   }
#'
#'	 #tau to max and min
#'   tau_max <- calibrate_lambda_max (function_to_run=wrapper_tau,
#'     var_to_check=paste0("gammas_hat_",p), lambda_start = 5, n_min_var=1)
#'   tau_min <- calibrate_lambda_min (function_to_run=wrapper_tau,
#'     var_to_check=paste0("gammas_hat_",p), lambda_start = 0.0001, max_var= 50)
#'
#'   #create a suited tau
#'   tau <- create_lambda(n=3, lmax=tau_max, lmin=tau_min)
#'   tau <- as.numeric(formatC(tau, format = "e", digits = 9))
#'
#'   #start cross-validation
#'   assign(name, cPYE_KS_compute_cv(penalty=p, df=df, X=X, y=y, C=C, trace=2,
#'     beta_start_default=param_start_default,
#'     gamma_start_default=param_start_default, n_folds=n_folds, lambda=lambda,
#'     tau=tau, alpha=alpha, regressors_betas=regressors_betas,
#'     regressors_gammas=regressors_gammas, used_cores=used_cores,
#'     max_iter=max_iter))
#'
#'   #take the best lambda per measures (with the same measure for diff lambdas,
#'   #we take the one associated to less betas)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_yi"), get(name)$lambda_hat_cPYE_KS_yi)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_auc"), get(name)$lambda_hat_cPYE_KS_auc)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_aauc"), get(name)$lambda_hat_cPYE_KS_aauc)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_aYI"), get(name)$lambda_hat_cPYE_KS_aYI)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_ccr"), get(name)$lambda_hat_cPYE_KS_ccr)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_gm"), get(name)$lambda_hat_cPYE_KS_gm)
#'   assign(paste0("lambda_hat_cPYE_KS_", p ,"_cPYE"), get(name)$lambda_hat_cPYE_KS_cPYE)
#'
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_yi"), get(name)$tau_hat_cPYE_KS_yi)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_auc"), get(name)$tau_hat_cPYE_KS_auc)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_aauc"), get(name)$tau_hat_cPYE_KS_aauc)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_aYI"), get(name)$tau_hat_cPYE_KS_aYI)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_ccr"), get(name)$tau_hat_cPYE_KS_ccr)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_gm"), get(name)$tau_hat_cPYE_KS_gm)
#'   assign(paste0("tau_hat_cPYE_KS_", p ,"_cPYE"), get(name)$tau_hat_cPYE_KS_cPYE)
#'
#' }
#'
#' @export
cPYE_KS_compute_cv <- function (n_folds, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C, lambda, tau,
                                trace = 2, alpha=0.5, a1=3.7, a2=3, penalty="L1", regressors_betas=NULL,
                                regressors_gammas=NULL, seed=1, used_cores=1, trend = "monotone",
                                delta = 1e-5,
                                max_iter=10000,
                                max_alpha=10000,
                                min_alpha=1e-10,
                                convergence_error=1e-7,
                                stepsizeShrink=0.8,
                                kernel="gaussian",
								beta_start_input=NULL,
								beta_start_default="zeros",
                                gamma_start_input=NULL,
                                gamma_start_default="zeros", scaling=FALSE,
								run_aauc=FALSE) {

  start_time <- Sys.time()
  
  #check df:
  if (nrow(df)==0){stop("df has no rows")}
  
  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(C) == "data.frame"){
    C <- names(C)
  } else if (class(C) != "character"){stop("C can only be of class character or data.frame.")}
  if (class(X) == "data.frame"){
    X <- names(X)
  } else if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) == "data.frame"){
    y <- names(y)
  } else if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  #check if const already exists in the dataset
  if("const" %in% colnames(df)){stop("const already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}
  if("const" %in% colnames(C)){stop("const already exists as column in C! Please delete or rename this column since I need this name to set the internal const")}

  ID <- rownames(df)
  const <- rep(1, length(ID)) #the constant for the coefficients gammas
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X))], const, df[,(names(df) %in% C)])
  #put "const" in C
  C1 <- c("const", C)

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}

  set.seed(seed)

  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye (df=df1, X=colnames(df1[, names(df1) %in% c(X,C)]), y="y")$df_scaled
  }

  #check if df is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1$y[df1$y==1])<2*n_folds)){stop("df contains too few 1s for this number of folds")
  } else if ((length(df1$y[df1$y==0])<2*n_folds)){stop("df contains too few 0s for this number of folds")}

  #divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1,y)),1:2]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==0,])), replace=FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==1,])), replace=FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0,fold_i_1))
  folds_i <- merge(df1[,1:2], df_sort[,1:3], by='ID', all = FALSE, sort = FALSE)[,4]

  #mix lambda and tau in only one vector
  lambda_tau <- c(lapply(lambda, function (x) c(x, mean(tau))), lapply(tau, function (x) c(mean(lambda), x)))

  #names of the columns
  names <- c(unlist(lapply(lambda, function (x) paste0(paste("lambda", x, sep="="),"_", paste("tau", mean(tau), sep="=")))), unlist(lapply(tau, function (x) paste0(paste("lambda", mean(lambda), sep="="),"_", paste("tau", x, sep="=")))))
  names(lambda_tau) <- names

  foldnames <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))

  #measures (train and test)
  cPYE_KS_L12 <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  cPYE_KS_L1 <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  cPYE_KS_EN <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  cPYE_KS_SCAD <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  cPYE_KS_MCP <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  auc <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  aauc <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  aYI <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  youden_index <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  sensitivity <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  specificity <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  geometric_mean <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  fdr <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  mcc <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  corrclass <- list (train = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)), test = matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames,names)))
  n_betas <- matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames, names))
  n_gammas <- matrix(NA, nrow = n_folds, ncol = (length(lambda_tau)), dimnames= list(foldnames, names))

  betas <- vector(mode="list", length=n_folds)
  gammas <- vector(mode="list", length=n_folds)

  results = list("cPYE_KS_L12"=cPYE_KS_L12, "cPYE_KS_L1"=cPYE_KS_L1, "cPYE_KS_EN" = cPYE_KS_EN,
                 "cPYE_KS_SCAD" = cPYE_KS_SCAD, "cPYE_KS_MCP" = cPYE_KS_MCP, "auc"=auc,
                 "aauc"=aauc, "aYI"=aYI, "youden_index"=youden_index,
                 "sensitivity"=sensitivity, "specificity"=specificity,
                 "geometric_mean"=geometric_mean, "fdr"=fdr,
                 "mcc"=mcc, "corrclass"=corrclass, "n_betas" = n_betas,
                 "betas"= betas, "n_gammas"=n_gammas, "gammas"=gammas)

  cat("Starting CV with the following estimation/convergence parameters: \n")
  cat("max_iter:", max_iter, "\n")
  cat("min_alpha:", min_alpha, "\n")
  cat("convergence_error:", convergence_error, "\n")
  cat("stepsizeShrink:", stepsizeShrink, "\n")

  #fill the matrices
  results <- mapply(function(k) cPYE_KS.cv(df=df1[, !(names(df1) %in% c("ID"))], X=X, y=y, C=C, lambda_tau=lambda_tau, trace=trace, alpha=alpha,
                                          a1=a1, a2=a2, penalty=penalty,
                                          folds_i, k, regressors_betas, regressors_gammas, cPYE_KS_L12=cPYE_KS_L12,
                                          cPYE_KS_L1=cPYE_KS_L1,
                                          cPYE_KS_EN=cPYE_KS_EN, cPYE_KS_SCAD=cPYE_KS_SCAD, cPYE_KS_MCP=cPYE_KS_MCP,
                                          auc=auc, aauc=aauc, aYI=aYI,
                                          youden_index=youden_index,
                                          sensitivity=sensitivity, specificity=specificity,
                                          geometric_mean=geometric_mean, fdr=fdr,
                                          mcc=mcc, corrclass=corrclass,
                                          max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                          stepsizeShrink=stepsizeShrink,
                                          n_betas=n_betas, n_gammas=n_gammas, used_cores=used_cores, kernel=kernel, trend=trend,
                                          delta=delta, max_alpha=max_alpha, beta_start_input=beta_start_input,
										  beta_start_default=beta_start_default,
                                          gamma_start_input=gamma_start_input,
                                          gamma_start_default=gamma_start_default,
										  run_aauc=run_aauc),
                                          seq(1:n_folds))

  betas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas)))
  gammas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

  temp_cPYE <- get(paste0("cPYE_KS_", penalty))
  temp_cPYE[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]], paste0("cPYE_KS_", penalty))$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]],paste0("cPYE_KS_", penalty))$test[k,]))))
  assign(paste0("cPYE_KS_", penalty) , temp_cPYE)
  auc[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@auc$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@auc$test[k,]))))
  aauc[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@aauc$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@aauc$test[k,]))))
  aYI[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@aYI$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@aYI$test[k,]))))
  youden_index[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@youden_index$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@youden_index$test[k,]))))
  sensitivity[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@sensitivity$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@sensitivity$test[k,]))))
  specificity[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@specificity$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@specificity$test[k,]))))
  geometric_mean[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@geometric_mean$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@geometric_mean$test[k,]))))
  fdr[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@fdr$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@fdr$test[k,]))))
  mcc[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@mcc$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@mcc$test[k,]))))
  corrclass[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@corrclass$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@corrclass$test[k,]))))
  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))
  n_gammas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[k,])))

  if (length(tau)==1){
    tau_hat_cPYE_KS_yi <- tau
    tau_hat_cPYE_KS_auc <- tau
    tau_hat_cPYE_KS_aauc <- tau
    tau_hat_cPYE_KS_aYI <- tau
    tau_hat_cPYE_KS_ccr <- tau
    tau_hat_cPYE_KS_sen <- tau
    tau_hat_cPYE_KS_spc <- tau
    tau_hat_cPYE_KS_gm <- tau
    tau_hat_cPYE_KS_cPYE <- tau
  } else {

    tau_hat_cPYE_KS_yi <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(youden_index$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(youden_index$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_auc <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(auc$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(auc$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_aauc <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(aauc$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(aauc$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_aYI <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(aYI$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(aYI$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_ccr <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(corrclass$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(corrclass$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_sen <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(sensitivity$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(sensitivity$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_spc <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(specificity$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(specificity$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_gm <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(geometric_mean$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(geometric_mean$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
    tau_hat_cPYE_KS_cPYE <- lambda_tau[which(colnames(n_gammas) == names(which.min(colMeans(n_gammas)[colMeans(get(paste0("cPYE_KS_", penalty))$test)[(length(lambda)+1):length(lambda_tau)] == max(colMeans(get(paste0("cPYE_KS_", penalty))$test)[(length(lambda)+1):length(lambda_tau)])])))][[1]][2]
  }

  if (length(lambda)==1){
    lambda_hat_cPYE_KS_yi <- lambda
    lambda_hat_cPYE_KS_auc <- lambda
    lambda_hat_cPYE_KS_aauc <- lambda
    lambda_hat_cPYE_KS_aYI <- lambda
    lambda_hat_cPYE_KS_ccr <- lambda
    lambda_hat_cPYE_KS_sen <- lambda
    lambda_hat_cPYE_KS_spc <- lambda
    lambda_hat_cPYE_KS_gm <- lambda
    lambda_hat_cPYE_KS_cPYE <- lambda
  } else {
    lambda_hat_cPYE_KS_yi <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(youden_index$test)[1:length(lambda)] == max(colMeans(youden_index$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_auc <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(auc$test)[1:length(lambda)] == max(colMeans(auc$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_aauc <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(aauc$test)[1:length(lambda)] == max(colMeans(aauc$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_aYI <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(aYI$test)[1:length(lambda)] == max(colMeans(aYI$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_ccr <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(corrclass$test)[1:length(lambda)] == max(colMeans(corrclass$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_sen <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(sensitivity$test)[1:length(lambda)] == max(colMeans(sensitivity$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_spc <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(specificity$test)[1:length(lambda)] == max(colMeans(specificity$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_gm <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(geometric_mean$test)[1:length(lambda)] == max(colMeans(geometric_mean$test)[1:length(lambda)])])))][[1]][1]
    lambda_hat_cPYE_KS_cPYE <- lambda_tau[which(colnames(n_betas) == names(which.min(colMeans(n_betas)[colMeans(get(paste0("cPYE_KS_", penalty))$test)[1:length(lambda)] == max(colMeans(get(paste0("cPYE_KS_", penalty))$test)[1:length(lambda)])])))][[1]][1]
  }

  cv_time = difftime(Sys.time(), start_time , units = "mins")

  if (trace %in% c(1,2)) {
	cat("----------------> END OF THE CROSS-VALIDATION OF THE cPYE METHODS <----------------- \n")
    cat("----------------> For the whoole Cross-validation it took:", cv_time, "minutes \n")
  }

  return(list(penalty=penalty, kernel=kernel, cv_time=cv_time, cPYE_KS_L12=cPYE_KS_L12, cPYE_KS_L1=cPYE_KS_L1,
              cPYE_KS_EN=cPYE_KS_EN, cPYE_KS_SCAD=cPYE_KS_SCAD, cPYE_KS_MCP=cPYE_KS_MCP, auc=auc,
              aauc=aauc, aYI=aYI, youden_index=youden_index,
              sensitivity=sensitivity, specificity=specificity, geometric_mean=geometric_mean,
              fdr=fdr, mcc=mcc, corrclass=corrclass, n_betas=n_betas, betas=betas,
              n_gammas=n_gammas, gammas=gammas,
              lambda_hat_cPYE_KS_yi=lambda_hat_cPYE_KS_yi,
              lambda_hat_cPYE_KS_auc=lambda_hat_cPYE_KS_auc, lambda_hat_cPYE_KS_aauc=lambda_hat_cPYE_KS_aauc,
              lambda_hat_cPYE_KS_aYI=lambda_hat_cPYE_KS_aYI, lambda_hat_cPYE_KS_ccr=lambda_hat_cPYE_KS_ccr,
              lambda_hat_cPYE_KS_sen=lambda_hat_cPYE_KS_sen, lambda_hat_cPYE_KS_gm=lambda_hat_cPYE_KS_gm,
              lambda_hat_cPYE_KS_cPYE=lambda_hat_cPYE_KS_cPYE,
              tau_hat_cPYE_KS_yi=tau_hat_cPYE_KS_yi,
              tau_hat_cPYE_KS_auc=tau_hat_cPYE_KS_auc, tau_hat_cPYE_KS_aauc=tau_hat_cPYE_KS_aauc,
              tau_hat_cPYE_KS_aYI=tau_hat_cPYE_KS_aYI, tau_hat_cPYE_KS_ccr=tau_hat_cPYE_KS_ccr,
              tau_hat_cPYE_KS_sen=tau_hat_cPYE_KS_sen, tau_hat_cPYE_KS_gm=tau_hat_cPYE_KS_gm,
              tau_hat_cPYE_KS_cPYE=tau_hat_cPYE_KS_cPYE))
}

