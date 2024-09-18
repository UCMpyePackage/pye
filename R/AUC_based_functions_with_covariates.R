#Bioinformatics, 2014, titled "AucPR: An AUC-based approach
#through penalized regression for disease prediction using
#high-dimensional omics data".


# calculate coefficient from penalized regression
# for elastic lasso regression
# alpha=1 indicates lasso
# --> used later
#' @importFrom glmnet glmnet
Coefs_Elas <- function(cvm, mu, lambda = 0.1, alpha = 1,
                       family = "gaussian", intercept = FALSE){
  fitm = glmnet::glmnet(x = cvm, y = mu, alpha = alpha, intercept = intercept,
                        lambda = lambda, family = family)
  coefs = fitm$beta
  return(as.numeric(coefs))
}

#' @title AucPR_estimation
#'
#' @description function to estimate the optimal value of betas using the
#' AucPR method.
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param lambda the penalization parameter of the regressors X
#' @param alpha parameter for the Elastic-Net penalization term. If 1 the method
#' uses the L1 penalization. If (0,1) it uses Elastic-Net. Default is 1
#' @param model_type model to use in the estimation. The user can choose
#' among "AucPR_L1","AucPR_EN". Default is NULL, i.e. if not specified, it is
#' defined based on alpha, and since the default of alpha is 1, implicitly the
#' default is "AucPR_L1".
#' @param fold number of the fold if the function is called in cross-validation.
#' Just for visualization purposes. Default is NULL
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param max.print number of elements to show if printing the results. Default
#' is 10
#' @param ... possibility to insert additional variables if necessary, not
#' affecting the result
#'
#' @return a list containing the optimal value of betas and the value
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
#' regressors_betas<-simMicroarrayData_cov02_dim50_covariates$nregressors
#' model_type <- "AucPR_L1"
#' lambda <- 0.1
#'
#' z_hat <- AucPR_estimation(df=df, X=X, y=y, model_type=model_type, lambda=lambda,
#'       	regressors_betas=regressors_betas, trace=1)
#'
#' print(z_hat)
#'
#' @importFrom glmnet glmnet
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom stats cov
#' @export
# AucL1 and AucEN using CV, alpha=1--AucL1; alpha=0.5--AucEN
AucPR_estimation <- function (df, X=names(df[,!(names(df) == y)]), y="y", lambda, alpha = 1, model_type=NULL, fold=NULL,
                             regressors_betas=NULL, trace = 1, max.print=10,...){


  start_time <- Sys.time()

  options(max.print = max.print)

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X))])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("glmnet_estimation stopped because lambda needs to have length 1")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}

  if (length(model_type)!=0){
    if (model_type == "AucPR_L1") {
	  alpha=1
    } else if (model_type == "AucPR_EN") {
	  alpha=0.5
    } else {
	  stop("The parameter model_type is wrongly specified")
    }
  } else {
	if (alpha == 1){
	  model_type = "AucPR_L1"
	} else if ((alpha < 1) & (alpha > 0)){
	  model_type = "AucPR_EN"
	} else {
	  stop("The parameter alpha needs to be a number between 0 (not included) and 1 (included)")
	}
  }

  #divide cases and control
  cases <- as.matrix(df1[which(df1["y"]==1), (names(df1) %in% X)])
  controls <- as.matrix(df1[which(df1["y"]==0), (names(df1) %in% X)])

  # split the data and calculate auc for each lambda
  n = nrow(cases)
  m = nrow(controls)
  cvmd = stats::cov(cases)
  cvmn = stats::cov(controls)
  dd = cvmd + cvmn
  len = length(lambda)
  auc = 0

  mu1 = colMeans(cases) - colMeans(controls)
  mu = mu1/max(abs(mu1))

  betas_hat = Coefs_Elas(dd, mu, lambda, alpha = alpha)
  names(betas_hat)= names(df1[,-(1:2)])

  #standardize coefficients
  maxcoef = max(abs(betas_hat))
  if(maxcoef > 0) betas_hat = betas_hat/maxcoef

  #predict
  z_hat = as.matrix(df1[, names(df1) %in% X])%*%betas_hat
  names(z_hat) ="z_hat"
  df_hat = cbind(df1, z_hat)

  #optimal cutpoint
  opt = OptimalCutpoints::optimal.cutpoints(data= df_hat, X = "z_hat", status = "y", methods = "Youden", tag.healthy =0)
  c <- if(length(opt$Youden$Global$optimal.cutoff$cutoff)==1){opt$Youden$Global$optimal.cutoff$cutoff} else {opt$Youden$Global$optimal.cutoff$cutoff[round(length(opt$Youden$Global$optimal.cutoff$cutoff)/2)]}

  #y_hat
  y_hat = ifelse(df_hat["z_hat"] < c, 0, 1)
  colnames(y_hat)="y_hat"
  df_hat = cbind(df_hat, y_hat)

  #measures
  #AUC e YI
  auc = opt$Youden$Global$measures.acc$AUC[1]
  youden_index <- opt$Youden$Global$optimal.criterion #I decide to re-calculate it

  #confusion matix
  TP <- sum(ifelse(df_hat["y"] == 1 & df_hat["z_hat"] >=  c ,1,0))
  TN <- sum(ifelse(df_hat["y"] == 0 & df_hat["z_hat"] < c ,1,0))
  FP <- sum(ifelse(df_hat["y"] == 0 & df_hat["z_hat"] >=  c ,1,0))
  FN <- sum(ifelse(df_hat["y"] == 1 & df_hat["z_hat"] < c ,1,0))

  #compute the Youden index
  spec <- TN/(TN+FP)
  fnr <- FN/(FN+TP)
  #youden_index <- spec - fnr #I comment this since I recalculate it above in another way

  #sensitivity
  sensitivity = 1-fnr

  #geometric mean
  gm <- sqrt(spec*sensitivity)

  #compute other measures: FDR:
  fdr <- FP/(FP + TP)

  #compute other measures: FDR:
  mcc <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  #compute the correct classification:
  corrclass <- (TP + TN)/nrow(df_hat)

  #if all the betas are zero, the measures are zeros
  if(sum(betas_hat)==0){
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    spec <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  #prepare the output
  model = model_type
  betas_hat= betas_hat
  auc = auc
  corrclass = corrclass
  lambda=lambda
  c_hat=c
  z_hat=df_hat[, c("ID","z_hat")]
  y_hat=df_hat[, c("ID","y_hat")]

  names(betas_hat)=c(names(df1[,-(1:2)]))

  if (trace %in% c(1,2)) {
    cat("Estimation done using the Penalized AUC method type:", model_type, "; \n")
    cat("->")
    if (length(fold)!=0) {cat("fold =", fold, "; ")}
    if (!is.na(lambda)) {cat("lambda:", lambda, ";")}
    visualize_betas <- c(betas_hat[which(betas_hat!=0)], stats::setNames(c_hat, "c"))
    cat("youden_index:", youden_index, "; sensitivity:", sensitivity, "; specificity:", spec, "; geometric_mean:", gm, "; fdr:", fdr, "; mcc:", mcc, "; auc:", auc, "; corrclass:", corrclass, "; \n")
    cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, "; betas_hat: \n")
    print(visualize_betas)
  }

  #total number of variables
  n_total_var <- length(betas_hat)

  #n of predicted zeros
  n_predicted_zeros <- sum(betas_hat==0)

  #n of predicted betas different from 0
  n_predicted_non_zeros <- sum(betas_hat!=0)

  #compute other measures
  if (length(regressors_betas)!=0){

    #n of beta catch
    n_caught_betas <- sum(betas_hat[regressors_betas!=0]!=0)

    #n of beta not catch
    n_non_caught_betas <- sum(betas_hat[regressors_betas!=0]==0)

    #n? of zeros catch
    n_caught_zero <- sum(betas_hat[regressors_betas==0]==0)

    #n of zeros not catch
    n_zero_not_caught <- sum(betas_hat[regressors_betas==0]!=0)
  } else {
    n_caught_betas<-NA
    n_non_caught_betas<-NA
    n_caught_zero<-NA
    n_zero_not_caught<-NA
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
    cat("\n")
  }

  return(list(model=model,
              model_type=model_type,
              betas_hat= betas_hat,
              youden_index=youden_index,
              sensitivity=sensitivity,
              specificity=spec,
              geometric_mean=gm,
              fdr=fdr, mcc=mcc,
              corrclass=corrclass,
              auc=auc, lambda=lambda, alpha=alpha,
              c_hat=c_hat, z_hat=z_hat, y_hat=y_hat,
              n_betas=length(betas_hat[which(betas_hat!=0)]),
              n_total_var=n_total_var,
              n_predicted_zeros=n_predicted_zeros,
              n_predicted_non_zeros=n_predicted_non_zeros,
              n_caught_betas=n_caught_betas,
              n_non_caught_betas=n_non_caught_betas,
              n_caught_zero=n_caught_zero,
              n_zero_not_caught=n_zero_not_caught,
              estimation_time=estimation_time))

}






#' @title AucPR_predict
#'
#' @description function to apply the prediction new data using the estimated
#' betas coming from the AucPR_estimation function using the AucPR method.
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param model_to_use list of parameters result of the function AucPR_estimation
#' containing all the information necessiry for the prediction
#' @param fold number of the fold if the function is called in cross-validation.
#' Just for visualization purposes. Default is NULL
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param c_function_of_covariates if TRUE, covYI is used to estimate the
#' cut-off point as function of the convariate information. If FALSE, the
#' covariate information is ignored. In this function is used only to define
#' if we have to print or not the result. If TRUE the final print of the result
#' is not done, since it should be done in covYI estimation. Default is FALSE
#'
#' @return a list containing the prediction using AucPR and the value
#' of the main accuracy measure.
#'
#' @examples
#' library(pye)
#' cols <- 2000
#' cols_cov <- 20
#' seed=1
#' simMicroarrayData_cov02_dim50_covariates <- create_sample_with_covariates(
#' 		rows_train=50, cols=cols, cols_cov=cols_cov, covar=0.2, seed=seed)
#' train_df <- simMicroarrayData_cov02_dim50_covariates$train_df_scaled
#' test_df <- simMicroarrayData_cov02_dim50_covariates$test_df_scaled
#' X <- simMicroarrayData_cov02_dim50_covariates$X
#' y <- simMicroarrayData_cov02_dim50_covariates$y
#' regressors_betas<-simMicroarrayData_cov02_dim50_covariates$nregressors
#' model_type <- "AucPR_L1"
#' lambda <- 0.1
#' c_function_of_covariates <- FALSE
#'
#' model <- AucPR_estimation(df=train_df, X=X, y=y, lambda=lambda,
#'		regressors_betas=regressors_betas, model_type=model_type,
#' 		trace=1)
#'
#' predictions <- AucPR_predict(df=test_df, X=X, model_to_use = model,
#'		trace=1, regressors_betas = regressors_betas,
#'		c_function_of_covariates=c_function_of_covariates)
#'
#' print(predictions)
#'
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @export
AucPR_predict <- function (df, X=names(df[,!(names(df) == y)]), y="y", model_to_use, fold=NULL,
                            regressors_betas=NULL, trace=1, c_function_of_covariates=FALSE){

  start_time <- Sys.time()

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID" %in% colnames(df)){stop("ID already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID <- rownames(df)
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X))])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}

  betas_hat <- model_to_use$betas_hat

  #predict
  z_hat <-  as.matrix(df1[,(names(df1) %in% X)])%*%betas_hat
  names(z_hat) <- "z_hat"
  df_hat <- cbind(df1, z_hat)

  #c_hat
  c_hat <- model_to_use$c_hat

  #y_hat
  y_hat = ifelse(df_hat["z_hat"] < c_hat, 0, 1)
  colnames(y_hat)="y_hat"
  df_hat = cbind(df_hat, y_hat)

  #measures
  #AUC e YI
  opt <- OptimalCutpoints::optimal.cutpoints(data= df_hat, X = "z_hat", status = "y", methods = "Youden", tag.healthy =0)
  auc <-  opt$Youden$Global$measures.acc$AUC[1]
  youden_index <- opt$Youden$Global$optimal.criterion[1]

  #confusion matix
  TP <- sum(ifelse(df_hat["y"] == 1 & df_hat["z_hat"] >=  c_hat ,1,0))
  TN <- sum(ifelse(df_hat["y"] == 0 & df_hat["z_hat"] < c_hat ,1,0))
  FP <- sum(ifelse(df_hat["y"] == 0 & df_hat["z_hat"] >=  c_hat ,1,0))
  FN <- sum(ifelse(df_hat["y"] == 1 & df_hat["z_hat"] < c_hat ,1,0))

  #compute the Youden index
  spec <- TN/(TN+FP)
  fnr <- FN/(FN+TP)
  #youden_index <- spec - fnr

  #sensitivity
  sensitivity = 1-fnr

  #geometric mean
  gm <- sqrt(spec*sensitivity)

  #compute other measures: FDR:
  fdr <- FP/(FP + TP)

  #compute other measures: FDR:
  mcc <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  #compute the correct classiication:
  corrclass <- (TP + TN)/nrow(df_hat)

  #if all the betas are zero, the measures are zeros
  if(sum(betas_hat)==0){
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    spec <- 0
    youden_index <- 0
    fdr <- 0
    mcc <- 0
    corrclass <- 0
  }

  #prepare the output
  lambda=model_to_use$lambda
  alpha=model_to_use$alpha
  z_hat=df_hat[, c("ID","z_hat")]
  y_hat=df_hat[, c("ID","y_hat")]

  #standardize coefficients
  #maxcoef = max(abs(coefs))
  #if(maxcoef > 0) coefs = coefs/maxcoef

  if (trace %in% c(1,2)) {
    #print the results only if c_function_of_covariates = FALSE
    if(c_function_of_covariates==FALSE){
      cat("Prediction executed using the AUC model type:", model_to_use$model_type, "; \n")
      cat("->")
      if (length(fold)!=0) {cat("fold =", fold, "; ")}
      if (!is.na(lambda)) {cat("lambda:", lambda, ";")}
      visualize_betas <- betas_hat[which(betas_hat!=0)]
      cat("youden_index:", youden_index, "; sensitivity:", sensitivity, "; specificity:", spec, "; geometric_mean:", gm, "; fdr:", fdr, "; mcc:", mcc, "; auc:", auc, "; corrclass:", corrclass, "; \n")
      cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, "; betas_hat: \n")
      print(visualize_betas)
    }
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

    #n? of zeros catch
    n_caught_zero = sum(betas_hat[regressors_betas==0]==0)

    #n of zeros not catch
    n_zero_not_caught = sum(betas_hat[regressors_betas==0]!=0)
  } else {
    n_caught_betas=NA
    n_non_caught_betas=NA
    n_caught_zero=NA
    n_zero_not_caught=NA
  }

  estimation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(estimation_time)
    cat("\n")
  }

  return(list(model=model_to_use$model, model_type=model_to_use$model_type, betas_hat= model_to_use$betas_hat,
              youden_index=youden_index, sensitivity=sensitivity, specificity=spec, geometric_mean=gm, fdr=fdr,
              mcc=mcc, corrclass=corrclass, auc=auc, lambda=lambda, alpha=alpha, c_hat=c_hat, z_hat=z_hat, y_hat=y_hat,
              n_total_var=n_total_var, n_predicted_zeros=n_predicted_zeros, n_predicted_non_zeros=n_predicted_non_zeros,
              n_caught_betas= n_caught_betas, n_non_caught_betas=n_non_caught_betas, n_caught_zero=n_caught_zero,
              n_zero_not_caught=n_zero_not_caught, TP=TP, TN=TN, FP=FP, FN=FN,
              estimation_time=estimation_time))
}


#create the output class of the AucPR.cv function
setClass(Class="AucPR_cross_validation_output",
         representation(
           cv_time = "ANY",
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

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom methods new
AucPR.cv <- function (df, X, y, C, lambda, tau, alpha_g, penalty_g, folds_i, k, regressors_betas,
                       trace, model_type, kernel_g, a1_g, a2_g, trend_g, gamma_start_input, gamma_start_default,
                       regressors_gammas=NULL, max_iter_g, delta_g, max_alpha_g, stepsizeShrink_g, min_alpha_g,
                       convergence_error_g,
                       auc, aauc, aYI, youden_index, sensitivity, specificity,
                       geometric_mean, fdr, mcc, corrclass, n_betas, n_gammas, used_cores, c_function_of_covariates,
                       simultaneous, run_aauc){

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  if(c_function_of_covariates==TRUE){greek<-"tau  "} else {greek<-"lambda"}
  if(trace %in% c(1,2)){
    cat(" --------------------------------------------------------\n |  starting with the", k ,"-th fold for the CV of ", greek,"   |\n --------------------------------------------------------\n")
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
           max.cores, ")! but maybe you are distributing so: going on...")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_AUC_models.txt", setup_strategy = "sequential")
    if (!("cluster" %in% class(cl)))
      stop("cluster is not of class 'cl'; see ?makeCluster")

    if(trace %in% c(1,2)){
      cat("Start parallel computing for cross-validation, I am using:", length(cl),"cores \n")
    }
    parallel::clusterExport(cl, c("train_df", "X", "y", "model_type", "regressors_betas", "k", "trace", "lambda"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))
    fitted_models <- parallel::parLapply(cl, lambda, function(x) AucPR_estimation(df=train_df1, X=X, y=y,
                                                                                   lambda=x,
                                                                                   regressors_betas=regressors_betas,
                                                                                   model_type=model_type, fold=k,
                                                                                   trace=trace))

    on.exit(parallel::stopCluster(cl))

  } else {

    fitted_models <- lapply(lambda, function(x) AucPR_estimation(df=train_df1, X=X, y=y,
                                                                 lambda=x,
                                                                 regressors_betas=regressors_betas,
                                                                 model_type=model_type, fold=k,
                                                                 trace=trace))

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

    for (i in 1:length(lambda)){
      #measures on the train set
      #these are multiple tables based on the number of considered TAUs
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
      #cl <- parallel::makeCluster(used_cores, outfile="log_AUC_models.txt", setup_strategy = "sequential")
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

    for (i in 1:length(lambda)){


      #measures on the train set
      #these are multiple tables based on the number of considered TAUs
      auc[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "auc")))
      aauc[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "aauc")))
      aYI[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "aYI")))
      youden_index[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "youden_index")))
      sensitivity[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "sensitivity")))
      specificity[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "specificity")))
      geometric_mean[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "geometric_mean")))
      fdr[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "fdr")))
      mcc[[i]]$train[k, ]  <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "mcc")))
      corrclass[[i]]$train[k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "corrclass")))

      n_gammas[[i]][k, ] <- unlist(lapply(c(1:length(tau)), function(yy) getElement(fitted_gammas[[i]][[yy]], "n_gammas")))
      gammas_hat[[i]] <- lapply(c(1:(length(tau))), function(yy) getElement(fitted_gammas[[i]][[yy]], paste0("gammas_hat_", penalty_g)))
      names(gammas_hat[[i]])<- lapply(tau, function (xx) paste("tau", xx, sep="="))
    }

  }

  #n_betas is only based on lambda, not tau!
  n_betas[k, ] <- unlist(lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "n_betas")))

  #create a list of the results (NB: betas_hat includes c_hat) - it is only based on lambda, not tau!
  betas_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "betas_hat"))
  names(betas_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  #c_hat is the fixed value of c in case we don't use the covariates to estimate a patient's specific cut-off value
  c_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "c_hat"))
  names(c_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  betas <- betas_hat
  names(betas) <- lambda
  gammas <- gammas_hat


  all_measures_test <- mapply(function(z) AucPR_predict(df=test_df, X=X,
                                                         model_to_use = fitted_models[[z]], fold=k, trace=trace,
                                                         c_function_of_covariates=c_function_of_covariates),
                              1:length(lambda))

  #then I have to evaluate the different taus with AucPR_predict and choose the best tau
  #based on the performance on the test set.

  #organize the results
  if(c_function_of_covariates == FALSE) {

    #measures on the test set - if we don't compute c with covariates, we are only dependent of lambda:
    for (i in 1:length(lambda)){
      #create a list of the results
      auc[[i]]$test[k, ] <- all_measures_test["auc",][[i]]
      aauc[[i]]$test[k, ] <- NA
      aYI[[i]]$test[k, ] <- NA
      youden_index[[i]]$test[k, ] <- all_measures_test["youden_index",][[i]]
      sensitivity[[i]]$test[k, ] <- all_measures_test["sensitivity",][[i]]
      specificity[[i]]$test[k, ] <- all_measures_test["specificity",][[i]]
      geometric_mean[[i]]$test[k, ] <- all_measures_test["geometric_mean",][[i]]
      fdr[[i]]$test[k, ] <- all_measures_test["fdr",][[i]]
      mcc[[i]]$test[k, ] <- all_measures_test["mcc",][[i]]
      corrclass[[i]]$test[k, ] <- all_measures_test["corrclass",][[i]]
    }

  } else {

    #put "const" in C if not already there - it is needed in covYI function
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
                                                                                                             z_hat=z_hat[[xx]]),
                                                                                                    z="z_hat", y=y, C=C1,
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
        cat("With lambda:", lambda[i], " and method:", model_type,"\n")
        visualize_betas <- c(betas_hat[[i]][which(betas_hat[[i]]!=0)], setNames(c_hat[[i]], "c"))
        cat(" youden_index:", all_measures_test["youden_index",][[i]], "; sensitivity:", all_measures_test["sensitivity",][[i]], "; specificity:", all_measures_test["specificity",][[i]], "; geometric_mean:", all_measures_test["geometric_mean",][[i]], "; fdr:", all_measures_test["fdr",][[i]], "; mcc:", all_measures_test["mcc",][[i]], "; auc:", all_measures_test["auc",][[i]], "; corrclass:", all_measures_test["corrclass",][[i]], "; \n")
        cat("TP:", all_measures_test["TP",][[i]], "; TN:", all_measures_test["TN",][[i]], "; FP:", all_measures_test["FP",][[i]], "; FN:", all_measures_test["FN",][[i]], ";  betas: \n")
        print(visualize_betas)
		cat("Cross-validation time:", fitted_models[[i]]$estimation_time, "\n")
        cat("\n")

        for (ii in 1:length(cov_results[[i]])){
          cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend_g,"version. ; \n")
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

  return(methods::new("AucPR_cross_validation_output", auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index,
             sensitivity=sensitivity, specificity=specificity, geometric_mean=geometric_mean, fdr=fdr,
             mcc=mcc, corrclass=corrclass, n_betas=n_betas, n_gammas=n_gammas, betas=betas, gammas=gammas))
}







#' @title AucPR_compute_cv
#'
#' @description function to perform the cross-validation to select the best
#' value of lambda (and possibly tau) for the estimation of betas and c (and
#' possibly gammas) using AucPR (and possibly covYI).
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
#' same covariates included in df, or just a vector of character. Default is NULL
#' @param lambda the penalization parameter of the regressors X
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param model_type model to use in the estimation. The user can choose
#' among "AucPR_L1","AucPR_EN".
#' @param seed fix the seed. Default is 1
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param scaling if TRUE, the dataset is scaled. FALSE otherwise, Default is
#' FALSE.
#' @param c_function_of_covariates if TRUE, covYI is used to estimate the
#' cut-off point as function of the convariate information. If FALSE, the
#' covariate information is ignored. Default is FALSE
#' @param simultaneous in case of c_function_of_covariates=TRUE, it defines if
#' gammas needs to be estimated simultaneously or as a second step with respect
#' of betas. Default is FALSE, i.e. not simultaneously
#' @param measure_to_select_lambda the measure used to select lambda if
#' simultaneous=FALSE, i.e. when the cross-validation proecess select fist
#' lambda and then tau. Default is "ccr", i.e. the correct classification rate
#' @param alpha_g parameter for the Elastic-Net penalization term in covYI.
#' Default is 0.5
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
#' value of the correlation of every regressor with the target variable.
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
#' @return a list containing the optimal values of lambda (and possibly tau) to
#' estimate betas (and possibly gammas) and the value of the main accuracy measure
#' for all the folds.
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
#' c <- 0
#' prox_penalty = get(paste0("proximal_operator_", penalty))
#' trend = "monotone" #or "nonmonotone"
#' pye_starting_point = "zeros" #c("zeros", "corr")
#' alpha = 0.5
#' c_zero_fixed <- FALSE
#' c_function_of_covariates <- TRUE
#' used_cores <- 1
#' used_penalty_pye <- c("L1", "MCP") #c("L12", "L1", "EN", "SCAD", "MCP")
#' max_iter <- 10
#' n_folds <- 5
#' simultaneous <- FALSE #with false, covYI is applied as a second step, once we
#' #have already estimated lambda
#'
#' #measure_to_select_lambda (in case simultaneous = FALSE, this is the measure
#' #we  use to select lambda_star)
#' measure_to_select_lambda = "ccr" #c("auc", "aauc", "aYI", "ccr", "yi", "gm")
#'
#' if (c_function_of_covariates==TRUE){
#'   used_penalty_covYI = c("L1") #c("L12", "L1", "EN", "SCAD", "MCP")
#' } else {
#'   used_penalty_covYI = "NO"
#' }
#'
#' #see printed warnings:
#' options(warn=1)
#'
#' for (p in used_penalty_covYI){
#'
#'   #competitors AUC based
#'   AUC_based_models = "AucPR_L1" #c("AucPR_L1", "AucPR_EN")
#'
#'   for (m in AUC_based_models){
#'
#'     name <- paste0("param_estimate_", m, "_", p)
#'
#'     #wrapper of the function
#'     wrapper_lambda <- function(lambda){
#'       return(AucPR_estimation(df=df, X=X, y=y, model_type=m, lambda=lambda,
#'	     	regressors_betas=regressors_betas, trace=1))
#'     }
#'     #lambda to max and min
#'     lambda_max <- calibrate_lambda_max (function_to_run=wrapper_lambda,
#'     var_to_check="betas_hat", lambda_start = 5)
#'     lambda_min <- calibrate_lambda_min (function_to_run=wrapper_lambda,
#'     	var_to_check="betas_hat", lambda_start = 0.0001, max_var = 100)
#'
#'     cat("\n Method:", m, "; lambda_max=", lambda_max, "; lambda_min=", lambda_min,
#'     "\n")
#'
#'     #create a suited lambda
#'     lambda = create_lambda(n=3, lmax=lambda_max, lmin=lambda_min)
#'     lambda = as.numeric(formatC(lambda, format = "e", digits = 9))
#'
#'     if(c_function_of_covariates == TRUE){
#'       #wrapper of the function for tau
#'       lambda_star <- (lambda_max+lambda_min)/2
#'       z_hat <- AucPR_estimation(df=df, X=X, y=y, model_type=m, lambda=lambda_star,
#'       	regressors_betas=regressors_betas, trace=1)$z_hat$z_hat
#'       wrapper_tau <- function(tau){
#'         return(covYI_KS_estimation(df=cbind(df, z_hat=z_hat), z="z_hat", y=y, C=C,
#'         	tau=tau, gamma_start_input=NULL, gamma_start_default=pye_starting_point,
#'         	trace=2, alpha=alpha, penalty=p, regressors_gammas=regressors_gammas,
#'         	trend=trend, max_iter=max_iter)) #max_iter=50
#'       }
#'
#'       #tau to max and min
#'       tau_max <- calibrate_lambda_max (function_to_run=wrapper_tau,
#'       	var_to_check=paste0("gammas_hat_",p), lambda_start = 5, n_min_var=1)
#'       tau_min <- calibrate_lambda_min (function_to_run=wrapper_tau,
#'       	var_to_check=paste0("gammas_hat_",p), lambda_start = 0.0001,
#'			max_var=(length(C)-1))
#'
#'       cat("\n tau_max=", tau_max, "; tau_min=", tau_min, "\n")
#'
#'       #create a suited tau
#'       tau <- create_lambda(n=3, lmax=tau_max, lmin=tau_min)
#'       tau <- as.numeric(formatC(tau, format = "e", digits = 9))
#'     }
#'
#'     #start cross-validation
#'     assign(name, AucPR_compute_cv(model_type=m, trace=1, n_folds=n_folds,
#'                                  df=df, X=X, y=y, C=C, lambda=lambda, tau=tau,
#'                                  used_cores=used_cores,
#'                                  regressors_betas=regressors_betas,
#'                                  regressors_gammas=regressors_gammas,
#'                                  c_function_of_covariates=c_function_of_covariates,
#'                                  simultaneous=simultaneous,
#'                                  measure_to_select_lambda=measure_to_select_lambda,
#'                                  penalty_g=p, max_iter_g=max_iter, run_aauc=FALSE))
#'
#'     #take the best lambda per measures (if we have the same measure for diff
#'     #lambdas, we take the one associated to less betas)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_yi"), get(name)$lambda_hat_yi)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_auc"), get(name)$lambda_hat_auc)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_aauc"), get(name)$lambda_hat_aauc)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_aYI"), get(name)$lambda_hat_aYI)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_ccr"), get(name)$lambda_hat_ccr)
#'     assign(paste0("lambda_hat_", m , "_", p ,"_gm"), get(name)$lambda_hat_gm)
#'
#'     assign(paste0("tau_hat_", m , "_", p ,"_yi"), get(name)$tau_hat_yi)
#'     assign(paste0("tau_hat_", m , "_", p ,"_auc"), get(name)$tau_hat_auc)
#'     assign(paste0("tau_hat_", m , "_", p ,"_aauc"), get(name)$tau_hat_aauc)
#'     assign(paste0("tau_hat_", m , "_", p ,"_aYI"), get(name)$tau_hat_aYI)
#'     assign(paste0("tau_hat_", m , "_", p ,"_ccr"), get(name)$tau_hat_ccr)
#'     assign(paste0("tau_hat_", m , "_", p ,"_gm"), get(name)$tau_hat_gm)
#'   }
#' }
#'
#' @export
AucPR_compute_cv <- function (n_folds, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C=NULL, lambda,
                               tau=0, regressors_betas=NULL, trace=1,
                               model_type=c("AucPR_L1","AucPR_EN"),
                               seed=1, used_cores, scaling=FALSE,
                               c_function_of_covariates=FALSE,
                               simultaneous=FALSE,
                               measure_to_select_lambda="ccr",
                               alpha_g=0.5, penalty_g="L1", kernel_g="gaussian", a1_g=3.7, a2_g=3,
                               trend_g="monotone", gamma_start_input=NULL, gamma_start_default="zeros",
                               regressors_gammas=NULL, max_iter_g=10000,
                               delta_g=1e-5, max_alpha_g=10000,
                               stepsizeShrink_g=0.8, min_alpha_g=1e-12, convergence_error_g=1e-7,
							   run_aauc=FALSE){

  start_time <- Sys.time()

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}
  if(c_function_of_covariates==TRUE){
    if (class(C) != "character"){stop("C can only be of class character or data.frame.")}
  }

  tau2 <- 0

  #check if tau exists when c_function_of_covariates=TRUE
  if(c_function_of_covariates==TRUE){
    if(is.null(tau)){stop("tau cannot be NULL if c_function_of_covariates=TRUE.")}
    if(length(tau)==1){if(tau==0){stop("tau cannot be 0 if c_function_of_covariates=TRUE.")}}
	if(sum(tau)==0){stop("tau cannot be a vector of 0s if c_function_of_covariates=TRUE.")}
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
  df1 <- cbind(ID, df[,(names(df) %in% c(y,X,C))])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (!(penalty_g %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(measure_to_select_lambda %in% c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr"))){stop("A wrong value has been assigned to the parameter measure_to_select_lambda")}


  if (model_type == "AucPR_L1") {
    alpha=1
  } else if (model_type == "AucPR_EN") {
    alpha=0.5
  } else {
    stop("The parameter model_type is wrongly specified")
  }

  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye (df=df1, X=colnames(df1[, names(df1) %in% c(X,C)]), y="y")$df_scaled
  }

  set.seed(seed)
  #check if df1 is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if ((length(df1$y[df1$y==1])<2*n_folds)){stop("df1 contains too few 1s for this number of folds")
  } else if ((length(df1$y[df1$y==0])<2*n_folds)){stop("df1 contains too few 0s for this number of folds")}

  #divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- df1[order(getElement(df1,y)),1:2]
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==0,])), replace=FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==1,])), replace=FALSE)
  df_sort <- cbind(df_sort, c(fold_i_0,fold_i_1))
  folds_i <- merge(df1[,1:2], df_sort[,1:3], by='ID', all = FALSE, sort = FALSE)[,4]

  #names of the columns
  lambdanames <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
  taunames <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
  foldnames <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
  #measures (train and test)
  list_of_measures <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                        "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

  auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures){
    #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)))
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

  results <- list(auc=auc, aauc=aauc, aYI=aYI, youden_index=youden_index, sensitivity=sensitivity,
                 specificity=specificity, geometric_mean=geometric_mean, fdr=fdr, mcc=mcc,
                 corrclass=corrclass, n_betas=n_betas, n_gammas=n_gammas, betas=betas, gammas=gammas)

  cat("Starting CV with the following method: \n")
  cat("model_type:", model_type, "\n")

  #fill the matrices
  results <- mapply(function(k) AucPR.cv(df=df1[,names(df1)!="ID"], X=X, y=y, C=C, lambda=lambda, tau=tau,
                                          alpha_g=alpha_g, penalty_g=penalty_g, folds_i=folds_i, k=k,
                                          regressors_betas=regressors_betas,
                                          trace=trace, model_type=model_type,
                                          auc=auc, aauc=aauc, aYI=aYI,
                                          youden_index=youden_index,
                                          sensitivity=sensitivity,
                                          specificity=specificity,
                                          geometric_mean=geometric_mean, fdr=fdr,
                                          mcc=mcc, corrclass=corrclass,
                                          used_cores=used_cores,
                                          n_betas=n_betas, n_gammas=n_gammas,
                                          c_function_of_covariates=c_function_of_covariates,
                                          kernel_g=kernel_g, a1_g=a1_g, a2_g=a2_g, trend_g=trend_g,
                                          gamma_start_input=gamma_start_input,
                                          gamma_start_default=gamma_start_default,
                                          regressors_gammas=regressors_gammas, max_iter_g=max_iter_g, delta_g=delta_g,
                                          max_alpha_g=max_alpha_g, stepsizeShrink_g=stepsizeShrink_g,
                                          min_alpha_g=min_alpha_g, convergence_error_g=convergence_error_g,
										  run_aauc=run_aauc),
                    seq(1:n_folds))

  #prepare the results
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

  for (mes in list_of_measures){
    #aauc[] <- lapply(lambda, function(i) list(train = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$test[k,])))))
    #names(aauc) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
    assign(mes[], lapply(1:length(lambda), function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
    eval(substitute(names(x)<- unlist(lapply(lambda, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
  }

  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))
  n_gammas[] <- lapply(1:length(lambda), function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k,]))))

  measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
  list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                         "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

  tau_hat_auc <- tau_hat_aauc <- tau_hat_aYI <- tau_hat_yi <- tau_hat_sen <- NULL
  tau_hat_spc <- tau_hat_gm <- tau_hat_fdr <- tau_hat_mcc <- tau_hat_ccr <- NULL
  lambda_hat_auc <- lambda_hat_aauc <- lambda_hat_aYI <- lambda_hat_yi <- lambda_hat_sen <- NULL
  lambda_hat_spc <- lambda_hat_gm <- lambda_hat_fdr <- lambda_hat_mcc <- lambda_hat_ccr <- NULL
  for (i in 1:length(measures)){
    #in general, create tau_hat equal to NA
    assign(paste0("tau_hat_", measures[i]), NA)

    if(length(tau) != 1){
      measures_matrix <- t(sapply(1:length(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    } else {
      measures_matrix <- as.matrix(sapply(1:length(lambda), function(ii) colMeans(get(list_of_measures2[i])[[ii]]$test)))
    }

    rownames(measures_matrix) <- lambdanames
    colnames(measures_matrix) <- taunames

    if (!is.na(max(measures_matrix))){

      max_measures <- which(measures_matrix == max(measures_matrix), arr.ind = TRUE)[1,]
      assign(paste0("lambda_hat_", measures[i]), lambda[max_measures[1]])

      if (c_function_of_covariates==TRUE){
        if (length(tau)>1){
          assign(paste0("tau_hat_", measures[i]), tau[max_measures[2]])
        } else {
          assign(paste0("tau_hat_", measures[i]), tau)
        }
      }
    } else {
      assign(paste0("lambda_hat_", measures[i]), NA)
    }
  }

  #if tau2 is not zero, but the condition is valud for vectors as well
  if(length(tau2) > 1 || sum(tau2) > 0) {
    c_function_of_covariates <- TRUE
    tau <- tau2

	#re-create the output tables
    #names of the columns
    taunames <- unlist(lapply(tau, function (x) paste("tau", x, sep="=")))
    #measures (train and test)
    list_of_measures <- c("auc", "aauc", "aYI", "youden_index", "sensitivity", "specificity", "geometric_mean",
                          "fdr", "mcc", "corrclass")
    auc <- aauc <- aYI <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
    for (mes in list_of_measures){
      #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)))
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

    #now we execute AucPR.cv using the best lambda with respect of the measure in variable measure_to_select_lambda
    lambda_star <- get(paste0("lambda_hat_", measure_to_select_lambda))
    #re-fill the matrices
    results <- mapply(function(k) AucPR.cv(df=df1[,names(df1)!="ID"], X=X, y=y, C=C, lambda=lambda_star, tau=tau,
                                           alpha_g=alpha_g, penalty_g=penalty_g, folds_i=folds_i, k=k,
                                           regressors_betas=regressors_betas,
                                           trace=trace, model_type=model_type,
                                           auc=auc, aauc=aauc, aYI=aYI,
                                           youden_index=youden_index,
                                           sensitivity=sensitivity,
                                           specificity=specificity,
                                           geometric_mean=geometric_mean, fdr=fdr,
                                           mcc=mcc, corrclass=corrclass,
                                           used_cores=used_cores,
                                           n_betas=n_betas, n_gammas=n_gammas,
                                           c_function_of_covariates=c_function_of_covariates,
                                           kernel_g=kernel_g, a1_g=a1_g, a2_g=a2_g, trend_g=trend_g,
                                           gamma_start_input=gamma_start_input,
                                           gamma_start_default=gamma_start_default,
                                           regressors_gammas=regressors_gammas, max_iter_g=max_iter_g, delta_g=delta_g,
                                           max_alpha_g=max_alpha_g, stepsizeShrink_g=stepsizeShrink_g,
                                           min_alpha_g=min_alpha_g, convergence_error_g=convergence_error_g,
										   run_aauc=run_aauc),
                      seq(1:n_folds))

    #prepare the results
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

    for (mes in list_of_measures){
      #aauc[] <- lapply(lambda, function(i) list(train = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$test[k,])))))
      #names(aauc) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
      assign(mes[], lapply(1, function(i) list(train = wrapper(results, mes, i, "train", tau), test = wrapper(results, mes, i, "test", tau))))
      eval(substitute(names(x)<- unlist(lapply(lambda_star, function (xx) paste("lambda", xx, sep="="))), list(x=as.symbol(mes))))
    }

    n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))
    n_gammas[] <- lapply(1, function(i) t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_gammas[[i]][k,]))))

    betas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@betas)))
    gammas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

    measures <- c("auc", "aauc", "aYI", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr")
    list_of_measures2 <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                           "specificity", "geometric_mean", "fdr", "mcc", "corrclass")

    for (i in 1:length(measures)){
      #in general, create tau_hat equal to NA
      assign(paste0("tau_hat_", measures[i]), NA)

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
            assign(paste0("tau_hat_", measures[i]), tau[max_measures[2]])
          } else {
            assign(paste0("tau_hat_", measures[i]), tau)
          }
        }
      }
    }
  } else {
    lambda_star <- NA
  }

  cv_time = difftime(Sys.time(), start_time , units = "mins")

  if (trace %in% c(1,2)) {
	cat("----------------> END OF THE CROSS-VALIDATION OF THE PENALIZED AUC METHODS <----------------- \n")
    cat("----------------> For the whoole Cross-validation it took:", cv_time, "minutes \n")
  }

  return(list(model_type=model_type, penalty_g=penalty_g, cv_time=cv_time,cv_time=cv_time, auc=auc,
              aauc=aauc, aYI=aYI, youden_index=youden_index,
              sensitivity=sensitivity, specificity=specificity,
              geometric_mean=geometric_mean, fdr=fdr, mcc=mcc,
              corrclass=corrclass, betas=betas, lambda_hat_yi=lambda_hat_yi,
              lambda_hat_auc=lambda_hat_auc, lambda_hat_aauc=lambda_hat_aauc, lambda_hat_aYI=lambda_hat_aYI,
              lambda_hat_ccr=lambda_hat_ccr, lambda_hat_sen=lambda_hat_sen,
              lambda_hat_spc=lambda_hat_spc, lambda_hat_gm=lambda_hat_gm,
              tau_hat_yi=tau_hat_yi, tau_hat_auc=tau_hat_auc,
              tau_hat_aauc=tau_hat_aauc, tau_hat_aYI=tau_hat_aYI,
              tau_hat_ccr=tau_hat_ccr, tau_hat_sen=tau_hat_sen,
              tau_hat_spc=tau_hat_spc, tau_hat_gm=tau_hat_gm,
              c_function_of_covariates=c_function_of_covariates,
              simultaneous=simultaneous, measure_to_select_lambda=measure_to_select_lambda,
              lambda_star=lambda_star, n_betas=n_betas, n_gammas=n_gammas, betas=betas, gammas=gammas))
}
