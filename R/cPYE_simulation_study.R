#PYE simulation study function code
#n = the number of conducted simulations
#df = the dataset to use
#y = the name of the target var

#to extract just part of the result of cPYE
cPYE_KS_with_print <- function(df, X, y, C, betas, gammas, lambda, tau, sim_n, n, alpha, a1, a2,
                               penalty, est_time, niter, kernel, trace, run_aauc){
  cPYE_KS_result <- cPYE_KS(df=df, X=X, y=y, C=C, betas=betas, gammas=gammas, lambda=lambda,
                            tau=tau, alpha=alpha, a1=a1, a2=a2, penalty=penalty, prediction=TRUE, 
							kernel=kernel, run_aauc=run_aauc)

  if (trace %in% c(1,2)){
    cat("-> Results on the TEST SET \n")
    cat("-> algorithm: cPYE_KS_proximal_gradient_method ; ")
    cat("simulation n.", sim_n, "of", n ,"; ")
    visualize_betas <- c(betas[which(betas!=0)])
    visualize_gammas <- c(gammas[which(gammas!=0)])
    cat("lambda=", lambda, "tau=", tau, "; penalty:", penalty, "; cPYE_KS:", getElement(cPYE_KS_result,  paste0("cPYE_",penalty)), "; youden_index:", cPYE_KS_result$youden_index, "; sensitivity:", cPYE_KS_result$sensitivity, "; fdr:", cPYE_KS_result$fdr, "; mcc:", cPYE_KS_result$mcc, "; auc:", cPYE_KS_result$auc, "; corrclass:", cPYE_KS_result$corrclass, " \n")
    cat("TP:", cPYE_KS_result$TP, "; TN:", cPYE_KS_result$TN, "; FP:", cPYE_KS_result$FP, "; FN:", cPYE_KS_result$FN, ";  betas: \n")
    print(visualize_betas)
    print(visualize_gammas)
    cat("Estimation time:", est_time, "; Number of iterations:", niter)
    cat("\n")
  }

  return(cPYE_KS_result)
}


#' @title cPYE_simulation_study_real_data
#'
#' @description function to perform the simulation on n different 
#' experiments, splitting the dataset randomly in train and test,
#' estimating cPYE in the train and evaluating the classification
#' measures on the test set. 
#'
#' @param n number of esperiments. Default is 1000
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
#' @param tau the penalization parameter of the covariates C
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point of gamma.
#' If "zeros", it starts with all zero values, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable. Default
#' is "zeros"
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
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param scaling if TRUE, the dataset is scaled. FALSE otherwise, Default is
#' FALSE.
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
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
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
#'
#' @return a list containing the classification measure related to every simulation.
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
#' pye_starting_point <- "zeros" #c("zeros", "corr")
#' alpha <- 0.5
#' used_cores <- 1
#' used_penalty_cPYE <- c("L1") #c("L12", "L1", "EN", "SCAD", "MCP")
#' measures_cPYE <- c("ccr") #c("auc", "aauc", "aYI", "ccr", "yi", "gm", "pye")
#' n <- 5 #number of simulations
#' trend <- "monotone" #or "nonmonotone"
#' max_iter <- 10
#' 
#' cat("This simulation is on with algorithm: ", trend, "\n")
#' 
#' #--------------------------------------------------------------------------
#' #                                  START                                  |
#' #--------------------------------------------------------------------------
#' #pye Gaussian Kernel Smooth + covYI
#' for (p in used_penalty_cPYE){
#'   for (mm in 1:length(measures_cPYE)){
#'     cat("---> Starting with measure", measures_cPYE[mm], ", 
#'			for the penalty", p, "<---\n")
#'     name <- paste0("cPYE_KS_", p ,"_sim_study_", measures_cPYE[mm])
#'     lambda_to_use <- 0.1
#'     tau_to_use <- 0.05
#'     assign(name, cPYE_simulation_study_real_data(n=n, df=df, X=X, y=y, C=C, 
#'                              lambda=lambda_to_use, 
#'                              tau=tau_to_use, 
#'                              trend=trend,
#'                              beta_start_default=pye_starting_point, 
#'								gamma_start_default=pye_starting_point, 
#'                              trace=1,
#'                              alpha=alpha,
#'                              a1=3.7, a2=3.7,
#'                              penalty=p,
#'								max_iter=max_iter,
#'                              kernel="gaussian", used_cores=used_cores))
#'   }
#' } 
#' cat("Computation finished! Start saving data.")
#' print(name)
#'
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom pROC coords
#' @importFrom pROC roc
#' @export
cPYE_simulation_study_real_data <- function(n=1000, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C, lambda, tau,
                                            beta_start_input=NULL, beta_start_default="zeros",
											gamma_start_input=NULL, gamma_start_default="zeros", 
											trace=1, alpha, a1=3.7, a2=3, penalty="L1",
                                            regressors_betas=NULL, regressors_gammas=NULL,
                                            used_cores=1, scaling=FALSE,
                                            max_iter=10000, trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                                            min_alpha=1e-10, convergence_error=1e-7, kernel="gaussian", run_aauc=FALSE){

  #start calc sim. time
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
  if (n < 2){stop("n needs to be at least 2.")}


  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye(df=df1, X=colnames(df1[, names(df1) %in% c(X,C)]), y="y")$df_scaled
  }

  #generate seeds
  seeds <- c(1:n)

  #names of the columns
  names <- unlist(lapply(seeds, function (x) paste("seed", x, sep="=")))

  #function to be computed
  func <- function(seed, n, df, X, y, C, lambda, tau, beta_start_input=NULL, beta_start_default, 
					gamma_start_input=NULL, gamma_start_default, trace, alpha, a1=a1, a2=a2,
					max_iter=max_iter, penalty, regressors_betas=regressors_betas, regressors_gammas=regressors_gammas, 
					trend=trend, stepsizeShrink=stepsizeShrink, delta=delta,
					max_alpha=max_alpha, min_alpha=min_alpha, convergence_error=convergence_error, 
					kernel=kernel, run_aauc=run_aauc){

    if (trace %in% c(1,2)){
      cat("--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda =", lambda, "; tau =", tau, "; penalty:", penalty, "\n")
    }
    set.seed(seed)
    #I have to stratify. If not may happen that some test-set have less elements then required in 1 of the 2 cat
    #split <- sample(rep(1:nrow(df1)), size=(round(nrow(df1)*0.7,0)), replace=FALSE)
    df_sort <- df[order(getElement(df, y)),1:2]
    split_0 <- sample(rep(1:nrow(df_sort[df_sort$y==0,])), size=(round(nrow(df_sort[df_sort$y==0,])*0.7,0)), replace=FALSE)
    split_1 <- nrow(df_sort[df_sort$y==0,]) + sample(rep(1:nrow(df_sort[df_sort$y==1,])), size=(round(nrow(df_sort[df_sort$y==1,])*0.7,0)), replace=FALSE)
    split <- sort(c(split_0, split_1))
    train_df <- df[which(df$ID %in% df_sort[split, "ID"]),]
    test_df <- df[which(df$ID %in% df_sort[-split, "ID"]),]

    #train PYE KS
    train_solution <- cPYE_KS_estimation(df=train_df[,!(names(train_df) %in% c("const","ID"))], X=X, y=y, C=C, lambda=lambda, tau=tau,
                                         beta_start_input=beta_start_input, beta_start_default=beta_start_default,
										 gamma_start_input=gamma_start_input, gamma_start_default=gamma_start_default, 
										 trace=trace, alpha=alpha,
                                         a1=a1, a2=a2, max_iter=max_iter, penalty=penalty, regressors_betas=regressors_betas,
                                         regressors_gammas=regressors_gammas, trend=trend, stepsizeShrink=stepsizeShrink, delta=delta, 
										 max_alpha=max_alpha, min_alpha=min_alpha, convergence_error=convergence_error, 
										 kernel=kernel, run_aauc=run_aauc)

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(train_df, y)), train_solution$z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    #identify the estimated betas
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    gammas <- getElement(train_solution, paste0("gammas_hat_", penalty))
    estimation_time <- train_solution$estimation_time
    niter <- train_solution$niter

    #test the results
    test_solution <- cPYE_KS_with_print(df=test_df[,names(test_df)!="ID"], X=X, y=y, C=C1, betas=betas,
                                        gammas=gammas, lambda=lambda, tau=tau, sim_n=seed, n=n,
                                        alpha=alpha, a1=a1, a2=a2, penalty=penalty, est_time=estimation_time,
                                        niter=niter, kernel=kernel, trace=trace, run_aauc=run_aauc)

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(test_df, y)), test_solution$z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_test <- pROC::coords(est_roc, 1-roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    return(list(estimation_time=estimation_time, seed=seed, train_solution=train_solution,
                test_solution=test_solution, roc_spec_points = roc_spec_points,
                roc_sens_points_train=roc_sens_points_train, roc_sens_points_test = roc_sens_points_test))
  }

  #measures (train and test)
  cPYE_L12 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_L1 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_EN <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_SCAD <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_MCP <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aauc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  n_total_var_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_betas"))
  n_total_var_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_gammas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_betas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_betas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_gammas"))
  n_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_betas_betas"))
  n_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_gammas"))
  n_non_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_betas_betas"))
  n_non_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_gammas"))
  n_caught_zero_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_betas"))
  n_caught_zero_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_betas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_gammas"))
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(rep(0, length(X)), nrow =  length(X), ncol = 1, dimnames= list(colnames(df1[,(names(df1) %in% X)]), "n_times_beta_diff_zero"))
  gammas_times_selected <- matrix(rep(0, length(C1)), nrow =  length(C1), ncol = 1, dimnames= list(colnames(df1[,(names(df1) %in% C1)]), "n_times_gamma_diff_zero"))

  if (trace %in% c(1,2)){
    cat("------------------------------------------------------------------\n|                   starting simulation study                    |\n------------------------------------------------------------------\n")
    cat("lambda =", lambda, " tau =", tau, "; penalty:", penalty, "\n")
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
           ") is larger than the number of avaiable cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_sim_cPYE_real.txt")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for the simulation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("n", "df1", "X", "y", "C", "lambda", "tau", "beta_start_input", "beta_start_default", "gamma_start_input",
                                  "gamma_start_default", "trace", "a1", "a2","alpha", "penalty", "regressors_betas",
                                  "regressors_gammas", "trend", "stepsizeShrink", "delta", "max_alpha",
								  "min_alpha", "convergence_error", "kernel", "max_iter", "run_aauc"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))

    #run
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed=x, n=n, df=df1, X=X, y=y, C=C, lambda=lambda, tau=tau,
                                                              beta_start_input=beta_start_input, beta_start_default=beta_start_default,
															  gamma_start_input=gamma_start_input, gamma_start_default=gamma_start_default,
                                                              trace=trace, alpha=alpha, a1=a1, a2=a2,
                                                              penalty=penalty, regressors_betas=regressors_betas,
                                                              regressors_gammas=regressors_gammas,
                                                              trend=trend, stepsizeShrink=stepsizeShrink, delta=delta,
															  max_alpha=max_alpha, min_alpha=min_alpha, 
															  convergence_error=convergence_error, 
															  kernel=kernel, max_iter=max_iter, run_aauc=run_aauc))

    on.exit(parallel::stopCluster(cl))

  } else {

    #run
    simulation_study <- lapply(seeds, function(x) func(seed=x, n=n, df=df1, X=X, y=y, C=C, lambda=lambda, tau=tau,
                                                       beta_start_input=beta_start_input,
													   beta_start_default=beta_start_default,
                                                       gamma_start_input=gamma_start_input,
                                                       gamma_start_default=gamma_start_default,
                                                       trace=trace, alpha=alpha, a1=a1, a2=a2,
                                                       penalty=penalty, regressors_betas=regressors_betas,
                                                       regressors_gammas=regressors_gammas,
                                                       trend=trend, stepsizeShrink=stepsizeShrink, delta=delta,
													   max_alpha=max_alpha, min_alpha=min_alpha,
													   convergence_error=convergence_error, 
													   kernel=kernel, max_iter=max_iter, run_aauc=run_aauc))
  }

  estimation_time <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time"))))

  #fill the matrices
  #measures on the train set
  temp_cPYE <- get(paste0("cPYE_", penalty))
  temp_cPYE[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("cPYE_KS_", penalty))))
  auc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
  aauc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "aauc")))
  aYI[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "aYI")))
  youden_index[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "youden_index")))
  sensitivity[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "sensitivity")))
  specificity[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "specificity")))
  geometric_mean[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "geometric_mean")))
  fdr[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "fdr")))
  mcc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "mcc")))
  corrclass[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "corrclass")))
  n_total_var_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var_betas")))
  n_total_var_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var_gammas")))
  n_predicted_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros_betas")))
  n_predicted_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros_gammas")))
  n_predicted_non_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros_betas")))
  n_predicted_non_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros_gammas")))
  n_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_gammas")))
  n_non_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_non_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_gammas")))
  n_caught_zero_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero_betas")))
  n_caught_zero_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero_gammas")))
  n_zero_not_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught_betas")))
  n_zero_not_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught_gammas")))
  betas_times_selected[,1] <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty)))!=0)
  gammas_times_selected[,1] <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("gammas_hat_", penalty)))!=0)

  #measures on test
  temp_cPYE[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, paste0("cPYE_", penalty))))
  assign(paste0("cPYE_", penalty), temp_cPYE)
  auc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "auc")))
  aauc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "aauc")))
  aYI[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "aYI")))
  youden_index[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "youden_index")))
  sensitivity[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "sensitivity")))
  specificity[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "specificity")))
  geometric_mean[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "geometric_mean")))
  fdr[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "fdr")))
  mcc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "mcc")))
  corrclass[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "corrclass")))

  #roc curve
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_test)

  betas_start <- getElement(simulation_study[[1]]$train_solution, "betas_start")
  gammas_start <- getElement(simulation_study[[1]]$train_solution, "gammas_start")
  #betas[] <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))

  #end computing sim. time
  simulation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(simulation_time)
  }

  results <- list(simulation_time=simulation_time, estimation_time=estimation_time, used_cores=used_cores, n=n, lambda=lambda, tau=tau,
                  cPYE_L12=cPYE_L12, cPYE_L1=cPYE_L1,
                  cPYE_EN = cPYE_EN, cPYE_SCAD = cPYE_SCAD, cPYE_MCP = cPYE_MCP, auc=auc,
                  youden_index=youden_index, sensitivity=sensitivity,
                  specificity=specificity, geometric_mean=geometric_mean,
                  fdr=fdr, mcc=mcc, corrclass=corrclass,
                  n_total_var_betas=n_total_var_betas, n_total_var_gammas=n_total_var_gammas,
                  n_predicted_zeros_betas=n_predicted_zeros_betas, n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                  n_predicted_non_zeros_betas=n_predicted_non_zeros_betas, n_predicted_non_zeros_gammas=n_predicted_non_zeros_gammas,
                  n_caught_betas=n_caught_betas, n_caught_gammas=n_caught_gammas,
                  n_non_caught_betas=n_non_caught_betas, n_non_caught_gammas=n_non_caught_gammas,
                  n_caught_zero_betas=n_caught_zero_betas, n_caught_zero_gammas=n_caught_zero_gammas,
                  n_zero_not_caught_betas=n_zero_not_caught_betas, n_zero_not_caught_gammas=n_zero_not_caught_gammas,
                  betas_times_selected= betas_times_selected, gammas_times_selected= gammas_times_selected,
                  betas_start=betas_start, gammas_start=gammas_start,
                  roc_spec_points=roc_spec_points, roc_sens_points_train=roc_sens_points_train, roc_sens_points_test=roc_sens_points_test,
				  run_aauc=run_aauc)#, betas= betas, gammas=gammas)

  return(results)
}

















#' @title cPYE_simulation_study_synthetic_data
#'
#' @description function to perform the simulation on n different 
#' experiments, for all the considered penalties of cPYE for simulated data.
#' In every iteration we start creating a new synthetic dataset and 
#' then we randomly split it in train and test, estimating pye in the 
#' train and evaluating the classification measures on the test set. 
#'
#' @param n number of esperiments. Default is 1000
#' @param rows_train number of rows of the training sample. Default is 50,
#' for an high-dimensional setting
#' @param rows_test number of rows of the test sample. Default is 1000.
#' We suggest to create a test sample much bigger than the training sample
#' to test your method/model on more data
#' @param cols number of regressor variables of both the training and the 
#' test samples. Default is 2000, for an high-dimensional setting 
#' @param cols_cov number of covariate variables of both the training and the 
#' test samples. Default is 20, increase it for an high-dimensional setting
#' @param covar covariance in the covariace matrix of the normal distribution.
#' Increasing it, the created features are more correlated. Default is 0.5
#' @param mu mean of the (multivariate) normal distribution of the regressors. 
#' Default is 0
#' @param mu_cov mean of the (multivariate) normal distribution of the 
#' covariates. Default is 0
#' @param lambda the penalization parameter of the regressors X
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
#' @param gamma_start_input vector of a specific starting point for gammas.
#' Default is NULL, i.e. no input vector
#' @param gamma_start_default set the default starting point of gamma.
#' If "zeros", it starts with all zero values, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable. Default
#' is "zeros"
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param max_iter maximum number of iterations in the algorithms mmAPG and
#' mnmAPG. Default is 500
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
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param run_aauc if FALSE the aAUC and aYI are not computed, to save stimation
#' time if not requested. Default is FALSE
#'
#' @return a list containing the classification measure related to every simulation.
#'
#' @examples
#' library(pye)
#' rows_train <- 50
#' rows_test <- 1000 
#' cols <- 2000
#' cols_cov <- 20
#' covar <- 0.5
#' mu <- rep(0,cols)
#' mu_cov <- rep(0,cols_cov)
#' pye_starting_point <- "zeros" #c("zeros", "corr")
#' alpha <- 0.5
#' used_cores <- 1
#' used_penalty_pye <- c("L1") #c("L12", "L1", "EN", "SCAD", "MCP")
#' measures_pye <- c("ccr") #c("auc", "aauc", "aYI", "ccr", "yi", "gm", "pye")
#' n <- 5 #number of simulations
#' trend <- "monotone" #or "nonmonotone"
#' max_iter <- 10
#' 
#' cat("This simulation is on with algorithm: ", trend, "\n")
#' 
#' #--------------------------------------------------------------------------
#' #                                  START                                  |
#' #--------------------------------------------------------------------------
#' #pye Gaussian Kernel Smooth + covYI
#' for (p in used_penalty_pye){
#'   for (mm in 1:length(measures_pye)){
#'     cat("---> Starting with measure", measures_pye[mm], ", for the penalty", p, "<---\n")
#'     name <- paste0("cPYE_KS_covYI_", p ,"_sim_study_", measures_pye[mm])
#'     lambda_to_use <- 0.1
#'     tau_to_use <- 0.05
#'     assign(name, cPYE_simulation_study_synthetic_data(n=n, rows_train=rows_train, 
#'							rows_test=rows_test,  cols=cols, 
#'							cols_cov=cols_cov, covar=covar, 
#'							mu=mu, mu_cov=mu_cov,
#'                          lambda=lambda_to_use, 
#'                          tau=tau_to_use, 
#'                          trend=trend,
#'                          beta_start_default=pye_starting_point,
#'                          gamma_start_default=pye_starting_point,
#'                          trace=1,
#'                          alpha=alpha,
#'                          a1=3.7, a2=3.7,
#'                          penalty=p,
#'                          kernel="gaussian",
#'							max_iter=max_iter,
#'							used_cores=used_cores,
#' 							run_aauc=FALSE))
#'   }
#' } 
#' cat("Computation finished! Start saving data.")
#' print(name)
#'
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom pROC coords
#' @importFrom pROC roc
#' @export
cPYE_simulation_study_synthetic_data <- function(n=1000, rows_train=50, rows_test=1000, cols=2000, cols_cov=20,
                                                 covar=0.5, mu=rep(0,cols), mu_cov=rep(0,cols_cov), lambda, tau,
                                                 beta_start_input=NULL, beta_start_default="zeros",
												 gamma_start_input=NULL, gamma_start_default="zeros", 
												 trace=1, alpha, a1=3.7, a2=3, penalty="L1", max_iter=10000, 
												 trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                                                 min_alpha=1e-10, convergence_error=1e-7, kernel="gaussian", used_cores=1,
												 run_aauc=FALSE){
  #start calc sim. time
  start_time <- Sys.time()

  #starting controls
  if (length(tau)>1){stop("cPYE_simulation_study_synthetic_data stopped because tau needs to have length 1")}
  if (length(lambda)>1){stop("cPYE_simulation_study_synthetic_data stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}
  if (n < 2){stop("n needs to be at least 2.")}

  #generate seeds
  seeds <- c(1:n)

  #names of the columns
  names <- unlist(lapply(seeds, function (x) paste("seed", x, sep="=")))

  #function to be computed
  func <- function(seed, n, rows_train, rows_test, cols, cols_cov, covar, mu, mu_cov,
				   lambda, tau,
                   beta_start_input=NULL, beta_start_default,
				   gamma_start_input=NULL, gamma_start_default, 
				   trace, alpha, a1=a1, a2=a2,
                   penalty, trend=trend, stepsizeShrink=stepsizeShrink, delta=delta,
				   max_alpha=max_alpha, min_alpha=min_alpha, convergence_error=convergence_error, 
				   kernel=kernel, max_iter=max_iter, run_aauc=run_aauc){

    if (trace %in% c(1,2)){
      cat("--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda =", lambda, "tau =", tau, "; penalty:", penalty, "\n")
    }

    #create the dataframe
    set.seed(1)

    #Vars 200, 400, 600, 800, 1200, 1400, 1600, 1800
    #Covars 3, 7, 12, 18
	
	df <- NULL
    df <- create_sample_with_covariates(rows_train=rows_train, cols=cols, cols_cov=cols_cov, covar=covar, mu=mu,
                                        mu_cov=mu_cov, rows_test=rows_test, seed=seed, varsN=c(200,400),
                                        varsB=c(600,800), varsE=c(1200,1400), varsP=c(1600,1800),
                                        varN_cov=3, varB_cov=7, varE_cov=12, varP_cov=18)

    train_df <- df$train_df_scaled
    test_df <- df$test_df_scaled
	
    X <- df$X
    y <- df$y
    C <- df$C
	#put "const" in C
    C1 <- c("const", C)
    regressors_betas <- df$nregressors
    regressors_gammas <- df$ncovariates

    #on the train set
    #ID <- rownames(train_df)
    #const <- rep(1, length(ID)) #the constant for the coefficients gammas
    #train_df1 <- cbind(ID, train_df[,(names(train_df) %in% c(y,X))], const, train_df[,(names(train_df) %in% C)])
    train_df1 <- train_df

    #on the test set
    ID <- rownames(test_df)
    const <- rep(1, length(ID)) #the constant for the coefficients gammas
    test_df1 <- cbind(test_df[,(names(test_df) %in% c(y,X))], const, test_df[,(names(test_df) %in% C)])

    #train PYE KS
    train_solution <- cPYE_KS_estimation(df=train_df1, X=names(train_df1[,!(names(train_df1) %in% c(y,C))]), y=y,
                                         C=C, lambda=lambda, tau=tau,
                                         beta_start_input=beta_start_input, beta_start_default=beta_start_default,
										 gamma_start_input=gamma_start_input, gamma_start_default=gamma_start_default,
                                         trace=trace, alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                         regressors_betas=regressors_betas, regressors_gammas=regressors_gammas,
                                         trend=trend, stepsizeShrink=stepsizeShrink, delta=delta, 
										 max_alpha=max_alpha, min_alpha=min_alpha, convergence_error=convergence_error, 
										 kernel=kernel, max_iter=max_iter, run_aauc=run_aauc)

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(train_df1, y)), train_solution$z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    #identify the estimated betas
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    gammas <- getElement(train_solution, paste0("gammas_hat_", penalty))
    #c <- getElement(train_solution, "c_hat")
    estimation_time <- train_solution$estimation_time
    niter <- train_solution$niter

    #test the results
    test_solution <- cPYE_KS_with_print(df=test_df1[,names(test_df1)!="ID"], X=X, y=y, C=C1, betas=betas,
                                        gammas=gammas, lambda=lambda, tau=tau, sim_n=seed, n=n,
                                        alpha=alpha, a1=a1, a2=a2, penalty=penalty, est_time=estimation_time,
                                        niter=niter, kernel=kernel, trace=trace, run_aauc=run_aauc)

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(test_df, y)), test_solution$z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_test <- pROC::coords(est_roc, 1-roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    return(list(estimation_time=estimation_time, seed=seed, train_solution=train_solution, test_solution=test_solution,
                roc_spec_points = roc_spec_points, roc_sens_points_train = roc_sens_points_train,
                roc_sens_points_test = roc_sens_points_test, regressors_betas=regressors_betas,
                regressors_gammas=regressors_gammas))
  }

  #simulate one dataset to take the names of the variables:
  #create the dataframe
  set.seed(1)
  #compute the df (the seed makes different all the simulations)
  df_for_the_names<- create_sample_with_covariates(rows_train=rows_train, cols=cols, cols_cov=cols_cov, covar=covar, mu=mu,
                                                   mu_cov=mu_cov, rows_test=rows_test, seed=1, varsN=c(200,400),
                                                   varsB=c(600,800), varsE=c(1200,1400), varsP=c(1600,1800),
                                                   varN_cov=3, varB_cov=7, varE_cov=12, varP_cov=18)

  if (trace %in% c(1,2)){
    cat("Real regressors are the variables", df_for_the_names$regressors, "\n")
    cat("Real covariates are the variables", df_for_the_names$covariates, "\n")
  }
  
  X <- df_for_the_names$X
  y <- df_for_the_names$y
  C <- df_for_the_names$C
 
  #measures (train and test)
  cPYE_L12 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_L1 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_EN <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_SCAD <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  cPYE_MCP <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aauc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  n_total_var_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_betas"))
  n_total_var_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_gammas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_betas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_betas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_gammas"))
  n_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_betas_betas"))
  n_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_gammas"))
  n_non_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_betas_betas"))
  n_non_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_gammas"))
  n_caught_zero_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_betas"))
  n_caught_zero_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_betas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_gammas"))
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))])), nrow =  ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))]), ncol = 1, dimnames= list(colnames(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))]), "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))])), nrow =  ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))]), ncol = 1, dimnames= list(colnames(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))]), "n_times_gamma_diff_zero"))
  #betas <- vector(mode="list", length=length(seeds))
  #names(betas) <- names

  #betas <- vector(mode="list", length=length(seeds))
  #names(betas) <- names

  if (trace %in% c(1,2)){
    cat("------------------------------------------------------------------\n|                   starting simulation study                    |\n------------------------------------------------------------------\n")
    cat("lambda =", lambda, "tau =", tau, "; penalty:", penalty, "\n")
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
           ") is larger than the number of avaiable cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_sim_cPYE_synthetic.txt")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for the simulation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("n", "rows_train", "rows_test", "cols", "cols_cov", "covar", "mu", "mu_cov",
								  "lambda", "tau", "beta_start_input", "beta_start_default", "gamma_start_input",
                                  "gamma_start_default", "trace", "a1", "a2","alpha",
                                  "penalty", "trend", "stepsizeShrink", "delta", "max_alpha",
								  "min_alpha", "convergence_error", "kernel", "max_iter", "run_aauc"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))

    #run
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed=x, n=n, rows_train=rows_train, rows_test=rows_test,
                                                                        cols=cols, cols_cov=cols_cov, covar=covar,
																		mu=mu, mu_cov=mu_cov,
                                                                        lambda=lambda, tau=tau,
                                                                        beta_start_input=beta_start_input,
																		beta_start_default=beta_start_default,
                                                                        gamma_start_input=gamma_start_input,
                                                                        gamma_start_default=gamma_start_default,
                                                                        trace=trace, alpha=alpha, a1=a1, a2=a2,
                                                                        penalty=penalty, trend=trend,
																		stepsizeShrink=stepsizeShrink, delta=delta,
																		max_alpha=max_alpha, min_alpha=min_alpha,
																		convergence_error=convergence_error, 
																		kernel=kernel, max_iter=max_iter, run_aauc=run_aauc))

    on.exit(parallel::stopCluster(cl))

  } else {

    #run
    simulation_study <- lapply(seeds, function(x) func(seed=x, n=n, rows_train=rows_train, rows_test=rows_test,
                                                       cols=cols, cols_cov=cols_cov, covar=covar,
													   mu=mu, mu_cov=mu_cov,
                                                       lambda=lambda, tau=tau,
                                                       beta_start_input=beta_start_input,
													   beta_start_default=beta_start_default,
                                                       gamma_start_input=gamma_start_input,
                                                       gamma_start_default=gamma_start_default,
                                                       trace=trace, alpha=alpha, a1=a1, a2=a2,
                                                       penalty=penalty, trend=trend,
													   stepsizeShrink=stepsizeShrink, delta=delta,
													   max_alpha=max_alpha, min_alpha=min_alpha,
													   convergence_error=convergence_error,
                                                       kernel=kernel, max_iter=max_iter, run_aauc=run_aauc))
  }

  estimation_time <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time"))))

  #fill the matrices
  #measures on the train set
  temp_cPYE <- get(paste0("cPYE_", penalty))
  temp_cPYE[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("cPYE_KS_", penalty))))
  auc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
  aauc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "aauc")))
  aYI[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "aYI")))
  youden_index[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "youden_index")))
  sensitivity[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "sensitivity")))
  specificity[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "specificity")))
  geometric_mean[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "geometric_mean")))
  fdr[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "fdr")))
  mcc[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "mcc")))
  corrclass[,1]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "corrclass")))
  n_total_var_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var_betas")))
  n_total_var_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var_gammas")))
  n_predicted_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros_betas")))
  n_predicted_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros_gammas")))
  n_predicted_non_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros_betas")))
  n_predicted_non_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros_gammas")))
  n_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_gammas")))
  n_non_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_non_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_gammas")))
  n_caught_zero_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero_betas")))
  n_caught_zero_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero_gammas")))
  n_zero_not_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught_betas")))
  n_zero_not_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught_gammas")))
  betas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty)))!=0)
  gammas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("gammas_hat_", penalty)))!=0)

  #measures on test
  temp_cPYE[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, paste0("cPYE_", penalty))))
  assign(paste0("cPYE_", penalty), temp_cPYE)
  auc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "auc")))
  aauc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "aauc")))
  aYI[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "aYI")))
  youden_index[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "youden_index")))
  sensitivity[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "sensitivity")))
  specificity[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "specificity")))
  geometric_mean[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "geometric_mean")))
  fdr[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "fdr")))
  mcc[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "mcc")))
  corrclass[,2]  <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, "corrclass")))

  #roc curve
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_test)

  betas_start <- getElement(simulation_study[[1]]$train_solution, "betas_start")
  gammas_start <- getElement(simulation_study[[1]]$train_solution, "gammas_start")

  #betas[] <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))
  #gammas[] <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_solution, paste0("gammas_hat_", penalty)))

  #regressors <- simulation_study[[1]]$regressors
  regressors_betas <- simulation_study[[1]]$nregressors
  regressors_gammas <- simulation_study[[1]]$ncovariates

  #end computing sim. time
  simulation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(simulation_time)
  }

  results <- list(simulation_time=simulation_time, estimation_time=estimation_time, used_cores=used_cores, n=n, lambda=lambda, tau=tau,
                  cPYE_L12=cPYE_L12, cPYE_L1=cPYE_L1,
                  cPYE_EN = cPYE_EN, cPYE_SCAD = cPYE_SCAD, cPYE_MCP = cPYE_MCP, auc=auc,
                  youden_index=youden_index, sensitivity=sensitivity,
                  specificity=specificity, geometric_mean=geometric_mean,
                  fdr=fdr, mcc=mcc, corrclass=corrclass,
                  n_total_var_betas=n_total_var_betas, n_total_var_gammas=n_total_var_gammas,
                  n_predicted_zeros_betas=n_predicted_zeros_betas, n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                  n_predicted_non_zeros_betas=n_predicted_non_zeros_betas, n_predicted_non_zeros_gammas=n_predicted_non_zeros_gammas,
                  n_caught_betas=n_caught_betas, n_caught_gammas=n_caught_gammas,
                  n_non_caught_betas=n_non_caught_betas, n_non_caught_gammas=n_non_caught_gammas,
                  n_caught_zero_betas=n_caught_zero_betas, n_caught_zero_gammas=n_caught_zero_gammas,
                  n_zero_not_caught_betas=n_zero_not_caught_betas, n_zero_not_caught_gammas=n_zero_not_caught_gammas,
                  betas_times_selected= betas_times_selected, gammas_times_selected= gammas_times_selected,
                  betas_start=betas_start, gammas_start=gammas_start,
                  roc_spec_points=roc_spec_points, roc_sens_points_train=roc_sens_points_train,
                  roc_sens_points_test=roc_sens_points_test, regressors_betas=regressors_betas,
                  regressors_gammas=regressors_gammas, run_aauc=run_aauc)#, betas= betas, gammas=gammas)
  return(results)
}
