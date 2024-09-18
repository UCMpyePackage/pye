#to extract just part of the result of pye
pye_KS_with_print <- function(df, X, y, betas, lambda, c, sim_n, n, alpha, a1, a2, penalty, est_time, niter, kernel, trace){
  pye_KS_result <- pye_KS(df=df, X=X, y=y, betas=betas, lambda=lambda, c=c, alpha=alpha, a1=a1, a2=a2, penalty=penalty, prediction=TRUE, kernel=kernel)

  if (trace %in% c(1,2)){
    cat("-> Results on the TEST SET \n")
    cat("-> algorithm: pye_KS_proximal_gradient_method ; ")
    cat("simulation n.", sim_n, "of", n ,"; ")
    visualize_betas <- c(betas[which(betas!=0)], c)
    cat("lambda:", lambda, "; penalty:", penalty, "; pye_KS:", getElement(pye_KS_result,  paste0("pye_",penalty)), "; youden_index:", pye_KS_result$youden_index, "; sensitivity:", pye_KS_result$sensitivity, "; fdr:", pye_KS_result$fdr, "; mcc:", pye_KS_result$mcc, "; auc:", pye_KS_result$auc, "; corrclass:", pye_KS_result$corrclass, " \n")
    cat("TP:", pye_KS_result$TP, "; TN:", pye_KS_result$TN, "; FP:", pye_KS_result$FP, "; FN:", pye_KS_result$FN, ";  betas: \n")
    print(visualize_betas)
    cat("Estimation time:", est_time, "; Number of iterations:", niter)
    cat("\n")
  }

  return(pye_KS_result)
}

#' @title PYE_simulation_study_real_data
#'
#' @description function to perform the simulation on n different
#' experiments, splitting the dataset randomly in train and test,
#' estimating pye in the train and evaluating the classification
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
#' @param tau the penalization parameter of the covariates C, in covYI. Default 0,
#' i.e. no penalization term
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
#' @param regressors_betas a vector containing the real betas (if known). Default
#' is NULL, i.e. we do not know the real regressors
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param scaling if TRUE, the dataset is scaled. FALSE otherwise, Default is
#' FALSE.
#' @param c_zero_fixed if TRUE the estimation process considers c, the cut-off
#' point, as fixed and equal to zero, to reduce the complexity of the estimation.
#' If FALSE c can vary and be different from zero, estimated by pye. Default
#' is FALSE
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
#' @param c_function_of_covariates if TRUE, covYI is used to estimate the
#' cut-off point as function of the convariate information. If FALSE, the
#' covariate information is ignored. Default is FALSE
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
#' value of the correlation of every regressor with the target variable. Default
#' is "zeros"
#' @param regressors_gammas a vector containing the real gammas (if known).
#' Default is NULL
#' @param max_iter_g maximum number of iterations in the algorithms mmAPG and
#' mnmAPG in covYI. Default is 500
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
#' used_penalty_pye <- c("L1") #c("L12", "L1", "EN", "SCAD", "MCP")
#' measures_pye <- c("ccr") #c("auc", "aauc", "aYI", "ccr", "yi", "gm", "pye")
#' n<- 5 #number of simulations
#' c_function_of_covariates <- TRUE
#' trend <- "monotone" #or "nonmonotone"
#' c_zero_fixed <- FALSE
#'
#' cat("This simulation is on with algorithm: ", trend, "and is c fixed: ",
#'   c_zero_fixed, "\n")
#'
#' #--------------------------------------------------------------------------
#' #                                  START                                  |
#' #--------------------------------------------------------------------------
#' #pye Gaussian Kernel Smooth + covYI
#' for (p in used_penalty_pye){
#'   for (mm in 1:length(measures_pye)){
#'     cat("---> Starting with measure", measures_pye[mm], ", for the penalty", p, "<---\n")
#'     name <- paste0("PYE_KS_covYI_", p ,"_sim_study_", measures_pye[mm])
#'     lambda_to_use <- 0.1
#'     tau_to_use <- 0.05
#'     assign(name, PYE_simulation_study_real_data(n=n, df=df, X=X, y=y, C=C,
#'                             lambda=lambda_to_use,
#'                             tau=tau_to_use,
#'                             trend=trend,
#'                             beta_start_default=pye_starting_point,
#'                             gamma_start_default=pye_starting_point,
#'                             trace=1,
#'                             alpha=alpha, alpha_g=alpha,
#'                             a1=3.7, a2=3.7,
#'                             a1_g=3.7, a2_g=3.7,
#'                             penalty=p, penalty_g=p,
#'                             kernel="gaussian", used_cores=used_cores,
#'                             c_function_of_covariates=c_function_of_covariates,
#'							   c_zero_fixed=c_zero_fixed,
#' 							   run_aauc=FALSE))
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
PYE_simulation_study_real_data <- function(n=1000, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y", C=NULL,
                                            lambda, tau=0,
                                            beta_start_input=NULL, beta_start_default="zeros", trace=1,
                                            alpha, a1=3.7, a2=3, penalty, regressors_betas=NULL,
                                            used_cores=1, scaling=FALSE, c_zero_fixed=FALSE,
                                            max_iter=10000, trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.8,
                                            min_alpha=1e-10, convergence_error=1e-7, kernel="gaussian",
                                            c_function_of_covariates=FALSE,
                                            alpha_g=0.5, penalty_g="L1", kernel_g="gaussian", a1_g=3.7, a2_g=3,
                                            trend_g="monotone", gamma_start_input=NULL, gamma_start_default="zeros",
                                            regressors_gammas=NULL, max_iter_g=10000,
                                            delta_g=1e-5, max_alpha_g=10000,
                                            stepsizeShrink_g=0.8, min_alpha_g=1e-12, convergence_error_g=1e-7,
                                            run_aauc=FALSE){

  #start calc est. time
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

  #check if tau exists when c_function_of_covariates=TRUE
  if(c_function_of_covariates==TRUE){
    if(is.null(tau)){stop("tau cannot be NULL if c_function_of_covariates=TRUE.")}
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
  if (length(tau)>1){stop("PYE_simulation_study_real_data stopped because tau needs to have length 1")}
  if (length(lambda)>1){stop("PYE_simulation_study_real_data stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(penalty_g %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}
  if (n < 2){stop("n needs to be at least 2.")}

  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye (df=df1, X=colnames(df1[, names(df1) %in% c(X,C)]), y="y")$df_scaled
  }

  #generate seeds
  seeds <- c(1:n)

  #names of the columns
  names <- unlist(lapply(seeds, function (x) paste("seed", x, sep="=")))

  #function to be computed
  func <- function(seed, n, df, X, y, C, lambda, tau,
                   beta_start_input,
                   beta_start_default,
                   gamma_start_input,
                   gamma_start_default,
                   trace, alpha,
                   alpha_g,
                   penalty_g,
                   a1, a2,
                   a1_g, a2_g,
                   penalty,
                   regressors_betas,
                   regressors_gammas,
                   trend, kernel,
                   c_zero_fixed,
                   c_function_of_covariates,
                   run_aauc){

    if (trace %in% c(1,2)){
      cat("--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda =", lambda, "; penalty:", penalty, "\n")
      if (c_function_of_covariates==TRUE){
        cat("tau =", tau, "; penalty_g =", penalty_g, "\n")
      }
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
    train_solution <- pye_KS_estimation(df=train_df[,names(train_df)!="ID"], X=X, y=y,
                                        lambda=lambda,
                                        beta_start_input=beta_start_input,
                                        beta_start_default=beta_start_default,
                                        trace=trace, alpha=alpha,
                                        a1=a1, a2=a2, max_iter=max_iter, penalty=penalty,
                                        regressors_betas=regressors_betas, trend=trend,
                                        stepsizeShrink=stepsizeShrink, delta=delta,
                                        max_alpha=max_alpha, min_alpha=min_alpha,
                                        convergence_error=convergence_error,
                                        kernel=kernel, c_zero_fixed=c_zero_fixed)

    estimation_time_original_method <- train_solution$estimation_time
    z_hat <- train_solution$z_hat

    if(c_function_of_covariates == TRUE) {
      #if c_function_of_covariates=TRUE, we compute covYI
      if(length(gamma_start_input) == 0){
        #if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
        gamma_start_input1 <- c(train_solution$c_hat, rep(0, length(C)))
        names(gamma_start_input1) <- c("const", C)
      }

      train_covYI_solution <- covYI_KS_estimation(df=cbind(train_df[,names(train_df)!="ID"], z_hat=train_solution$z_hat[,"z_hat"]),
                                                  z="z_hat", y=y, C=C, tau=tau,
                                                  gamma_start_input=gamma_start_input1,
                                                  gamma_start_default=gamma_start_default, trace=trace,
                                                  alpha=alpha_g, a1=a1_g, a2=a2_g, penalty=penalty_g,
                                                  max_iter=max_iter_g,
                                                  min_alpha=min_alpha_g,
                                                  convergence_error=convergence_error_g,
                                                  regressors_gammas=regressors_gammas,
                                                  trend=trend_g,
                                                  stepsizeShrink=stepsizeShrink_g,
                                                  delta=delta_g, max_alpha=max_alpha_g, kernel=kernel_g,
                                                  run_aauc=run_aauc)

      estimation_time_covYI <- train_covYI_solution$estimation_time
      z_hat <- train_covYI_solution$z_hat

    } else {
      train_covYI_solution <- NULL
    }

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(train_df, y)), z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    #identify the estimated betas
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    c <- getElement(train_solution, "c_hat")
    niter <- train_solution$niter

    #test the results
    test_solution <- pye_KS_with_print(df=test_df[,names(test_df)!="ID"], X=X, y=y, betas=betas, lambda=lambda, c=c,
                                       sim_n=seed, n=n, alpha=alpha, a1=a1, a2=a2, penalty=penalty, est_time=estimation_time_original_method,
                                       niter=niter, kernel=kernel, trace=trace)

    z_hat <- test_solution$z_hat

    if(c_function_of_covariates == TRUE) {
      #put "const" in C
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df=cbind(test_df[,names(train_df)!="ID"], z_hat=test_solution$z_hat[,"z_hat"]),
                                      z="z_hat", y=y, C=C1, gammas=train_covYI_solution$gammas_hat,
                                      tau=tau, kernel=kernel_g, alpha=alpha_g, a1=a1_g, a2=a2_g,
                                      penalty=penalty_g, prediction=TRUE, run_aauc=run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1,2)){
        cat("-> Results on the TEST SET \n")
        cat("-> algorithm: covYI_KS_proximal_gradient_method ; ")
        visualize_gammas = train_covYI_solution$gammas_hat[which(train_covYI_solution$gammas_hat!=0)]
        cat("tau:", tau, "; penalty:", penalty_g, "; covYI_KS:", getElement(test_covYI_solution,  paste0("covYI_KS_",penalty_g)), "; youden_index:", test_covYI_solution$youden_index, "; aYI:", test_covYI_solution$aYI, "; sensitivity:", test_covYI_solution$sensitivity, "; specificity:", test_covYI_solution$specificity, "; geometric_mean:", test_covYI_solution$geometric_mean, "; fdr:", test_covYI_solution$fdr, "; mcc:", test_covYI_solution$mcc, "; auc:", test_covYI_solution$auc, "; aauc:", test_covYI_solution$aauc, "; corrclass:", test_covYI_solution$corrclass, " \n")
        cat("TP:", test_covYI_solution$TP, "; TN:", test_covYI_solution$TN, "; FP:", test_covYI_solution$FP, "; FN:", test_covYI_solution$FN, "; gammas_hat: \n")
        cat("\n")
        print(visualize_gammas)
        cat("\n")
      }
    } else {
      test_covYI_solution <- NULL
    }

    #ROC curve @ certain levels
    est_roc <-  pROC::roc(as.numeric(getElement(test_df, y)), z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_test <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    return(list(estimation_time_original_method=estimation_time_original_method,
                estimation_time_covYI=estimation_time_covYI, seed=seed,
                train_solution=train_solution, train_covYI_solution=train_covYI_solution,
                test_solution=test_solution, test_covYI_solution=test_covYI_solution,
                roc_spec_points = roc_spec_points, roc_sens_points_train = roc_sens_points_train,
                roc_sens_points_test = roc_sens_points_test))
  }

  #measures (train and test)
  pye_L12 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_L1 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_EN <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_SCAD <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_MCP <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aauc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aYI_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  n_total_var_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_betas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros"))
  n_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_betas"))
  n_non_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_betas"))
  n_caught_zero_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_betas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_betas"))
  n_total_var_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_gammas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_gammas"))
  n_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_gammas"))
  n_non_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_gammas"))
  n_caught_zero_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_gammas"))
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(rep(0, ncol(df1[,(names(df1) %in% c(X))])), nrow =  ncol(df1[,(names(df1) %in% c(X))]), ncol = 1, dimnames= list(colnames(df1[,(names(df1) %in% c(X))]), "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(rep(0, ncol(df1[,(names(df1) %in% c(C))])), nrow =  ncol(df1[,(names(df1) %in% c(C))]), ncol = 1, dimnames= list(colnames(df1[,(names(df1) %in% c(C))]), "n_times_gamma_diff_zero"))
  #betas <- vector(mode="list", length=length(seeds))
  #names(betas) <- names

  if (trace %in% c(1,2)){
    cat("------------------------------------------------------------------\n|                   starting simulation study                    |\n------------------------------------------------------------------\n")
  }

  #start cores
  #if no cores are defined lets use the 70% of the available cores
  if (length(used_cores)==0){
    if (parallel::detectCores()>1){used_cores=round(parallel::detectCores()*0.7,0)} else {used_cores=1}
  }

  if (used_cores > 1) {
    cat("I am using this number of cores: used_cores =", used_cores, "\n")
    max.cores <- parallel::detectCores()
    if (used_cores > max.cores) {
      #stop("The number of cores specified (", used_cores,
      warning("The number of cores specified (", used_cores,
           ") is larger than the number of avaiable cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_sim_PYE_real.txt")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for the simulation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("n", "df1", "X", "y", "C", "lambda", "tau",
                                  "beta_start_input", "beta_start_default", "gamma_start_input", "gamma_start_default",
                                  "trace", "a1", "a2", "a1_g", "a2_g",
                                  "alpha", "alpha_g", "penalty", "penalty_g", "trend", "kernel", "c_zero_fixed",
                                  "regressors_betas", "regressors_gammas", "c_function_of_covariates", "run_aauc"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))

    #run
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed=x, n=n, df=df1, X=X, y=y, C=C,
                                                                        lambda=lambda, tau=tau,
                                                                        beta_start_input=beta_start_input,
                                                                        beta_start_default=beta_start_default,
                                                                        gamma_start_input=gamma_start_input,
                                                                        gamma_start_default=gamma_start_default,
                                                                        c_zero_fixed=c_zero_fixed,
                                                                        trace=trace, alpha=alpha,
                                                                        alpha_g=alpha_g,
                                                                        a1=a1, a2=a2,
                                                                        a1_g=a1_g, a2_g=a2_g,
                                                                        penalty=penalty,
                                                                        penalty_g=penalty_g,
                                                                        regressors_betas=regressors_betas,
                                                                        regressors_gammas=regressors_gammas,
                                                                        trend=trend, kernel=kernel,
                                                                        c_function_of_covariates=c_function_of_covariates,
                                                                        run_aauc=run_aauc))

    on.exit(parallel::stopCluster(cl))

  } else {

    #run
    simulation_study <- lapply(seeds, function(x) func(seed=x, n=n, df=df1, X=X, y=y, C=C,
                                                       lambda=lambda, tau=tau,
                                                       beta_start_input=beta_start_input,
                                                       beta_start_default=beta_start_default,
                                                       gamma_start_input=gamma_start_input,
                                                       gamma_start_default=gamma_start_default,
                                                       c_zero_fixed=c_zero_fixed,
                                                       trace=trace, alpha=alpha,
                                                       alpha_g=alpha_g,
                                                       a1=a1, a2=a2,
                                                       a1_g=a1_g, a2_g=a2_g,
                                                       penalty=penalty,
                                                       penalty_g=penalty_g,
                                                       regressors_betas=regressors_betas,
                                                       regressors_gammas=regressors_gammas,
                                                       trend=trend, kernel=kernel,
                                                       c_function_of_covariates=c_function_of_covariates,
                                                       run_aauc=run_aauc))
  }

  estimation_time_original_method <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  #fill the matrices
  #measures on the train set
  temp_pye <- get(paste0("pye_", penalty))
  temp_pye[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("pye_KS_", penalty))))


  list_of_measures <- c("auc", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
  for (i in 1:length(list_of_measures)){
    mes <- get(list_of_measures[i])
    #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][,1], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i]))))
    mes[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if(c_function_of_covariates==TRUE){
    list_of_measures_covYI <- c("auc", "aauc", "aYI", "youden_index", "sensitivity",
                          "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
    for (i in 1:length(list_of_measures_covYI)){
      mes <- get(paste0(list_of_measures_covYI[i],"_covYI"))
      mes[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, list_of_measures_covYI[i])))
      #assign(list_of_measures[i][,1], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      assign(paste0(list_of_measures_covYI[i],"_covYI"), mes)
    }

    n_total_var_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_total_var_gammas")))
    n_predicted_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_zeros_gammas")))
    n_predicted_non_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_non_zeros_gammas")))
    n_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_gammas")))
    n_non_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_non_caught_gammas")))
    n_caught_zero_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_zero_gammas")))
    n_zero_not_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_zero_not_caught_gammas")))
    gammas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))!=0)

  }

  n_total_var_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var")))
  n_predicted_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros")))
  n_predicted_non_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros")))
  n_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_non_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_caught_zero_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero")))
  n_zero_not_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught")))
  betas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty)))!=0)

  #measures on test
  temp_pye[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, paste0("pye_", penalty))))
  assign(paste0("pye_", penalty), temp_pye)
  #measures on test

  for (i in 1:length(list_of_measures)){
    mes <- get(list_of_measures[i])
    #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][,2], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i]))))
    mes[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if(c_function_of_covariates==TRUE){
    for (i in 1:length(list_of_measures_covYI)){
      mes <- get(paste0(list_of_measures_covYI[i],"_covYI"))
      #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      #assign(list_of_measures[i][,2], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      mes[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures_covYI[i])))
      assign(paste0(list_of_measures_covYI[i],"_covYI"), mes)
    }
  }

  #roc curve
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_test)

  betas_start <- getElement(simulation_study[[1]]$train_solution, "betas_start")
  betas <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))
  gammas <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))

  #end computing est. time
  simulation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(simulation_time)
  }

  results <- list(simulation_time=simulation_time,
                  estimation_time_original_method=estimation_time_original_method,
                  estimation_time_covYI=estimation_time_covYI,
                  used_cores=used_cores, n=n, lambda=lambda, tau=tau,
                  pye_L12=pye_L12, pye_L1=pye_L1, pye_EN= pye_EN, pye_SCAD=pye_SCAD, pye_MCP=pye_MCP,
                  auc=auc, youden_index=youden_index, sensitivity=sensitivity,
                  specificity=specificity, geometric_mean=geometric_mean,
                  fdr=fdr, mcc=mcc, corrclass=corrclass,
                  auc_covYI=auc_covYI, aauc_covYI=aauc_covYI, aYI_covYI=aYI_covYI, youden_index_covYI=youden_index_covYI,
                  sensitivity_covYI=sensitivity_covYI, specificity_covYI=specificity_covYI,
                  geometric_mean_covYI=geometric_mean_covYI, fdr_covYI=fdr_covYI, mcc_covYI=mcc_covYI,
                  corrclass_covYI=corrclass_covYI,
                  n_total_var_betas=n_total_var_betas,
                  n_predicted_zeros_betas=n_predicted_zeros_betas,
                  n_predicted_non_zeros_betas=n_predicted_non_zeros_betas,
                  n_caught_betas=n_caught_betas,
                  n_non_caught_betas=n_non_caught_betas,
                  n_caught_zero_betas=n_caught_zero_betas,
                  n_zero_not_caught_betas=n_zero_not_caught_betas,
                  n_total_var_gammas=n_total_var_gammas,
                  n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                  n_predicted_non_zeros_gammas=n_predicted_non_zeros_gammas,
                  n_caught_gammas=n_caught_gammas,
                  n_non_caught_gammas=n_non_caught_gammas,
                  n_caught_zero_gammas=n_caught_zero_gammas,
                  n_zero_not_caught_gammas=n_zero_not_caught_gammas,
                  betas_times_selected=betas_times_selected,
                  gammas_times_selected=gammas_times_selected,
                  betas_start=betas_start, c_zero_fixed=c_zero_fixed,
                  roc_spec_points=roc_spec_points,
                  roc_sens_points_train=roc_sens_points_train,
                  roc_sens_points_test=roc_sens_points_test,
                  regressors_betas=regressors_betas,
                  regressors_gammas=regressors_gammas,
                  c_function_of_covariates=c_function_of_covariates,
                  run_aauc=run_aauc, betas=betas, gammas=gammas)

  return(results)
}


















#' @title PYE_simulation_study_synthetic_data
#'
#' @description function to perform the simulation on n different
#' experiments, for all the considered penalties of pye for simulated data.
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
#' @param trace 2:visualize all the steps, 1:visualize just the result,
#' 0:visualize nothing. Default is 1
#' @param used_cores number of core used for the parallelization of the
#' process. if equal to 1, then no parallelization is adopted. Default is 1.
#' @param c_function_of_covariates if TRUE, covYI is used to estimate the
#' cut-off point as function of the convariate information. If FALSE, the
#' covariate information is ignored. Default is FALSE
#' @param c_zero_fixed if TRUE the estimation process considers c, the cut-off
#' point, as fixed and equal to zero, to reduce the complexity of the estimation.
#' If FALSE c can vary and be different from zero, estimated by pye. Default
#' is FALSE
#' @param beta_start_input vector of a specific starting point for betas.
#' Default is NULL, i.e. no input vector
#' @param beta_start_default set the default starting point of betas. If
#' "zeros", it starts with a vector of all zeros, if "corr" it starts with the
#' value of the correlation of every regressor with the target variable y.
#' Default is "zeros"
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
#' @param alpha parameter for the Elastic-Net penalization term. Default is 0.5
#' @param alpha_g parameter for the Elastic-Net penalization term in covYI.
#' Default is 0.5
#' @param a1 parameter for the SCAD and MCP penalization term. Default is 3.7
#' @param a2 parameter for the MCP penalization term. Default is 3.0
#' @param kernel the kernel type to use for the estimation of the density
#' function (tested only for "gaussian").  Default is "gaussian"
#' @param penalty the considered penalty. To be chosen among L12, L1, EN, SCAD
#' and MCP. Default is "L1"
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
#' @param max_iter_g maximum number of iterations in the algorithms mmAPG and
#' mnmAPG in covYI. Default is 500
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
#' n<- 5 #number of simulations
#' c_function_of_covariates <- TRUE
#' trend <- "monotone" #or "nonmonotone"
#' c_zero_fixed <- FALSE
#'
#' cat("This simulation is on with algorithm: ", trend, "and is c fixed: ", c_zero_fixed, "\n")
#'
#' #--------------------------------------------------------------------------
#' #                                  START                                  |
#' #--------------------------------------------------------------------------
#' #pye Gaussian Kernel Smooth + covYI
#' for (p in used_penalty_pye){
#'   for (mm in 1:length(measures_pye)){
#'     cat("---> Starting with measure", measures_pye[mm], ", for the penalty", p, "<---\n")
#'     name <- paste0("PYE_KS_covYI_", p ,"_sim_study_", measures_pye[mm])
#'     lambda_to_use <- 0.1
#'     tau_to_use <- 0.05
#'     assign(name, PYE_simulation_study_synthetic_data(n=n, rows_train=rows_train,
#'								rows_test=rows_test,  cols=cols,
#'								cols_cov=cols_cov, covar=covar,
#'								mu=mu, mu_cov=mu_cov,
#'                              lambda=lambda_to_use,
#'                              tau=tau_to_use,
#'                              trend=trend,
#'                              beta_start_default=pye_starting_point,
#'                              gamma_start_default=pye_starting_point,
#'                              trace=1,
#'                              alpha=alpha, alpha_g=alpha,
#'                              a1=3.7, a2=3.7,
#'                              a1_g=3.7, a2_g=3.7,
#'                              penalty=p, penalty_g=p,
#'                              kernel="gaussian", used_cores=used_cores,
#'                              c_function_of_covariates=c_function_of_covariates,
#'								c_zero_fixed=c_zero_fixed,
#' 								run_aauc=FALSE))
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
PYE_simulation_study_synthetic_data <- function(n=1000, rows_train=50, rows_test=1000, cols=2000, cols_cov=20,
                                                 covar=0.5, mu=rep(0,cols), mu_cov=rep(0,cols_cov),
                                                 lambda, tau=0, trace=1, used_cores=1,
                                                 c_function_of_covariates=FALSE, c_zero_fixed=FALSE,
                                                 beta_start_input=NULL, beta_start_default="zeros",
                                                 max_iter=10000, trend = "monotone", delta = 1e-5, max_alpha=10000,
                                                 stepsizeShrink=0.8, min_alpha=1e-10, convergence_error=1e-7,
                                                 alpha=0.5, alpha_g=0.5, a1=3.7, a2=3, kernel="gaussian",
                                                 penalty="L12", penalty_g="L1", kernel_g="gaussian", a1_g=3.7, a2_g=3,
                                                 trend_g="monotone", gamma_start_input=NULL, gamma_start_default="zeros",
                                                 max_iter_g=10000, delta_g=1e-5, max_alpha_g=10000,
                                                 stepsizeShrink_g=0.8, min_alpha_g=1e-12, convergence_error_g=1e-7,
												 run_aauc=FALSE){

  #start calc sim. time
  start_time <- Sys.time()

  #starting controls
  if (length(tau)>1){stop("PYE_simulation_study_synthetic_data stopped because tau needs to have length 1")}
  if (length(lambda)>1){stop("PYE_simulation_study_synthetic_data stopped because lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(penalty_g %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}
  if (n < 2){stop("n needs to be at least 2.")}

  #check if tau exists when c_function_of_covariates=TRUE
  if(c_function_of_covariates==TRUE){
    if(is.null(tau)){stop("tau cannot be NULL if c_function_of_covariates=TRUE.")}
  } else if(c_function_of_covariates ==FALSE) {
    tau <- 0
  } else {stop("c_function_of_covariates shoould be equal to TRUE or FALSE.")}

  #generate seeds
  seeds <- c(1:n)

  #names of the columns
  names <- unlist(lapply(seeds, function (x) paste("seed", x, sep="=")))

  #function to be computed
  func <- function(seed, n, rows_train, rows_test, cols, cols_cov, covar,
				   mu, mu_cov,
				   lambda, tau,
                   beta_start_input,
                   beta_start_default,
                   gamma_start_input,
                   gamma_start_default,
                   trace,
                   alpha, alpha_g,
                   a1, a2, a1_g, a2_g,
                   penalty, penalty_g, trend, kernel,
                   max_iter,
                   c_zero_fixed,
                   c_function_of_covariates,
				   run_aauc){

    if (trace %in% c(1,2)){
      cat("--------------------> experiment n",  seed, "of", n, "<---------------------- \n")
      cat("lambda =", lambda, "; penalty:", penalty, "\n")
      if (c_function_of_covariates==TRUE){
        cat("tau =", tau, "; penalty_g =", penalty_g, "\n")
      }
    }

    #create the dataframe
    set.seed(1)
    #choose the regressors. I keep it commented since I keep default ones (below the default code)
    #varsN <- c(1,2)
    #varsB <- c((cols/4+1),(cols/4+2))
    #varsE <- c((cols/2+1),(cols/2+2))
    #varsP <- c((cols/4*3+1),(cols/4*3+2))
    #varsN_cov <- 1
    #varsB_cov <- (cols_cov/4+1)
    #varsE_cov <- (cols_cov/2+1)
    #varsP_cov <- (cols_cov/4*3+1)

    #compute the df (the seed makes different all the simulations)
	  df <- NULL
    df <- create_sample_with_covariates(rows_train=rows_train, cols=cols, cols_cov=cols_cov, covar=covar, mu=mu,
                                        mu_cov=mu_cov, rows_test=rows_test, seed=seed)
                                        #varsN=varsN, varsB=varsB, varsE=varsE, varsP=varsP,
                                        #varN_cov=varN_cov, varB_cov=varB_cov, varE_cov=varE_cov, varP_cov=varP_cov)
    train_df <- df$train_df_scaled
    test_df <- df$test_df_scaled

    X <- df$X
    y <- df$y
    C <- df$C
    regressors_betas <- df$nregressors
    regressors_gammas <- df$ncovariates

    #train PYE KS
    train_solution <- pye_KS_estimation(df=train_df, X=X, y=y, lambda=lambda,
                                        beta_start_input=beta_start_input, beta_start_default=beta_start_default,
                                        trace=trace, alpha=alpha, a1=a1, a2=a2, max_iter=max_iter, penalty=penalty,
                                        regressors_betas=regressors_betas, trend=trend,
                                        stepsizeShrink=stepsizeShrink, delta=delta,
                                        max_alpha=max_alpha, min_alpha=min_alpha,
                                        convergence_error=convergence_error,
                                        kernel=kernel, c_zero_fixed=c_zero_fixed)

    estimation_time_original_method <- train_solution$estimation_time
    z_hat <- train_solution$z_hat

    #here I have to put the covYI estimation
    #computing c
    if(c_function_of_covariates == TRUE) {
      #if c_function_of_covariates=TRUE, we compute covYI
      if(length(gamma_start_input) == 0){
        #if gamma_start_input is not present, we use the optimal c of the betas estimation as the starting point of the constant
        gamma_start_input1 <- c(train_solution$c_hat, rep(0, length(C)))
        names(gamma_start_input1) <- c("const", C)
      }

      train_covYI_solution <- covYI_KS_estimation(df=cbind(train_df[,names(train_df)!="ID"], z_hat=train_solution$z_hat[,"z_hat"]),
                                                  z="z_hat", y=y, C=C, tau=tau,
                                                  gamma_start_input=gamma_start_input1,
                                                  gamma_start_default=gamma_start_default, trace=trace,
                                                  alpha=alpha_g, a1=a1_g, a2=a2_g, penalty=penalty_g,
                                                  max_iter=max_iter_g,
                                                  min_alpha=min_alpha_g,
                                                  convergence_error=convergence_error_g,
                                                  regressors_gammas=regressors_gammas,
                                                  trend=trend_g,
                                                  stepsizeShrink=stepsizeShrink_g,
                                                  delta=delta_g, max_alpha=max_alpha_g, kernel=kernel_g,
                                                  run_aauc=run_aauc)

      estimation_time_covYI <- train_covYI_solution$estimation_time

      z_hat <- train_covYI_solution$z_hat
    } else {
      train_covYI_solution <- NULL
    }

    #ROC curve @ certain levels
    est_roc <- pROC::roc(as.numeric(getElement(train_df, y)), z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_train <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    #identify the estimated betas
    betas <- getElement(train_solution, paste0("betas_hat_", penalty))
    c <- getElement(train_solution, "c_hat")
    niter <- train_solution$niter

    #test the results
    test_solution <- pye_KS_with_print(df=test_df, X=X, y=y, betas=betas,
                                       lambda=lambda, c=c, sim_n=seed, n=n, alpha=alpha, a1=a1, a2=a2,
                                       penalty=penalty, est_time=estimation_time_original_method, niter=niter, kernel=kernel, trace=trace)

    z_hat <- test_solution$z_hat

    if(c_function_of_covariates == TRUE) {
      #put "const" in C
      C1 <- c("const", C)
      test_covYI_solution <- covYI_KS(df=cbind(test_df[,names(train_df)!="ID"], z_hat=test_solution$z_hat[,"z_hat"]),
                                      z="z_hat", y=y, C=C1, gammas=train_covYI_solution$gammas_hat,
                                      tau=tau, kernel=kernel_g, alpha=alpha_g, a1=a1_g, a2=a2_g,
                                      penalty=penalty_g, prediction=TRUE, run_aauc=run_aauc)

      z_hat <- test_covYI_solution$z_hat

      if (trace %in% c(1,2)){
        cat("-> Results on the TEST SET \n")
        cat("-> algorithm: covYI_KS_proximal_gradient_method ; ")
        visualize_gammas = train_covYI_solution$gammas_hat[which(train_covYI_solution$gammas_hat!=0)]
        cat("tau:", tau, "; penalty:", penalty_g, "; covYI_KS:", getElement(test_covYI_solution,  paste0("covYI_KS_",penalty_g)), "; youden_index:", test_covYI_solution$youden_index, "; aYI:", test_covYI_solution$aYI, "; sensitivity:", test_covYI_solution$sensitivity, "; specificity:", test_covYI_solution$specificity, "; geometric_mean:", test_covYI_solution$geometric_mean, "; fdr:", test_covYI_solution$fdr, "; mcc:", test_covYI_solution$mcc, "; auc:", test_covYI_solution$auc, "; aauc:", test_covYI_solution$aauc, "; corrclass:", test_covYI_solution$corrclass, " \n")
        cat("TP:", test_covYI_solution$TP, "; TN:", test_covYI_solution$TN, "; FP:", test_covYI_solution$FP, "; FN:", test_covYI_solution$FN, "; gammas_hat: \n")
        cat("\n")
        print(visualize_gammas)
        cat("\n")
      }
    } else {
      test_covYI_solution <- NULL
    }

    #ROC curve @ certain levels
    est_roc <-  pROC::roc(as.numeric(getElement(test_df, y)), z_hat[,"z_hat"], levels=c(0, 1), direction="<")
    roc_spec_points <- seq(0, 1, by = 0.05)
    roc_sens_points_test <- pROC::coords(est_roc, 1 - roc_spec_points, input="specificity", ret="sensitivity", transpose=TRUE)

    return(list(estimation_time_original_method=estimation_time_original_method,
                estimation_time_covYI=estimation_time_covYI, seed=seed,
                train_solution=train_solution, train_covYI_solution=train_covYI_solution,
                test_solution=test_solution, test_covYI_solution=test_covYI_solution,
                roc_spec_points = roc_spec_points, roc_sens_points_train = roc_sens_points_train,
                roc_sens_points_test = roc_sens_points_test,
                regressors_betas=regressors_betas, regressors_gammas=regressors_gammas))
  }

  #simulate one dataset to take the names of the variables:
  #create the dataframe
  set.seed(1)
  #choose the regressors. I keep it commented since I keep default ones (below the default code)
  #varsN <- c(1,2)
  #varsB <- c((cols/4+1),(cols/4+2))
  #varsE <- c((cols/2+1),(cols/2+2))
  #varsP <- c((cols/4*3+1),(cols/4*3+2))
  #varsN_cov <- 1
  #varsB_cov <- (cols_cov/4+1)
  #varsE_cov <- (cols_cov/2+1)
  #varsP_cov <- (cols_cov/4*3+1)
  #compute the df (the seed makes different all the simulations)
  df_for_the_names <- create_sample_with_covariates(rows_train=rows_train, cols=cols, cols_cov=cols_cov, covar=covar, mu=mu,
                                                    mu_cov=mu_cov, rows_test=rows_test, seed=1)
                                                    #varsN=varsN, varsB=varsB, varsE=varsE, varsP=varsP,
                                                    #varN_cov=varN_cov, varB_cov=varB_cov, varE_cov=varE_cov, varP_cov=varP_cov)

  if (trace %in% c(1,2)){
    cat("Real regressors of betas are: ", df_for_the_names$regressors, "\n")
    if(c_function_of_covariates==TRUE){
      cat("Real regressors of gammas are: ", df_for_the_names$covariates, "\n")
    }
  }

  X <- df_for_the_names$X
  y <- df_for_the_names$y
  C <- df_for_the_names$C

  #measures (train and test)
  pye_L12 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_L1 <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_EN <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_SCAD <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  pye_MCP <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  auc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aauc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  aYI_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  youden_index_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  sensitivity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  specificity_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  geometric_mean_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  fdr_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  mcc_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  corrclass_covYI <- matrix(NA, nrow = length(seeds), ncol = 2, dimnames= list(names,c("train", "test")))
  n_total_var_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_betas"))
  n_predicted_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros"))
  n_predicted_non_zeros_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros"))
  n_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_betas"))
  n_non_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_betas"))
  n_caught_zero_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_betas"))
  n_zero_not_caught_betas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_betas"))
  n_total_var_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_total_var_gammas"))
  n_predicted_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_zeros_gammas"))
  n_predicted_non_zeros_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_predicted_non_zeros_gammas"))
  n_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_gammas"))
  n_non_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_non_caught_gammas"))
  n_caught_zero_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_caught_zero_gammas"))
  n_zero_not_caught_gammas <- matrix(NA, nrow = length(seeds), ncol = 1, dimnames= list(names, "n_zero_not_caught_gammas"))
  #betas_times_selected: how many times the single betas have been selected in the simulation
  betas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))])), nrow =  ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))]), ncol = 1, dimnames= list(colnames(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(X))]), "n_times_beta_diff_zero"))
  #gammas_times_selected: how many times the single betas have been selected in the simulation
  gammas_times_selected <- matrix(rep(0, ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))])), nrow =  ncol(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))]), ncol = 1, dimnames= list(colnames(df_for_the_names$train_df_scaled[,(names(df_for_the_names$train_df_scaled) %in% c(C))]), "n_times_gamma_diff_zero"))
  #betas <- vector(mode="list", length=length(seeds))
  #names(betas) <- names

  if (trace %in% c(1,2)){
    cat("------------------------------------------------------------------\n|                   starting simulation study                    |\n------------------------------------------------------------------\n")
  }

  #start cores
  #if no cores are defined lets use the 70% of the available cores
  if (length(used_cores)==0){
    if (parallel::detectCores()>1){used_cores=round(parallel::detectCores()*0.7,0)} else {used_cores=1}
  }

  if (used_cores > 1) {
    cat("I am using this number of cores: used_cores =", used_cores, "\n")
    max.cores <- parallel::detectCores()
    if (used_cores > max.cores) {
      #stop("The number of cores specified (", used_cores,
      warning("The number of cores specified (", used_cores,
           ") is larger than the number of avaiable cores (",
           max.cores, ")!")
    }
    cl <- parallel::makeCluster(used_cores, outfile="log_sim_PYE_synthetic.txt")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for the simulation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("n", "rows_train", "rows_test", "cols", "cols_cov", "covar", "mu", "mu_cov",
								  "lambda", "tau",
                                  "beta_start_input", "beta_start_default",
                                  "gamma_start_input", "gamma_start_default",
                                  "trace", "max_iter",
                                  "a1", "a2", "a1_g", "a2_g",
                                  "alpha", "alpha_g",
                                  "penalty", "penalty_g", "trend", "kernel", "c_zero_fixed",
                                  "c_function_of_covariates", "run_aauc"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))

    #run
    simulation_study <- parallel::parLapply(cl, seeds, function(x) func(seed=x,  n=n, rows_train=rows_train, rows_test=rows_test,
                                                                        cols=cols, cols_cov=cols_cov, covar=covar,
																		                                    mu=mu, mu_cov=mu_cov,
                                                                        lambda=lambda, tau=tau,
                                                                        beta_start_input=beta_start_input,
                                                                        beta_start_default=beta_start_default,
                                                                        gamma_start_input=gamma_start_input,
                                                                        gamma_start_default=gamma_start_default,
                                                                        max_iter=max_iter,
                                                                        trace=trace, alpha=alpha, alpha_g=alpha_g,
                                                                        a1=a1, a2=a2,
                                                                        a1_g=a1_g, a2_g=a2_g,
                                                                        penalty=penalty, penalty_g=penalty_g, trend=trend,
                                                                        kernel=kernel, c_zero_fixed=c_zero_fixed,
                                                                        c_function_of_covariates=c_function_of_covariates,
																		run_aauc=run_aauc))

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
                                                       max_iter=max_iter,
                                                       trace=trace, alpha=alpha, alpha_g=alpha_g,
                                                       a1=a1, a2=a2,
                                                       a1_g=a1_g, a2_g=a2_g,
                                                       penalty=penalty, penalty_g=penalty_g, trend=trend,
                                                       kernel=kernel, c_zero_fixed=c_zero_fixed,
                                                       c_function_of_covariates=c_function_of_covariates,
													   run_aauc=run_aauc))
  }

  estimation_time_original_method <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time_original_method"))))
  estimation_time_covYI <- mean(unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]], "estimation_time_covYI"))))

  #fill the matrices
  #measures on the train set
  temp_pye <- get(paste0("pye_", penalty))
  temp_pye[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("pye_KS_", penalty))))


  #measures on the train set
  list_of_measures <- c("auc", "youden_index", "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
  for (i in 1:length(list_of_measures)){
    mes <- get(list_of_measures[i])
    #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][,1], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i]))))
    mes[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if(c_function_of_covariates==TRUE){
    #measures on the train set
    list_of_measures_covYI <- c("auc", "aauc", "aYI", "youden_index",
                          "sensitivity", "specificity", "geometric_mean", "fdr", "mcc", "corrclass")
    for (i in 1:length(list_of_measures_covYI)){
      mes <- get(paste0(list_of_measures_covYI[i],"_covYI"))
      mes[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, list_of_measures_covYI[i])))
      #assign(list_of_measures[i][,1], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      assign(paste0(list_of_measures_covYI[i],"_covYI"), mes)
    }

    n_total_var_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_total_var_gammas")))
    n_predicted_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_zeros_gammas")))
    n_predicted_non_zeros_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_predicted_non_zeros_gammas")))
    n_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_gammas")))
    n_non_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_non_caught_gammas")))
    n_caught_zero_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_caught_zero_gammas")))
    n_zero_not_caught_gammas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, "n_zero_not_caught_gammas")))
    gammas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))!=0)

  }

  n_total_var_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_total_var")))
  n_predicted_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_zeros")))
  n_predicted_non_zeros_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_predicted_non_zeros")))
  n_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_betas")))
  n_non_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_non_caught_betas")))
  n_caught_zero_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_caught_zero")))
  n_zero_not_caught_betas[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "n_zero_not_caught")))
  betas_times_selected <- rowSums(sapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, paste0("betas_hat_", penalty)))!=0)

  #measures on test
  temp_pye[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, paste0("pye_", penalty))))
  assign(paste0("pye_", penalty), temp_pye)
  #measures on test

  for (i in 1:length(list_of_measures)){
    mes <- get(list_of_measures[i])
    #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
    #assign(list_of_measures[i][,2], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i]))))
    mes[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_solution, list_of_measures[i])))
    assign(list_of_measures[i], mes)
  }
  if(c_function_of_covariates==TRUE){
    for (i in 1:length(list_of_measures_covYI)){
      mes <- get(paste0(list_of_measures_covYI[i],"_covYI"))
      #auc[,1] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$train_solution, "auc")))
      #assign(list_of_measures[i][,2], unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures[i]))))
      mes[,2] <- unlist(lapply(c(1:length(seeds)), function(x) getElement(simulation_study[[x]]$test_covYI_solution, list_of_measures_covYI[i])))
      assign(paste0(list_of_measures_covYI[i],"_covYI"), mes)
    }
  }

  #roc curve
  roc_spec_points <- simulation_study[[1]]$roc_spec_points
  roc_sens_points_train <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_train)
  roc_sens_points_test <- lapply(c(1:length(seeds)), function(x) simulation_study[[x]]$roc_sens_points_test)

  betas <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_solution, paste0("betas_hat_", penalty)))
  gammas <- lapply(c(1:length(seeds)), function(k) getElement(simulation_study[[k]]$train_covYI_solution, paste0("gammas_hat_", penalty_g)))

  regressors_betas <- simulation_study[[1]]$regressors_betas
  regressors_gammas <- simulation_study[[1]]$regressors_gammas

  #end computing sim. time
  simulation_time = difftime(Sys.time(), start_time , units = "mins")
  if (trace %in% c(1,2)){
    print(simulation_time)
  }

  results <- list(simulation_time=simulation_time,
                  estimation_time_original_method=estimation_time_original_method,
                  estimation_time_covYI=estimation_time_covYI,
                  used_cores=used_cores, n=n, lambda=lambda, tau=tau,
                  pye_L12=pye_L12, pye_L1=pye_L1, pye_EN=pye_EN, pye_SCAD= pye_SCAD, pye_MCP=pye_MCP,
                  auc=auc, youden_index=youden_index,
                  sensitivity=sensitivity, specificity=specificity,
                  geometric_mean=geometric_mean, fdr=fdr,
                  mcc=mcc, corrclass=corrclass,
                  auc_covYI=auc_covYI, aauc_covYI=aauc_covYI, aYI_covYI=aYI_covYI, youden_index_covYI=youden_index_covYI,
                  sensitivity_covYI=sensitivity_covYI, specificity_covYI=specificity_covYI,
                  geometric_mean_covYI=geometric_mean_covYI, fdr_covYI=fdr_covYI, mcc_covYI=mcc_covYI,
                  corrclass_covYI=corrclass_covYI,
                  n_total_var_betas=n_total_var_betas,
                  n_predicted_zeros_betas=n_predicted_zeros_betas,
                  n_predicted_non_zeros_betas=n_predicted_non_zeros_betas,
                  n_caught_betas=n_caught_betas,
                  n_non_caught_betas=n_non_caught_betas,
                  n_caught_zero_betas=n_caught_zero_betas,
                  n_zero_not_caught_betas=n_zero_not_caught_betas,
                  n_total_var_gammas=n_total_var_gammas,
                  n_predicted_zeros_gammas=n_predicted_zeros_gammas,
                  n_predicted_non_zeros_gammas=n_predicted_non_zeros_gammas,
                  n_caught_gammas=n_caught_gammas,
                  n_non_caught_gammas=n_non_caught_gammas,
                  n_caught_zero_gammas=n_caught_zero_gammas,
                  n_zero_not_caught_gammas=n_zero_not_caught_gammas,
                  betas_times_selected=betas_times_selected,
                  gammas_times_selected=gammas_times_selected,
                  roc_spec_points=roc_spec_points,
                  roc_sens_points_train=roc_sens_points_train,
                  roc_sens_points_test=roc_sens_points_test,
                  regressors_betas=regressors_betas,
                  regressors_gammas=regressors_gammas,
                  c_function_of_covariates=c_function_of_covariates,
				          run_aauc=run_aauc, betas=betas, gammas=gammas)

  return(results)
}
