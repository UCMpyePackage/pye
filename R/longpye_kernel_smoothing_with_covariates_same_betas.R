#' @title longpye_KS_same_betas
#'
#' @description The Penalizad Youden Index (PYE) function for longitudinal data
#' based on the Kernel Smooth density estimator. It not only evaluate PYE in
#' logitudinal case, but it returns all the necessary for the estimation
#' process, like measure of fit and derivatives. It works for all the
#' considered penalties (L12, L1, EN, SCAD and MCP)
#'
#' @param df the input dataset
#' @param X regressors to consider in the estimation. It can be of type
#' dataframe, containing also the same name of the regressors included in df,
#' of just a vector of character. Default is all not present in y
#' @param y the target variable. It can be only binomial 0,1. It can be of type
#' dataframe, containing also the same name of the same target variable included
#' in df, or just a character. Default is "y".
#' @param t is the name of the variable that refers to the time in the input
#' dataset.
#' @param betas the coefficients of the biomarker combination used for the
#' evaluation of PYE
#' @param lambda the penalization parameter of the regressors X
#' @param c the cut-off points, a vector having the length of the times
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
#' @param print.CDF.plot if TRUE it prints also the plot of the CDFs of cases
#' and controls for the compination of regressors Z for each dimension. Default
#' is FALSE
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
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
#'
#' PYE_result <- pye_KS(df=df[,names(df) %in% c(X,y)], X=X, y=y, betas=betas,
#'   lambda=lambda, c=c, alpha=0.5, a1=3.7, a2=3, penalty=penalty)
#' print(PYE_result)
#'
#' @importFrom evmix kdz
#' @importFrom OptimalCutpoints optimal.cutpoints
#' @importFrom plyr join
#' @export
longpye_KS_same_betas <- function(df_train, df_test=NULL, id="id", X=names(df_train[,!(names(df_train) %in% c(y,id,t))]),
                       y="y", t="t", betas, lambda,
                       c=rep(0,length(unique(df_train[,t]))),
                       kernel="gaussian", alpha=0.5, a1=3.7, a2=3,
                       penalty="L1", h_exponent=0.2, prediction=FALSE,
                       print.CDF.plot=FALSE, c_hat_discriminant=NULL){

  #check df_train:
  if (nrow(df_train)==0){stop("df_train has no rows")}

  #checks
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}
  if (class(t) != "character"){stop("t can only be of class character or data.frame.")}
  if (class(c) == "character"){ c <- df_train[,c]
  } else if ((length(c) == 1) & (sum(nrow(c))<2)){ #it is just a single element
    c <- rep(c, length(unique(df_train[,t])))
    names(c) <- order(unique(df_train[,t]))
  } else if (length(c) == length(unique(df_train[,t]))){
    c <- c
    names(c) <- order(unique(df_train[,t]))
  } else {stop("c can only be of class numeric on length 1 (same value repeated) or a matrix with t columns and the same number rows as df_train")}

  #check if ID already exists in the dataset
  if("ID_rows" %in% colnames(df_train)){stop("ID_rows already exists as column in df_train! Please delete or rename this column since I need this name to set the internal ID")}
  if("c" %in% colnames(df_train)){stop("c already exists as column in df_train! Please delete or rename this column since I need this name to set the internal column")}

  #Create a new df considering only the columns included in X (the regressors to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  ID_rows <- rownames(df_train)
  df1 <- cbind(ID_rows, id=df_train[,id], t=df_train[,t], df_train[,(names(df_train) %in% c(y)), drop=FALSE], df_train[,(names(df_train) %in% c(X)), drop=FALSE])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("longpye_KS_same_betas stopped becuase lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}

  if (!(kernel %in% c("gaussian", "normal", "uniform", "rectangular", "triangular",
                      "epanechnikov", "biweight", "triweight", "tricube", "parzen",
                      "cosine", "optcosine"))){
    stop("kernel parameter is not in the available options. Options are: gaussian,
          normal, uniform, rectangular, triangular, epanechnikov, biweight,
          triweight, tricube, parzen, cosine, optcosine")
  }

  if (length(betas) != length(X)){
    stop("The number of element of betas is different then the number of columns in df_train")
  }

  #divide control and diseased and multiply by betas
  n_times <- length(unique(df1$t))
  z_y0_i <- list()
  z_y1_i <- list()
  z_i <- list()

  # Vectorized calculations
  z_y0_i <- lapply(1:n_times, function (i) stats::setNames(transform(merge(df1["id"],
    as.matrix(df1[df1$t == i & df1$y == 0, names(df1) %in% X, drop=FALSE])%*%betas, by=0), row.names=Row.names, Row.names=NULL), c("id", paste0("z_", i, "_hat"))))
  z_y1_i <- lapply(1:n_times, function (i) stats::setNames(transform(merge(df1["id"],
    as.matrix(df1[df1$t == i & df1$y == 1, names(df1) %in% X, drop=FALSE])%*%betas, by=0), row.names=Row.names, Row.names=NULL), c("id", paste0("z_", i, "_hat"))))

  #append the z(t)
  z_i <- lapply(1:n_times, function (i) stats::setNames(transform(merge(df1["id"],
              as.matrix(df1[df1$t == i, names(df1) %in% X, drop=FALSE])%*%betas, by=0), row.names=Row.names), c("ID_rows", "id", "z_i_hat")))

  #add z_hat and c to the original dataframe
  df2 <- merge(df1, do.call(rbind, z_i)[,c("ID_rows", "z_i_hat")], by = 'ID_rows')
  df2["c"] <- NULL
  df2 <- merge(df2, base::data.frame("t" = as.integer(names(c)), c), by = "t")[, union(names(df2), names(base::data.frame("t" = as.integer(names(c)), c)))]
  rownames(df2) <- df2$ID_rows

  z_y0_matrixNA <- Reduce(function(x, y) merge(x, y, by="id", all.x=TRUE, all.y=TRUE), z_y0_i)
  z_y0_matrix <- z_y0_matrixNA
  z_y0_matrix[is.na(z_y0_matrix)] <- 0
  z_y1_matrixNA <- Reduce(function(x, y) merge(x, y, by="id", all.x=TRUE, all.y=TRUE), z_y1_i)
  z_y1_matrix <- z_y1_matrixNA
  z_y1_matrix[is.na(z_y1_matrix)] <- 0





  #frames for the derivatives
  #OLD:
  #list_of_features_y0 <- lapply(1:length(X), function(x){
  #                                #working frames for the derivatives
  #                                as.data.frame(sapply(1:n_times, function (i) {
  #                                  temp <- df2[((df2$y == 0) & (df2$t == 1)),id, drop=FALSE] #keep the same number of element of time 1
  #                                  temp <- plyr::join(temp, df2[((df2$y == 0) & (df2$t == i)), names(df2) %in% c(id, X[x]), drop=FALSE], by=id)
  #                                  names(temp)[2] <- paste0(X[x], "_", i)
  #                                  temp[is.na(temp)] <- 0
  #                                  if(i == 1){
  #                                    temp
  #                                  } else {
  #                                    temp[paste0(X[x], "_", i)]
  #                                  }
  #                                }))
  #                              })
  #optimized
  list_of_features_y0 <- lapply(X, function(x) {
                                  # Subset data for y == 0 and t == 1
                                  temp_base <- df2[df2$y == 0 & df2$t == 1, c("id", x), drop = FALSE]

                                  # Initialize an empty data frame
                                  result <- data.frame(id = temp_base$id)

                                  # Loop over time points
                                  for (i in 1:n_times) {
                                    temp_i <- df2[df2$y == 0 & df2$t == i, c("id", x), drop = FALSE]
                                    col_name <- paste0(x, "_", i)

                                    # Left join with previous result
                                    result <- merge(result, temp_i, by = "id", all.x = TRUE)
                                    names(result)[names(result) == x] <- col_name

                                    # Replace NA values with 0
                                    result[is.na(result)] <- 0
                                  }

                                  result
                              })
  #OLD:
  #list_of_features_y1 <- lapply(1:length(X), function(x){
  #                                #working frames for the derivatives
  #                                as.data.frame(sapply(1:n_times, function (i) {
  #                                  temp <- df2[((df2$y == 1) & (df2$t == 1)),id, drop=FALSE] #keep the same number of element of time 1
  #                                  temp <- plyr::join(temp, df2[((df2$y == 1) & (df2$t == i)), names(df2) %in% c(id, X[x]), drop=FALSE], by=id)
  #                                  names(temp)[2] <- paste0(X[x], "_", i)
  #                                  temp[is.na(temp)] <- 0
  #                                  if(i == 1){
  #                                    temp
  #                                  } else {
  #                                    temp[paste0(X[x], "_", i)]
  #                                  }
  #                                }))
  #                              })
  #optimized
  list_of_features_y1 <- lapply(X, function(x) {
                                    # Subset data for y == 1 and t == 1
                                    temp_base <- df2[df2$y == 1 & df2$t == 1, c("id", x), drop = FALSE]

                                    # Initialize an empty data frame
                                    result <- data.frame(id = temp_base$id)

                                    # Loop over time points
                                    for (i in 1:n_times) {
                                      temp_i <- df2[df2$y == 1 & df2$t == i, c("id", x), drop = FALSE]
                                      col_name <- paste0(x, "_", i)

                                      # Left join with previous result
                                      result <- merge(result, temp_i, by = "id", all.x = TRUE)
                                      names(result)[names(result) == x] <- col_name

                                      # Replace NA values with 0
                                      result[is.na(result)] <- 0
                                    }

                                    result
                                })


  Gamma0 <- cov(z_y0_matrixNA[-1], use = 'pairwise.complete.obs', method="pearson")
  n0 <- colSums(!is.na(z_y0_matrixNA[-1]))

  #if Gamma0 is a zero matrix we equal it to a diag(1) matrix
  if (all(Gamma0==0)){
    Gamma0 <- diag(1, n_times)
  }

  #if not positive-definite
  if(!matrixcalc::is.positive.definite(Gamma0, tol=1e-8)){
    #tranform the matrix in positive-definte:
    #method 1)
    if (!require(Matrix)) {
      install.packages("Matrix")
    }

    # Find the nearest positive definite matrix
    Gamma0 <- as.matrix(Matrix::nearPD(Gamma0)$mat)

    #method 2)
    # Compute the eigenvalues and eigenvectors
    #eigen_decomp <- eigen(Gamma0)
    # Make a diagonal matrix with the eigenvalues, replacing negative ones with a small positive number
    #eigen_values <- pmax(eigen_decomp$values, .Machine$double.eps)
    # Recreate the covariance matrix
    #pd_matrix2 <- eigen_decomp$vectors %*% diag(eigen_values) %*% t(eigen_decomp$vectors)
  }

  Gamma0_sqrt <- expm::sqrtm(Gamma0)
  H0 <- n0^(-1/(ncol(Gamma0)+4))*Gamma0_sqrt*diag(1, ncol(Gamma0_sqrt)) #+ sign(Gamma0) #<<- we added a regularization term equal to the sign of the element to avoind having an H0_inv matrix with values too high
  if(sum(abs(H0))<0.1){H0 <- sign(Gamma0_sqrt)*diag(0.5, ncol(Gamma0_sqrt))} #if the elements of H0 are too small, we increase their size
  #H0 <- diag(1, ncol(Gamma0_sqrt)) #<-test
  H0_inv <- solve(H0) #inverse of H0
  H0_det <- det(H0)

  #if(h0<0.1) {h0=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  #x0 <- as.matrix(-sweep(z_y0_matrix[!(names(z_y0_matrix) %in% id)], 2, c))

  t0 <- as.matrix(-sweep(z_y0_matrix[!(names(z_y0_matrix) %in% "id")], 2, c)) %*% H0_inv

  #multivariate normal distribution used as kernel
  #f0 <- sum(mnormt::pmnorm(x=t0, mean=rep(0,ncol(t0)), varcov=diag(1, ncol(t0)), maxpts=10000, abseps=1e-10)/determinant(H0)$modulus)/nrow(t0)
  f0 <- sum(mnormt::pmnorm(x=t0, mean=rep(0,ncol(t0)), varcov=diag(1, ncol(t0)), maxpts=10000, abseps=1e-10))/nrow(t0)

  #diseased
  Gamma1 <- cov(z_y1_matrixNA[-1], use = 'pairwise.complete.obs', method="pearson")
  n1 <- Matrix::colSums(!is.na(z_y1_matrixNA[-1]))

  #if Gamma0 is a zero matrix we equal it to a diag(1) matrix
  if (all(Gamma1==0)){
    Gamma1 <- diag(1, n_times)
  }

  #if not positive-definite
  if(!matrixcalc::is.positive.definite(Gamma1, tol=1e-8)){
    #tranform the matrix in positive-definte:
    #method 1)
    if (!require(Matrix)) {
      install.packages("Matrix")
    }
    # Find the nearest positive definite matrix
    Gamma1 <- as.matrix(Matrix::nearPD(Gamma1)$mat)

    #method 2)
    # Compute the eigenvalues and eigenvectors
    #eigen_decomp <- eigen(Gamma0)
    # Make a diagonal matrix with the eigenvalues, replacing negative ones with a small positive number
    #eigen_values <- pmax(eigen_decomp$values, .Machine$double.eps)
    # Recreate the covariance matrix
    #pd_matrix2 <- eigen_decomp$vectors %*% diag(eigen_values) %*% t(eigen_decomp$vectors)
  }

  Gamma1_sqrt <- expm::sqrtm(Gamma1)
  H1 <- n1^(-1/(ncol(Gamma1)+4))*Gamma1_sqrt*diag(1, ncol(Gamma1_sqrt)) #+ sign(Gamma1) #<<- we added a regularization term equal to the sign of the element to avoind having an H1_inv matrix with values too high
  if(sum(abs(H1))<0.1){H1 <- sign(Gamma1_sqrt)*diag(0.5,ncol(Gamma1_sqrt))} #if the elements of H0 are too small, we increase their size
  #H1 <- diag(1, ncol(Gamma1_sqrt)) #<-test
  H1_inv <- solve(H1) #inverse of H1
  H1_det <- det(H1)

  #if(h0<0.1) {h0=0.1} #cannot divide for 0 or a number too small (does not makes sense)
  #x1 <- as.matrix(-sweep(z_y1_matrix[!(names(z_y1_matrix) %in% id)], 2, c))
  t1 <- as.matrix(-sweep(z_y1_matrix[!(names(z_y1_matrix) %in% "id")], 2, c)) %*% H1_inv

  #multivariate normal distribution used as kernel
  #f1 <- sum(mnormt::pmnorm(x=t1, mean=rep(0,ncol(t1)), varcov=diag(1, ncol(t1)), maxpts=10000, abseps=1e-10)/determinant(H1)$modulus)/nrow(t1)
  f1 <- sum(mnormt::pmnorm(x=t1, mean=rep(0,ncol(t1)), varcov=diag(1, ncol(t1)), maxpts=10000, abseps=1e-10))/nrow(t1)

  #youden index
  yi <- round(f0 - f1, 15)

  #evaluate the gradient
  #control
  #old:
  #gr0_betas_i <- sapply(list_of_features_y0, function (x) mean((t(mnormt::dmnorm(t0, rep(0,ncol(t0)), varcov=diag(1, ncol(t0)))) %*% as.matrix(-x[!(names(x) %in% id)]) %*% H0_inv)/(nrow(t0)*H0_det)))
  #gr0_c_i <- (t(mnormt::dmnorm(x=t0, rep(0,ncol(t0)), varcov=diag(1, ncol(t0)))) %*% matrix(1, nrow = nrow(t0), ncol = ncol(t0)) %*% H0_inv)/(nrow(t0)*H0_det)
  #colnames(gr0_c_i)=names(c)
  #new:
  #test:
  #numDeriv::grad(func0, t0[1,])
  #a1 <- mnormt::dmnorm(t0[1,1], rep(0,1), varcov=diag(1, 1))
  #a2 <- mnormt::pmnorm(t0[1,2], rep(0,1), varcov=diag(1, 1))
  #a3 <- mnormt::pmnorm(t0[1,3], rep(0,1), varcov=diag(1, 1))
  #a4 <- mnormt::pmnorm(t0[1,4], rep(0,1), varcov=diag(1, 1))
  #a1*a2*a3*a4
  #it works!
  #new/old:
  #p_grad0_normal_betas <- base::Reduce(`*`,lapply(1:ncol(t0), function (i) mnormt::dmnorm(t0[,i], rep(0,1), varcov=diag(1, 1))))
  #gr0_betas_i <- sapply(list_of_features_y0, function (x) {
  #                      mean((p_grad0_normal_betas %*% -(as.matrix(x[!(names(x) %in% id)]) %*% H0_inv))/(nrow(t0)*H0_det))
  #                    })
  #
  #p_grad0_normal_c <- do.call(rbind,lapply(1:nrow(t0), function(x){
  #                      func0 <- function(z) {mnormt::pmnorm(z, rep(0,ncol(t0)), varcov=diag(1, ncol(t0)))}
  #                      numDeriv::grad(func0, t0[x,])
  #                    }))
  #gr0_c_i <- (colSums(p_grad0_normal_c) %*% H0_inv)/(nrow(t0)*H0_det)
  #colnames(gr0_c_i)=names(c)
  #New/New:
  func0 <- function(z) {mnormt::pmnorm(z, rep(0, ncol(t0)), varcov = diag(1, ncol(t0))) }
  p_grad0_normal_c <- t(apply(t0, 1, numDeriv::grad, func = func0))
  #p_grad0_normal_betas <- mnormt::dmnorm(t0, rep(0, ncol(t0)), varcov = diag(1, ncol(t0)))
  gr0_betas_i <- sapply(list_of_features_y0, function(x) {
                          mean(colSums(p_grad0_normal_c * -(as.matrix(x[!(names(x) %in% "id")]) %*% H0_inv)) / (nrow(t0) * H0_det))
                          #sum((p_grad0_normal_betas %*% -(as.matrix(x[!(names(x) %in% "id")]) %*% H0_inv)) / (nrow(t0) * H0_det))
                        })

  #func0 <- function(z) {mnormt::pmnorm(z, rep(0, ncol(t0)), varcov = diag(1, ncol(t0))) }
  #p_grad0_normal_c <- t(apply(t0, 1, numDeriv::grad, func = func0))
  gr0_c_i <- (colSums(p_grad0_normal_c) %*% H0_inv) / (nrow(t0) * H0_det)
  colnames(gr0_c_i) <- names(c)

  #cases
  #old:
  #gr1_betas_i <- sapply(list_of_features_y1, function (x) mean((t(mnormt::dmnorm(t1, rep(0,ncol(t1)), varcov=diag(1, ncol(t1)))) %*% as.matrix(-x[!(names(x) %in% id)]) %*% H1_inv)/(nrow(t1)*H1_det)))
  #gr1_c_i <- (t(mnormt::dmnorm(x=t1, rep(0,ncol(t1)), varcov=diag(1, ncol(t1)))) %*% matrix(1, nrow = nrow(t1), ncol = ncol(t1)) %*% H1_inv)/(nrow(t1)*H1_det)
  #colnames(gr1_c_i)=names(c)

  #new/old:
  #p_grad1_normal_betas <- base::Reduce(`*`,lapply(1:ncol(t1), function (i) mnormt::dmnorm(t1[,i], rep(0,1), varcov=diag(1, 1))))
  #gr1_betas_i <- sapply(list_of_features_y1, function (x) {
  #                      mean((p_grad1_normal_betas %*% -(as.matrix(x[!(names(x) %in% id)]) %*% H1_inv)))/(nrow(t1)*H1_det)
  #                     })
  #p_grad1_normal_c <- do.call(rbind,lapply(1:nrow(t1), function(x){
  #                      func1 <- function(z) {mnormt::pmnorm(z, rep(0,ncol(t1)), varcov=diag(1, ncol(t1)))}
  #                      numDeriv::grad(func1, t1[x,])
  #                    }))
  #gr1_c_i <- (colSums(p_grad1_normal_c) %*% H1_inv)/(nrow(t1)*H1_det)
  #colnames(gr1_c_i)=names(c)

  #New/New:
  func1 <- function(z) {mnormt::pmnorm(z, rep(0, ncol(t1)), varcov = diag(1, ncol(t1))) }
  p_grad1_normal_c <- t(apply(t1, 1, numDeriv::grad, func = func1))
  #p_grad1_normal_betas <- mnormt::dmnorm(t1, rep(0, ncol(t1)), varcov = diag(1, ncol(t1)))
  gr1_betas_i <- sapply(list_of_features_y1, function(x) {
                          mean(colSums(p_grad1_normal_c * -(as.matrix(x[!(names(x) %in% "id")]) %*% H1_inv)) / (nrow(t1) * H1_det))
                          #sum((p_grad1_normal_betas %*% -(as.matrix(x[!(names(x) %in% "id")]) %*% H1_inv)) / (nrow(t1) * H1_det))
                        })

  #func1 <- function(z) {mnormt::pmnorm(z, rep(0, ncol(t1)), varcov = diag(1, ncol(t1))) }
  #p_grad1_normal_c <- t(apply(t1, 1, numDeriv::grad, func = func1))
  gr1_c_i <- (colSums(p_grad1_normal_c) %*% H1_inv) / (nrow(t1) * H1_det)
  colnames(gr1_c_i) <- names(c)

  #gradient
  #gr_yi0 <- mapply('-', c(gr0_betas_i, list(gr0_c_i)), c(gr1_betas_i, list(gr1_c_i)), SIMPLIFY = FALSE)
  #gr_yi <- t(sapply(1:length(gr_yi0), function (x) gr_yi0[[x]]))
  gr_yi_betas <- gr0_betas_i - gr1_betas_i
  gr_yi_c <- as.numeric(gr0_c_i - gr1_c_i)
  names(gr_yi_betas) <- X
  names(gr_yi_c) <- names(c)

  df2["y_hat"] = ifelse(df2$z_i_hat >= df2$c,1,0)

  #prepare the solution in all the scenarios
  z_hat_test<-NA
  df_results_by_id<-NA
  TP_test<-NA
  TN_test<-NA
  FP_test<-NA
  FN_test<-NA
  spec_test<-NA
  fnr_test<-NA
  sensitivity_test<-NA
  gm_test<-NA
  fdr_test<-NA
  mcc_test<-NA
  auc_test<-NA
  yi_test<-NA
  corrclass_test<-NA
  youden_index_test<-NA
  z_matrix_test<-NA
  df_pred<-NA
  f0_test<-NA
  f1_test<-NA
  post_prob1_test<-NA
  post_prob0_test<-NA
  ratio_test<-NA

  #if we have a df_test we compute the probabilities of being healthy or diseased for each subject
  if (length(df_test) > 0) {
    ID_rows_test <- rownames(df_test)
    df1_test <- cbind(ID_rows_test, id=df_test[,id], t=df_test[,t], df_test[,(names(df_test) %in% c(y,X)), drop=FALSE])

    z_i_test <- list()

    # Vectorized calculations
    z_i_test <- lapply(1:n_times, function (i) stats::setNames(transform(merge(df1_test["id"],
                    as.matrix(df1_test[df1_test$t == i, names(df1_test) %in% X, drop=FALSE])%*%betas, by=0), row.names=Row.names, Row.names=NULL), c("id", paste0("z_", i, "_hat"))))

    z_matrixNA_test <- Reduce(function(x, y) merge(x, y, by="id", all.x=TRUE, all.y=TRUE), z_i_test)
    z_matrix_test <- z_matrixNA_test
    z_matrix_test[is.na(z_matrix_test)] <- 0

    #example
    #valuetest <- c(-1,-1,-1,-1)
    #test0 <- as.matrix(-sweep(z_y0_matrix[!(names(z_y0_matrix) %in% id)], 2, valuetest)) %*% H0_inv
    #f0_test <- sum(mnormt::dmnorm(x=test0, mean=rep(0,ncol(test0)), varcov=diag(1, ncol(test0))))/nrow(test0)

    f0_test <- apply(z_matrix_test[,-1], 1, function(x){
      t0_test <- as.matrix(-sweep(z_y0_matrix[!(names(z_y0_matrix) %in% "id")], 2, x)) %*% H0_inv
      sum(mnormt::dmnorm(x=t0_test, mean=rep(0,ncol(t0_test)), varcov=diag(1, ncol(t0_test))))/nrow(t0_test)
      })

    f1_test <- apply(z_matrix_test[,-1], 1, function(x){
      t1_test <- as.matrix(-sweep(z_y1_matrix[!(names(z_y1_matrix) %in% "id")], 2, x)) %*% H1_inv
      sum(mnormt::dmnorm(x=t1_test, mean=rep(0,ncol(t1_test)), varcov=diag(1, ncol(t1_test))))/nrow(t1_test)
      })

    #discriminant analysis
    p0 <- length(t0) / (length(t0) + length(t1))
    p1 <- length(t1) / (length(t0) + length(t1))
    post_prob1 <- f1_test*p1 / (f0_test*p0 + f1_test*p1)
    post_prob0 <- f0_test*p0 / (f0_test*p0 + f1_test*p1)
    ratio <- log(post_prob1/ post_prob0)
    #if y is available in df_test we compute the performance measures
    if (y %in% colnames(df1_test)){
      z_matrix_test <- merge(unique(df1_test[c("id", y)]), z_matrix_test, by="id")
      df_pred <- cbind(z_matrix_test, f0_test=f0_test, f1_test=f1_test, post_prob1=post_prob1, post_prob0=post_prob0, ratio=ratio)
      opt <- OptimalCutpoints::optimal.cutpoints(data=df_pred, X = "ratio", status = "y", methods = "Youden", tag.healthy =0)
      auc_test <- opt$Youden$Global$measures.acc$AUC[1]
      yi_test <- opt$Youden$Global$optimal.criterion[1]
      if(length(c_hat_discriminant)==0){
        c_hat_discriminant <- opt$Youden$Global$optimal.cutoff$cutoff[1] #if(length(opt$Youden$Global$optimal.cutoff$cutoff)==1){opt$Youden$Global$optimal.cutoff$cutoff} else {opt$Youden$Global$optimal.cutoff$cutoff[round(length(opt$Youden$Global$optimal.cutoff$cutoff)/2)]}
      }
      pred <- ifelse(ratio >= c_hat_discriminant, 1,0)
      df_pred <- cbind(df_pred, c_hat_discriminant =c_hat_discriminant, y_hat_test=pred)
      TP_test <- sum(ifelse(df_pred[df_pred[y] == 1, "ratio"] >= c_hat_discriminant ,1,0), na.rm=TRUE)
      TN_test <- sum(ifelse(df_pred[df_pred[y] == 0, "ratio"] < c_hat_discriminant ,1,0), na.rm=TRUE)
      FP_test <- sum(ifelse(df_pred[df_pred[y] == 0, "ratio"] >= c_hat_discriminant ,1,0), na.rm=TRUE)
      FN_test <- sum(ifelse(df_pred[df_pred[y] == 1, "ratio"] < c_hat_discriminant ,1,0), na.rm=TRUE)

      #compute the performance measures
      spec_test  <- TN_test/(TN_test+FP_test)
      fnr_test  <- FN_test/(FN_test+TP_test)

      #sensitivity
      sensitivity_test  = 1-fnr_test

      #geometric mean
      gm_test  <- sqrt(spec_test *sensitivity_test)

      #compute other measures: FDR:
      fdr_test  <- FP_test/(FP_test + TP_test)

      #compute other measures: FDR:
      mcc_test  <- ((TP_test*TN_test)-(FP_test*FN_test))/sqrt((TP_test+FP_test)*(TP_test+FN_test)*(TN_test+FP_test)*(TN_test+FN_test))

      #compute the correct classification:
      corrclass_test <- (TP_test+TN_test)/(TP_test+TN_test+FP_test+FN_test)

      if(print.CDF.plot == TRUE){
        #lets print the empirical CDFs to visualize the result
        z_y0 <- df_pred[df_pred[y]==0, ][,"ratio"]
        z_y1 <- df_pred[df_pred[y]==1, ][,"ratio"]
        z_y0_ord <- z_y0[order(z_y0)]
        z_y1_ord <- z_y1[order(z_y1)]
        #correction of the outliers just for plotting
        z_y0_ord[z_y0_ord %in% boxplot.stats(z_y0_ord)$out] <- min(boxplot.stats(z_y0_ord)$stats[5])
        z_y1_ord[z_y1_ord %in% boxplot.stats(z_y1_ord)$out] <- min(boxplot.stats(z_y1_ord)$stats[5])
        ecdf0 <- ecdf(z_y0_ord)
        ecdf1 <- ecdf(z_y1_ord)
        plot_data <- data.frame(x = c(z_y0_ord, z_y1_ord),
                                y = c(ecdf0(z_y0_ord), ecdf1(z_y1_ord)),
                                group = c(rep("CDF_Y0",length(z_y0_ord)), rep("CDF_Y1", length(z_y1_ord))))

        # Crea il plot
        print(ggplot2::ggplot(plot_data, ggplot2::aes(x=x, y=y, color=group)) +
                ggplot2::geom_line() +
                ggplot2::labs(x="log[Pr(Y=1)/Pr(Y=0)]", y = "CDF of the log-posterior probability ratiofor case and conntrol patients", color="Groups") +
                ggplot2::theme_minimal())
      }
    } else {
      df_pred <- cbind(z_matrix_test, f0=f0_test, f1=f1_test, post_prob1=post_prob1, post_prob0=post_prob0, ratio=ratio)
    }
    f0_test<-df_pred[,c("id", "f0_test")]
    f1_test<-df_pred[,c("id", "f1_test")]
    post_prob1_test<-df_pred[,c("id", "post_prob1")]
    post_prob0_test<-df_pred[,c("id", "post_prob0")]
    ratio_test<-df_pred[,c("id", "ratio")]
  }



  #confusion matrix
  TP <- as.numeric(sum(z_y1_matrixNA[-1] >= c, na.rm=TRUE))
  TN <- as.numeric(sum(z_y0_matrixNA[-1] < c, na.rm=TRUE))
  FP <- as.numeric(sum(z_y0_matrixNA[-1] >= c, na.rm=TRUE))
  FN <- as.numeric(sum(z_y1_matrixNA[-1] < c, na.rm=TRUE))

  TP_i <- colSums(z_y1_matrixNA[-1] >= c, na.rm=TRUE)
  TN_i <- colSums(z_y0_matrixNA[-1] < c, na.rm=TRUE)
  FP_i <- colSums(z_y0_matrixNA[-1] >= c, na.rm=TRUE)
  FN_i <- colSums(z_y1_matrixNA[-1] < c, na.rm=TRUE)
  Tot_i <- colSums(!is.na(z_y0_matrixNA[-1]), na.rm=TRUE) + colSums(!is.na(z_y1_matrixNA[-1]), na.rm=TRUE)

  #compute the Youden index
  spec <- TN/(TN+FP)
  fnr <- FN/(FN+TP)

  #sensitivity
  sensitivity = 1-fnr

  #geometric mean
  gm <- sqrt(spec*sensitivity)

  #compute other measures: FDR:
  fdr <- FP/(FP + TP)

  #compute other measures: MCC:
  mcc <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  #compute the correct classification:
  corrclass <- (TP+TN)/(TP+TN+FP+FN)
  corrclass_i <- (TP_i+TN_i)/Tot_i

  #AUC e YI
  #find the optimal cut-point
  opt <- OptimalCutpoints::optimal.cutpoints(data=df2, X="z_i_hat", status="y", methods="Youden", tag.healthy=0)
  auc <- opt$Youden$Global$measures.acc$AUC[1]

  #if all the betas are zero, the measures are zeros
  if(sum(betas)==0){
    spec <- 0
    fnr <- 0
    sensitivity <- 0
    gm <- 0
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
  longpye_L12  <- yi - phi_L12
  longpye_L1   <- yi - phi_L1
  longpye_EN   <- yi - phi_EN
  longpye_SCAD <- yi - phi_SCAD
  longpye_MCP  <- yi - phi_MCP

  #gr_pye_L12 <- c(gr_yi[-length(gr_yi)] - gr_phi_L12, gr_yi[length(gr_yi)])
  #gr_pye_L1 <- c(gr_yi[-length(gr_yi)] - gr_phi_L1, gr_yi[length(gr_yi)])
  #gr_pye_EN <- c(gr_yi[-length(gr_yi)] - gr_phi_EN, gr_yi[length(gr_yi)])
  #gr_pye_SCAD <- c(gr_yi[-length(gr_yi)] - gr_phi_SCAD, gr_yi[length(gr_yi)])
  #gr_pye_MCP <- c(gr_yi[-length(gr_yi)] - gr_phi_MCP, gr_yi[length(gr_yi)])


  return(list("longpye_L12"=longpye_L12, "longpye_L1"=longpye_L1,
              "longpye_EN"=longpye_EN, "longpye_SCAD"=longpye_SCAD, "longpye_MCP"=longpye_MCP,
              #"gr_pye_L12"=gr_pye_L12, "gr_pye_L1"=gr_pye_L1, "gr_pye_EN"=gr_pye_EN, "gr_pye_SCAD"=gr_pye_SCAD, "gr_pye_MCP"=gr_pye_MCP,
              #"gr_phi_L12"=gr_phi_L12, "gr_phi_L1"=gr_phi_L1, "gr_phi_EN"=gr_phi_EN, "gr_phi_SCAD"=gr_phi_SCAD, "gr_phi_MCP"=gr_phi_MCP,
              "gr_yi_betas"=gr_yi_betas, "gr_yi_c"=gr_yi_c,
              "TP"=TP, "TN"=TN, "FP"=FP, "FN"=FN,
              "TP_i"=TP_i, "TN_i"=TN_i, "FP_i"=FP_i, "FN_i"=FN_i,
              "youden_index"=yi,
              "sensitivity"=sensitivity, "specificity"=spec, "geometric_mean"=gm,
              "fdr"=fdr, "mcc"=mcc, "auc"=auc, "corrclass"=corrclass, "corrclass_i"=corrclass_i,
              z_hat_test=z_matrix_test,
              f0_test=f0_test,
              f1_test=f1_test,
              post_prob1_test=post_prob1_test,
              post_prob0_test=post_prob0_test,
              ratio_test=ratio_test,
              df_results_by_id=df_pred,
              c_hat_discriminant=c_hat_discriminant,
              TP_test=TP_test,
              TN_test=TN_test,
              FP_test=FP_test,
              FN_test=FN_test,
              specificity_test=spec_test,
              fnr_test=fnr_test,
              sensitivity_test=sensitivity_test,
              geometric_mean_test=gm_test,
              fdr_test=fdr_test,
              mcc_test=mcc_test,
              auc_test=auc_test,
              corrclass_test=corrclass_test,
              youden_index_test=yi_test,
              #"y_hat"=df2[,c("ID_rows", id,"y_hat")],
              "input_data" = list("beta"=betas, "c_hat"=c, "lambda"=lambda, "alpha"=alpha,
                                  "a1"=a1, "a2"=a2, "kernel"=kernel)
              ))
}





#' @title longpye_KS_estimation_same_betas
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
#' c <- 0
#' prox_penalty <- get(paste0("proximal_operator_", penalty))
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
longpye_KS_estimation_same_betas <- function(df, id="id", X=names(df[,!(names(df) %in% c(y,id,t))]), y="y",
                                   t="t", lambda, penalty="L1", beta_start_input=NULL,
                                   beta_start_default=NULL, max.print=20,
                                   alpha=0.5, a1=3.7, a2=3, regressors_betas=NULL, fold=NULL, max_iter=10,
                                   trend = "monotone", delta = 1e-5, max_alpha=10000, stepsizeShrink=0.5,
                                   min_alpha=1e-7, convergence_error=1e-7,
                                   trace=1, seed=1, c_zero_fixed = FALSE, kernel = "gaussian"){

  start_time <- Sys.time()

  options(max.print = max.print)

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors_betas to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}
  if (class(t) != "character"){stop("t can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID_rows" %in% colnames(df)){stop("ID_rows already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID_rows <- rownames(df)
  df1 <- cbind(ID_rows, id=df[,id], t=df[,t], df[,(names(df) %in% c(y)), drop=FALSE], df[,(names(df) %in% c(X)), drop=FALSE])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (length(lambda)>1){stop("longpye_KS_estimation_same_betas stopped becuase lambda needs to have length 1")}
  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}

  betas1_initial_zeros <- rep(0, length(X))
  betas1_initial_corr <- NULL

  #define c
  c <- rep(0, length(unique(df1[,"t"])))
  names(c) <- order(unique(df1[,"t"]))

  #if not assigned, if we have a low number of covariates it is better to start betas with cov,
  #otherwise it is better to start with zeros

  if(length(X) <20){
    zeros_stay_zeros_from_iteration <- 1
    if(length(beta_start_default)== 0 ){
      beta_start_default <- "zeros" #"corr"
    }
  } else if(length(X) <1000){
    zeros_stay_zeros_from_iteration <- 1
    if(length(beta_start_default)== 0 ){
      beta_start_default <- "zeros" #"corr"
    }
  } else {
    zeros_stay_zeros_from_iteration <- 20
    if(length(beta_start_default)== 0 ){
      beta_start_default <- "zeros"
    }
  }

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

  betas_start <- betas1
  c_start <- c
  #add c to the betas
  betas1 <- c(betas1[!(names(betas1) %in% names(c))], c) #merge betas1 and c

  if (penalty == "SCAD") {a=a1} else {a=a2}
  prox_penalty = get(paste0("proximal_operator_", penalty)) #proximal oper. to be used

  #wrappers
  delta_fx <- function(x){
    execution <- longpye_KS_same_betas(df_train=df1[,names(df1)!="ID_rows", drop=FALSE], id="id", X=X[X %in% names(x)], t="t",
                            y=y, betas=x[!(names(x) %in% names(c))], lambda=lambda, c=x[(names(x) %in% names(c))],
                            kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty)
    #result <- c(-execution$gr_yi_betas, rep(0, length(unique(df1[,"t"]))))
    #result <- c(-execution$gr_yi_betas, 0, -execution$gr_yi_c[-1])
    result <- c(-execution$gr_yi_betas, -execution$gr_yi_c)
    #print(result)
    return(result)
  }

  proxx <- function(x, eta){
    result <- c(prox_penalty(betas=x[!(names(x) %in% names(c))], lambda=eta*lambda, alpha=alpha, a=a), x[(names(x) %in% names(c))])
    return(result)
  }


  Fx <- function(x){
    result <- -getElement(longpye_KS_same_betas(df_train=df1[,names(df1)!="ID_rows", drop=FALSE], id="id",  X=X[X %in% names(x)], t="t",
                                     y=y, betas=x[!(names(x) %in% names(c))], lambda=lambda, c=x[(names(x) %in% names(c))],
                                     kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty),  paste0("longpye_", penalty)) #"youden_index")
    return(result)
  }


  if (trend == "monotone"){
    estim <- mmAPG(x0=betas1, c_pos=which(names(betas1) %in% names(c)), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda, penalty=penalty,
                                      fold=fold, stepsizeShrink=stepsizeShrink, delta=delta, max_alpha=max_alpha,
                                      max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                      trace=trace, seed=seed, max.print=max.print,
                                      zeros_stay_zeros_from_iteration=zeros_stay_zeros_from_iteration)
  } else if (trend == "nonmonotone"){
    estim <- mnmAPG(x0=betas1, c_pos=which(names(betas1) %in% names(c)), delta_fx=delta_fx, proxx=proxx, Fx=Fx, lambda=lambda, penalty=penalty,
                                      fold=fold, stepsizeShrink=stepsizeShrink, delta=delta, max_alpha=max_alpha,
                                      max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                      trace=trace, seed=seed, max.print=max.print,
                                      zeros_stay_zeros_from_iteration=zeros_stay_zeros_from_iteration)
  }

  #divide betas1 and c
  betas_hat <- estim$x1[!(names(estim$x1) %in% names(c))]
  c_hat <- estim$x1[(names(estim$x1) %in% names(c))]

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
    c_hat <- c_start
  }

  #test
  #longpye_KS_value1 <- longpye_KS_same_betas(df_train=df1[,names(df1)!="ID_rows", drop=FALSE],
  #                               df_test=df1[,names(df1)!="ID_rows", drop=FALSE],
  #                               id="id", X=X, y=y, betas=betas_start, lambda=lambda, t="t", c=c_start,
  #                               kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty,
  #                               print.CDF.plot=TRUE)

  #compute z_hat
  longpye_KS_value <- longpye_KS_same_betas(df_train=df1[,names(df1)!="ID_rows", drop=FALSE],
                                 df_test=df1[,names(df1)!="ID_rows", drop=FALSE],
                                 id="id", X=X, y=y, betas=betas_hat, lambda=lambda, t="t", c=c_hat,
                                 kernel=kernel, alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                 print.CDF.plot=TRUE)

  #z_hat <- longpye_KS_value$z_hat
  #df_hat <- merge(df1, z_hat[, c("ID_rows", "z_i_hat")], by = 'ID_rows')
  #rownames(df_hat) <- df_hat$ID_rows
  df_hat_by_id <- longpye_KS_value$df_results_by_id

  youden_index <- longpye_KS_value$youden_index
  sensitivity <- longpye_KS_value$sensitivity_test
  specificity <- longpye_KS_value$specificity_test
  geometric_mean <- longpye_KS_value$geometric_mean_test
  fdr <- longpye_KS_value$fdr_test
  mcc <- longpye_KS_value$mcc_test
  auc <- longpye_KS_value$auc_test
  corrclass <- longpye_KS_value$corrclass_test
  c_hat_discriminant <- longpye_KS_value$c_hat_discriminant
  #corrclass_i <- longpye_KS_value$corrclass_i

  if (trace %in% c(1,2)){
    #print the estimation
    cat("Final estimation: \n")
    cat("-> algorithm: Accelerated Proximal Gradient for Nonconvex Programming (APG), the", trend, "version. ; \n")
    if (length(fold)!=0) {cat("fold =", fold, "; \n")}
    visualize_betas <- c(betas_hat[which(betas_hat!=0)], c_hat)
    cat("total iters:", estim$tot_iters, "; Backtraking iters:" , estim$backtrack_iters , "; lambda:", lambda, "; penalty:", penalty, "; longpye_KS_same_betas:", getElement(longpye_KS_value,  paste0("longpye_",penalty)), "; youden_index:", longpye_KS_value$youden_index, "; sensitivity:", longpye_KS_value$sensitivity_test, "; specificity:", longpye_KS_value$specificity_test, "; geometric_mean:", longpye_KS_value$geometric_mean_test, "; fdr:", longpye_KS_value$fdr_test, "; mcc:", longpye_KS_value$mcc_test, "; auc:", longpye_KS_value$auc_test, "; corrclass:", longpye_KS_value$corrclass_test, "; \n")
    cat("TP:", longpye_KS_value$TP_test, "; TN:", longpye_KS_value$TN_test, "; FP:", longpye_KS_value$FP_test, "; FN:", longpye_KS_value$FN_test, ";  betas: \n")
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
  name2 = paste0("longpye_KS_", penalty)

  results <-list("t1" = betas_hat,
                 "t2" = getElement(longpye_KS_value, paste0("longpye_", penalty)),
                 lambda=lambda,
                 penalty = penalty,
                 betas_start = betas_start,
                 kernel=kernel,
                 c_hat=c_hat,
                 c_hat_discriminant=c_hat_discriminant,
                 z_hat=longpye_KS_value$z_hat,
                 df_results_by_id=longpye_KS_value$df_results_by_id,
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
                 input_parameters=c(X=X, y=y, id=id, t=t, alpha=alpha)
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
setClass(Class="longpye_KS_CV_output",
         representation(
           penalty="character",
           kernel="character",
           longpye_KS_L12= "ANY",
           longpye_KS_L1= "ANY",
           longpye_KS_EN= "ANY",
           longpye_KS_SCAD= "ANY",
           longpye_KS_MCP= "ANY",
           auc= "ANY",
           youden_index="ANY",
           sensitivity="ANY",
           specificity="ANY",
           geometric_mean="ANY",
           fdr="ANY",
           mcc="ANY",
           corrclass="ANY",
           n_betas="ANY",
           betas="list"
         )
)

#to extract just part of the estim of pye
subset_longpye_KS_same_betas <- function(df_train, df_test, X, y, id, t, betas, lambda, c, fold, alpha, trace, a1, a2, penalty,
                           cv_time, niter, kernel, c_hat_discriminant){

  longpye_KS_result <- longpye_KS_same_betas(df_train=df_train[,names(df_train)!="ID_rows"], df_test=df_test[,names(df_test)!="ID_rows"],
                                  X=X, y=y, id=id, t=t, betas=betas, lambda=lambda, c=c,
                                  alpha=alpha, a1=a1, a2=a2, penalty=penalty, prediction=TRUE, kernel=kernel,
                                  c_hat_discriminant=c_hat_discriminant)

  #z_hat <- longpye_KS_result$z_hat$z_hat

  #optimal cutpoint (the one estimated in the training)
  youden_index <- longpye_KS_result$youden_index_test
  sensitivity <- longpye_KS_result$sensitivity_test
  specificity <- longpye_KS_result$specificity_test
  geometric_mean <- longpye_KS_result$geometric_mean_test
  fdr <- longpye_KS_result$fdr_test
  mcc <- longpye_KS_result$mcc_test
  auc <- longpye_KS_result$auc_test
  corrclass <- longpye_KS_result$corrclass_test
  TP <- longpye_KS_result$TP_test
  TN <- longpye_KS_result$TN_test
  FP <- longpye_KS_result$FP_test
  FN <- longpye_KS_result$FN_test
  c_hat <- c

  if (trace %in% c(1,2)) {
    cat("-> Results on the TEST SET of longPYE\n")
    cat("-> algorithm: longpye_KS_proximal_gradient_method ; ")
    if (length(fold)!=0) {cat("fold =", fold, "; \n")}
    visualize_betas <- c(betas[which(betas!=0)], c_hat)
    cat("lambda:", lambda, "; penalty:", penalty, "; longpye_KS_same_betas:", getElement(longpye_KS_result,  paste0("longpye_",penalty)), "; youden_index:", longpye_KS_result$youden_index_test, "; sensitivity:", longpye_KS_result$sensitivity_test, "; specificity:", longpye_KS_result$specificity_test, "; geometric_mean:", longpye_KS_result$geometric_mean_test, "; fdr:", longpye_KS_result$fdr_test, "; mcc:", longpye_KS_result$mcc_test, "; auc:", longpye_KS_result$auc_test, "; corrclass:", longpye_KS_result$corrclass_test, " \n")
    cat("TP:", TP, "; TN:", TN, "; FP:", FP, "; FN:", FN, ";  betas: \n")
    print(visualize_betas)
    cat("Cross-validation time:", cv_time, "; Number of iterations:", niter, "\n")
    cat("\n")
  }

  return(list(longpye_KS_L12=longpye_KS_result$longpye_L12, longpye_KS_L1=longpye_KS_result$longpye_L1,
              longpye_KS_EN=longpye_KS_result$longpye_EN,
              longpye_KS_SCAD= longpye_KS_result$longpye_SCAD, longpye_KS_MCP=longpye_KS_result$longpye_MCP,
              youden_index=youden_index, sensitivity=sensitivity,
              specificity=specificity, geometric_mean=geometric_mean,
              fdr=fdr, mcc=mcc, auc=auc, corrclass=corrclass,
              TP=TP, TN=TN, FP=FP, FN=FN,
              z_hat=longpye_KS_result$z_hat))
}

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterCall
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom stats setNames
longpye_KS_same_betas.cv <- function (df, X, y, C, id, t, lambda, trace=1, alpha, a1, a2, penalty, folds_i, k,
                        regressors_betas=NULL,
                        longpye_KS_L12, longpye_KS_L1, longpye_KS_EN, longpye_KS_SCAD, longpye_KS_MCP,
                        auc, youden_index, sensitivity, specificity, geometric_mean, fdr,
                        mcc, corrclass, n_betas, n_gammas, used_cores, trend, delta, max_alpha, kernel,
                        beta_start_input, beta_start_default, c_zero_fixed, max_iter, min_alpha, convergence_error,
                        stepsizeShrink){

  test_i <- which(folds_i == k)
  train_df <- df[-test_i, ]
  test_df <- df[test_i, ]

  if(trace %in% c(1,2)){
    cat(" --------------------------------------------------------\n |  starting with the", k ,"-th fold for the CV of lambda   |\n --------------------------------------------------------\n")
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
    cl <- parallel::makeCluster(used_cores, outfile="log_longpye_ks_models.txt", setup_strategy = "sequential")
    if (!("cluster" %in% class(cl)))
      stop("cl is not of class 'cl'; see ?makeCluster")

    cat("Start parallel computing for cross-validation, I am using:", length(cl),"cores \n")
    parallel::clusterExport(cl, c("train_df", "X", "y", "id", "t", "trace", "alpha", "penalty",
                                  "regressors_betas", "k", "trend", "delta",
                                  "max_alpha", "kernel", "beta_start_input", "beta_start_default",
                                  "c_zero_fixed", "a2", "a1",
                                  "max_iter", "min_alpha", "convergence_error", "stepsizeShrink"), envir = environment())
    parallel::clusterCall(cl, function() require(pye))
    fitted_models <- parallel::parLapply(cl, lambda, function(x) longpye_KS_estimation_same_betas(df=train_df[,!(names(train_df) %in% "ID_rows")],
                                                                                    X=X, y=y, lambda=x, id="id", t="t",
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
                                                                                    c_zero_fixed=c_zero_fixed
                                                                                    ))

    on.exit(parallel::stopCluster(cl))

  } else {

    fitted_models <- lapply(lambda, function(x) longpye_KS_estimation_same_betas(df=train_df[,!(names(train_df) %in% "ID_rows")],
                                                                   X=X, y=y, lambda=x, id="id", t="t",
                                                                   beta_start_input=beta_start_input,
                                                                   beta_start_default=beta_start_default,
                                                                   trace=trace,
                                                                   alpha=alpha,
                                                                   a1=a1, a2=a2, penalty=penalty, max_iter=max_iter,
                                                                   min_alpha=min_alpha,
                                                                   convergence_error=convergence_error,
                                                                   regressors_betas=regressors_betas, fold=k,
                                                                   trend=trend,
                                                                   stepsizeShrink=stepsizeShrink,
                                                                   delta=delta, max_alpha=max_alpha,
                                                                   kernel=kernel,
                                                                   c_zero_fixed=c_zero_fixed
                                                                   ))
  }

  z_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "z_hat"))

  #computing c
  #measures on the train set - if we don't compute c with covariates, we are only dependent of lambda:
  temp_pye = get(paste0("longpye_KS_", penalty))

  temp_pyetemp_pye <<- temp_pye
  fitted_modelsfitted_models <<- fitted_models

  for (i in 1:length(lambda)){
    #measures on the train set
    temp_pye$train[k,i] <- getElement(fitted_models[[i]], paste0("longpye_KS_", penalty))
    auc$train[k,i]  <- getElement(fitted_models[[i]], "auc")
    youden_index$train[k,i] <- getElement(fitted_models[[i]], "youden_index")
    sensitivity$train[k,i] <- getElement(fitted_models[[i]], "sensitivity")
    specificity$train[k,i]  <- getElement(fitted_models[[i]], "specificity")
    geometric_mean$train[k,i]  <- getElement(fitted_models[[i]], "geometric_mean")
    fdr$train[k,i]  <- getElement(fitted_models[[i]], "fdr")
    mcc$train[k,i]  <- getElement(fitted_models[[i]], "mcc")
    corrclass$train[k,i] <- getElement(fitted_models[[i]], "corrclass")
  }

  n_betas[k, ] <- unlist(lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "n_betas")))

  #create a list of the results (NB: betas_hat includes c_hat) - it is only based on lambda
  betas_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], paste0("betas_hat_", penalty)))
  names(betas_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  #c_hat is the fixed value of c in case we don't use the covariates to estimate a patient's specific cut-off value
  c_hat <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "c_hat"))
  names(c_hat)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  c_hat_discriminant <- lapply(c(1:length(lambda)), function(x) getElement(fitted_models[[x]], "c_hat_discriminant"))
  names(c_hat_discriminant)<- lapply(lambda, function (xx) paste("lambda", xx, sep="="))

  betas <- betas_hat
  #names(betas) <- lambda

  #measures on the test set
  all_measures_test <- mapply(function(x, xx, xxx, z, zz) subset_longpye_KS_same_betas(df_train=train_df, df_test=test_df,
                                                           X=X, y=y, betas=x,lambda=z, id=id, t="t",
                                                           c=xx, fold = k, alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                                           cv_time=fitted_models[[zz]]$estimation_time,
                                                           niter=fitted_models[[zz]]$niter, kernel=kernel, trace=trace,
                                                           c_hat_discriminant=xxx),
                              betas_hat, c_hat, c_hat_discriminant, lambda, 1:length(lambda))


  #measures on the test set - if we don't compute c with covariates, we are only dependent of lambda:
  for (i in 1:length(lambda)){
    #measures on the train set
    temp_pye$test[k,i] <- all_measures_test[paste0("longpye_KS_", penalty),][[i]]
    auc$test[k,i]  <- all_measures_test["auc",][[i]]
    youden_index$test[k,i] <- all_measures_test["youden_index",][[i]]
    sensitivity$test[k,i] <- all_measures_test["sensitivity",][[i]]
    specificity$test[k,i]  <- all_measures_test["specificity",][[i]]
    geometric_mean$test[k,i]  <- all_measures_test["geometric_mean",][[i]]
    fdr$test[k,i]  <- all_measures_test["fdr",][[i]]
    mcc$test[k,i]  <- all_measures_test["mcc",][[i]]
    corrclass$test[k,i] <- all_measures_test["corrclass",][[i]]
  }

  assign(paste0("longpye_KS_", penalty), temp_pye)


  return(new("longpye_KS_CV_output", penalty=penalty, kernel=kernel,
             longpye_KS_L12=longpye_KS_L12, longpye_KS_L1=longpye_KS_L1, longpye_KS_EN=longpye_KS_EN,
             longpye_KS_SCAD=longpye_KS_SCAD, longpye_KS_MCP=longpye_KS_MCP,
             auc=auc, youden_index=youden_index, sensitivity=sensitivity,
             specificity=specificity, geometric_mean=geometric_mean, fdr=fdr, mcc=mcc, corrclass=corrclass,
             n_betas=n_betas, betas=betas))
}


#' @title longpye_KS_compute_cv_same_betas
#'
#' @description function to perform the cross-validation to select the best
#' value of lambda for the estimation of betas and c (and
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
#' @param measure_to_select_lambda the measure used to select lambda in the
#' cross-validation process. Default is "ccr", i.e. correct classification rate
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
#' c <- 0
#' prox_penalty = get(paste0("proximal_operator_", penalty))
#' trend = "monotone" #or "nonmonotone"
#' pye_starting_point = "zeros"
#' alpha = 0.5
#' c_zero_fixed <- FALSE
#' used_cores <- 1
#' used_penalty_pye <- c("L1", "MCP") #c("L12", "L1", "EN", "SCAD", "MCP")
#' max_iter <- 10
#' n_folds <- 5
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
#'   #start cross-validation
#'   assign(name, pye_KS_compute_cv(penalty=p, df=df,
#'     X=names(df[,!(names(df) %in% c(y,C))]), y=y, C=C, trace=1,
#'     beta_start_default=pye_starting_point, trend = trend,
#'     beta_start_input=NULL, n_folds=n_folds, lambda=lambda,
#'     alpha=alpha, a1=3.7, a2=3,
#'     regressors_betas=regressors_betas, regressors_gammas=regressors_gammas,
#'     used_cores=used_cores, kernel="gaussian", c_zero_fixed=c_zero_fixed,
#'     penalty_g=p, max_iter=max_iter, max_iter_g=max_iter))
#'
#'   #take the best lambda per measures (in case of same measure for diff lambdas,
#'   #we take the one associated to less betas)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_yi"), get(name)$lambda_hat_pye_KS_yi)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_auc"), get(name)$lambda_hat_pye_KS_auc)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_ccr"), get(name)$lambda_hat_pye_KS_ccr)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_gm"), get(name)$lambda_hat_pye_KS_gm)
#'   assign(paste0("lambda_hat_PYE_KS_covYI_", p ,"_pye"), get(name)$lambda_hat_pye_KS_pye)
#' }
#'
#'
#' @export
longpye_KS_compute_cv_same_betas <- function (n_folds, df, X=names(df[,!(names(df) %in% c(y,C))]), y="y",
                                id="id", t="t", lambda, trace=1,
                                alpha=0.5,
                                a1=3.7, a2=3, penalty="L1", regressors_betas=NULL,
                                seed=1, used_cores=1, trend = "monotone", delta = 1e-5, max_alpha=10000,
                                kernel="gaussian", beta_start_input=NULL,
                                max_iter=20, #reduced for the CV
                                min_alpha=1e-9,
                                convergence_error=1e-7,
                                stepsizeShrink=0.5,
                                beta_start_default=NULL, scaling=FALSE, c_zero_fixed=FALSE,
                                measure_to_select_lambda="ccr") {

  start_time <- Sys.time()

  #check df:
  if (nrow(df)==0){stop("df has no rows")}

  #Create a new df considering only the columns included in X (the regressors_betas to consider) and y, the target variable
  #Create also the variable ID and add it at the beginning (useful for merges)
  if (class(X) != "character"){stop("X can only be of class character or data.frame.")}
  if (class(y) != "character"){stop("y can only be of class character or data.frame.")}
  if (class(t) != "character"){stop("t can only be of class character or data.frame.")}

  #check if ID already exists in the dataset
  if("ID_rows" %in% colnames(df)){stop("ID_rows already exists as column in df! Please delete or rename this column since I need this name to set the internal ID")}

  ID_rows <- rownames(df)
  df1 <- cbind(ID_rows, id=df[,id], t=df[,t], df[,(names(df) %in% c(y)), drop=FALSE], df[,(names(df) %in% c(X)), drop=FALSE])

  #starting controls
  #control that the target variable y is (0,1)
  if (!all(levels(factor(df1$y)) == c("0", "1"))){stop("The target variable y contains values different from 0 and 1, this is not allowed.")}

  if (!(penalty %in% c("L12","L1","EN","SCAD","MCP"))){stop("A wrong value has been assigned to the parameter penalty.")}
  if (!(trace %in% c(0,1,2))){stop("The parameter trace has been wrongly assigned. It can be 0 (no print),1 (partial print) or 2 (full print).")}
  if (!(trend %in% c("monotone", "nonmonotone"))){stop("The parameter trend has been wrongly assigned. It can be monotone or nonmonotone")}
  if (!(measure_to_select_lambda %in% c("auc", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye"))){stop("A wrong value has been assigned to the parameter measure_to_select_lambda")}

  set.seed(seed)

  #standardize df1
  if (scaling==TRUE){
    df1 <- scaling_df_for_pye (df=df1, X=colnames(df1[, names(df1) %in% c(X)]), y="y")$df_scaled
  }

  #check if df is well populated for the variable y: we need at least 2 element of 1 and 0 per fold
  if (((length(df1$y[df1$y==1])/length(base::unique(df1$t)))<2*n_folds)){stop("df contains too few 1s for this number of folds")
  } else if (((length(df1$y[df1$y==0])/length(base::unique(df1$t)))<2*n_folds)){stop("df contains too few 0s for this number of folds")}

  #divide the dataset in folds: to equalize the number of 0 and 1 in each sample I stratify
  df_sort <- base::unique(df1[order(getElement(df1,y)), c("id", y)])
  fold_i_0 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==0,])), replace=FALSE)
  fold_i_1 <- sample(rep(1:n_folds, length.out = nrow(df_sort[df_sort$y==1,])), replace=FALSE)
  df_sort <- cbind(df_sort, fold=c(fold_i_0,fold_i_1))
  folds_i <- merge(df1[,c("ID_rows","id")], df_sort[,c("id","fold")], by="id", all = FALSE, sort = FALSE)[,3]

  #names of the columns
  lambdanames <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
  foldnames <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
  #measures (train and test)
  list_of_measures <- c("longpye_KS_L12", "longpye_KS_L1", "longpye_KS_EN", "longpye_KS_SCAD", "longpye_KS_MCP",
                        "auc", "youden_index", "sensitivity", "specificity", "geometric_mean",
                        "fdr", "mcc", "corrclass")
  longpye_KS_L12 <- longpye_KS_L1 <- longpye_KS_EN <- longpye_KS_SCAD <- longpye_KS_MCP <- NULL
  auc <- youden_index <- sensitivity <- specificity <- geometric_mean <- fdr <- mcc <- corrclass <- NULL
  for (mes in list_of_measures){
    #auc <- list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,taunames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(c(1:n_folds),taunames)))
    assign(mes, list (train = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,lambdanames)), test = matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames,lambdanames))))
  }

  n_betas <- matrix(NA, nrow = n_folds, ncol = length(lambda), dimnames= list(foldnames, unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))))

  results = list(longpye_KS_L12=longpye_KS_L12, longpye_KS_L1=longpye_KS_L1, longpye_KS_EN=longpye_KS_EN,
                 longpye_KS_SCAD=longpye_KS_SCAD, longpye_KS_MCP=longpye_KS_MCP,
                 auc=auc, youden_index=youden_index,
                 sensitivity=sensitivity, specificity=specificity,
                 geometric_mean=geometric_mean, fdr=fdr, mcc=mcc,
                 corrclass=corrclass, n_betas=n_betas)

  cat("Starting CV with the following estimation/convergence parameters: \n")
  cat("max_iter:", max_iter, "\n")
  cat("min_alpha:", min_alpha, "\n")
  cat("convergence_error:", convergence_error, "\n")
  cat("stepsizeShrink:", stepsizeShrink, "\n")

  #fill the matrices
  results <- mapply(function(k) longpye_KS_same_betas.cv(df=df1[,names(df1)!="ID_rows"], X=X, y=y, id="id", t="t",
                                          lambda=lambda, trace=trace,
                                          alpha=alpha, a1=a1, a2=a2, penalty=penalty,
                                          folds_i=folds_i, k=k, regressors_betas, longpye_KS_L12=longpye_KS_L12,
                                          longpye_KS_L1=longpye_KS_L1, longpye_KS_EN=longpye_KS_EN,
                                          longpye_KS_SCAD=longpye_KS_SCAD, longpye_KS_MCP=longpye_KS_MCP,
                                          auc=auc, youden_index=youden_index,
                                          sensitivity=sensitivity, specificity=specificity,
                                          geometric_mean=geometric_mean, fdr=fdr,
                                          mcc=mcc, corrclass=corrclass,
                                          n_betas=n_betas,
                                          max_iter=max_iter, min_alpha=min_alpha, convergence_error=convergence_error,
                                          stepsizeShrink=stepsizeShrink,
                                          used_cores=used_cores, kernel=kernel, trend=trend,
                                          delta=delta, max_alpha=max_alpha, beta_start_input=beta_start_input,
                                          beta_start_default=beta_start_default, c_zero_fixed=c_zero_fixed
                                          ),
                                        seq(1:n_folds))

  wrapper <- function(results, mes, dataset){
    m <- matrix(apply(expand.grid(c(1:n_folds), 1:length(lambda)), 1, function(v) getElement(getElement(results[[v[1]]],mes),dataset)[v[1],v[2]]), n_folds, length(lambda))
    colnames(m) <- lambdanames
    rownames(m) <- foldnames
    return(m)
  }

  #prepare the results
  temp_pye <- list(train = wrapper(results, paste0("longpye_KS_", penalty), "train") , test = wrapper(results, paste0("longpye_KS_", penalty), "test"))

  #temp_pye[] <- list(train = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]], paste0("longpye_KS_", penalty))$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) getElement(results[[k]],paste0("longpye_KS_", penalty))$test[k,]))))
  assign(paste0("longpye_KS_", penalty) , temp_pye)

  for (mes in list_of_measures){
    #aauc[] <- lapply(lambda, function(i) list(train = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$train[k,]))), test = t(as.matrix(sapply(c(1:n_folds), function(k) resultsresults[[k]]@aauc[[i]]$test[k,])))))
    #names(aauc) <- unlist(lapply(lambda, function (x) paste("lambda", x, sep="=")))
    assign(mes, list(train = wrapper(results, mes, "train"), test = wrapper(results, mes, "test")))
  }

  n_betas[] <- t(as.matrix(sapply(c(1:n_folds), function(k) results[[k]]@n_betas[k,])))

  betas <- t(sapply(c(1:n_folds), function(k) results[[k]]@betas))
  rownames(betas) <- unlist(lapply(1:n_folds, function (x) paste("fold", x, sep="=")))
  #gammas[] <- t(as.list(sapply(c(1:n_folds), function(k) results[[k]]@gammas)))

  measures <- c("auc", "yi", "sen", "spc", "gm", "fdr", "mcc", "ccr", "pye")
  list_of_measures2 <- c("auc", "youden_index", "sensitivity", "specificity",
                        "geometric_mean", "fdr", "mcc", "corrclass", paste0("longpye_KS_", penalty))

  lambda_hat_pye_KS_yi <- lambda_hat_pye_KS_auc <- lambda_hat_pye_KS_ccr <- NULL
  lambda_hat_pye_KS_sen <- lambda_hat_pye_KS_spc <- lambda_hat_pye_KS_gm <- lambda_hat_pye_KS_pye <- NULL
  for (i in 1:length(measures)){

    measures_matrix <- as.matrix(sapply(1:length(lambda), function(ii) colMeans(get(list_of_measures2[i])$test[,ii, drop=FALSE])))

    if (!is.na(max(measures_matrix))){
      max_measures <- which(measures_matrix == max(measures_matrix), arr.ind = TRUE)[1,]
      assign(paste0("lambda_hat_pye_KS_", measures[i]), lambda[max_measures[1]])
    } else {
      assign(paste0("lambda_hat_pye_KS_", measures[i]), NA)
    }
  }

  cv_time = difftime(Sys.time(), start_time , units = "mins")

  if (trace %in% c(1,2)) {
	cat("----------------> END OF THE CROSS-VALIDATION OF THE PYE METHODS <----------------- \n")
    cat("----------------> For the whoole Cross-validation it took:", cv_time, "minutes \n")
  }

  return(list(penalty=penalty, kernel=kernel, cv_time=cv_time, pye_KS_L12=pye_KS_L12,
              longpye_KS_L1=longpye_KS_L1, longpye_KS_EN=longpye_KS_EN, longpye_KS_SCAD=longpye_KS_SCAD,
              longpye_KS_MCP=longpye_KS_MCP,
              auc=auc, youden_index=youden_index, sensitivity=sensitivity,
              specificity=specificity, geometric_mean=geometric_mean,
              fdr=fdr, mcc=mcc, corrclass=corrclass,
              lambda_hat_pye_KS_yi=lambda_hat_pye_KS_yi, lambda_hat_pye_KS_auc=lambda_hat_pye_KS_auc,
              lambda_hat_pye_KS_ccr=lambda_hat_pye_KS_ccr, lambda_hat_pye_KS_sen=lambda_hat_pye_KS_sen,
              lambda_hat_pye_KS_spc=lambda_hat_pye_KS_spc,lambda_hat_pye_KS_gm=lambda_hat_pye_KS_gm,
              lambda_hat_pye_KS_pye=lambda_hat_pye_KS_pye, measure_to_select_lambda=measure_to_select_lambda,
              lambda_star=lambda_star, n_betas=n_betas, betas=betas))
}

