# Dataset Creation
# Generate:
# High-dimensional
# 1000 datasets 50x10.000
# 1000 datasets 100x10.000
# 1000 datasets 200x10.000
# Low-dimensional
# 1000 datasets 50x20
# 1000 datasets 100x20
# 1000 datasets 200x20
#
# Corr non-regressors variables
# 0.2   0.5  0.8
#
# function to generate the correlated regressors following different distributions
#' @importFrom Matrix bdiag
#' @importFrom MASS mvrnorm
create_variables = function(rows, cols, cov, mu=rep(0,cols), seed=1){
  set.seed(seed)
  #Var-cov matrix
  #library(Matrix)

  #block-diagonal matrix
  little_sigma <- matrix(cov, nrow=cols/4, ncol=cols/4) + diag(cols/4)*(1-cov)
  Sigma <- as.matrix(Matrix::bdiag(little_sigma,little_sigma,little_sigma,little_sigma))
  #diagonal matrix
  # Sigma <- matrix(0, nrow=cols, ncol=cols)
  # Sigma[1:(cols/4),1:(cols/4)] <- matrix(cov, nrow=cols/4, ncol=cols/4) + diag(cols/4)*(1-cov)
  # Sigma[((cols/4)+1):(cols/2),((cols/4)+1):(cols/2)] <- matrix(cov, nrow=cols/4, ncol=cols/4) + diag(cols/4)*(1-cov)
  # Sigma[((cols/2)+1):((cols/4)*3),((cols/2)+1):((cols/4)*3)] <- matrix(cov, nrow=cols/4, ncol=cols/4)+ diag(cols/4)*(1-cov)
  # Sigma[(((cols/4)*3)+1):cols,(((cols/4)*3)+1):cols] <- matrix(cov, nrow=cols/4, ncol=cols/4)+ diag(cols/4)*(1-cov)

  #Creation
  normal <- MASS::mvrnorm(n=rows, mu=mu, Sigma=Sigma)
  #CDF
  pvars <- stats::pnorm(normal)

  #normal variables
  #normvars <- as.data.frame(stats::qnorm(pvars[,1:(cols/4)], mean=3, sd= 2))
  #colnames(normvars)<- paste("norm", 1:(cols/4), sep="")
  #keep only the normal variables
  normvars <- as.data.frame(stats::qnorm(pvars[,1:cols], mean=0, sd=1))
  colnames(normvars)<- paste("norm", 1:(cols), sep="")
  #apply normal copula
  #install.packages("copula", repos="http://R-Forge.R-project.org")
  #involved packages to put in Description file: gsl, pcaPP, pspline, copula
  # myCop <- normalCopula(param=c(rep(cov, cols)), dim = cols, dispstr = "un")
  # myMvd <- mvdc(copula=myCop, margins=c(rep("normal", cols)),
  #               paramMargins=list(list(shape=2, scale=1),
  #                                 list(shape1=2, shape2=2),
  #                                list(df=5)) )

  #student-t variables
  #stuvars <- as.data.frame(qt(pvars[,1:(cols/4)], df=1))
  #colnames(stuvars)<- paste("stt", 1:(cols/4), sep="")
  #chi-squared variables
  #chivars <- as.data.frame(qchisq(pvars[,1:(cols/4)], df=7))
  #colnames(chivars)<- paste("chisq", 1:(cols/4), sep="")
  #gamma variables
  #gammavars <- as.data.frame(qgamma(pvars[,1:(cols/4)], shape=2, rate=2))
  #colnames(gammavars)<- paste("gamma", 1:(cols/4), sep="")
  #bernoulli variables
  #bervars <- as.data.frame(apply(pvars[,((cols/4)+1):(cols/2)], c(1,2) , function(x) if(x>0.5){1} else {0}))
  #colnames(bervars)<- paste("ber", 1:(cols/4), sep="")
  #exponential variables
  #expvars <- as.data.frame(qexp(pvars[,((cols/2)+1):((cols/4)*3)], rate = 0.5))
  #colnames(expvars)<- paste("exp", 1:(cols/4), sep="")
  #poisson variables
  #poisvars <- as.data.frame(qpois(pvars[,(((cols/4)*3)+1):cols], 5))
  #colnames(poisvars)<- paste("pois", 1:(cols/4), sep="")

  #merge
  #df4dist <- cbind(normvars, bervars, expvars, poisvars)
  df4dist <- cbind(normvars)

  return(df4dist)
}

# function to generate the correlated regressors following different distributions
#' @importFrom Matrix bdiag
#' @importFrom MASS mvrnorm
create_covariates <- function(rows, cols_cov, cov, mu_cov=rep(0,cols_cov), seed=1){
  set.seed(seed)
  #Var-cov matrix
  #library(Matrix)

  #block-diagonal matrix
  little_sigma <- matrix(cov, nrow=cols_cov/4, ncol=cols_cov/4) + diag(cols_cov/4)*(1-cov)
  Sigma <- as.matrix(Matrix::bdiag(little_sigma,little_sigma,little_sigma,little_sigma))

  #Creation
  normal <- MASS::mvrnorm(n=rows, mu=mu_cov, Sigma=Sigma)
  #CDF
  pvars <- stats::pnorm(normal)

  #normal variables
  normvars1 <- as.data.frame(stats::qnorm(pvars[,1:(cols_cov/4)], mean=0, sd=1))
  colnames(normvars1)<- paste("norm_cov", 1:(cols_cov/4), sep="")
  normvars2 <- as.data.frame(stats::qnorm(pvars[,((cols_cov/2)+1):((cols_cov/4)*3)], mean=0, sd=1))
  colnames(normvars2)<- paste("norm_cov", ((cols_cov/2)+1):((cols_cov/4)*3), sep="")
  #apply normal copula
  #install.packages("copula", repos="http://R-Forge.R-project.org")
  #involved packages to put in Description file: gsl, pcaPP, pspline, copula
  # myCop <- normalCopula(param=c(rep(cov, cols_cov)), dim = cols_cov, dispstr = "un")
  # myMvd <- mvdc(copula=myCop, margins=c(rep("normal", cols_cov)),
  #               paramMargins=list(list(shape=2, scale=1),
  #                                 list(shape1=2, shape2=2),
  #                                list(df=5)) )

  #student-t variables
  #stuvars <- as.data.frame(qt(pvars[,1:(cols_cov/4)], df=1))
  #colnames(stuvars)<- paste("stt", 1:(cols_cov/4), sep="")
  #chi-squared variables
  #chivars <- as.data.frame(qchisq(pvars[,1:(cols_cov/4)], df=7))
  #colnames(chivars)<- paste("chisq", 1:(cols_cov/4), sep="")
  #gamma variables
  #gammavars <- as.data.frame(qgamma(pvars[,1:(cols_cov/4)], shape=2, rate=2))
  #colnames(gammavars)<- paste("gamma", 1:(cols_cov/4), sep="")
  #bernoulli variables
  bervars1 <- as.data.frame(apply(pvars[,((cols_cov/4)+1):(cols_cov/2)], c(1,2) , function(x) if(x>0.5){1} else {0}))
  colnames(bervars1)<- paste("ber_cov", ((cols_cov/4)+1):(cols_cov/2), sep="")
  bervars2 <- as.data.frame(apply(pvars[,(((cols_cov/4)*3)+1):cols_cov], c(1,2) , function(x) if(x>0.5){1} else {0}))
  colnames(bervars2)<- paste("ber_cov", (((cols_cov/4)*3)+1):cols_cov, sep="")
  #exponential variables
  #expvars <- as.data.frame(qexp(pvars[,((cols_cov/2)+1):((cols_cov/4)*3)], rate = 0.5))
  #colnames(expvars)<- paste("exp", 1:(cols_cov/4), sep="")
  #poisson variables
  #poisvars <- as.data.frame(qpois(pvars[,(((cols_cov/4)*3)+1):cols_cov], 5))
  #colnames(poisvars)<- paste("pois", 1:(cols_cov/4), sep="")

  #merge
  #df4dist <- cbind(normvars, bervars, expvars, poisvars)
  df4dist <- cbind(normvars1, bervars1, normvars2, bervars2)

  return(df4dist)
}


#' @title create_sample
#'
#' @description function to create a synthetic sample with a target variable.
#' all the created variables follow a normal distribution
#'
#' @param rows_train number of rows of the training sample. Default is 50,
#' for an high-dimensional setting
#' @param cols number of variables of both the training and the test samples.
#' Default is 2000, for an high-dimensional setting
#' @param cov covariance in the covariace matrix of the normal distribution.
#' Increasing it, the created features are more correlated. Default is 0.5
#' @param mu mean of the (multivariate) normal distribution. Default is 0
#' @param rows_test number of rows of the test sample. Default is 1000.
#' We suggest to create a test sample much bigger than the training sample
#' to test your method/model on more data
#' @param seed fix the seed. Default is 1
#' @param varsN select 2 variables of the first cols/4 variables to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsB select 2 variables of the second cols/4 variables to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsE select 2 variables of the third cols/4 variables to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsP select 2 variables of the forth cols/4 variables to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#'
#' @return a list of elements containing the created dataframe in
#' different versions (scaled and not) and additional information
#'
#' @examples
#' library(pye)
#' df <- create_sample(cols = 20, rows_test = 10,
#' varsN = c(2, 4), varsB = c(6, 8), varsE = c(12, 14), varsP = c(16, 18))
#'
#' @importFrom stats quantile
#' @export
create_sample <- function (rows_train=50, cols=2000, cov=0.5, mu=rep(0, cols),
                           rows_test=1000, seed=1,
                           varsN=NULL, varsB=NULL, varsE=NULL, varsP=NULL){


  #set the seed
  set.seed(seed)
  if(length(varsN)==0){ varsN <- sort(sample(1:(cols/4), 2, replace=FALSE))} #var v167 (norm167) and v324 (norm324)
  if(length(varsB)==0){ varsB <- sort(sample((cols/4+1):(cols/2), 2, replace=FALSE))} #var v629 (ber129) and v918 (ber418)
  if(length(varsE)==0){ varsE <- sort(sample((cols/2+1):(cols/4*3), 2, replace=FALSE))} #var v1299 (exp299) and v1471 (exp471)
  if(length(varsP)==0){ varsP <- sort(sample((cols/4*3+1): cols, 2, replace=FALSE))} #var v1770 (pois270) and v1966 (pois466)

  rows <- rows_train + rows_test

  df <- create_variables(rows = rows, cols = cols, cov = cov, mu = mu, seed = seed)

  set.seed(seed)

  nregressors <- rep(0, cols)
  for (n in c(varsN, varsB, varsE, varsP)) {
    nregressors[n] = n
  }

  regressors <- names(df[, which(nregressors != 0)])

  betas <- c(-16, -4, 12, 8, -8, -12, 4, 16)
  coefficients <- rep(0, cols)

  for (i in 1:length(betas)) {
    coefficients[c(which(nregressors != 0))][i] <- betas[i]
  }

  df["z"] <- as.matrix(df[, nregressors]) %*% coefficients[nregressors]
  cutoff <- stats::quantile(df$z, 0.7)
  df["y"] <- ifelse(df$z > cutoff, 1, 0)
  z <- df$z

  df <- df[, c(ncol(df), 1:(ncol(df) - 2))]
  df_scaled <- scaling_df_for_pye(df=df, X=colnames(df[-1]), y="y")

  split <- sample(rep(1:rows), size = (rows - rows_test), replace = FALSE)

  train_df_scaled <- df_scaled$df_scaled[split, ]
  test_df_scaled <- df_scaled$df_scaled[-split, ]

  linearformula = "z <- as.matrix(df[,nregressors])%*%coefficients"

  return(list(df = df, df_scaled = df_scaled, train_df_scaled = train_df_scaled,
              test_df_scaled = test_df_scaled, nregressors = nregressors, regressors = regressors,
              coefficients = coefficients, linearformula = linearformula, z = z, cutoff = cutoff))
}


#' @title create_sample_with_covariates
#'
#' @description function to create a synthetic sample with a target variable.
#' This function creates 2 kind of variables: regressors and covariates,
#' setting useful to test, among others, pye, covYI and cPYE. The created
#' regressors follow a normal distribution while some of created covariates
#' follow a normal distribution and some others a bernoulli distribution.
#'
#' @param rows_train number of rows of the training sample. Default is 50,
#' for an high-dimensional setting
#' @param cols number of regressor variables of both the training and the
#' test samples. Default is 2000, for an high-dimensional setting
#' @param cols_cov number of covariate variables of both the training and the
#' test samples. Default is 20, increase it for an high-dimensional setting
#' @param covar covariance in the covariace matrix of the normal distribution
#' of the biomarkers. Increasing it, the created features are more correlated.
#' Default is 0.3
#' @param covar_cov covariance in the covariace matrix of the normal distribution
#' of the covariates. Increasing it, the created covariates are more correlated.
#' Default is 0.5
#' @param mu mean of the (multivariate) normal distribution of the regressors.
#' Default is 0
#' @param mu_cov mean of the (multivariate) normal distribution of the
#' covariates. Default is 0
#' @param rows_test number of rows of the test sample. Default is 1000.
#' We suggest to create a test sample much bigger than the training sample
#' to test your method/model on more data
#' @param seed fix the seed. Default is 1
#' @param varsN select 2 variables of the first cols/4 regressors to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsB select 2 variables of the second cols/4 regressors to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsE select 2 variables of the third cols/4 regressors to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varsP select 2 variables of the forth cols/4 regressors to be
#' real regressors, combined with the others for the creation of the
#' synthetic target variable. Default is NULL
#' @param varN_cov select 1 variable of the first cols_cov/4 covariates
#' to be real covariates, combined with the others for the creation of
#' the synthetic target cut-off. Default is NULL
#' @param varB_cov select 1 variable of the second cols_cov/4 covariates
#' to be real covariates, combined with the others for the creation of
#' the synthetic target cut-off. Default is NULL
#' @param varE_cov select 1 variable of the third cols_cov/4 covariates
#' to be real covariates, combined with the others for the creation of
#' the synthetic target cut-off. Default is NULL
#' @param varP_cov select 1 variable of the forth cols_cov/4 covariates
#' to be real covariates, combined with the others for the creation of
#' the synthetic target cut-off. Default is NULL
#'
#' @return a list of elements containing the created dataframe in
#' different versions (scaled and not) and additional information
#'
#' @examples
#' library(pye)
#' df <- create_sample_with_covariates(cols = 20, rows_test = 10,
#'   varsN = c(2, 4), varsB = c(6, 8), varsE = c(12, 14), varsP = c(16, 18))
#'
#' @export
create_sample_with_covariates <- function(rows_train=50, cols=2000, cols_cov=20, covar=0.3, covar_cov=0.5, mu=rep(0,cols),
                                  mu_cov=rep(0,cols_cov), rows_test=1000, seed=1,
                                  varsN=NULL, varsB=NULL, varsE=NULL, varsP=NULL,
                                  varN_cov=NULL, varB_cov=NULL, varE_cov=NULL, varP_cov=NULL){

  #set the seed
  set.seed(seed)
  varsN <- c(1,2)
  varsB <- c((cols/4+1),(cols/4+2))
  varsE <- c((cols/2+1),(cols/2+2))
  varsP <- c((cols/4*3+1),(cols/4*3+2))
  #if(length(varsN)==0){ varsN <- sort(sample(1:(cols/4), 2, replace=FALSE))} #var v167 (norm167) and v324 (norm324)
  #if(length(varsB)==0){ varsB <- sort(sample((cols/4+1):(cols/2), 2, replace=FALSE))} #var v629 (ber129) and v918 (ber418)
  #if(length(varsE)==0){ varsE <- sort(sample((cols/2+1):(cols/4*3), 2, replace=FALSE))} #var v1299 (exp299) and v1471 (exp471)
  #if(length(varsP)==0){ varsP <- sort(sample((cols/4*3+1): cols, 2, replace=FALSE))} #var v1770 (pois270) and v1966 (pois466)

  varsN_cov <- 1
  varsB_cov <- (cols_cov/4+1)
  varsE_cov <- (cols_cov/2+1)
  varsP_cov <- (cols_cov/4*3+1)
  #if(length(varsN_cov)==0){ varN_cov <- sort(sample(1:(cols_cov/4), 1, replace=FALSE))}
  #if(length(varsB_cov)==0){ varB_cov <- sort(sample((cols_cov/4+1):(cols_cov/2), 1, replace=FALSE))}
  #if(length(varsE_cov)==0){ varE_cov <- sort(sample((cols_cov/2+1):(cols_cov/4*3), 1, replace=FALSE))}
  #if(length(varsP_cov)==0){ varP_cov <- sort(sample((cols_cov/4*3+1): cols_cov, 1, replace=FALSE))}

  #increase the sample of other 1000 var (to test the model)
  rows<-rows_train+rows_test

  cv <- create_covariates(rows=rows, cols_cov=cols_cov, cov=covar_cov, mu_cov=mu_cov, seed=seed)
  addOn_for_mu <- (rowMeans(cv + sin(cv)) - mean(rowMeans(cv + sin(cv))))/ stats::sd(rowMeans(cv + sin(cv)))
  addOn_for_cov <- rowMeans(cv + sin(cv))

  df <- 0.6*addOn_for_mu + create_variables(rows=rows, cols=cols, cov=covar, mu=mu, seed=seed+1)

  #X and C
  X <- C <- NULL
  X <- names(df)
  C <- names(cv)

  #create the response variable as
  #select 10 random number between 1 and 2000
  set.seed(seed)

  nregressors <- rep(0,cols)
  for (n in c(varsN, varsB, varsE, varsP)){
    nregressors[n]=n
  }
  regressors <- names(df[, which(nregressors!= 0)])

  ncovariates <- rep(0,cols_cov)
  for (n in c(varsN_cov, varsB_cov, varsE_cov, varsP_cov)){
    ncovariates[n]=n
  }
  covariates <- names(cv[, which(ncovariates!= 0)])

  #intercept
  #b0=2

  #coefficients regressors
  betas <- c(1,-2,3,-4,4,-3,2,-1) #chosen ad hoc
  coefficients <- rep(0,cols)
  for (i in 1:length(betas)){
    coefficients[c(which(nregressors!= 0))][i] <- betas[i]
  }

  #coefficients covariates
  gammas_cv <- c(1, -2, 2, -1) #chosen ad hoc
  coefficients_cv <- rep(0,cols_cov)
  for (i in 1:length(gammas_cv)){
    coefficients_cv[c(which(ncovariates!= 0))][i] <- gammas_cv[i]
  }

  #linear formula
  df["z"] <- as.matrix(df[,nregressors])%*%coefficients[nregressors]
  #define c as the 75th percentile of z
  cutoff_cv <- as.matrix(cv[,ncovariates])%*%coefficients_cv[ncovariates]
  #create y
  df['y'] <- ifelse(df$z >= cutoff_cv, 1,0)
  z <- df$z
  target <- "y"

  #y
  y <- target

  #check correlation
  #corr.xy <- data.frame(matrix(ncol = cols+1, nrow = 0))#-4 since the const here is included. I want it always be the fist var to be computed!
  #names(corr.xy) <- names(train_df[,1:(cols+1)])
  #corr.xy[1:3,] <- t(sapply(c("pearson", "kendall", "spearman"), function (h) unlist(parallel::parLapply(cl, seq(1, cols+1, 1), function (x) stats::cor(train_df[,'y'], train_df[,x],  method = h)))))
  #row.names(corr.xy) <- c("pearson", "kendall", "spearman")
  #mean.corr.xy <- as.data.frame(t(colMeans(corr.xy)))
  #mean.corr.xy <- mean.corr.xy[,order(abs(mean.corr.xy),decreasing = TRUE)]

  #cor(df[,'z'], df[,1:(ncol(df)-2)])

  #Create an ID
  #df['ID']=as.numeric(rownames(df))
  #reorder
  df1 <- cbind(df[,c(ncol(df),1:(ncol(df)-2))], cv)
  df_scaled <- scaling_df_for_pye (df=df1, X=colnames(df1[-1]), y="y")

  #divide in training-set and test-set (rows-1000)
  split <- sample(rep(1:rows), size=(rows-rows_test), replace=FALSE)
  sum_unos <- sum(df_scaled$df_scaled[split,]$y)
  sum_zeros <- length(df_scaled[df_scaled$df_scaled[split,]$y==0])
  #lets put at least 20 1s and 20 0s (to ensure a good n-fold cross validation):
  while ((sum_unos < 20) | (sum_zeros < 20)) {
    split <- sample(rep(1:rows), size=(rows-rows_test), replace=FALSE)
    sum_unos <- sum(df_scaled$df_scaled[split,]$y)
    sum_zeros <- length(df_scaled[df_scaled$df_scaled[split,]$y==0])
  }

  train_df_scaled <- df_scaled$df_scaled[split,]
  test_df_scaled <- df_scaled$df_scaled[-split,]

  #preparing output
  linearformula_z = "z <- as.matrix(df[,nregressors])%*%coefficients"
  linearformula_cutoff = "as.matrix(cv[,ncovariates])%*%coefficients_cv[ncovariates]"

  return(list(df=df1, df_scaled=df_scaled, train_df_scaled=train_df_scaled, test_df_scaled=test_df_scaled,
              nregressors=nregressors, regressors=regressors, coefficients_regressors=coefficients,
              ncovariates=ncovariates, covariates=covariates, coefficients_covariates=coefficients_cv,
              target=target, linearformula_z=linearformula_z, linearformula_cutoff=linearformula_cutoff,
              X=X, y=y, C=C, z=z, cutoff=cutoff_cv))
}


#' @title PBC_data_load
#'
#' @description function to load the Mayo Clinic Primary Biliary Cholangitis
#' Data (PBC) ready to be used in longpye.
#'
#' @param T latest time to consider, i.e. from 1 to T. 4 is the default because
#' we already excluded the subjects with less than 4 visits.
#'
#' @return the Mayo Clinic Primary Biliary Cholangitis Data (PBC) ready to be
#' used in longpye
#'
#' @examples
#' library(pye)
#' df <- PBC_Mayo_Clinic_data_for_longpye()
#'
#' @export
PBC_Mayo_Clinic_data_for_longpye <- function(T=4){

  require("survival")
  require("pye")
  data(pbc, package="survival")
  df_PBC <- pbcseq #this is the dataset to be used, not pbc
  #df_PBC2 <- pbc #we could use pbc just tu add copper and trig variables at the starting time point
  df_PBC_ordered <- df_PBC[order(df_PBC$id, df_PBC$futime, df_PBC$day),]
  #unique(df_PBC_ordered$id)
  #19 variables, 312 subjects (IDs) with 1.945 clinical visits recorded

  #DETAILS:
  #Primary sclerosing cholangitis is an autoimmune disease leading to destruction
  #of the small bile ducts in the liver. Progression is slow but inexhortable,
  #eventually leading to cirrhosis and liver decompensation. The condition has been
  #recognized since at least 1851 and was named "primary biliary cirrhosis" in 1949.
  #Because cirrhosis is a feature only of advanced disease, a change of its name
  #to "primary biliary cholangitis" was proposed by patient advocacy groups in 2014.

  #This data is from the Mayo Clinic trial in PBC conducted between 1974 and 1984.
  #A total of 424 PBC patients, referred to Mayo Clinic during that ten-year
  #interval, met eligibility criteria for the randomized placebo controlled trial
  #of the drug D-penicillamine. The first 312 cases in the data set participated
  #in the randomized trial and contain largely complete data. The additional 112
  #cases did not participate in the clinical trial, but consented to have basic
  #measurements recorded and to be followed for survival. Six of those cases were
  #lost to follow-up shortly after diagnosis, so the data here are on an
  #additional 106 cases as well as the 312 randomized participants.


  #VARIABLE DESCRIPTION and DATA PREPARATION:
  #1) id: case number
  #2) futime: number of days between registration and the earlier of death, transplantion, or study analysis in July, 1986
  #3) status: status at endpoint - 0=alive, 1=transplanted, 2=dead
  #4) trt: drug - 1= D-penicillamine, 0=placebo
  #5) age in years, at registration
  #6) sex: m/f
  #df_PBC_ordered$sex
  sex <- as.data.frame(list("sex"=ifelse(df_PBC_ordered$sex=="m",1,0)))
  #7) day: number of days between enrollment and this visit date, remaining values on the line of data refer to this visit date.
  #8) ascites: presence of ascites: 0=no 1=yes and NAs
  #df_PBC_ordered$ascites
  ascites <- as.data.frame(list("ascites"=factor(df_PBC_ordered$ascites, exclude = NULL)))
  ascites <- model.matrix(~.-1, data = ascites, contrasts.arg = list(ascites = contrasts(ascites$ascites, contrasts=FALSE)))
  #9) hepato: presence of hepatomegaly or enlarged liver 0=no 1=yes
  #df_PBC_ordered$hepato
  hepato <- as.data.frame(list("hepato"=factor(df_PBC_ordered$hepato, exclude = NULL)))
  hepato <- model.matrix(~.-1, data = hepato, contrasts.arg = list(hepato = contrasts(hepato$hepato, contrasts=FALSE)))
  #10) spiders: presence of spiders 0=no 1=yes (blood vessel malformations in the skin)
  #df_PBC_ordered$spiders
  spiders <- as.data.frame(list("spiders"=factor(df_PBC_ordered$spiders, exclude = NULL)))
  spiders <- model.matrix(~.-1, data = spiders, contrasts.arg = list(spiders = contrasts(spiders$spiders, contrasts=FALSE)))
  #11) edema: presence of edema 0=no edema and no diuretic therapy for edema; 0.5 = edema present without diuretics (untreated), or edema resolved by diuretics (successfully treated); 1 = edema despite diuretic therapy
  #df_PBC_ordered$edema
  edema <- as.data.frame(list("edema"=factor(df_PBC_ordered$edema, exclude = NULL)))
  edema <- model.matrix(~.-1, data = edema, contrasts.arg = list(edema = contrasts(edema$edema, contrasts=FALSE)))
  #12) bili: serum bilirubin in mg/dl
  #13) chol: serum cholesterol in mg/dl
  #df_PBC_ordered$chol
  cholNA <- as.data.frame(list("cholNA"=ifelse(is.na(df_PBC_ordered$chol),1,0))) #create a flag variable when chol is NA
  chol <- df_PBC_ordered["chol"]
  chol[is.na(chol)] <- mean(chol[!is.na(chol)]) #input the mean value to chol
  #14) albumin: serum albumin in gm/dl
  #15) alk.phos: alkaline phosphatase in U/liter
  #df_PBC_ordered$alk.phos
  alk.phosNA <- as.data.frame(list("alk.phosNA"=ifelse(is.na(df_PBC_ordered$alk.phos),1,0))) #create a flag variable when alk.phos is NA
  alk.phos <- df_PBC_ordered["alk.phos"]
  alk.phos[is.na(alk.phos)] <- mean(chol[!is.na(alk.phos)]) #input the mean value to alk.phos
  #16) ast: aspartate aminotransferase, once called SGOT in U/ml (serum glutamic-oxaloacetic transaminase, the enzyme name has subsequently changed to "ALT" in the medical literature)
  #17) platelet: platelets per cubic ml / 1000
  #df_PBC_ordered$platelet
  plateletNA <- as.data.frame(list("plateletNA"=ifelse(is.na(df_PBC_ordered$platelet),1,0))) #create a flag variable when platelet is NA
  platelet <- df_PBC_ordered["platelet"]
  platelet[is.na(platelet)] <- mean(platelet[!is.na(platelet)]) #input the mean value to platelet
  #18) protime: prothrombin time in seconds (standardised blood clotting time)
  #19) stage: histologic stage of disease (needs biopsy)
  #df_PBC_ordered$stage
  stage <- as.data.frame(list("stage"=factor(df_PBC_ordered$stage, exclude = NULL)))
  stage <- model.matrix(~.-1, data = stage, contrasts.arg = list(stage = contrasts(stage$stage, contrasts=FALSE)))

  df_PBC_ordered_new <- cbind(df_PBC_ordered["id"], df_PBC_ordered["futime"], df_PBC_ordered["status"],
                              df_PBC_ordered["trt"], df_PBC_ordered["age"], sex, df_PBC_ordered["day"],
                              ascites, hepato, spiders, edema, df_PBC_ordered["bili"],
                              chol, cholNA, df_PBC_ordered["albumin"], alk.phos, alk.phosNA,
                              df_PBC_ordered["ast"], platelet, plateletNA, df_PBC_ordered["protime"],
                              stage)

  #check
  #sum(!complete.cases(df_PBC_ordered_new))
  #no more NAs

  #TARGET VARIABLE "status": status at endpoint, 0/1/2 for censored, transplant, dead
  #table(df_PBC_ordered_new$status)
  #  0    1    2
  #1073  147  725

  #let's delete the transplanted patients (they can bias the outcome)
  df_PBC_ordered2 <- df_PBC_ordered_new[!(df_PBC_ordered_new$status == 1),]
  #unique(df_PBC_ordered2$id)
  #283 subjects (IDs) with 1798 clinical visits recorded

  #TARGET VARIABLE "status": status at endpoint, 0/1/2 for censored, transplant, dead
  #table(df_PBC_ordered2$status)
  #  0    2
  #1073  725

  #We will consider the bivariate variable dead (1) / alive (0) at the endpoint (or censoring time) as target variable
  df_PBC_ordered2["y"] <- ifelse(df_PBC_ordered2$status < 2, 0, 1)

  #let's consider only the ones with no left censored
  #unique(df_PBC_ordered2[df_PBC_ordered2$day==0,][,"id"])
  #all the subjects had the starting visit at t=0

  #Let's consider only subjects with at least 4 visits
  df_PBC_ordered2 <- transform(df_PBC_ordered2, count=ave(id, id, FUN=seq_along))
  t1 <- df_PBC_ordered2[order(df_PBC_ordered2$id, -df_PBC_ordered2$count), ]
  t1 <- t1[!duplicated(t1[,c('id')]),][,]
  filtered_ids <- t1[t1$count > 3,]$id
  df_PBC_filtered <- df_PBC_ordered2[df_PBC_ordered2$id %in% filtered_ids,]
  #nrow(df_PBC_filtered)
  #length(unique(df_PBC_filtered$id))
  #207 subjects remained with 1645 visits
  #table(df_PBC_filtered$y)
  #  0    1
  #1011  634

  #distribution of count and y:
  t1 <- df_PBC_filtered[order(df_PBC_filtered$id, -df_PBC_filtered$count), ]
  t1 <- t1[!duplicated(t1[,c('id')]),][,]
  cbind(table(subset(t1, select = c('count', 'y'))), tot = table(t1$count))

  #times 14,15 and 16 has zero diseased. I have to exclude these times:
  t1 <- df_PBC_filtered[order(df_PBC_filtered$id, -df_PBC_filtered$count), ]
  t1 <- t1[!duplicated(t1[,c('id')]),][,]
  filtered_ids2 <- t1[t1$count < 14,]$id
  df_PBC_filtered2 <- df_PBC_filtered[df_PBC_filtered$id %in% filtered_ids2,]
  #nrow(df_PBC_filtered2)
  #length(unique(df_PBC_filtered2$id))
  #193 subjects remained with 1437 visits
  table(df_PBC_filtered2$y)
  #  0    1
  # 803  634

  #distribution of count and y:
  t1 <- df_PBC_filtered2[order(df_PBC_filtered2$id, -df_PBC_filtered2$count), ]
  t1 <- t1[!duplicated(t1[,c('id')]),][,]
  cbind(table(subset(t1, select = c('count', 'y'))), tot = table(t1$count))

  #scaled version of the dataset
  df_PBC_scaled_all <- scaling_df_for_pye (df=df_PBC_filtered2, X=colnames(df_PBC_filtered2[,!colnames(df_PBC_filtered2) %in% c("id", "futime", "status", "y", "count")]), y="y")

  #train-test split (70-30)
  split_df_PBC_all <- sample(rep(1:nrow(df_PBC_filtered2)), size=(round(nrow(df_PBC_filtered2)*0.7,0)), replace=FALSE)
  #standardize df
  train_df_PBC_scaled <- df_PBC_scaled_all$df_scaled[split_df_PBC_all,]
  test_df_PBC_scaled <- df_PBC_scaled_all$df_scaled[-split_df_PBC_all,]

  #TEST: consider only the full dataframe, i.e. only the t=1-4 periods
  #df <- df_PBC_scaled_all$df_scaled #<- the full dataset
  #df <- df_PBC_scaled_all$df_scaled[df_PBC_scaled_all$df_scaled$count %in% c(1,2,3,4),] #<- we consider only t=4 (193 complete subjects)
  subjects_to_select <- df_PBC_scaled_all$df_scaled$id %in% df_PBC_scaled_all$df_scaled$id[df_PBC_scaled_all$df_scaled$count == T]
  df <- df_PBC_scaled_all$df_scaled[subjects_to_select & df_PBC_scaled_all$df_scaled$count %in% c(1,2,3,4,5,6),] #<- we consider only t=6 (129 complete subjects)
  X <- colnames(df[,!colnames(df) %in% c("id", "futime", "status", "y", "count")])
  y <- "y"
  t <- "count"
  id <- "id"

  return(list(df=df, X=X, y=y, t=t, id=id,
              train_df_PBC_scaled=train_df_PBC_scaled, test_df_PBC_scaled=test_df_PBC_scaled))
}
