# ------------------------------------------------------------------------------
# Variational Bayes Patient Phenotyping from Electronic Health Records
#
# References:
#
# Hubbard, Rebecca A., et al. "A Bayesian latent class approach for EHR‐based
# phenotyping." Statistics in Medicine 38.1 (2019): 74-87.
#
# Brian Buckley December 2024
# ------------------------------------------------------------------------------

VBphenoR.env <- new.env(parent = emptyenv())

#' Run the Variational Bayes patient Phenotyping model
#'
#' @param ehr_data The EHR data.
#' @param biomarkers The EHR variables that are biomarkers. This is a vector of data column names corresponding to the biomarker variables.
#' @param gmm_X The input design matrix. Note the intercept column vector is assumed included.
#' @param gmm_k The maximum guess for how many disease components exist.  Minimum is 2 (binary response).
#' @param gmm_delta Change in ELBO that triggers algorithm stopping.
#' @param gmm_init Initialize the clusters c("random", "kmeans", proportion R{ > 0 and < 1}).
#' @param gmm_initParams Parameters for an initialiser requiring its own parameters e.g. dbscan.
#' @param gmm_maxiters The maximum iterations for VB GMM.
#' @param gmm_prior An informative prior for the GMM
#' @param gmm_stopIfELBOReverse Stop the VB iterations if the ELBO reverses direction.
#' @param gmm_verbose Print out information per iteration to track progress in case of long-running experiments.
#' @param logit_X The input design matrix. Note the intercept column vector is assumed included.
#' @param logit_prior An informative prior for the logit
#' @param logit_tol Change in ELBO that triggers algorithm stopping.
#' @param logit_maxiter The maximum iterations for VB logit.
#' @param logit_verbose Print out information per iteration to track progress in case of long-running experiments.
#'
#' @return A list containing:
#' * prevalence - The mean probability of latent phenotype given the data and priors.
#' * biomarker_shift - A data frame containing the biomarker shifts from normal for the phenotype.
#' * gmm - The VB GMM results. For details see help(vb_gmm_cavi).
#' * logit - The VB Logit results.  For details see help(logit_CAVI).
#' @export
#'
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @import knitr
#'
runModel <- function(ehr_data, biomarkers,
                     gmm_X, gmm_k, gmm_delta=1e-6,
                            gmm_init="kmeans", gmm_initParams=NULL,
                            gmm_maxiters=200, gmm_prior=NULL,
                            gmm_stopIfELBOReverse=FALSE,
                            gmm_verbose=FALSE,
                     logit_X, logit_prior,
                            logit_tol=1e-16, logit_maxiter=10000,
                            logit_verbose=FALSE) {

  gmm_result <- vb_gmm_cavi(X=gmm_X, k=gmm_k, delta=gmm_delta,
                            init=gmm_init, initParams=gmm_initParams,
                            maxiters = gmm_maxiters,
                            prior=gmm_prior,
                            stopIfELBOReverse=gmm_stopIfELBOReverse,
                            verbose=gmm_verbose)

  # Set 1,2 to 0,1 where 0 is the main class
  z <- gmm_result$z_post
  ztab <- as.data.frame(table(z))
  ztab$z <- as.character(ztab$z)
  ztab <- ztab[order(ztab$Freq),]
  if(ztab$z[1] == "2") {
    z[z==1] <- 0
    z[z==2] <- 1
  } else {
    z[z==2] <- 0
    z[z==1] <- 1
  }

  # Now Run regression model for shift in biomarkers using the GMM cluster
  # assignments as the response
  y <- z
  y <- as.numeric(y)

  logit_X <- as.matrix(logit_X)
  y <- as.matrix(y)

  logit_result <- logit_CAVI(X=logit_X, y=y, prior=logit_prior,
                             tol=logit_tol, maxiter=logit_maxiter,
                             verbose=logit_verbose)
  coeff <- logit_result$mu

  # Add the log odds and probability of latent phenotype
  . <- log_odds <- NULL
  . <- prob <- NULL
  logit_dt <- as.data.table(logit_X)
  logit_dt[,log_odds:=sum(coeff * .SD), by = 1:nrow(logit_dt)]
  logit_dt[,prob:=exp(log_odds)/(1 + exp(log_odds))]

  # Mean probability of latent phenotype in the EHR cohort
  prevalence <- mean(logit_dt$prob) * 100

  # Shift in Biomarkers for latent phenotype
  bio_shift <- vector(mode="list", length=length(biomarkers))
  df <- as.data.frame(logit_X)
  for (i in 1:length(biomarkers)) {
    idx <- grep(biomarkers[i], colnames(df))
    if(sign(coeff[idx]) < 0) {
      bio_shift[[i]] <- colMeans(df[,idx,drop = FALSE])*abs(coeff[idx]) - colMeans(df[,idx,drop = FALSE])
    } else {
      bio_shift[[i]] <- colMeans(df[,idx,drop = FALSE])*coeff[idx]
    }
  }
  biomarker_shift <- data.frame(biomarker=biomarkers, shift=unlist(bio_shift))

  return(list(prevalence=prevalence,
              biomarker_shift=biomarker_shift,
              gmm=gmm_result,
              logit=logit_result))
}
