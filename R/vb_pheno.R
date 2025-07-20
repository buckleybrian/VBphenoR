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

vbPhenoR.env <- new.env(parent = emptyenv())

#' Load Patient Data
#'
#' @param path The path to the patient data.
#'
#' @return loaded data.
#' @export
#'
#' @importFrom utils head
#'
loadData <- function(path) {

  # Dummy content to be filled
  data <- data.frame(x=seq(from=1, to=1000, by=1), y=rep(from=1, to=1000, by=1))
  vbPhenoR.env$data <- data

  class(data) <- c("vbPheno", class(data))
  return(data)
}

#' Run the Variational Bayes patient Phenotyping model
#'
#' @param ehr_data The EHR data.
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
#' * prob - The probability of phenotype given the data and priors.
#' * biomarker_shift - A data frame containing the biomarker shifts from normal for the phenotype.
#' * sensitivity - The sensitivity and specificity for the dichotomous clinical codes.
#' * gmm - The VB GMM results. For details see help(vb_gmm_cavi).
#' * logit - The VB Logit results.  For details see help(logit_CAVI).
#' @export
#'
runModel <- function(ehr_data,
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

  # Use proper label-switching - see label.switching R package
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

  # Now Run regression model for shift in biomarkers and sensitivity analysis
  # using the GMM cluster assignments as the response

  y <- z
  y <- as.numeric(y)

  logit_X <- as.matrix(logit_X)
  y <- as.matrix(y)

  logit_result <- logit_CAVI(X=logit_X, y=y, prior=logit_prior,
                             tol=logit_tol, maxiter=logit_maxiter,
                             verbose=logit_verbose)

  coeff <- logit_result$mu

  # [1,]  6.400845359 Intercept
  # [2,]  0.004672085 age
  # [3,]  0.151123043 highrisk
  # [4,] -1.560289320 CBC
  # [5,]  4.661329073 RC

  log_odds_SCD <- coeff[1] +(coeff[2]*mean(logit_X[,'age'])) +
                  (coeff[3]*mean(logit_X[,'highrisk'])) +
                  (coeff[4]*mean(logit_X[,'CBC'])) +
                  (coeff[5]*mean(logit_X[,'RC']))
  prob_SCD <- exp(log_odds_SCD)/(1 + exp(log_odds_SCD))

  options(scipen = 999)
  prob_SCD*100 # % probability of SCD

  # Shift in Biomarkers for SCD phenotype
  CBC_shift <- (mean(logit_X[,'CBC'])*abs(coeff[4])) - mean(logit_X[,'CBC'])   # 7.925
  RC_shift <- mean(logit_X[,'RC'])*abs(coeff[5]) # 3.67

  return(list(prob=prob_SCD,
              biomarker_shift=data.frame(CBC_shift=CBC_shift,RC_shift=RC_shift),
              sensitivity="TBD",
              gmm=gmm_result,
              logit=logit_result))
}

#' Variational Bayes Patient Phenotyping.
#'
#' Print a vbPheno object
#' @param x An object of class 'vbPheno'
#' @param ... Additional arguments to be passed onto lower-level functions at a later stage of development.
#' @export
print.vbPheno <- function(x, ...) {
  # TODO - A formatted table of mean biomarker shifts and clinical code sensitivity/specificity
  knitr::kable(head(as.data.frame(x),10), format = "markdown", row.names = FALSE)
}

#' Variational Bayes Patient Phenotyping.
#'
#' Summary method for class 'vbPheno'
#'
#' @param object An object of class 'vbPheno'.
#' @param ... Additional arguments to be passed onto lower-level functions at a later stage of development.
#' @export
summary.vbPheno <- function(object, ...) {
  # TODO - a formatted list of posterior coefficients with credibility intervals
  summary.data.frame(object, ...)
}

#' Variational Bayes Patient Phenotyping.
#'
#' Plot Model Summaries, Biomarker shifts and Model Diagnostics for Variational Bayes Patient Phenotyping.
#'
#' @param x An object of class 'vbPheno'.
#' @param ... Additional arguments passed to the plot method
#' @export
plot.vbPheno <- function(x, ...) {
  # TODO - a grid of plots - ROC for sensitivity analysis and forest plot for biomarkers
  plot(x[,1], x[,-1], pch = 19,
       xlab = names(x)[8], ylab = names(x)[7],
       main = "Example VB Phenotype Plot")
}
