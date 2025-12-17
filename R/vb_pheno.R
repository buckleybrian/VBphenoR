# ------------------------------------------------------------------------------
# Variational Bayes Patient Phenotyping from Electronic Health Records
#
# References:
#
# Hubbard, Rebecca A., et al. "A Bayesian latent class approach for EHR‐based
# phenotyping." Statistics in Medicine 38.1 (2019): 74-87.
#
# Brian Buckley August 2025
# ------------------------------------------------------------------------------

VBphenoR.env <- new.env(parent = emptyenv())

#' Run the Variational Bayes patient phenotyping model
#'
#' @param biomarkers The EHR variables that are biomarkers. This is a vector of data column names corresponding to the biomarker variables.
#' @param gmm_X n x p data matrix (or data frame that will be converted to a matrix).
#' @param logit_X The input design matrix. Note the intercept column vector is assumed included.
#' @param gmm_delta Change in ELBO that triggers algorithm stopping.
#' @param logit_delta Change in ELBO that triggers algorithm stopping.
#' @param gmm_maxiters The maximum iterations for VB GMM.
#' @param logit_maxiters The maximum iterations for VB logit.
#' @param gmm_init Initialize the clusters c("random", "kmeans", "dbscan").
#' @param gmm_initParams Parameters for an initialiser requiring its own parameters e.g. dbscan requires 'eps' and 'minPts'.
#' @param gmm_prior An informative prior for the GMM.
#' @param logit_prior An informative prior for the logit.
#' @param gmm_stopIfELBOReverse Stop the VB iterations if the ELBO reverses direction (TRUE or FALSE).
#' @param gmm_verbose Print out information per iteration to track progress in case of long-running experiments.
#' @param logit_verbose Print out information per iteration to track progress in case of long-running experiments.
#' @param progressbar Show a progressbar driven by the GMM & logit variational iterations.
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
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' ##Example 1: Use the internal Sickle Cell Disease data to find the rare
#' ##           phenotype.  SCD is extremely rare so we use DBSCAN to initialise
#' ##           the VB GMM. We also use an informative prior for the mixing
#' ##           coefficient and stop iterations when the ELBO starts to reverse
#' ##           so that we stop when the minor (SCD) component is reached.
#'
#' library(data.table)
#'
#' # Load the SCD example data supplied with the VBphenoR package
#' data(scd_cohort)
#'
#' # We will use the SCD biomarkers to discover the SCD latent class.
#' # X1 is the data matrix for the VB GMM.
#' X1 <- scd_cohort[,.(CBC,RC)]
#'
#' # We need to supply DBSCAN hyper-parameters as we will initialise VBphenoR
#' # with DBSCAN. See help(DBSCAN) for details of these parameters.
#' initParams <- c(0.15, 5)
#' names(initParams) <- c('eps','minPts')
#'
#' # Set an informative prior for the VB GMM mixing coefficient alpha
#' # hyper-parameter
#' prior_gmm <- list(
#'   alpha = 0.001
#' )
#'
#' # Set informative priors for the beta coefficients of the VB logit
#' prior_logit <- list(mu=c(1,
#'                    mean(scd_cohort$age),
#'                    mean(scd_cohort$highrisk),
#'                    mean(scd_cohort$CBC),
#'                    mean(scd_cohort$RC)),
#'               Sigma=diag(1,5))           # Simplest isotropic case
#'
#' # X2 is the design matrix for the VB logit
#' X2 <- scd_cohort[,.(age,highrisk,CBC,RC)]
#' X2[,age:=as.numeric(age)]
#' X2[,highrisk:=as.numeric(highrisk)]
#' X2[,Intercept:=1]
#' setcolorder(X2, c("Intercept","age","highrisk","CBC","RC"))
#'
#' # Run the patient phenotyping model
#'
#' # Need to state what columns are the biomarkers
#' biomarkers <- c('CBC', 'RC')
#' set.seed(123)
#'
#' pheno_result <- run_Model(biomarkers,
#'                         gmm_X=X1, gmm_init="dbscan",
#'                         gmm_initParams=initParams,
#'                         gmm_maxiters=20, gmm_prior=prior_gmm,
#'                         gmm_stopIfELBOReverse=TRUE,
#'                         logit_X=X2, logit_prior=prior_logit
#')
#'
#' # S3 print method
#' print(pheno_result)
#' }
#'
run_Model <- function(biomarkers, gmm_X, logit_X, gmm_delta=1e-6, logit_delta=1e-16,
                    gmm_maxiters=200, logit_maxiters=10000,
                    gmm_init="kmeans", gmm_initParams=NULL,
                    gmm_prior=NULL, logit_prior=NULL,
                    gmm_stopIfELBOReverse=FALSE,
                    gmm_verbose=FALSE, logit_verbose=FALSE,
                    progressbar=FALSE) {

  gmm_result <- vb_gmm_cavi(X=gmm_X, k=2, delta=gmm_delta,
                            init=gmm_init, initParams=gmm_initParams,
                            maxiters=gmm_maxiters, prior=gmm_prior,
                            stopIfELBOReverse=gmm_stopIfELBOReverse,
                            verbose=gmm_verbose, progressbar=progressbar)

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
  # assignments as the response (either using binary z or soft probabilities)
  y <- z
  y <- as.numeric(y)

  logit_X <- as.matrix(logit_X)
  y <- as.matrix(y)

  logit_result <- logit_CAVI(X=logit_X, y=y, prior=logit_prior,
                             delta=logit_delta, maxiters=logit_maxiters,
                             verbose=logit_verbose, progressbar=progressbar)
  coeff <- logit_result$mu

  # Add the log odds and probability of latent phenotype
  . <- log_odds <- NULL
  . <- prob <- NULL
  logit_dt <- as.data.table(logit_X)
  logit_dt[,log_odds:=sum(coeff * .SD), by = 1:nrow(logit_dt)]
  logit_dt[,prob:=exp(log_odds)/(1 + exp(log_odds))]

  # Mean probability of latent phenotype in the EHR cohort (as a %)
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

  structure(
    list(
      prevalence=prevalence,
      biomarker_shift=biomarker_shift,
      gmm=gmm_result,
      logit=logit_result
    ),
    class = "vbphenor"
  )
}

#' @export
print.vbphenor <- function(x, ...) {
  writeLines(c(
    paste0("VBphenoR phenotyping results."),
    paste0("-----------------------------")
  ))

  cat("\n")

  writeLines(c(
    paste0("EHR data with ",formatC(length(x$gmm$z_post), big.mark=',')," rows")
  ))

  cat("\n")

  writeLines(c(
    paste0("Prevalence = ",round(x$prevalence,3),"% in these data"),
    paste0("Biomarker shift for patients with the phenotype of interest:")
  ))

  bio <- data.frame(lapply(x$biomarker_shift, function(y) if(is.numeric(y)) round(y, 2) else y))
  print(bio,row.names=FALSE)
}
