#' The vbPhenoR package
#'
#' vbPhenoR is a library and R package for a variational Bayes approach to clinical patient phenotyping using EHR data.
#'
#' @rdname aaa-vbPhenoR-package
#' @name vbPheno-package
#' @keywords internal
#' @aliases vbPhenoR-package vbPhenoR
#'
#' @section Introduction:
#' The main goal of the vbPhenoR library is to provide a comprehensive variational
#' Bayes approach to clinical patient phenotyping using Electronic Health Records (EHR)
#' data.  The phenotyping model is adapted from Hubbard et al (2019). The motivation
#' for using variational Bayes rather than the gold-standard Monte Carlo Bayes
#' approach is computational performance.  Monte Carlo is unable to cope with many
#' EHR clinical studies due to the number of observations and variables. We explain in
#' more detail in our paper, Buckley et al. (2024).
#'
#' @section VB Phenotyping algorithm:
#' The implementation of vbPhenoR performs the following steps:
#'
#' 1. Run a variational Gaussian Mixture Model using EHR-derived patient characteristics
#'    to discover the latent variable \eqn{D_i} indicating the phenotype of interest for
#'    the \eqn{i}th patient. Patient characteristics can be any patient variables typically
#'    found in EHR data e.g.
#'    - Gender
#'    - Age
#'    - Ethnicity (for disease conditions with ethnicity-related increased risk)
#'    - Physical e.g. BMI for Type 2 Diabetes
#'    - Clinical codes (diagnosis, observations, specialist visits, etc.)
#'    - Prescription medications related to the disease condition
#'    - Biomarkers (usually by laboratory tests)
#'
#' 2. Run a variational Regression model using the latent variable \eqn{D_i} derived
#'    in step 1 as an interaction effect to determine the shift in biomarker levels
#'    from baseline for patients with the phenotype versus those without. Appropriately
#'    informative priors are used to set the biomarker baseline.
#'
#' 3. Run a variational Regression model using the latent variable \eqn{D_i} derived
#'    in step 1 as an interaction effect to determine the sensitivity and specificity
#'    of binary indicators for clinical codes, medications and availability of biomarkers
#'    (since biomarker laboratory tests will include a level of missingness).
#'
#' Further details about the model can be found in the package vignette.
#'
#' @references
#' Hubbard RA, Huang J, Harton J, Oganisian A, Choi G, Utidjian L, et al. A
#' Bayesian latent class approach for EHR-based phenotyping. StatMed. (2019) 38:74–87.
#' doi:10.1002/sim.7953
#'
#' Buckley, Brian, Adrian O'Hagan, and Marie Galligan. Variational Bayes latent
#' class analysis for EHR-based phenotyping with large real-world data.
#' Frontiers in Applied Mathematics and Statistics 10 (2024): 1302825.
#' doi:10.3389/fams.2024.1302825
#'
#' @examples
#' \dontrun{
#' ##Example 1: Use the internal Sickle Cell Disease data to find the rare
#' ##           phenotype.  SCD is extremely rare so we use DBSCAN to initialise
#' ##           the VB GMM. We also use an informative prior for the mixing
#' ##           coefficient and stop iterations when the ELBO starts to reverse
#' ##           so that we stop when the minor (SCD) component is reached.
#'
#' library(data.table)
#'
#' # Load the SCD example data supplied with the vbPhenoR package
#' data(scd_cohort)
#'
#' # We will use the SCD biomarkers to discover the SCD latent class
#' X1 <- scd_cohort[,.(CBC,RC)]
#'
#' # We need to supply DBSCAN hyper-parameters as we will initialise vbPhenoR
#' # with DBSCAN. See help(DBSCAN) for details of these parameters.
#' initParams <- c(0.15, 5)
#' names(initParams) <- c('eps','minPts')
#'
#' # X1 is the data matrix for the VB GMM
#' X1 <- t(X1)
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
#' set.seed(123)
#' pheno_result <- runModel(scd_cohort,
#'                         gmm_X=X1, gmm_k=k, gmm_init="dbscan",
#'                         gmm_initParams=initParams,
#'                         gmm_maxiters=20, gmm_prior=prior_gmm,
#'                         gmm_stopIfELBOReverse=TRUE,
#'                         logit_X=X2, logit_prior=prior_logit
#')
#'
#' # Biomarker shifts for phenotype of interest
#' pheno_result
#' }
#'
"_PACKAGE"
