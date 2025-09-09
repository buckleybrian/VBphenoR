#' The VBphenoR package
#'
#' VBphenoR is a library and R package for a variational Bayes approach to clinical patient phenotyping using EHR data.
#'
#' @rdname VBphenoR
#' @name VBphenoR-package
#' @keywords internal
#' @aliases VBphenoR VBphenoR
#'
#' @section Introduction:
#' The main goal of the VBphenoR package is to provide a comprehensive variational
#' Bayes approach to clinical patient phenotyping using Electronic Health Records (EHR)
#' data.  The phenotyping model is adapted from Hubbard et al (2019). The motivation
#' for using variational Bayes rather than the gold-standard Monte Carlo Bayes
#' approach is computational performance.  Monte Carlo is unable to cope with many
#' EHR clinical studies due to the number of observations and variables. We explain in
#' more detail in our paper, Buckley et al. (2024).
#'
#' @section VB Phenotype algorithm:
#' The implementation of VBphenoR performs the following steps:
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
"_PACKAGE"
