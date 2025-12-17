# ------------------------------------------------------------------------------
# Variational Bayes Gaussian Mixture Model (VB-GMM) implemented by Coordinate
# Ascent Variational Inference (CAVI)
#
# References:
#
# Pattern Recognition and Machine Learning by Christopher M. Bishop,
# 2006 (Section 10.2)
# Adapted from Matlab code for above by Mo Chen (2016),
# Variational Bayesian Inference for Gaussian Mixture Model
#
# Brian Buckley December 2024
# ------------------------------------------------------------------------------

# A better version of bsxfun from the Github of the user 'shouldsee':
# https://github.com/shouldsee/Rutil/blob/master/R/bsxfun.R
#source("bsxfun.R")

#' Main algorithm function for the VB CAVI GMM
#'
#' @param X n x p data matrix (or data frame that will be converted to a matrix).
#' @param k guess for the number of mixture components.
#' @param prior Prior for the GMM parameters.
#' @param delta change in ELBO that triggers algorithm stopping.
#' @param maxiters maximum iterations to run if delta does not stop the algorithm already.
#' @param init initialize the clusters c("random", "kmeans", "dbscan").
#' @param initParams initialization parameters for dbscan.  NULL if dbscan not selected for init.
#' @param stopIfELBOReverse stop the run if the ELBO at iteration t is detected to have reversed from iteration t-1.
#' @param verbose print out information per iteration to track progress in case of long-running experiments.
#' @param logDiagnostics log detailed diagnostics.  If TRUE, a diagnostics RDS file will be created using the path specified in logFilename.
#' @param logFilename the filename of the diagnostics log.
#' @param progressbar A visual progress bar to indicate iterations (on by default).
#'
#' @return A list containing:
#' * error - An error message if convergence failed or the number of iterations to achieve convergence.
#' * z_post - A vector of posterior cluster mappings.
#' * q_post - A list of the fitted posterior Q family. q_post includes:
#'   * alpha - The Dirichlet prior for the mixing coefficient
#'   * lambda - The proportionality of the precision for the Normal-Wishart prior (called 'beta' in Bishop)
#'   * m - The mean vector for the Normal part of the Normal-Wishart prior
#'   * v - The degrees of freedom of the Wishart gamma ensuring it is well defined
#'   * U - The W hyperparameter of the Wishart conjugate prior for the precision of the mixtures. k sets of D x D symmetric, positive definite matrices.
#'   * logW - The logW term used to calculate the expectation of the mixture component precision. A vector of k.
#'   * R - The responsibilities. An n x k matrix.
#'   * logR - The log of responsibilities. An n x k matrix.
#' * elbo - A vector of the ELBO at each iteration.
#' * elbo_d - A vector of ELBO differences between each iteration.
#' @export
#'
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @importFrom stats model.matrix
#' @importFrom stats kmeans
#' @importFrom utils modifyList
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom dbscan dbscan
#' @importFrom CholWishart lmvgamma
#' @importFrom pracma bsxfun
#'
#' @example man/examples/vb_gmm_example_priors.R
#'
vb_gmm_cavi <- function(X, k,
                   prior=NULL,
                   delta=1e-6,
                   maxiters=5000,
                   init="kmeans",
                   initParams=NULL,
                   stopIfELBOReverse=FALSE,
                   verbose=FALSE,
                   logDiagnostics=FALSE,
                   logFilename="vb_gmm_log.txt",
                   progressbar=TRUE) {

  # Cannot have an empty data
  if (is.null(n <- nrow(X)))
    stop("X must be a data frame or matrix")

  # Cannot have missing values
  if(anyNA(X))
    stop("The VB GMM cannot run with missing values in X")

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X <- t(X)

  if(init=="dbscan" & is.null(initParams)) {
    stop("dbscan initialisation requires initParams - c(ets, minPts)")
  }

  if(logDiagnostics) {
    stage <- c("Base", "Init", "Expectation", "Maximisation", "ELBO")
    logs <- c()
    log_q_post <- vector(mode="list")
  }

  # Initialise with uninformative priors.  Overwrite hyper-parameters where set
  this_prior <- list(
    alpha = 1,
    lambda = 1,
    m = matrix(rowMeans(X),ncol=1),
    v = p+1,
    W = diag(1,p),
    logW = -2*sum(log(diag(chol(diag(1,p)))))
  )

  if(!missing(prior) | !is.null(prior)) {
    this_prior <- modifyList(this_prior, prior)
  }

  if(verbose) {
    print("Using Prior:\n")
    print(this_prior)
  }

  if(logDiagnostics) {
    logs <- c(logs, paste0("stage: ", "Base", " itr: ", 0, " q_post: ", ""))
    log_q_post[[1]] <- list("Prior", this_prior)
  }

  q_post <- VB_GMM_Init(X, k, n, this_prior, init, initParams)

  if(logDiagnostics) {
    logs <- c(logs, paste0("stage: ", "Init", " itr: ", 0, " q_post: ", q_post[1:6]) )
    I <- list("Init", q_post)
    log_q_post <- c(log_q_post, I)
  }

  lb <- rep(-Inf, maxiters)                # This is a vector of the ELBO lower bound calculations
  lbd <- rep(0, maxiters)                  # This is a vector of the ELBO differences per iteration
  converge <- FALSE                        # Indicator when converged
  itr <- 1                                 # Iteration counter

  if(progressbar==TRUE) {
    # Initialise a progress bar
    pb <- txtProgressBar(title = "VB GMM", min = 1, max = maxiters, style = 3)
  }

  while(!converge & itr <= maxiters) {

    itr <- itr + 1                         # Increment here so that we can compare with first -Inf

    if(progressbar==TRUE) setTxtProgressBar(pb, itr)

    q_post  <- VB_GMM_Expectation(X, n, q_post)
    if(logDiagnostics) {
      logs <- c(logs, paste0("stage: ", "Expectation", " itr: ", itr, " q_post: ", q_post[1:6]) )
      E <- list("Expectation", q_post)
      log_q_post <- c(log_q_post, E)
    }

    q_post  <- VB_GMM_Maximisation(X, q_post, this_prior)
    if(logDiagnostics) {
      logs <- c(logs, paste0("stage: ", "Maximisation", " itr: ", itr, " q_post: ", q_post[1:6]) )
      M <- list("Maximisation", q_post)
      log_q_post <- c(log_q_post, M)
    }

    lb[itr] <- VB_GMM_ELBO(X, p, n, q_post, this_prior)/n
    if(logDiagnostics) {
      logs <- c(logs, paste0("stage: ", "ELBO", " itr: ", itr, " q_post: ", q_post[1:6]) )
      EB <- list("ELBO", q_post)
      log_q_post <- c(log_q_post, EB)
    }

    lbd[itr] <- lb[itr] - lb[itr-1]

    # Decide what to do if the ELBO starts reversing direction (check after initial settle)
    elbo_reversed <- FALSE
    if(stopIfELBOReverse & itr > 3) {
      if(sign(lbd[itr]) != sign(lb[itr-1])) {
        elbo_reversed <- TRUE
        break
      }
    }

    # Stopping rule for convergence of the ELBO
    converge <- abs(lbd[itr]) < delta * abs(lb[itr])

    if(verbose) {
      if(sign(lbd[itr]) != sign(lb[itr-1])) {
        print(paste0("VB iteration ",itr,": ELBO = ",lb[itr]," : ELBO Diff = ", lbd[itr], " ELBO is reversing!!"))
      } else {
        print(paste0("VB iteration ",itr,": ELBO = ",lb[itr]," : ELBO Diff = ", lbd[itr]))
      }
    }
  }
  if(progressbar==TRUE) close(pb)

  if(logDiagnostics) {
    logDiagnostics(filename=logFilename, logs=logs)
    saveRDS(log_q_post, paste0(gsub("(^.*\\.)(.*$)", "\\1", logFilename),"RDS"))
  }

  elbo   <- lb[c(2,2:itr)]                   # Ignore the first -Inf element so we have a numeric vector
  elbo_d <- lbd[c(2,2:itr)]
  z_post <- max.col(q_post$R,ties.method="random")

  conv <- ifelse(!converge,
                  ifelse(elbo_reversed,
                         paste0("The algorithm has not reached convergence after ",itr," iterations.  The ELBO reversed direction."),
                         paste0("The algorithm has not reached convergence after ",itr," iterations.")
                  ),
                  paste0("The algorithm converged after ",itr," iterations.")
  )

  return(list(conv=conv,z_post=z_post, q_post=q_post, elbo=elbo, elbo_d=elbo_d))
}


#' Initialise the variational parameters and the hyper parameters
#'
#' @param X n x p data matrix (or data frame that will be converted to a matrix).
#' @param k guess for the number of mixture components.
#' @param n number of rows in X.
#' @param prior Prior for the GMM parameters.
#' @param init initialize the clusters c("random", "kmeans", "dbscan")
#' @param initParams initialization parameters for dbscan.  NULL if dbscan not selected for init.
#'
#' @return A list of the initially fitted posterior Q family
#'
VB_GMM_Init <- function(X, k, n, prior, init, initParams) {

  q_post <- list()

  # BB to account for constraining of the mixing weights via the prior
  # The else choice is an informative prior

  # Dirichlet for q(pi)
  if(length(prior$alpha) == 1) {
    q_post$alpha <- rep(prior$alpha, k)
  } else {
    q_post$alpha <- prior$alpha
  }

  # Governs the proportionality of the precision, Λ and influences the variance of μ.
  if(length(prior$lambda) == 1) {
    q_post$lambda <- rep(prior$lambda, k)
  } else {
    q_post$lambda <- prior$lambda
  }

  q_post$m     <- rep(prior$m, k)
  q_post$v     <- prior$v
  q_post$U     <- 0.0
  q_post$logW  <- prior$logW
  q_post$R     <- rep(0, n)
  q_post$logR  <- 0.0

  if(init == "kmeans") {
    # try k-means to initialize the clusters - this will provide a good guess for assignment of
    # observations to the k clusters
    km <- kmeans(t(X), k, nstart = 25)
    z <- km$cluster
  } else if(init == "dbscan") {
    dbs <- dbscan(t(X), eps = initParams['eps'], minPts = initParams['minPts'])
    z <- dbs$cluster

    # Align DBSCAN clusters with k by assigning smallest clusters to highest cluster
    # - replace if/when HDSCAN solution completes in reasonable time
    if(length(unique(z)) > k) {
      ztab <- as.data.frame(table(z))
      ztab$z <- as.character(ztab$z)
      ztab <- ztab[order(ztab$Freq),]
      dropz <- ztab$z[1:(nrow(ztab)-k)]
      z[z %in% dropz] <- ztab[nrow(ztab),]$z
    }
  } else if(init == "random") {
    # Random assignment of clusters - this will evenly distribute k clusters across the observations
    z <- round(as.integer(runif(n,min=1,max=k+1)))
  } else {
    # Should be a proportion between 0 and 1 - this will randomly distribute k clusters across the observations subject to the proportion selected
    z <- rbinom(n, 1, as.numeric(init))
  }

  q_post$R <- matrix(as.integer(model.matrix(~factor(z) - 1)),ncol=k)

  # Initial call to maximization
  q_post <- VB_GMM_Maximisation(X, q_post, prior)

  return(q_post)
}


#' Variational Bayes Expectation step
#'
#' @param X n x p data matrix (or data frame that will be converted to a matrix).
#' @param n number of rows in X.
#' @param q_post A list of the fitted posterior Q family at iteration t-1.
#'
#' @return A list of the fitted posterior Q family after the E step at iteration t.
#'
VB_GMM_Expectation <- function(X, n, q_post) {

  alpha <- q_post$alpha
  lambda <- q_post$lambda
  m     <- q_post$m
  v     <- q_post$v
  U     <- q_post$U
  logW  <- q_post$logW

  p <- nrow(m)
  k <- ncol(m)

  EQ <- array(0, dim=c(n,k))
  for(i in 1:k) {
    Q      <- solve(t(U[,,i]), m_bsxfun("-", X, m[,i]))       # Equations in Bishop, 2006 (Section 10.2)
    EQ[,i] <- p/lambda[i] + v[i] * pracma::dot(Q,Q)                                  # 10.64
  }
  vp         <- m_bsxfun("-", matrix(rep(v+1, p),nrow=p,byrow=TRUE), as.matrix(1:p))
  ElogLambda <- colSums(digamma(vp)) + p*log(2) + logW                             # 10.65
  Elogpi     <- digamma(alpha) - digamma(sum(alpha))                               # 10.66
  logRho     <- (-0.5) * (m_bsxfun("-", EQ, ElogLambda - p*log(2*pi)))             # 10.46
  logRho     <- m_bsxfun("+", logRho, Elogpi)                                      # 10.46
  logR       <- m_bsxfun("-", logRho, logsumexp(logRho, 1))                        # 10.49
  R          <- exp(logR)                                                          # 10.67

  q_post$logR <- logR
  q_post$R <- R

  return(q_post)
}


#' Variational Bayes Maximisation step
#'
#' @param X n x p data matrix (or data frame that will be converted to a matrix).
#' @param q_post A list of the fitted posterior Q family at iteration t-1.
#' @param prior Prior for the GMM parameters.
#'
#' @return A list of the fitted posterior Q family after the M step at iteration t.
#'
VB_GMM_Maximisation <- function(X, q_post, prior) {

  alpha0 <- prior$alpha
  # lambda0  <- prior$lambda
  lambda0  <- prior$lambda[1]
  m0     <- prior$m
  v0     <- prior$v
  W0     <- prior$W
  R      <- q_post$R
                                                          # Equations in Bishop, 2006 (Section 10.2)
  nk <- colSums(R)                                        # 10.51
  alpha <- alpha0 + nk                                    # 10.58
  lambda <- lambda0 + nk                                  # 10.60

  m <- m_bsxfun("+", X %*% R, lambda0 * m0)
  m <- m_bsxfun("*", m, 1/lambda)                         # 10.61
  v <- v0 + nk                                            # 10.63

  k <- ncol(m)
  p <- nrow(m)
  U <- array(0, c(p, p, k))
  logW <- array(0, dim=c(1,k))
  r <- sqrt(t(R))

  for(i in 1:k) {
    Xm <- m_bsxfun("-", X, m[,i])
    Xm <- m_bsxfun("*", Xm, r[i,])
    m0m <- m0-m[,i]
    W <- W0+Xm %*% t(Xm) + lambda0*(m0m %*% t(m0m))       # 10.62 (equivalent)
    U[,,i] <- chol(W)
    logW[1,i] = -2*sum(log(diag(U[,,i])))
  }

  q_post$alpha <- alpha
  q_post$lambda <- lambda
  q_post$m     <- m
  q_post$v     <- v
  q_post$U     <- U
  q_post$logW  <- logW

  return(q_post)
}


#' Calculate the Evidence Lower Bound (ELBO)
#'
#' @param X n x p data matrix (or data frame that will be converted to a matrix).
#' @param p number of parameters
#' @param n number of rows
#' @param q_post A list of the fitted posterior Q family at iteration t.
#' @param prior Prior for the GMM parameters.
#'
#' @return elbo the Evidence Lower Bound at iteration t following E and M steps.
#'
VB_GMM_ELBO <- function(X, p, n, q_post, prior) {

  alpha0 <- prior$alpha
  lambda0  <- prior$lambda
  v0     <- prior$v
  logW0  <- prior$logW
  alpha  <- q_post$alpha
  lambda   <- q_post$lambda
  v      <- q_post$v
  logW   <- q_post$logW
  R      <- q_post$R
  logR   <- q_post$logR

  k <-  ncol(R)

  # Posterior Expectations
  Eqz  <- as.numeric(pracma::dot(R, logR))
  Eqz  <- mean(Eqz) # There is one Eqz value per k
  Eqpi <- as.numeric(lgamma(sum(alpha)) - sum(lgamma(alpha)))
  Eqmu <- as.numeric(0.5 * p * sum(log(lambda)))
  Eqlambda <- as.numeric(sum(-0.5 * v * (logW + p * log(2)) -
                            matrix(CholWishart::lmvgamma(0.5 * v, p),nrow=1) ))

  Eppi <- as.numeric(lgamma(k*alpha0) - k*(lgamma(alpha0)))
  Epmu <- as.numeric(0.5 * p * k * log(lambda0))
  Eplambda <- as.numeric(k * (-0.5 * v0 * (logW0 + p * log(2)) -
                            matrix(CholWishart::lmvgamma(0.5 * v0, p),nrow=1) ))

  Epx <- as.numeric(-0.5 * p * n * log(2 * pi))

  elbo <- mean(Eqz + Eppi - Eqpi + Epmu - Eqmu + Eplambda - Eqlambda + Epx)

  return(elbo)
}

#' An internal numerically stable implementation of matrixStats::logSumExp.
#'
#' @description
#' The implementation of logSumExp from the matrixStats package is not numerically
#' stable for the variational updates. This implementation from the Matlab code by Mo Chen (see ref)
#' effectively handles numerical underflow.
#'
#' @param x lhs matrix.
#' @param margin margin size.
#'
#' @return calculated logSumExp.
#'
#' @noRd
#'
logsumexp <- function(x, margin=1) {

  if ( ! is.matrix(x) ) {
    x <- as.matrix(x)
  }

  # subtract the largest in each column
  y <- apply(x, margin, max)

  if (nrow(x) == ncol(x)) {
    x <- m_bsxfun("-", x, y, expandByRow = FALSE)
  } else {
    x <- m_bsxfun("-", x, y)
  }

  s <- y + log(apply(exp(x), margin, sum))
  i <- which(!is.finite(s))
  if(length(i) > 0) s[i] <- y[i]

  return(s)
}

#' An internal implementation of the MATLAB bsxfun() function.
#'
#' @description
#' R implementation of Matlab bsxfun() that fixes limitations in the version from the pracma R package.
#' Specifically, the Matlab version allows for implicit expansion but the R version assumes equal sized
#' data structures.
#'
#' @param func Algebraic operator (as a string) passed to pracma::bsxfun.
#' @param x lhs matrix.
#' @param y rhs matrix.
#' @param expandByRow expand by row if true, by column if false.
#'
#' @return calculated bsxfun.
#'
#' @noRd
#'
m_bsxfun <- function(func, x, y, expandByRow=TRUE) {
  if(length(y) == 1) return(pracma::arrayfun(func, x, y))

  expandCol <- nrow(x) == length(y)
  expandRow <- ncol(x) == length(y)

  if(expandCol & expandRow & expandByRow) expandCol <- FALSE
  if(expandCol & expandRow & !expandByRow) expandRow <- FALSE
  if(expandRow) y.repmat <- matrix(rep(as.numeric(y), each=nrow(x)), nrow=nrow(x))
  if(expandCol) y.repmat <- matrix(rep(as.numeric(y), ncol(x)), ncol=ncol(x))

  return(pracma::bsxfun(func, x, y.repmat))
}

#' An internal diagnostics logger
#'
#' @description
#' A simple function that logs the diagnostics of each iteration in a text file.
#'
#' @param filename  the log file.
#' @param logs      the vector of all logged iterations.
#'
#' @noRd
#'
logDiagnostics <- function(filename, logs) {
  writeLines(logs, filename)
}
