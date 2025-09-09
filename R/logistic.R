# ----------------------------------------------------------------------------
# Variational Bayes logistic regression.
#
# This code is from Durante & Rigon (https://github.com/tommasorigon/logisticVB)
# and is subject to the following licence:
#
# MIT License
#
# Copyright (c) 2017 Tommaso Rigon
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# We have made some changes and amendments to the original software.
#
# Brian Buckley December 2024
# ----------------------------------------------------------------------------


#' Variational inference for Bayesian logistic regression using CAVI algorithm
#'
#' @param X The input design matrix. Note the intercept column vector is assumed included.
#' @param y The binary response.
#' @param prior Prior for the logistic parameters.
#' @param delta The ELBO difference tolerance for conversion.
#' @param maxiters The maximum iterations to run if convergence is not achieved.
#' @param verbose A diagnostics flag (off by default).
#' @param progressbar A visual progress bar to indicate iterations (on by default).
#'
#' @return A list containing:
#' * error - An error message if convergence failed or the number of iterations to achieve convergence.
#' * mu - A vector of posterior means.
#' * Sigma - A vector of posterior variances.
#' * Convergence - A vector of the ELBO at each iteration.
#' * LBDifference - A vector of ELBO differences between each iteration.
#' * xi - A vector of log-odds per X row.
#' @export
#'
#' @importFrom stats plogis
#'
#' @example man/examples/vb_gmm_example_priors.R
#'
logit_CAVI <- function(X,
                       y,
                       prior,
                       delta=1e-16,
                       maxiters=10000,
                       verbose=FALSE,
                       progressbar=TRUE){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  # Internal class function to compute the log-determinant of a matrix
  ldet <- function(X) {
    if(!is.matrix(X)) return(log(X))
    determinant(X,logarithm = TRUE)$modulus
  }

  lowerbound <- numeric(maxiters) # vector to store ELBO iterations
  lbdiff <- numeric(maxiters)     # vector to store difference in ELBO iterations
  p    <- ncol(X)
  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)
  Pdet <- ldet(P)

  # Initialization for omega equal to 0.25
  X1 <- X*rep(1/4,n)
  P_vb       <- crossprod(X1,X) + P
  Sigma_vb   <- solve(P_vb)
  mu_vb      <- Sigma_vb %*% (crossprod(X,y - 0.5) + Pmu)
  eta        <- c(X%*%mu_vb)
  xi         <- sqrt(eta^2 + rowSums(X %*% Sigma_vb * X))
  omega      <- tanh(xi/2)/(2*xi);
  omega[is.nan(omega)] <- 0.25

  lowerbound[1]  <- 0.5*p +
                    0.5*ldet(Sigma_vb) +
                    0.5*Pdet -
                    0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) +
                    sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) -
                    0.5*sum(diag(P %*% Sigma_vb))

  if(progressbar==TRUE) {
    # Initialise a progress bar
    pb <- txtProgressBar(min = 1, max = maxiters, style = 3)
  }

  for(t in 2:maxiters){

    if(progressbar==TRUE) setTxtProgressBar(pb, t)

    P_vb       <- crossprod(X*omega,X) + P
    Sigma_vb   <- solve(P_vb)
    mu_vb      <- Sigma_vb %*% (crossprod(X,y-0.5) + Pmu)
    eta        <- c(X%*%mu_vb)
    xi         <- sqrt(eta^2 +  rowSums(X %*% Sigma_vb * X))
    omega      <- tanh(xi/2)/(2*xi);
    omega[is.nan(omega)] <- 0.25

    lowerbound[t]  <- 0.5*p +
                      0.5*ldet(Sigma_vb) +
                      0.5*Pdet -
                      0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) +
                      sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) -
                      0.5*sum(diag(P %*% Sigma_vb))

    lbdiff[t] <- lowerbound[t] - lowerbound[t-1]

    if(verbose) print(paste0("[",t,"]: ", lowerbound[t], " : ", lbdiff[t]))

    if(abs(lbdiff[t]) < delta) {
      if(progressbar==TRUE) close(pb)
      if(verbose) print(paste0("The algorithm has converged in ", t, " iterations"))
      return(list(mu=matrix(mu_vb,p,1),
                  Sigma=matrix(Sigma_vb,p,p),
                  Convergence=cbind(Iteration=(1:t)-1, Lowerbound=lowerbound[1:t]),
                  LBDifference=cbind(Iteration=(1:t)-1, LBDiff=lbdiff[1:t]),
                  xi=xi))
    }
  }
  if(progressbar==TRUE) close(pb)
  if(verbose) print("The algorithm has not reached convergence")

  return(list(error="The algorithm has not reached convergence",
              mu=matrix(mu_vb,p,1),
              Sigma=matrix(Sigma_vb,p,p),
              Convergence=cbind(Iteration=(1:t)-1, Lowerbound=lowerbound[1:t]),
              LBDifference=cbind(Iteration=(1:t)-1, LBDiff=lbdiff[1:t]),
              xi=xi))
}

