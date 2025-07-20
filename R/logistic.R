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
#' @param tol The ELBO difference tolerance for conversion.
#' @param maxiter The maximum iterations if convergence is not achieved.
#' @param verbose A diagnostics flag added by Buckley et al.
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
logit_CAVI <- function(X, y, prior, tol=1e-16, maxiter=10000, verbose=FALSE){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  # Internal class function to compute the log-determinant of a matrix
  ldet <- function(X) {
    if(!is.matrix(X)) return(log(X))
    determinant(X,logarithm = TRUE)$modulus
  }

  lowerbound <- numeric(maxiter) # vector to store ELBO iterations
  lbdiff <- numeric(maxiter)     # vector to store difference in ELBO iterations
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

  # Initialise a progress bar
  pb <- txtProgressBar(min = 1, max = maxiter, style = 3)

  for(t in 2:maxiter){
    # Print progress
    setTxtProgressBar(pb, t)

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

    # BB added the difference for diagnostics
    lbdiff[t] <- lowerbound[t] - lowerbound[t-1]

    if(verbose) print(paste0("[",t,"]: ", lowerbound[t], " : ", lbdiff[t]))


    if(abs(lbdiff[t]) < tol) {

      close(pb)
      return(list(mu=matrix(mu_vb,p,1),
                  Sigma=matrix(Sigma_vb,p,p),
                  Convergence=cbind(Iteration=(1:t)-1, Lowerbound=lowerbound[1:t]),
                  LBDifference=cbind(Iteration=(1:t)-1, LBDiff=lbdiff[1:t]),
                  xi=xi))
    }
  }
  close(pb)

  print("The algorithm has not reached convergence")

  return(list(error="The algorithm has not reached convergence",
              mu=matrix(mu_vb,p,1),
              Sigma=matrix(Sigma_vb,p,p),
              Convergence=cbind(Iteration=(1:t)-1, Lowerbound=lowerbound[1:t]),
              LBDifference=cbind(Iteration=(1:t)-1, LBDiff=lbdiff[1:t]),
              xi=xi))
}


#' Variational inference for Bayesian logistic regression using SVI algorithm
#'
#' @param X The input design matrix. Note the intercept column vector is assumed included.
#' @param y The binary response.
#' @param prior Prior for the logistic parameters.
#' @param iter The maximum number of iterations if convergence is not achieved.
#' @param tau The delay down-weighting early iterations if >= 0.
#' @param kappa The forgetting rate.  Must be between 0.5 and 1.0
#'
#' @return A list containing:
#' * mu - A vector of posterior means.
#' * Sigma - A vector of posterior variances.
#' @export
#'
logit_SVI <- function(X, y, prior, iter, tau, kappa){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  p    <- ncol(X)
  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)

  # Initialization of Eta1 and Eta2
  Eta1     <- Pmu
  Eta2     <- P

  Eta1_out <- Eta1
  Eta2_out <- Eta2

  # Iterative procedure
  for(t in 1:iter){

    # Sample the observation
    id  <- sample.int(n,1)
    x_i <- X[id,]
    y_i <- y[id]

    # Update the local parameter
    Sigma_vb  <- solve(Eta2_out)
    mu_vb     <- Sigma_vb%*%Eta1_out

    eta_i   <- c(crossprod(x_i, mu_vb))
    xi_i    <- sqrt(eta_i^2 +  rowSums(x_i %*% Sigma_vb * x_i))
    omega_i <- tanh(xi_i/2)/(2*xi_i);

    Eta1    <- n*x_i*(y_i-0.5) + Pmu
    Eta2    <- n*tcrossprod(x_i*omega_i,x_i) + P

    # Update the final estimates
    rho      <- 1/(t + tau)^kappa
    Eta1_out <- (1 - rho)*Eta1_out + rho*Eta1
    Eta2_out <- (1 - rho)*Eta2_out + rho*Eta2
  }

  # Output
  Sigma_vb   <- matrix(solve(Eta2_out),p,p)
  mu_vb      <- matrix(Sigma_vb%*%Eta1_out,p,1)

  return(list(mu=mu_vb,Sigma=Sigma_vb))
}


#' Maximum likelihood for the logistic regression using Newton-Raphson
#'
#' @param X The input design matrix. Note the intercept column vector is assumed included.
#' @param y The binary response.
#' @param tol The Newton-Raphson difference tolerance for conversion.
#' @param beta_start The starting guess for beta.
#' @param maxiter The maximum number of iterations if convergence is not achieved.
#'
#' @return A list containing:
#' * beta - A vector of posterior means.
#' * Convergence - A vector of iterations up to convergence.
#' * Loglikelihood - A vector of log likelihoods for all iterations.
#' @export
#'
logit_NR <- function(X, y,  tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitly initialized at beta=0)
  if(is.null(beta_start)) {
    beta <- solve(crossprod(X/4,X),crossprod(X,y-0.5))
  } else {
    beta <- beta_start
  }
  # Initialization
  eta       <- c(X%*%beta)
  prob      <- 1/(1+exp(-eta))
  w         <- prob*(1-prob)

  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- beta + solve(qr(crossprod(X*w,X)),crossprod(X,y-prob))
    eta        <- c(X%*%beta)
    prob       <- 1/(1+exp(-eta))
    w          <- prob*(1-prob)
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol)
      return(list(beta=beta, Convergence=cbind(Iteration=(1:t)-1,
                  Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}


#' Maximum likelihood for the logistic regression using Polya-gamma EM
#'
#' @param X The input design matrix. Note the intercept column vector is assumed included.
#' @param y The binary response.
#' @param tol The Polya-gamma EM difference tolerance for conversion.
#' @param beta_start The starting guess for beta.
#' @param maxiter The maximum number of iterations if convergence is not achieved.
#'
#' @return A list containing:
#' * beta - A vector of posterior means.
#' * Convergence - A vector of iterations up to convergence.
#' * Loglikelihood - A vector of log likelihoods for all iterations.
#' @export
#'
logit_EM <- function(X, y,  tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  loglik <- numeric(maxiter)
  Xy     <- crossprod(X,y-0.5)

  # Initialization (If null, implicitly initialized at beta=0)
  if(is.null(beta_start)) {
    beta <- solve(crossprod(X/4,X),crossprod(X,y-0.5))
  } else {
    beta <- beta_start
  }

  # Initialization
  eta        <- c(X%*%beta)
  w          <- tanh(eta/2)/(2*eta);
  w[is.nan(w)] <- 0.25


  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- solve(qr(crossprod(X*w,X)),Xy)
    eta        <- c(X%*%beta)
    w          <- tanh(eta/2)/(2*eta);
    w[is.nan(w)] <- 0.25

    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol)
      return(list(beta=beta, Convergence=cbind(Iteration=(1:t)-1,Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}


#' Maximum likelihood for the logistic regression using Bohning bound
#'
#' @param X The input design matrix. Note the intercept column vector is assumed included.
#' @param y The binary response.
#' @param tol The Bohning bound difference tolerance for conversion.
#' @param beta_start The starting guess for beta.
#' @param maxiter The maximum number of iterations if convergence is not achieved.
#'
#' @return A list containing:
#' * beta - A vector of posterior means.
#' * Convergence - A vector of iterations up to convergence.
#' * Loglikelihood - A vector of log likelihoods for all iterations.
#' @export
#'
logit_MM <- function(X, y,  tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitly initialized at beta=0)
  B <- qr(crossprod(X/4,X)) # Bohning matrix (QR version)
  if(is.null(beta_start)) {
    beta <- solve(B,crossprod(X,y-0.5))
  } else {
    beta <- beta_start
  }

  # Initialization
  eta        <- c(X%*%beta)
  prob       <- 1/(1+exp(-eta))
  w          <- prob*(1-prob)

  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- beta + solve(B,crossprod(X,y-prob))
    eta        <- c(X%*%beta)
    prob       <- 1/(1+exp(-eta))
    w          <- prob*(1-prob)
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol)
      return(list(beta=beta,Convergence=cbind(Iteration=(1:t)-1,Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}
