\donttest{
  # Use rare Sickle Cell Disease data to show the effect of VB GMM initialisation
  # method
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Plotting utility
  # ------------------------------------------------------------------------------
  #' Plots the GMM components with centroids
  #'
  #' @param scd_cohort The SCD data set
  #' @param z The cluster allocation vector from the VB GMM
  #'
  #' @returns The ggplot figure (p)
  #' @export
  #'
  do_init_plot <- function(scd_cohort, gmm_result) {
    dd <- as.data.frame(cbind(scd_cohort, cluster=gmm_result$z_post))
    dd$cluster <- as.factor(dd$cluster)
    dd$highrisk <- as.factor(dd$highrisk)
    mu <- as.data.frame( t(gmm_result$q_post$m) )

    cols <- c("#1170AA", "#55AD89", "#EF6F6A", "#D3A333", "#5FEFE8", "#11F444")
    p <- ggplot() +
      geom_point(dd, mapping=aes(x=CBC, y=RC, color=cluster)) +
      scale_color_discrete(cols, guide = 'none') +
      geom_point(mu, mapping=aes(x = CBC, y = RC), color="black", pch=7, size=2) +
      stat_ellipse(dd, geom="polygon",
                   mapping=aes(x=CBC, y=RC, fill=cluster),
                   alpha=0.25)

    return(p)
  }

  # ------------------------------------------------------------------------------
  # Run the VB GMM with SCD data
  # ------------------------------------------------------------------------------

  #' Runs the VB GMM with selectable k and initialiser
  #'
  #' @param k The number of centroids required.
  #' @param init The initialiser ('kmeans' or 'dbscan').
  #' @param initParams The initialiser parameters if 'dbscan' is selected for init. NULL otherwise.
  #'
  #' @returns A plot of the cluster assignments with centroids on the faithful data.
  #'
  run_gmm <- function(k, init, initParams) {
    data(scd_cohort)
    x <- scd_cohort[,.(CBC,RC)]
    p <- ncol(x)
    prior <- list(
      alpha = 0.001,                              # Dirichlet for q(pi)                             eq 10.57
      beta = 1,                                   # Gaussian in Normal-Wishart for q(mu,lambda)     eq 10.59
      m = matrix(rowMeans(t(x)),ncol=1),          # Gaussian in Normal-Wishart for q(mu,lambda)     eq 10.59
      v = p+1,                                    # Wishart (nu) in Normal-Wishart for q(mu,lambda) eq 10.59
      W = diag(1,p),                              # Wishart (W) in Normal-Wishart for q(mu,lambda)  eq 10.59
      logW = -2*sum(log(diag(chol(diag(1,p)))))
    )

    # With informative prior and stop iterations when ELBO starts to reverse
    gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 5000,
                              init=init, initParams=initParams,
                              stopIfELBOReverse=TRUE, verbose=FALSE,
                              progressbar=TRUE)

    do_init_plot(scd_cohort, gmm_result)
  }

  # ------------------------------------------------------------------------------
  # init = "kmeans" and k=2.                                            Figure 5a
  # ------------------------------------------------------------------------------
  run_gmm(k = 2, init = "kmeans", initParams = NULL)

  # ------------------------------------------------------------------------------
  # init = "kmeans" and k=3.                                            Figure 5b
  # ------------------------------------------------------------------------------
  run_gmm(k = 3, init = "kmeans", initParams = NULL)

  # ------------------------------------------------------------------------------
  # init = "kmeans" and k=10.                                           Figure 5c
  # ------------------------------------------------------------------------------
  # Note this may take a few minutes runtime
  run_gmm(k = 10, init = "kmeans", initParams = NULL)

  # ------------------------------------------------------------------------------
  # init = "dbscan" and k=2.                                            Figure 6
  # ------------------------------------------------------------------------------
  initParams <- c(0.15, 5)
  names(initParams) <- c('eps','minPts')
  run_gmm(k = 2, init = "dbscan", initParams = initParams)
}
