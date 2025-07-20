\dontrun{
  # Use Old Faithful data to show the effect of VB GMM Priors along with choice
  # of initialiser (kmeans vs dbscan), stopping on ELBO reverse parameter and delta
  # threshold
  # ------------------------------------------------------------------------------

  require(ggplot2)

  # gen_path <- tempdir()
  gen_path <- "c:/Users/buckl/Documents/vbPhenoR/man/examples/plots"
  data("faithful")
  X <- faithful
  P <- ncol(X)

  # Run with 4 presumed components for demonstration purposes
  k = 4

  # ------------------------------------------------------------------------------
  # Plotting
  # ------------------------------------------------------------------------------

  #' Plots the GMM components with centroids
  #'
  #' @param i List index to place the plot
  #' @param gmm_result Results from the VB GMM run
  #' @param var_name Variable to hold the GMM hyperparameter name
  #' @param grid Grid element used in the plot file name
  #' @param fig_path Path to the directory where the plots should be stored
  #'
  #' @returns
  #' @importFrom ggplot2 ggplot
  #' @export

  do_plots <- function(i, gmm_result, var_name, grid, fig_path) {
    dd <- as.data.frame(cbind(X, cluster = gmm_result$z_post))
    dd$cluster <- as.factor(dd$cluster)

    # The group means
    # ---------------------------------------------------------------------------
    mu <- as.data.frame( t(gmm_result$q_post$m) )

    # Plot the posterior mixture groups
    # ---------------------------------------------------------------------------
    cols <- c("#1170AA", "#55AD89", "#EF6F6A", "#D3A333", "#5FEFE8", "#11F444")
    p <- ggplot() +
      geom_point(dd, mapping=aes(x=eruptions, y=waiting, color=cluster)) +
      scale_color_discrete(cols, guide = 'none') +
      geom_point(mu, mapping=aes(x = eruptions, y = waiting), color="black",
                 pch=7, size=2) +
      stat_ellipse(dd, geom="polygon",
                   mapping=aes(x=eruptions, y=waiting, fill=cluster),
                   alpha=0.25)

    plots[[i]] <- p
    rm(gmm_result, dd, mu)
    gc()

    grids <- paste((grid[i,]), collapse = "_")
    ggsave(filename=paste0(var_name,"_",grids,".png"), plot=p, path=fig_path,
           width=12, height=12, units="cm", dpi=600, create.dir = TRUE)
  }

  # ------------------------------------------------------------------------------
  # Dirichlet alpha
  # ------------------------------------------------------------------------------
  alpha_grid <- data.frame(c1=c(1,1,183),
                           c2=c(1,92,92),
                           c3=c(1,183,198),
                           c4=c(1,198,50))
  init <- "kmeans"

  z <- vector(mode="list", length=nrow(alpha_grid))
  plots <- vector(mode="list", length=nrow(alpha_grid))

  # Just plot the first 5 grids to save time
  for (i in 1:nrow(alpha_grid)) {
    prior <- list(
      alpha = as.integer(alpha_grid[i,]) # set most of the weight on one component
    )

    gmm_result <- vb_gmm_cavi(X=X, k=k, prior=prior, delta=1e-9, init=init,
                              verbose=FALSE, logDiagnostics=FALSE)
    z[[i]] <- table(gmm_result$z_post)
    do_plots(i, gmm_result, "alpha", alpha_grid, gen_path)
  }

  # ------------------------------------------------------------------------------
  # Normal-Wishart beta for precision proportionality
  # ------------------------------------------------------------------------------
  beta_grid <- data.frame(c1=c(0.1,0.9),
                          c2=c(0.1,0.9),
                          c3=c(0.1,0.9),
                          c4=c(0.1,0.9))
  init <- "kmeans"

  z <- vector(mode="list", length=nrow(beta_grid))
  plots <- vector(mode="list", length=nrow(beta_grid))

  # Just plot the first 5 grids to save time
  for (i in 1:nrow(beta_grid)) {
    prior <- list(
      beta = as.numeric(beta_grid[i,])
    )

    gmm_result <- vb_gmm_cavi(X=X, k=k, prior=prior, delta=1e-9, init=init,
                              verbose=FALSE, logDiagnostics=FALSE)
    z[[i]] <- table(gmm_result$z_post)
    do_plots(i, gmm_result, "beta", beta_grid, gen_path)
  }

  # ------------------------------------------------------------------------------
  # Normal-Wishart W0 (assuming variance matrix) & logW
  # ------------------------------------------------------------------------------

  w_grid <- data.frame(w1=c(0.001,2.001),
                       w2=c(0.001,2.001))
  init <- "kmeans"

  z <- vector(mode="list", length=nrow(w_grid))
  plots <- vector(mode="list", length=nrow(w_grid))

  # Just plot the first 5 grids to save time
  for (i in 1:nrow(w_grid)) {
    w0 = diag(w_grid[i,],P)
    prior <- list(
      W = w0,
      logW = -2*sum(log(diag(chol(w0))))
    )

    gmm_result <- vb_gmm_cavi(X=X, k=k, prior=prior, delta=1e-9, init=init,
                              verbose=FALSE, logDiagnostics=FALSE)
    z[[i]] <- table(gmm_result$z_post)
    do_plots(i, gmm_result, "w", w_grid, gen_path)
  }
}


