
require(ggplot2)
require(dbscan)
require(data.table)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

plot_multi_histogram <- function(df, feature, label_column, density=TRUE) {
  if(density==TRUE) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)),
                          fill=eval(parse(text=label_column)))) +
      geom_histogram(alpha=0.7,
                     position="identity",
                     aes(y = after_stat(density)),
                     color="black", bins=60) +
      geom_density(alpha=0.7) +
      geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", linewidth=1) +
      labs(x=feature, y = "Density")
  } else {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)),
                          fill=eval(parse(text=label_column)))) +
      geom_histogram(alpha=0.7,
                     position="identity",
                     #aes(y = after_stat(density)),
                     color="black", bins=60) +
      geom_density(alpha=0.7) +
      geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", linewidth=1) +
      labs(x=feature, y = "Density")
  }
  plt + guides(fill=guide_legend(title=label_column))
}

# -----------------------------------------------------------------------------
# Synth-MD Sickle-Cell Disease
# -----------------------------------------------------------------------------

data("scd_cohort")

# ----------------------------------------------------------------------------
# Data transformation to match VB patient phenotyping model variables
# ----------------------------------------------------------------------------

# Introduce missing biomarkers
missingCBC <- scd_cohort[, .SD[sample(.N, floor(.5 * .N))], by = CBC]$idx
missingRC <- scd_cohort[, .SD[sample(.N, floor(.5 * .N))], by = RC]$idx

scd_cohort_miss <- copy(scd_cohort)
scd_cohort_miss[idx %in% missingCBC,CBC:=NA]
scd_cohort_miss[idx %in% missingRC,RC:=NA]
scd_cohort_miss[is.na(RC), SCD:=FALSE]
scd_cohort_miss[is.na(CBC), SCD:=FALSE]

# -----------------------------------------------------------------------------
# Test model
# -----------------------------------------------------------------------------

x <- scd_cohort[,.(CBC,RC)]
k <- 3

# For DBSCAN
initParams <- c(0.15, 5)
names(initParams) <- c('eps','minPts')

# For HDBSCAN
#names(initParams) <- c('minPts')

p <- ncol(x)
X <- t(x)
prior <- list(
  alpha = 0.001,                              # Dirichlet for q(pi)                             eq 10.57
  beta = 1,                                   # Gaussian in Normal-Wishart for q(mu,lambda)     eq 10.59
  m = matrix(rowMeans(X),ncol=1),             # Gaussian in Normal-Wishart for q(mu,lambda)     eq 10.59
  v = p+1,                                    # Wishart (nu) in Normal-Wishart for q(mu,lambda) eq 10.59
  W = diag(1,p),                              # Wishart (W) in Normal-Wishart for q(mu,lambda)  eq 10.59
  logW = -2*sum(log(diag(chol(diag(1,p)))))
)

# With informative prior
gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-7, maxiters = 5000,
                          init="dbscan", initParams=initParams, verbose=FALSE)

# With informative prior and stop iterations when ELBO starts to reverse
gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 5000,
                          init="dbscan", initParams=initParams,
                          stopIfELBOReverse=TRUE, verbose=FALSE)

table(gmm_result$z_post)

eb <- data.frame(lb=gmm_result$elbo)
ggplot(eb, aes(x=1:nrow(eb), y=lb)) + geom_line() + xlab("Iteration") + ylab("ELBO")

dd <- as.data.frame(cbind(scd_cohort, cluster=gmm_result$z_post))
dd$cluster <- as.factor(dd$cluster)
dd$race <- as.factor(dd$race)
mu <- as.data.frame( t(gmm_result$q_post$m) )

plot_multi_histogram(dd, 'RC', 'cluster')
plot_multi_histogram(dd, 'CBC', 'cluster')

cols <- c("#1170AA", "#55AD89", "#EF6F6A", "#D3A333", "#5FEFE8", "#11F444")
ggplot() +
  geom_point(dd, mapping=aes(x=CBC, y=RC, color=cluster)) +
  scale_color_discrete(cols, guide = 'none') +
  geom_point(mu, mapping=aes(x = CBC, y = RC), color="black", pch=7, size=2) +
  stat_ellipse(dd, geom="polygon",
               mapping=aes(x=CBC, y=RC, fill=cluster),
               alpha=0.25)

# -----------------------------------------------------------------------------
# Try to merge clusters using Hierarchical DBSCAN
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# DBSCAN comparison
# -----------------------------------------------------------------------------
dbscan_res <- dbscan(x, eps = 0.15, minPts = 5, borderPoints = TRUE)

dd <- as.data.frame(cbind(scd_cohort, cluster=dbscan_res$cluster))
dd$cluster <- as.factor(dd$cluster)
dd$race <- as.factor(dd$race)

#plot_multi_histogram(dd, 'RC', 'cluster')
#plot_multi_histogram(dd, 'CBC', 'cluster')

ggplot() +
  geom_point(dd, mapping=aes(x=CBC, y=RC, color=cluster)) +
  scale_color_discrete(cols, guide = 'none') +
  stat_ellipse(dd, geom="polygon",
               mapping=aes(x=CBC, y=RC, fill=cluster),
               alpha=0.25)

# -----------------------------------------------------------------------------
# HDBSCAN comparison
# -----------------------------------------------------------------------------
hdbscan_res <- hdbscan(x, minPts = 289, cluster_selection_epsilon = 0.15)

dd <- as.data.frame(cbind(scd_cohort, cluster=dbscan_res$cluster))
dd$cluster <- as.factor(dd$cluster)
dd$race <- as.factor(dd$race)

#plot_multi_histogram(dd, 'RC', 'cluster')
#plot_multi_histogram(dd, 'CBC', 'cluster')

ggplot() +
  geom_point(dd, mapping=aes(x=CBC, y=RC, color=cluster)) +
  scale_color_discrete(cols, guide = 'none') +
  stat_ellipse(dd, geom="polygon",
               mapping=aes(x=CBC, y=RC, fill=cluster),
               alpha=0.25)



