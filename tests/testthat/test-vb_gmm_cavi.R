
test_that("A prior can be omitted ", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  gmm_result <- vb_gmm_cavi(X=x, k=k, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("NA values in data are rejected ", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  x$CBC[3] <- NA
  k <- 2
  p <- ncol(x)
  X <- t(x)

  expect_error(gmm_result <- vb_gmm_cavi(X=x, k=k, delta=1e-6, maxiters = 50, init="kmeans"))
})

test_that("A partial prior can be submitted (alpha) ", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(alpha = 0.001)

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (beta)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(beta = 1)

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (m)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(m = matrix(rowMeans(X),ncol=1))

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (v)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(v = p+1)

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (W)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(W = diag(1,p))

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (logW)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(logW = -2*sum(log(diag(chol(diag(1,p))))))

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})

test_that("A partial prior can be submitted (three hyperparameters out of six)", {
  x <- read.csv(test_path("testdata", "gmm_dat.csv"))
  k <- 2
  p <- ncol(x)
  X <- t(x)

  prior <- list(beta = 1,
                v = p+1,
                logW = -2*sum(log(diag(chol(diag(1,p))))))

  gmm_result <- vb_gmm_cavi(X=x, k=k, prior=prior, delta=1e-6, maxiters = 50, init="kmeans")
  expect_equal(length(gmm_result$z_post), 1000)
})
