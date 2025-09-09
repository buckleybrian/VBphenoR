
test_that("Simple logit call with flat prior works ", {
  X <- read.csv(test_path("testdata", "X.csv"))
  y <- read.csv(test_path("testdata", "y.csv"))
  X$age - as.numeric(X$age)
  X$highrisk <- as.numeric(X$highrisk)
  X$Intercept <- 1
  X <- X[, c(5,1,2,3,4)]

  X <- as.matrix(X)
  y <- as.matrix(y)

  prior <- list(mu=c(1,
                      mean(scd_cohort$age),
                      mean(scd_cohort$highrisk),
                      mean(scd_cohort$CBC),
                      mean(scd_cohort$RC)),
                 Sigma=diag(1,5))

  results <- logit_CAVI(X, y, prior, delta=1e-16, maxiters=100)
  expect_equal(length(results$mu), 5)
  expect_equal(length(results$xi), 1000)
})

test_that("NA values in data are rejected ", {
  X <- read.csv(test_path("testdata", "X.csv"))
  y <- read.csv(test_path("testdata", "y.csv"))
  X$age - as.numeric(X$age)
  X$highrisk <- as.numeric(X$highrisk)
  X$Intercept <- 1
  X <- X[, c(5,1,2,3,4)]

  X$age[3] <- NA

  X <- as.matrix(X)
  y <- as.matrix(y)

  prior <- list(mu=c(1,
                     mean(scd_cohort$age),
                     mean(scd_cohort$highrisk),
                     mean(scd_cohort$CBC),
                     mean(scd_cohort$RC)),
                Sigma=diag(1,5))

  expect_error(results <- logit_CAVI(X, y, prior, delta=1e-16, maxiters=100))
})

