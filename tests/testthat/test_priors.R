# Title     : General tests
# Objective : Test the code
# Created by: jonlachmann
# Created on: 2021-02-25

context("Priors")

test_that("Test various priors through the fbms function", {
  set.seed(123)
  x <- matrix(rnorm(300), 100)
  y <- rnorm(100, 0, 0.5) + rowSums(x[, 1:2])

  validate.model <- function (model, x, y) {
    expect_true(all(c(model$marg.probs[1:2] > 0.9, model$marg.probs[3] < 0.9)))
    summary <- summary(model, labels = c("a", "b", "c"), tol = -1)
    expect_true(all(c(summary$marg.probs[1:2] > 0.9, summary$marg.probs[3] < 0.9)))
    plot(model)
    pred <- predict(model, x)
    # Handle paralell runs
    if (!is.null(pred$aggr)) {
      pred <- pred$aggr
    }
    rmse <- sqrt(mean((pred$mean - y)^2))
    expect_true(rmse < 0.6)
    best_model <- get.best.model(model)
    mpm_model <- get.mpm.model(model, y, x)
  }

  # No intercept
  data <- as.data.frame(cbind(y, x))
  colnames(data) <- c("y", "a", "b", "c")
  family <- "gaussian"
  beta_prior <- list(type = "g-prior", g = 5, a = 3, b = 1, s = 1, rho = 0, v = 1, k = 1)

  gaussian_priors <- c(
    "CH", "tCCH", "TG","beta.prime", "intrinsic", "ZS-adapted", "uniform","Jeffreys", "benchmark", "robust",
    "g-prior", "hyper-g", "EB-local", "ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS",
    "Jeffreys-BIC", "g-prior"
  )

  for (prior in gaussian_priors) {
    beta_prior$type <- prior
    mod1 <- fbms(y ~ ., family = family, beta_prior = beta_prior, method = "mjmcmc", data = data, verbose = FALSE)
    if (!(prior %in% c("hyper-g", "EB-local", "uniform", "ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS"))) {
      validate.model(mod1, x, y)
    }
  }
})

