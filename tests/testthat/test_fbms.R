# Title     : General tests
# Objective : Test the code
# Created by: jonlachmann
# Created on: 2021-02-25

context("MJMCMC")

test_that("Test (G)MJMCMC through the fbms function", {
  set.seed(123)
  x <- matrix(rnorm(300), 100)
  y <- rnorm(100, 0, 0.5) + rowSums(x[, 1:2])
  y_shift <- y + 10
  y_sin <- rnorm(100, 0, 0.05) + sin(x[, 1]) * 3 + sin(x[, 2]) * 2
  y_sin_shift <- y_sin + 10

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

  validate.gmodel <- function (model, x, y) {
    suppressMessages(summary <- summary(model, labels = c("a", "b", "c"), tol = -1))
    expect_true(all(c(summary$marg.probs[1:2] > 0.9, summary$marg.probs[-(1:2)] < 0.9)))
    expect_true(all(summary$feats.strings[1:2] %in% c("sin(a)", "sin(b)")))
    plot(model)
    pred <- predict(model, x)
    # Handle paralell runs
    if (!is.null(pred$aggr)) {
      pred <- pred$aggr
    }
    rmse <- sqrt(mean((pred$mean - y)^2))
    expect_true(rmse < 0.2)
    best_model <- get.best.model(model)
    mpm_model <- get.mpm.model(model, y, x)
  }

  params <- gen.params.gmjmcmc(list(x = x, y = y_sin))
  probs <- gen.probs.gmjmcmc("sin")
  probs$gen <- c(0, 1, 0, 0)
  params$feat$D <- 1
  params$feat$L <- 2

  # No intercept
  data <- as.data.frame(cbind(y, x))
  colnames(data) <- c("y", "a", "b", "c")
  mod1 <- fbms(y ~ . - 1, family = "gaussian", method = "mjmcmc", data = data, verbose = FALSE)
  mod1p <- fbms(y ~ . - 1, family = "gaussian", method = "mjmcmc.parallel", data = data, verbose = FALSE)
  validate.model(mod1, x, y)
  validate.model(mod1p, x, y)

  data$y <- y_sin
  gmod1 <- fbms(y ~ . - 1, family = "gaussian", method = "gmjmcmc", data = data, transforms = "sin", params = params, probs = probs, P = 20, verbose = FALSE)
  gmod1p <- fbms(y ~ . - 1, family = "gaussian", method = "gmjmcmc.parallel", data = data, transforms = "sin", params = params, probs = probs, P = 20, verbose = FALSE)
  validate.gmodel(gmod1, x, y_sin)
  validate.gmodel(gmod1p, x, y_sin)

  # Model defined intercept
  data$y <- y_shift
  mod2 <- fbms(y ~ ., family = "gaussian", method = "mjmcmc", data = data, verbose = FALSE)
  mod2p <- fbms(y ~ ., family = "gaussian", method = "mjmcmc.parallel", data = data, verbose = FALSE)
  validate.model(mod2, x, y_shift)
  validate.model(mod2p, x, y_shift)

  data$y <- y_sin_shift
  gmod2 <- fbms(y ~ ., family = "gaussian", method = "gmjmcmc", data = data, transforms = "sin", params = params, probs = probs, P = 20, verbose = FALSE)
  gmod2p <- fbms(y ~ ., family = "gaussian", method = "gmjmcmc.parallel", data = data, transforms = "sin", params = params, probs = probs, P = 20, verbose = FALSE)
  validate.gmodel(gmod2, x, y_sin_shift)
  validate.gmodel(gmod2p, x, y_sin_shift)

  # User defined intercept
  data <- cbind(data[, 1], 1, data[, -1])
  colnames(data) <- c("y", "const", "a", "b", "c")
  data$y <- y_shift
  mod3 <- fbms(y ~ . - 1, family = "gaussian", method = "mjmcmc", data = data, fixed = 1, verbose = FALSE)
  mod3p <- fbms(y ~ . - 1, family = "gaussian", method = "mjmcmc.parallel", data = data, fixed = 1, verbose = FALSE)
  validate.model(mod3, cbind(1, x), y_shift)
  validate.model(mod3p, cbind(1, x), y_shift)

  data$y <- y_sin_shift
  gmod3 <- fbms(y ~ . - 1, family = "gaussian", method = "gmjmcmc", data = data, transforms = "sin", params = params, probs = probs, P = 20, fixed = 1, verbose = FALSE)
  gmod3p <- fbms(y ~ . - 1, family = "gaussian", method = "gmjmcmc.parallel", data = data, transforms = "sin", params = params, probs = probs, P = 20, fixed = 1, verbose = FALSE)
  validate.gmodel(gmod3, cbind(1, x), y_sin_shift)
  validate.gmodel(gmod3p, cbind(1, x), y_sin_shift)
})

