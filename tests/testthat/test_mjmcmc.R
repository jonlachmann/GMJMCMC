# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-25

context("MJMCMC")

test_that("Testing MJMCMC algorithm", {
  # Dummy test likelihood function
  loglik.tester <- function (y, x, model, complex, params) {
    lmmod <- lm.fit(x[, model, drop = FALSE], y)
    n <- length(y)
    k <- length(lmmod$coefficients)
    rss <- sum(resid(lmmod)^2)
    aic <- n * log(rss / n) + 2 * k
    return(list(crit = aic, coefs = lmmod$coefficients))
  }

  data <- matrix(rnorm(600), 100)
  resm <- mjmcmc(data, loglik.tester)
  summary(resm)
  plot(resm)
  predm <- predict(resm, cbind(1, data[, -1, drop = FALSE]))

  resg <- gmjmcmc(data, loglik.tester, NULL, c("p0", "exp.dbl"))
  summary(resg)
  plot(resg)
  prediction <- predict(resg, cbind(1, data[, -1, drop = FALSE]))

  respm <- mjmcmc.parallel(2, 2, data, loglik.tester)
  summary(respm)
  plot(respm)
  pred_pm <- predict(respm, cbind(1, data[, -1, drop = FALSE]))

  respg <- gmjmcmc.parallel(2, 2, NULL, data, loglik.tester, NULL, c("p0", "exp.dbl"))
  summary(respg)
  plot(respg)
  pred_pg <- predict(respg, cbind(1, data[, -1, drop = FALSE]))
})
