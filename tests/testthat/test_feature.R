# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-11

context("Feature class")
library(GMJMCMC)

test_that("Feature is created", {
  transforms <- list("sin", "cos")
  expect_equal(print(GMJMCMC::create.feature(1, 1), transforms), "sin(x1)")
})