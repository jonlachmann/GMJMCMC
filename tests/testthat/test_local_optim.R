# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-25

context("Local optimizers")

test_that("Testing Greedy algorithm", {
  greedy_kern_test <- list(probs = c(0,0,0,0,1,0),
                      neigh.size = 1, neigh.min = 1, neigh.max = 2)
  greedy_params_test <- list(steps = 8, kern = greedy_kern_test, tries = 10)
  # Dummy test likelihood function
  loglik.tester <- function (x, y, model, formula, complex, params) {
    return(list(crit = sum(model)))
  }

  # Optimize empty model but dont allow all indices, should set all to true except disallowed
  optmod <- greedy.optim(
    c(F, F, F, F, F, F, F, F, F, F),
    list(fixed = 0),
    loglik.tester,
    c(F, F, T, T, T, T, T, T, T, T),
    NULL,
    greedy_params_test
  )$model

  expect_equal(sum(optmod), 8)
  expect_equal(optmod, c(F, F, T, T, T, T, T, T, T, T))
})
