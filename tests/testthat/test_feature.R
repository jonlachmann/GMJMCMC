# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-11

context("Feature class")

test_that("Simple features are created", {
  transforms <- list("sin", "cos")
  set.transforms(transforms)
  feature1 <- create.feature(1, 1, 1) # sin(x1)
  feature2 <- create.feature(1, list(1,2), 1) # sin(1*x1+1*x2)
  feature3 <- create.feature(0, list(feature1,feature2), 1) # (sin(x1)*sin(1*x1+1*x2))
  expect_equal(print(feature1), "sin(x1)")
  expect_equal(print(feature2), "sin(1*x1+1*x2)")
  expect_equal(print(feature3), "(sin(x1)*sin(1*x1+1*x2))")
})

test_that("Projection feature is created", {
  transforms <- list("sin", "cos")
  set.transforms(transforms)
  feature1 <- create.feature(1, 1, 1) # sin(x1)
  feature2 <- create.feature(1, list(1,2), 1) # sin(x1+x2)
  feature3 <- create.feature(2, list(feature1,feature2,3), c(1, 1), c(0.2,-1.5,7,9)) # cos(0.2+-1.5*sin(x1)+7*sin(1*x1+1*x2)+9*x3)
  expect_equal(print(feature3), "cos(0.2+-1.5*sin(x1)+7*sin(1*x1+1*x2)+9*x3)")
})