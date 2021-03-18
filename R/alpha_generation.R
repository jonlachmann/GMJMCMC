# Title     : Alpha generators
# Objective : Generation of alpha parameters for features that are kept static throughout a population
# Created by: jonlachmann
# Created on: 2021-03-16

#' Alpha generator using strategy 1 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#' @param transforms The transforms used
alpha_1 <- function (feature, transforms) {
  return(feature)
}

#' Alpha generator using strategy 1 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#' @param transforms The transforms used
alpha_2 <- function (feature, transforms) {
  return(feature)
}

#' Alpha generator using strategy 3 as per Hubin et. al.
#'
#' @param feature The feature to generate alphas for
#' @param transforms The transforms used
#' @param link The link function to use
#' @param loglik log likelihood function to use
alpha_3 <- function (feature, transforms, link, loglik) {
  modfun <- model.function(T, list(feature), transforms, link)

  # Set initial range for Simulated Annealing
  range <- 10
  done <- FALSE
  while(!done) {
    # Run simulated annealing on current range
    sares <- GenSA(rep(0,modfun$count), loglik,
                      rep(-range/2,modfun$count), rep(range/2,modfun$count),
                      control=list(max.call=1e4), modfun$formula)
    # Check if any estimate is on the edge of the range, if so, extend the range and run again
    if (sum((sares$par==(-range/2))+(sares$par==(range/2))) != 0) range <- range*2
    else done <- TRUE
  }
  feature <- update.alphas(feature, sares$par)
  return(feature)
}