# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for logistic regression with a prior p(m)=sum(total_width)
#'
#' @param data Data to be used for the estimation
#' @param model The model as a logical vector to estimate
#' @param formula The formula to be used for estimation
#' @param complex A list of complexity measures for the features
#'
#' @export logistic.loglik
logistic.loglik <- function (data, model, formula, complex) {
  r <- 20/223
  suppressWarnings({model <- glm(formula = formula, data=data, family = "binomial", maxit=100)})
  ret <- (-(model$deviance -2*log(r)*sum(complex$width)))/2
  return(ret)
}