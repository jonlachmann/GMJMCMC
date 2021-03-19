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
logistic.loglik <- function (y, x, model, complex) {
  r <- 20/223
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial(), method=2)})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$width)))/2
  return(ret)
}

logistic.loglik.alpha <- function (a, mu_func) {
  m <- eval(parse(text=mu_func))
  -sum((y2 * log(m) + (1-y2) * log(1 - m)))
}