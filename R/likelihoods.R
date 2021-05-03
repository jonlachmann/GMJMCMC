# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for logistic regression with a prior p(m)=sum(total_width)
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param data Data to be used for the estimation
#' @param model The model as a logical vector to estimate
#' @param formula The formula to be used for estimation
#' @param complex A list of complexity measures for the features
#'
#' @export logistic.loglik
logistic.loglik <- function (y, x, model, complex, params) {
  r <- 20/223
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial())})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$width)))/2
  return(ret)
}

#' Log likelihood function for logistic regression for alpha calculation
#' This function is just the bare likelihood function
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a[i].
#'
#' @export logistic.loglik.alpha
logistic.loglik.alpha <- function (a, data, mu_func) {
  m <- 1/(1+exp(-eval(parse(text=mu_func))))
  -sum((data[,1] * log(m) + (1-data[,1]) * log(1 - m)))
}

#' Log likelihood function for gaussian regression with a prior p(m)=sum(total_width).
#'
#' @param data Data to be used for the estimation
#' @param model The model as a logical vector to estimate
#' @param formula The formula to be used for estimation
#' @param complex A list of complexity measures for the features
#'
#' @export gaussian.loglik
gaussian.loglik <- function (y, x, model, complex, params) {
  r <- 250/223
  if (sum(model) > 7) return(-1000000)
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=gaussian())})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$width)))/2
  return(ret)
}

#' Log likelihood function for gaussian regression for alpha calculation
#' This function is just the bare likelihood function
#' Note that it only gives a proportional value and is equivalent to least squares
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a[i].
#'
#' @export gaussian.loglik.alpha
gaussian.loglik.alpha <- function (a, data, mu_func) {
  m <- eval(parse(text=mu_func))
  sum((data[,1]-m)^2)
}

#' Log likelihood function for linear regression using Zellners g-prior
#'
#' @param data Data to be used for the estimation
#' @param model The model as a logical vector to estimate
#' @param formula The formula to be used for estimation
#' @param complex A list of complexity measures for the features
#'
#' @export linear.g.prior.loglik
linear.g.prior.loglik <- function (y, x, model, complex, params) {
  out <- lm.fit(as.matrix(x[,model]),y)
  rsquared <- 1-sum(var(out$residuals))/sum(var(y))
  p <- out$rank
  n <- nrow(x)
  logmarglik <- 0.5*(log(1+params$g)*(n-p) - log(1+params$g*(1-rsquared))*(n-1))*(p!=1)
  return(logmarglik)
}