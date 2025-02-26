# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for logistic regression with a prior p(m)=sum(total_width)
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' logistic.loglik(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#' 
#'
#' @export logistic.loglik
logistic.loglik <- function (y, x, model, complex, params = list(r = exp(-0.5))) {
  if (length(params) == 0)
    params <- list(r = 1/dim(x)[1])
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial())})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(params$r) * sum(complex$oc))) / 2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for logistic regression with an approximate Laplace approximations used
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' logistic.loglik.ala(as.integer(rnorm(100) > 0), matrix(rnorm(100)), TRUE, list(oc = 1))
#' 
#'
#' @export logistic.loglik.ala
logistic.loglik.ala <- function (y, x, model, complex, params = list(r = exp(-0.5))) {
  if (length(params) == 0)
    params <- list(r = 1/dim(x)[1])
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial(),maxit = 1)})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) -2 * log(params$r) * sum(complex$oc))) / 2
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log model prior function 
#' @param params list of passed parameters of the likelihood in GMJMCMC
#' @param complex list of complexity measures of the features included into the model 
#' 
#' @return A numeric with the log model prior.
#' 
#' @examples
#' log_prior(params = list(r=2), complex = list(oc = 2))
#' 
#' @export log_prior
log_prior <- function (params, complex) {
  pl <- log(params$r) * (sum(complex$oc))
  return(pl)
}

#' Log likelihood function for logistic regression for alpha calculation
#' This function is just the bare likelihood function
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a\[i\].
#'
#' @return A numeric with the log likelihood.
#'
#' @export logistic.loglik.alpha
logistic.loglik.alpha <- function (a, data, mu_func) {
  m <- 1 / (1 + exp(-eval(parse(text = mu_func))))
  -sum((data[,1] * log(m) + (1 - data[, 1]) * log(1 - m)))
}

#' Log likelihood function for gaussian regression with a prior p(m)=r*sum(total_width).
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' gaussian.loglik(rnorm(100), matrix(rnorm(100)), TRUE, list(oc = 1), NULL)
#'   
#'
#' @export gaussian.loglik
gaussian.loglik <- function (y, x, model, complex, params) {
  if(length(params)==0)
    params <- list()
  if (length(params$r) == 0)
    params$r <- 1/dim(x)[1]
  if(length(params$var) == 0)
    params$var <- 1
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())})
  
  if(params$var == "unknown")
    ret <- (-(mod$aic + (log(length(y))-2) * (mod$rank) - 2 * log(params$r) * (sum(complex$oc)))) / 2
  else
    ret <- (-(mod$deviance/params$var + log(length(y)) * (mod$rank - 1) - 2 * log(params$r) * (sum(complex$oc)))) / 2
  
  return(list(crit=ret, coefs=mod$coefficients))
}

#' Log likelihood function for gaussian regression for alpha calculation
#' This function is just the bare likelihood function
#' Note that it only gives a proportional value and is equivalent to least squares
#'
#' @param a A vector of the alphas to be used
#' @param data The data to be used for calculation
#' @param mu_func The function linking the mean to the covariates,
#' as a string with the alphas as a\[i\].
#'
#' @return A numeric with the log likelihood.
#' @examples 
#'\dontrun{
#'gaussian.loglik.alpha(my_alpha,my_data,my_mu)
#'}
#' @export gaussian.loglik.alpha
gaussian.loglik.alpha <- function (a, data, mu_func) {
  m <- eval(parse(text=mu_func))
  sum((data[,1]-m)^2)
}

#' Log likelihood function for linear regression using Zellners g-prior
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' linear.g.prior.loglik(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @export linear.g.prior.loglik
linear.g.prior.loglik <- function (y, x, model, complex, params = list(g = 4)) {
  out <- lm.fit(as.matrix(x[, model]), y)
  rsquared <- 1 - sum(var(out$residuals)) / sum(var(y))
  p <- out$rank
  n <- nrow(x)
  logmarglik <- 0.5 * (log(1 + params$g) * (n - p) - log(1 + params$g * (1 - rsquared)) * (n - 1)) * (p != 1)
  return(list(crit=logmarglik, coefs=out$coefficients))
}


