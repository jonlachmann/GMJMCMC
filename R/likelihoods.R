# Title     : Log likelihood functions
# Objective : Log likelihood functions with priors to be used as templates or directly in GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-24

#' Log likelihood function for glm regression with parameter priors from BAS package
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user, important to specify the tuning parameters of beta priors and family that BAS uses in glm models
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' glm.logpost.bas(as.integer(rnorm(100) > 0), cbind(1, matrix(rnorm(100))), c(TRUE, TRUE), list(oc = 1))
#' 
#' @importFrom BAS uniform Jeffreys g.prior
#' @importFrom stats poisson Gamma glm.control
#' @export glm.logpost.bas
glm.logpost.bas <- function (y, x, model, complex, params = list(r = exp(-0.5), family = "binomial", prior_beta = Jeffreys(), laplace = FALSE)) {
  if (length(params) == 0)
    params <- list(r = 1 / dim(x)[1], family = "binomial", prior_beta = g.prior(max(dim(x)[1], sum(model) - 1)), laplace = FALSE)
  p <- sum(model) - 1 
  if (p == 0) {
    probinit <- as.numeric(c(1, 0.99))
    model[2] <- T
  } else {
    probinit <- as.numeric(c(1, rep(0.99, p)))
  }
  
  mod <- NULL

  if (params$family == "binomial")
    family_use <- binomial()
  else if (params$family == "poisson")
    family_use <- poisson()
  else
    family_use <- Gamma()
  
  tryCatch({ suppressWarnings({
      mod <- .Call(
        BAS:::C_glm_deterministic,
        y = as.numeric(y),
        X = as.matrix(x[, model]),
        Roffset = as.numeric(rep(0, length(y))),
        Rweights = as.numeric(rep(1, length(y))),
        Rprobinit = probinit,
        Rmodeldim = as.integer(rep(0, ifelse(p == 0,2,1))),
        modelprior = uniform(),
        betaprior = params$prior_beta,
        family = family_use,
        Rcontrol = glm.control(),
        Rlaplace =  as.integer(params$laplace)
      )
    })
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })

  if (length(mod) == 0) {
    return(list(crit = -.Machine$double.xmax + log(params$r * sum(complex$oc)), coefs = rep(0, p + 1)))
  }

  if (p == 0) {
    ret <- mod$logmarg[2] + log(params$r) * sum(complex$oc)
    return(list(crit = ret, coefs = mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(params$r) * sum(complex$oc)
  return(list(crit = ret, coefs = mod$mle[[1]]))
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user, important to specify the tuning parameters of beta priors where the corresponding integers as prior_beta must be provided "g-prior" = 0, "hyper-g" = 1, "EB-local" = 2, "BIC" = 3, "ZS-null" = 4, "ZS-full" = 5, "hyper-g-laplace" = 6, "AIC" = 7, "EB-global" = 2, "hyper-g-n" = 8, "JZS" = 9 and in Gaussian models
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' lm.logpost.bas(rnorm(100), cbind(1,matrix(rnorm(100))), c(TRUE,TRUE), list(oc = 1))
#'
#'
#' @export lm.logpost.bas
lm.logpost.bas <- function (y, x, model, complex, params = list(r = exp(-0.5), beta_prior = list(method = 1))) {
  if (length(params) == 0)
    params <- list(
      r = 1/dim(x)[1],
      beta_prior = list(method = 0, alpha = max(dim(x)[1], sum(model)^2))
    )

  p <- sum(model) - 1
  if (p == 0) {
    probinit <- as.numeric(c(1, 0.99))
    model[2] <- T
  } else {
    probinit <- as.numeric(c(1, rep(0.99, p)))
  }

  mod <- NULL

  tryCatch({
      suppressWarnings({
        mod <- .Call(BAS:::C_deterministic,
                     y = y, X = as.matrix(x[, model]),
                     as.numeric(rep(1, length(y))),
                     probinit,
                     as.integer(rep(0, ifelse(p == 0,2,1))),
                     incint = as.integer(F),
                     alpha = ifelse(length(params$beta_prior$alpha) > 0, as.numeric(params$beta_prior$alpha),NULL),
                     method = as.integer(params$beta_prior$method),
                     modelprior = uniform(),
                     Rpivot = TRUE,
                     Rtol = 1e-7)
      })
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })

  if (length(mod) == 0) {
    return(list(crit = -.Machine$double.xmax + log(params$r * sum(complex$oc)), coefs = rep(0, p + 1)))
  }

  if (p == 0) {
    ret <- mod$logmarg[2] + log(params$r) * sum(complex$oc)
    return(list(crit = ret, coefs = mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(params$r) * sum(complex$oc)
  return(list(crit = ret, coefs = mod$mle[[1]]))
}


#' Log likelihood function for logistic regression with a Jeffreys parameter prior and BIC approximations of the posterior
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
    params <- list(r = 1 / dim(x)[1])
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = binomial())})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(params$r) * sum(complex$oc))) / 2
  return(list(crit = ret, coefs = mod$coefficients))
}

#' Log likelihood function for glm regression with a Jeffreys parameter prior and BIC approximations of the posterior
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user, family must be specified
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' glm.loglik(abs(rnorm(100))+1, matrix(rnorm(100)), TRUE, list(oc = 1))
#'
#'
#' @export glm.loglik
glm.loglik <- function (y, x, model, complex, params = list(r = exp(-0.5), family = "Gamma")) {
  if (length(params) == 0)
    params <- list(r = 1 / dim(x)[1])

  if (params$family == "binomial") {
    fam = binomial()
  } else if (params$family == "poisson") {
    fam = poisson()
  } else {
    fam = Gamma()
  }

  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = fam)})
  ret <- (-(mod$deviance + log(length(y)) * (mod$rank - 1) - 2 * log(params$r) * sum(complex$oc))) / 2
  return(list(crit = ret, coefs = mod$coefficients))
}


#' Log likelihood function for gaussian regression with a Jeffreys prior and BIC approximation of MLIK with both known and unknown variance of the responses
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
  if (sum(model) == 0)
    return(list(crit = -Inf, coefs = numeric()))
  if (length(params) == 0)
    params <- list()
  if (length(params$r) == 0)
    params$r <- 1/dim(x)[1]
  if (length(params$var) == 0)
    params$var <- 1
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())})

  if (params$var == "unknown")
    ret <- (-(mod$aic + (log(length(y))-2) * (mod$rank) - 2 * log(params$r) * (sum(complex$oc)))) / 2
  else
    ret <- (-(mod$deviance/params$var + log(length(y)) * (mod$rank - 1) - 2 * log_prior(params, complex))) / 2

  return(list(crit = ret, coefs = mod$coefficients))
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
#' gaussian.loglik.g(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @export gaussian.loglik.g
gaussian.loglik.g <- function (y, x, model, complex, params = NULL) {
  if (sum(model) == 0)
    return(list(crit = -Inf, coefs = numeric()))
  if (length(params) == 0)
    params <- list()
  if (length(params$r) == 0)
    params$r <- 1/dim(x)[1]
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  # Calculate R-squared
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  Rsquare <- 1 - (RSS / TSS)

  if (length(params$r) == 0 || length(params$g) == 0)
  {
    params$r <- 1/dim(x)[1]
    params$g <- max(mod$rank^2,length(y))
  }

  # logarithm of marginal likelihood
  mloglik <- 0.5*(log(1.0 + params$g) * (dim(x)[1] - mod$rank)  - log(1.0 + params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1))*(mod$rank!=1)

  # logarithm of model prior
   # default value or parameter r
  lp <- log_prior(params, complex)

  return(list(crit = mloglik + lp, coefs = mod$coefficients))
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#'
#' This function computes the marginal likelihood of a Gaussian regression model under different priors.
#'
#' @param y A numeric vector containing the dependent variable.
#' @param x A matrix containing the independent variables, including an intercept column.
#' @param model A logical vector indicating which variables to include in the model.
#' @param complex A list containing complexity measures for the features.
#' @param params A list of parameters for the log likelihood, specifying the tuning parameters of beta priors.
#'
#' @return A list with elements:
#'   \item{crit}{Log marginal likelihood combined with the log prior.}
#'   \item{coefs}{Posterior mode of the coefficients.}
#'
#' @examples
#' gaussian_tcch_log_likelihood(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @importFrom BAS phi1 hypergeometric1F1 hypergeometric2F1
#' @importFrom tolerance F1
#' @export
gaussian_tcch_log_likelihood <- function(y, x, model, complex, params = list(r = exp(-0.5), beta_prior = "intrinsic")) {
  # Fit the linear model using fastglm
  fitted_model <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  log_likelihood <- -(fitted_model$aic  -2 * (fitted_model$rank))/2
  # Compute R-squared manually
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(fitted_model$residuals^2)
  R2_M <- 1 - (RSS / TSS)

  p_M <- fitted_model$rank
  n <- length(y)

  # Switch-like structure to assign hyperparameters based on prior
  hyper <- params$beta_prior$hyper.parameters
  if (params$beta_prior$type == "CH") {
    # CH prior: b and s should be user-specified, with defaults if not provided
    a <- ifelse(!is.null(hyper$a), hyper$a, 1)  # Default to 1 if not specified
    b <- ifelse(!is.null(hyper$b), hyper$b, 2)  # Default to 1 if not specified
    r <- 0
    s <- ifelse(!is.null(hyper$s), hyper$s, 1)  # Default to 1 if not specified
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "hyper-g") {
    a <- 1
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "uniform") {
    a <- 2
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "Jeffreys") {
    a <- 0.0001
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "beta.prime") {
    a <- 1/2
    b <- n - p_M - 1.5
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "benchmark") {
    a <- 0.02
    b <- 0.02 * max(n, p_M^2)
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "TG") {
    a <- 2 * ifelse(!is.null(hyper$a), hyper$a, 1)
    b <- 2
    r <- 0
    s <- 2 * ifelse(!is.null(hyper$s), hyper$s, 1)
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "ZS-adapted") {
    a <- 1
    b <- 2
    r <- 0
    s <- n + 3
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "robust") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- (n + 1) / (p_M + 1)
    k <- 1
  } else if (params$beta_prior$type == "hyper-g-n") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- 1
    k <- 1
  } else if (params$beta_prior$type == "intrinsic") {
    a <- 1
    b <- 1
    r <- 1
    s <- 0
    v <- (n + p_M + 1) / (p_M + 1)
    k <- (n + p_M + 1) / n
  } else if (params$beta_prior$type == "tCCH") {
    a <- hyper$a
    b <- hyper$b
    r <- hyper$r
    s <- hyper$s
    v <- hyper$v
    k <- hyper$k
  } else {
    stop("Unknown prior name: ", params$beta_prior$type)
  }

  if (!is.null(r) & r == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- phi1(b / 2, (n - 1) / 2, (a + b + p_M) / 2, s / (2 * v), min(0.8, R2_M / (v - (v - 1) * R2_M), log = TRUE))

    if (R2_M / (v - (v - 1) * R2_M) > 0.8) {
      warning("Infinite marginal log likelihood! phi1 last argument reduced to 0.8. Use a different prior_beta (Robust, Hyper-g/n, Intrinsic, or g-prior)")
    }

    term3 <- lbeta(a / 2, b / 2)
    term4 <- hypergeometric1F1(b / 2, (a + b) / 2, s / (2 * v), log = TRUE)
    marginal_likelihood <- log_likelihood + (term1) + (term2) - (p_M / 2) * log(v) - ((n - 1) / 2) * log(1 - (1 - 1 / v) * R2_M) - (term3) - (term4)
  } else if (!is.null(s) & s == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- hypergeometric2F1(r, b / 2, (a + b) / 2, 1 - k, log = TRUE)
    term3 <- F1((a + p_M) / 2, (a + b + p_M + 1 - n - 2 * r) / 2, (n - 1) / 2, (a + b + p_M) / 2, 1 - k, 1 - k - (R2_M^2 * k) / ((1 - R2_M) * v))
    marginal_likelihood <- log_likelihood + (a + p_M - 2 * r) / 2 * log(k) + (term1) - (term2) - (term3) - (p_M / 2) * log(v) - log(1 - R2_M) * ((n - 1) / 2) - lbeta(a / 2, b / 2)

  } else {
    stop("Invalid inputs: either r = 0 or s = 0 must be specified.")
  }

  if (length(params$r) == 0) params$r <- 1 / dim(x)[1]  # default value or parameter r

  lp <- log_prior(params, complex)

  return(list(crit = marginal_likelihood + lp, coefs = fitted_model$coefficients))
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


#' Master Log Marginal Likelihood Function
#'
#' This function serves as a unified interface to compute the log marginal likelihood
#' for different regression models and priors by calling specific log likelihood functions.
#'
#' @param y A numeric vector containing the dependent variable.
#' @param x A matrix containing the precalculated features (independent variables).
#' @param model A logical vector indicating which variables to include in the model.
#' @param complex A list of complexity measures for the features.
#' @param params A list of parameters controlling the model family, prior, and tuning parameters.
#'   Key elements include:
#'   - family: "binomial", "poisson", "gamma" (all three referred to as GLM below), or "gaussian" (default: "gaussian")
#'   - prior_beta: Type of prior as a string (default: "g-prior"). Possible values include:
#'     - "beta.prime": Beta-prime prior (GLM/Gaussian, no additional args)
#'     - "CH": Compound Hypergeometric prior (GLM/Gaussian, requires `a`, `b`, optionally `s`)
#'     - "EB-local": Empirical Bayes local prior (GLM/Gaussian, requires `a` for Gaussian)
#'     - "EB-global": Empirical Bayes local prior (Gaussian, requires `a` for Gaussian)
#'     - "g-prior": Zellner's g-prior (GLM/Gaussian, requires `g`)
#'     - "hyper-g": Hyper-g prior (GLM/Gaussian, requires `a`)
#'     - "hyper-g-n": Hyper-g/n prior (GLM/Gaussian, requires `a`)
#'     - "tCCH": Truncated Compound Hypergeometric prior (GLM/Gaussian, requires `a`, `b`, `s`, `rho`, `v`, `k`)
#'     - "intrinsic": Intrinsic prior (GLM/Gaussian, no additional args)
#'     - "TG": Truncated Gamma prior (GLM/Gamma, requires `a`, `s`)
#'     - "Jeffreys": Jeffreys prior (GLM/Gaussian, no additional args)
#'     - "uniform": Uniform prior (GLM/Gaussian, no additional args)
#'     - "benchmark": Benchmark prior (Gaussian/GLM, no additional args)
#'     - "ZS-adapted": Zellner-Siow adapted prior (Gaussian TCCH, no additional args)
#'     - "robust": Robust prior (Gaussian/GLM, no additional args)
#'     - "Jeffreys-BIC": Jeffreys prior with BIC approximation of marginal likelihood (Gaussian/GLM)
#'     - "ZS-null": Zellner-Siow null prior (Gaussian, requires `a`)
#'     - "ZS-full": Zellner-Siow full prior (Gaussian, requires `a`)
#'     - "hyper-g-laplace": Hyper-g Laplace prior (Gaussian, requires `a`)
#'     - "AIC": AIC prior from BAS (Gaussian, requires penalty `a`)
#'     - "BIC": BIC prior from BAS (Gaussian/GLM)
#'     - "JZS": Jeffreys-Zellner-Siow prior (Gaussian, requires `a`)
#'   - r: Model complexity penalty (default: 1/n)
#'   - g: Tuning parameter for g-prior (default: max(n, p^2))
#'   - a, b, s, v, rho, k: Hyperparameters for various priors
#'   - n: Sample size for some priors (default: length(y))
#'   - var: Variance assumption for Gaussian models ("known" or "unknown", default: "unknown")
#'   - laplace: Logical for Laplace approximation in GLM only (default: FALSE)
#'
#' @return A list with elements:
#'   \item{crit}{Log marginal likelihood combined with the log prior.}
#'   \item{coefs}{Posterior mode of the coefficients.}
#'
#' @examples
#' fbms.mlik.master(rnorm(100), matrix(rnorm(100)), TRUE, list(oc = 1), list(family = "gaussian", prior_beta = "g-prior"))
#'
#' @importFrom BAS beta.prime bic.prior CCH EB.local g.prior hyper.g hyper.g.n tCCH intrinsic TG Jeffreys uniform
#' @export
fbms.mlik.master <- function(y, x, model, complex, params = list(family = "gaussian", prior_beta = "g-prior", r = exp(-0.5))) {
  # Extract dimensions
  n <- length(y)
  p <- sum(model) - 1  # Number of predictors excluding intercept

  # Set default parameters if not fully specified
  if (is.null(params$family)) params$family <- "gaussian"
  if (is.null(params$prior_beta)) params$prior_beta <- "g-prior"
  if (is.null(params$g)) params$g <- max(p^2, n)
  if (is.null(params$n)) params$n <- n
  if (is.null(params$r)) params$r <- 1/n

  # Ensure complex has oc if not provided, ignore by default
  if (is.null(complex$oc)) complex$oc <- 0

  # Homogenize and prepare params for nested calls
  params_master <- params
  params_nested <- list(r = params_master$r)

  # Define valid priors for each family
  #glm_only_priors <- c("CCH", "tCCH", "TG")
  glm_and_gaussian_priors <- c("CH", "tCCH", "TG","beta.prime", "EB-local", "g-prior", "hyper-g", "hyper-g-n",
                               "intrinsic", "ZS-adapted", "Jeffreys", "uniform", "benchmark", "robust", "Jeffreys-BIC")
  gaussian_only_priors <- c("ZS-null", "ZS-full", "hyper-g-laplace","BIC", "AIC", "JZS","EB-global")

  #review a bit
  gaussian_not_robust <-  c("CH", "tCCH", "ZS-adapted", "TG","beta.prime", "benchmark","Jeffreys")
  gaussian_robust <- c("g-prior", "hyper-g", "EB-local","BIC", "Jeffreys-BIC", "ZS-null", "ZS-full", "hyper-g-laplace",
                       "AIC",  "hyper-g-n",  "JZS")
  gaussian_tcch <- c("CH", "tCCH", "TG","beta.prime", "intrinsic", "ZS-adapted", "uniform","Jeffreys", "benchmark", "robust")
  gaussian_bas <- c("g-prior", "hyper-g", "EB-local","ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS")

  all_priors <- c(glm_and_gaussian_priors, gaussian_only_priors)
  #browser()
  # Validate prior_beta
  if (!params_master$prior_beta %in% all_priors) {
    stop(sprintf("Prior '%s' is not supported. Supported priors: %s",
                 params_master$prior_beta, paste(all_priors, collapse = ", ")))
  }

  # Decision logic based on family and prior_beta
  if (params_master$family %in% c("binomial", "poisson", "gamma")) {
    if (params_master$prior_beta %in% gaussian_only_priors) {
      stop(sprintf("Prior '%s' is not supported for GLM family '%s'. Supported GLM priors: %s",
                   params_master$prior_beta, params_master$family,
                   paste(c(glm_only_priors, glm_and_gaussian_priors), collapse = ", ")))
    }

    params_nested$family <- params_master$family
    if (is.null(params_master$laplace)) params_nested$laplace <- FALSE else params_nested$laplace <- params_master$laplace

    #if(params_nested$laplace)
    #  print("Laplace approximations will be used")

    if (params_master$prior_beta == "Jeffreys-BIC") {
      if(params_nested$family == "binomial")
        result <- logistic.loglik(y, x, model, complex, params_nested)
      else if(params_nested$family%in% c("poisson", "gamma"))
        result <- glm.loglik(y, x, model, complex, params_nested)

    } else {
      params_nested$prior_beta <- switch(
        params_master$prior_beta,
        "beta.prime" = beta.prime(n = n),
        "CH" = CCH(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    beta = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "EB-local" = EB.local(),
        "g-prior" = g.prior(g = params_master$g),
        "hyper-g" = hyper.g(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "tCCH" = tCCH(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                      beta = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                      s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s,
                      r = if (is.null(params_master$rho)) stop("rho must be provided") else params_master$rho,
                      v = if (is.null(params_master$v)) stop("v must be provided") else params_master$v,
                      theta = if (is.null(params_master$k)) stop("k must be provided") else params_master$k),
        "intrinsic" = intrinsic(n = params_master$n),
        "TG" = TG(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "Jeffreys" = Jeffreys(),
        "uniform" = tCCH(alpha = 2,
                         beta = 2,
                         s = 0,
                         r = 0,
                         v = 1,
                         theta = 1),
        "benchmark" =  tCCH(alpha = 0.02,
                            beta = 0.02*max(n,p^2),
                            s = 0,
                            r = 0,
                            v = 1,
                            theta = 1),
        "ZS-adapted" =  tCCH(alpha = 1,
                            beta = 2,
                            s = n + 3,
                            r = 0,
                            v = 1,
                            theta = 1),
        "TG" = TG(alpha = if (is.null(params_master$a)) stop("a must be provided") else params_master$a),
        "robust" = robust(n = if (is.null(params_master$n)) as.numeric(n) else as.numeric(params_master$n)),
        "hyper-g-n" = hyper.g.n(alpha = if (is.null(params_master$a)) 3 else params_master$a,
                                n = params_master$n),
        "BIC" = bic.prior(n = if (is.null(params_master$n)) n else params_master$n),
        stop("Unrecognized prior_beta for GLM: ", params_master$prior_beta)
      )
      result <- glm.logpost.bas(y, x, model, complex, params_nested)
    }
  } else if (params_master$family == "gaussian") {

    if (params_master$prior_beta %in% gaussian_not_robust) {
      warning(sprintf("Prior '%s' is not reccomended supported for Gaussian family '%s'. Can be unstable for strong signals (R^2 > 0.9). Reccomended priors under Gaussian family: %s",
                   params_master$prior_beta, params_master$family,
                   paste(gaussian_robust, collapse = ", ")))
    }

    params_nested$r <- params_master$r

    if (params_master$prior_beta %in% gaussian_tcch) {

      params_nested$prior_beta <- switch(
        params_master$prior_beta,
        "beta.prime" = list("beta.prime"),
        "CH" = list("CH",a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    b = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "tCCH" = list("tCCH", a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                      b = if (is.null(params_master$b)) stop("b must be provided") else params_master$b,
                      s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s,
                      rho = if (is.null(params_master$rho)) stop("rho must be provided") else params_master$rho,
                      v = if (is.null(params_master$v)) stop("v must be provided") else params_master$v,
                      k = if (is.null(params_master$k)) stop("k must be provided") else params_master$k),
        "intrinsic" = list("intrinsic"),
        "TG" = list("TG",a = if (is.null(params_master$a)) stop("a must be provided") else params_master$a,
                    s = if (is.null(params_master$s)) stop("s must be provided") else params_master$s),
        "Jeffreys" = list("Jeffreys"),
        "ZS-adapted" = list("ZS-adapted"),
        "benchmark" = list("benchmark"),
        "robust" = list("robust"),
        "uniform" = list("uniform"),
        stop("Unrecognized prior_beta for Gaussian GLM: ", params_master$prior_beta)
      )
      result <- gaussian_tcch_log_likelihood(y, x, model, complex, params_nested)

    }else if (params_master$prior_beta == "Jeffreys-BIC") {
      if (is.null(params_master$var)) params_nested$var <- "unknown" else params_nested$var <- params_master$var
      result <- gaussian.loglik(y, x, model, complex, params_nested)
    } else if (params_master$prior_beta %in% gaussian_bas) {

      params_nested$prior_beta  <- switch(
        params_master$prior_beta,
        "g-prior" = 0,
        "hyper-g" = 1,
        "EB-local" = 2,
        "BIC" = 3,
        "ZS-null" = 4,
        "ZS-full" = 5,
        "hyper-g-laplace" = 6,
        "AIC" = 7,
        "EB-global" = 2,
        "hyper-g-n" = 8,
        "JZS" = 9
      )
      if(params_master$prior_beta == "g-prior")
      {
        if (!is.null(params_master$g)) params_nested$g <- params_master$g else stop("g must be provided")
        result <- gaussian.loglik.g(y, x, model, complex, params_nested)
      }
      else{
        if (!is.null(params_master$a)) params_nested$alpha <- params_master$a else params_nested$alpha = -1
        result <- lm.logpost.bas(y, x, model, complex, params_nested)
      }

    } else {
      stop("Unexpected error in prior_beta logic for Gaussian.")
    }
  } else {
    stop("Unsupported family: ", params_master$family, ". Supported families are 'binomial', 'poisson', 'gamma', or 'gaussian'.")
  }
  
  # Ensure consistent return structure
  if (is.null(result$crit) || is.null(result$coefs)) {
    stop("Error in computation: Returned result does not contain 'crit' and 'coefs'.")
  }
  
  return(list(crit = result$crit, coefs = result$coefs))
}