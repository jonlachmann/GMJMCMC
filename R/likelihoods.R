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
glm.logpost.bas <- function (y, x, model, complex, params = list(r = NULL, family = "binomial", beta_prior = Jeffreys(), laplace = FALSE)) {
  if (length(params) == 0)
    params <- list(r = 1 / dim(x)[1], family = "binomial", beta_prior = g.prior(max(dim(x)[1], sum(model) - 1)), laplace = FALSE)
  else if(length(params$r) == 0)
    params$r = 1 / dim(x)[1]
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
        betaprior = params$beta_prior,
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
      r = 1 / dim(x)[1],
      beta_prior = list(method = 0, alpha = max(dim(x)[1], sum(model)^2))
    ) else if(length(params$r) == 0) params$r = 1 / dim(x)[1]

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
                     alpha = ifelse(length(params$beta_prior$a) > 0, as.numeric(params$beta_prior$a),NULL),
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
  #browser()
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
  else if(length(params$r) == 0)
    params$r = 1 / dim(x)[1]
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
  else if(length(params$r) == 0)
    params$r = 1 / dim(x)[1]

  if (params$family == "binomial") {
    fam = binomial()
  } else if (params$family == "poisson") {
    fam = poisson()
  } else {
    fam = Gamma()
  }

  #browser()
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
  if (length(params$beta_prior$var) == 0)
    params$beta_prior$var <- 1
  suppressWarnings({mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())})

  if (params$beta_prior$var == "unknown")
    ret <- (-(mod$aic + (log(length(y))-2) * (mod$rank) - 2 * log(params$r) * (sum(complex$oc)))) / 2
  else
    ret <- (-(mod$deviance/params$beta_prior$var + log(length(y)) * (mod$rank - 1) - 2 * log_prior(params, complex))) / 2

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

  if (length(params$r) == 0 || length(params$g) == 0) {
    params$r <- 1 / dim(x)[1]
    if (!is.null(params$beta_prior$g))
      params$g <- params$beta_prior$g
    else
      params$g <- max(mod$rank^2, length(y))
  }

  # logarithm of marginal likelihood
  mloglik <- 0.5 * (log(1.0 + params$g) * (dim(x)[1] - mod$rank) - log(1.0 + params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1)) * (mod$rank != 1)

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
  hyper <- params$beta_prior
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
    r <- hyper$rho
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
fbms.mlik.master <- function(y, x, model, complex, params = list(family = "gaussian", beta_prior = list(type = "g-prior"), r = exp(-0.5))) {
  # Extract dimensions
  n <- length(y)
  p <- sum(model) - 1  # Number of predictors excluding intercept
  params_use <- list()

  if(params$family == "gaussian")
    params_use$beta_prior <- gen.mlpost.params.lm(params$beta_prior$type, params$beta_prior, p+1, n)
  else
  {
    params_use$beta_prior <- gen.mlpost.params.glm(params$beta_prior$type, params$beta_prior, p+1, n)
    params_use$family <- params$family
  }
  
  loglik.pi <- select.mlpost.fun(params$beta_prior$type, params$family)
  
  result <- loglik.pi(y,x,model,complex,params_use)

  
  return(list(crit = result$crit, coefs = result$coefs))
}