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
#' glm.logpost.bas(as.integer(rnorm(100) > 0),cbind(1,matrix(rnorm(100))),c(TRUE,TRUE),list(oc = 1))
#' 
#' @importFrom BAS uniform Jeffreys g.prior
#' @importFrom stats poisson Gamma glm.control
#' @export glm.logpost.bas
glm.logpost.bas <- function (y, x, model, complex, params = list(r = exp(-0.5), family = "binomial", betaprior = Jeffreys(), laplace = FALSE)) {
  if (length(params) == 0)
    params <- list(r =  1/dim(x)[1], family = "binomial", betaprior = g.prior(max(dim(x)[1],sum(model)-1)), laplace = FALSE)
  p <- sum(model) - 1 
  if(p==0)
  { 
    probinit <- as.numeric(c(1,0.99))
    model[2] <- T
  }else{
    probinit <- as.numeric(c(1,rep(0.99,p)))
  }
  
  mod<-NULL
  
  tryCatch({
  if(params$family == "binomial")
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y), X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = params$betaprior,
                                                         family = binomial(), 
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(params$laplace))
    })
  else if(params$family == "poisson")
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y),  X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = params$betaprior,
                                                         family = poisson(), 
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(params$laplace))
    })
  else
    suppressWarnings({
      mod <- .Call(BAS:::C_glm_deterministic,
                   y = as.numeric(y), X = as.matrix(x[,model]),
                                  Roffset = as.numeric(rep(0, length(y))),
                                  Rweights = as.numeric(rep(1, length(y))),
                                  Rprobinit = probinit,
                                  Rmodeldim = as.integer(rep(0,ifelse(p==0,2,1))),
                                                         modelprior = uniform(),
                                                         betaprior = params$betaprior,
                                                         family = Gamma(), 
                                                         Rcontrol = glm.control(),
                                                         Rlaplace =  as.integer(params$laplace))
    })
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  if(length(mod)==0) {
    return(list(crit = -.Machine$double.xmax + log(params$r * sum(complex$oc)),coefs = rep(0,p+1)))
  }
  
  if(p == 0)
  {  
    ret <- mod$logmarg[2] + log(params$r) * sum(complex$oc) 
    return(list(crit=ret, coefs=mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(params$r) * sum(complex$oc)
  return(list(crit=ret, coefs=mod$mle[[1]]))
}


#' Log likelihood function for Gaussian regression with parameter priors from BAS package
#' This function is created as an example of how to create an estimator that is used
#' to calculate the marginal likelihood of a model.
#'
#' @param y A vector containing the dependent variable
#' @param x The matrix containing the precalculated features
#' @param model The model to estimate as a logical vector
#' @param complex A list of complexity measures for the features
#' @param params A list of parameters for the log likelihood, supplied by the user, important to specify the tuning parameters of beta priors and in Gaussian models
#'
#' @return A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
#'
#' @examples
#' lm.logpost.bas(rnorm(100), cbind(1,matrix(rnorm(100))), c(TRUE,TRUE), list(oc = 1))
#' 
#'
#' @export lm.logpost.bas
lm.logpost.bas <- function (y, x, model, complex, params = list(r = exp(-0.5),betaprior = "g-prior",alpha = 4)) {
  if (length(params) == 0)
    params <- list(r =  1/dim(x)[1], betaprior = "g-prior",alpha = max(dim(x)[1],sum(model)^2))
  
  method.num <- switch(
    params$betaprior,
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
  
  p <- sum(model) - 1 
  if(p==0)
  { 
    probinit <- as.numeric(c(1,0.99))
    model[2] <- T
  }else{
    probinit <- as.numeric(c(1,rep(0.99,p)))
  }
  
  mod<-NULL
  #browser()
  tryCatch({
      suppressWarnings({
        mod <- .Call(BAS:::C_deterministic,
                     y = y, X = as.matrix(x[,model]),
                     as.numeric(rep(1, length(y))),
                     probinit,
                     as.integer(rep(0,ifelse(p==0,2,1))),
                     incint = as.integer(F),
                     alpha = as.numeric(params$alpha),
                     method = as.integer(method.num),
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
  
  if(length(mod)==0) {
    return(list(crit = -.Machine$double.xmax + log(params$r * sum(complex$oc)),coefs = rep(0,p+1)))
  }
  
  if(p == 0)
  {  
    ret <- mod$logmarg[2] + log(params$r) * sum(complex$oc)  
    return(list(crit=ret, coefs=mod$mle[[2]]))
  }
  ret <- mod$logmarg + log(params$r) * sum(complex$oc) 
  return(list(crit=ret, coefs=mod$mle[[1]]))
}



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
#' linear.g.prior.loglik(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))
#'
#' @importFrom BAS phi1 hypergeometric1F1 hypergeometric2F1
#' @importFrom tolerance F1
#' @export 
gaussian_tcch_log_likelihood <- function(y, x, model, complex, params = list(r = exp(-0.5), prior_beta = "Intrinsic")) {
  
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
  if (params$prior_beta == "CH") {
    # CH prior: b and s should be user-specified, with defaults if not provided
    a <- ifelse(!is.null(params$a),params$a, 1)  # Default to 1 if not specified
    b <- ifelse(!is.null(params$b),params$b, 2)  # Default to 1 if not specified
    r <- 0
    s <- ifelse(!is.null(params$s), params$s, 1)  # Default to 1 if not specified
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "Hyper-g") {
    a <- 1
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "Uniform") {
    a <- 2
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "Jeffreys") {
    a <- 0.0001
    b <- 2
    r <- 0
    s <- 0
    v <- 1
    k <- 1
  } else if (params$prior_beta == "Beta-prime") {
    a <- 1/2
    b <- n - p_M - 1.5
    r <- 0
    s <- 0
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "Benchmark") {
    a <- 0.02
    b <- 0.02 * max(n, p_M^2)
    r <- 0
    s <- 0
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "TruncGamma") {
    
    a <- 2 * ifelse(!is.null(params$at),params$at, 1)
    b <- 2
    r <- 0
    s <- 2 * ifelse(!is.null(params$st),params$st, 1)
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "ZS adapted") {
    a <- 1
    b <- 2
    r <- 0
    s <- n + 3
    v <- 1
    k <- 1
  } else if (params$prior_beta == "Robust") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- (n + 1) / (p_M + 1)
    k <- 1
    
  } else if (params$prior_beta == "Hyper-g/n") {
    a <- 1
    b <- 2
    r <- 1.5
    s <- 0
    v <- 1
    k <- 1
    
  } else if (params$prior_beta == "Intrinsic") {
    a <- 1
    b <- 1
    r <- 1
    s <- 0
    v <- (n + p_M + 1) / (p_M + 1)
    k <- (n + p_M + 1) / n
    
  } else {
    stop("Unknown prior name: ", params$prior_beta)
  }
  
  #
  if (!is.null(r) & r == 0) {
    #browser()
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- phi1(b / 2, (n - 1) / 2, (a + b + p_M) / 2, s / (2 * v), min(0.8,R2_M/(v - (v - 1) * R2_M),log = T))
    
    if(R2_M/(v - (v - 1) * R2_M)>0.8)
    {
      warning("Infinite marginal log likelihood! phi1 last argument reduced to 0.8. Use a different prior_beta (Robust, Hyper-g/n, Intrinsic, or g-prior)")
    }
    
    term3 <- lbeta(a / 2, b / 2) 
    term4 <- hypergeometric1F1(b / 2, (a + b) / 2, s / (2 * v),log = T)
    marginal_likelihood <- log_likelihood + (term1) + (term2) - (p_M / 2) * log(v) - ((n - 1) / 2)*log(1 - (1 - 1 / v) * R2_M) - (term3) - (term4)
  } else if (!is.null(s) & s == 0) {
    term1 <- lbeta((a + p_M) / 2, b / 2)
    term2 <- hypergeometric2F1(r, b / 2, (a + b) / 2, 1 - k,log = T)
    term3 <- F1((a + p_M) / 2, (a + b + p_M + 1 - n - 2 * r) / 2, (n - 1) / 2, (a + b + p_M) / 2, 1 - k, 1 - k - (R2_M^2 * k) / ((1 - R2_M) * v))
    marginal_likelihood <- log_likelihood + (a+p_M-2*r)/2*log(k) + (term1) - (term2) - (term3) - (p_M / 2) * log(v) - log(1 - R2_M) * ((n - 1) / 2) - lbeta(a / 2, b / 2)
    
  } else {
    stop("Invalid inputs: either r = 0 or s = 0 must be specified.")
  }
  
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  
  lp <- log_prior(params, complex)
  
  return(list(crit = marginal_likelihood + lp, coefs = fitted_model$coefficients))
}
