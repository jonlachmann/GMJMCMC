#' Fit a BGNLM  model using Genetically Modified Mode Jumping Markov Chain Monte Carlo (MCMC) sampling.
#' Or Fit a BGLM model using Modified Mode Jumping Markov Chain Monte Carlo (MCMC) sampling.
#' 
#' This function fits a model using the relevant MCMC sampling. The user can specify the formula,
#' family, data, transforms, and other parameters to customize the model.
#'
#' @param formula A formula object specifying the model structure. Default is NULL.
#' @param family The distribution family of the response variable. Currently supports "gaussian", "binomial" and  "custom". Default is "gaussian".
#' @param loglik.pi Custom function to compute the logarithm of the posterior mode based on logarithm of marginal likelihood and logarithm of prior functions (needs specification only used if family = "custom")
#' @param data A data frame containing the variables in the model. If NULL, the variables are taken from the environment of the formula. Default is NULL.
#' @param method Which fitting algorithm should be used, currently implemented options include "gmjmcmc", "gmjmcmc.parallel", "mjmcmc" and "mjmcmc.parallel" with "mjmcmc" being the default and 'mjmcmc' means that only linear models will be estimated
#' @param verbose If TRUE, print detailed progress information during the fitting process. Default is TRUE.
#' @param ... Additional parameters to be passed to the underlying method.
#'
#' @return An object containing the results of the fitted model and MCMC sampling.
#'
#' @examples
#' # Fit a Gaussian multivariate time series model
#' fbms_result <- fbms(
#'  X1 ~ .,
#'  family = "gaussian",
#'  method = "gmjmcmc.parallel",
#'  data = data.frame(matrix(rnorm(600), 100)),
#'  transforms = c("sin","cos"),
#'  P = 10,
#'  runs = 1,
#'  cores = 1
#' )
#' summary(fbms_result)
#' plot(fbms_result)
#' 
#'
#' @seealso \code{\link{mjmcmc}}, \code{\link{gmjmcmc}}, \code{\link{gmjmcmc.parallel}}
#' @export
fbms <- function(formula = NULL, family = "gaussian", data = NULL, 
                 loglik.pi = gaussian.loglik,
                 method = "mjmcmc", verbose = TRUE, ...) {
  if (family == "gaussian")
    loglik.pi <- gaussian.loglik
  else if (family == "binomial")
    loglik.pi <- logistic.loglik
  else if (family == "custom")
    loglik.pi <- loglik.pi
  if (!is.null(formula)) {
    if (missing(data)) {
      data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "any")
    X <- model.matrix(formula, data = data)[, -1]
    df <- data.frame(Y, X)
  } else {
    df <- data
  }

  if (method == "mjmcmc.parallel")
    res <- mjmcmc.parallel(df, loglik.pi, verbose = verbose, ...)
  else if (method == "mjmcmc")
    res <- mjmcmc(df, loglik.pi, verbose = verbose, ...)
  else if (method == "gmjmcmc.parallel") {
    res <- gmjmcmc.parallel(data = df, loglik.pi = loglik.pi, verbose = verbose,...)
  }

  else if (method == "gmjmcmc")
    res <- gmjmcmc(df, loglik.pi, verbose = verbose, ...)
  else
    stop("Error: Method must be one of gmjmcmc, gmjmcmc.parallel,mjmcmc or mjmcmc.parallel!")
  
  return(res)
}