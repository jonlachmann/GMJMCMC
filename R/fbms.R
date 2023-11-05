#' Fit a BGNLM  model using Genetically Modified Mode Jumping Markov Chain Monte Carlo (MCMC) sampling.
#' Or Fit a BGLM model using Modified Mode Jumping Markov Chain Monte Carlo (MCMC) sampling.
#' 
#' This function fits a model using the relevant MCMC sampling. The user can specify the formula,
#' family, data, transforms, and other parameters to customize the model.
#'
#' @param formula A formula object specifying the model structure. Default is NULL.
#' @param family The distribution family of the response variable. Currently supports "gaussian" and "binomial". Default is "gaussian".
#' @param data A data frame containing the variables in the model. If NULL, the variables are taken from the environment of the formula. Default is NULL.
#' @param transforms A list of transformations for BGNLM model. Default is NULL.
#' @param loglik.pi The log-likelihood function for estimating the marginal likelihood and posterior modes (only used if family = "custom")
#' @param loglik.alpha The log-likelihood function for the alpha parameter in the model. Default is gaussian.loglik.alpha.
#' @param P The number of GMJMCMC generations. Default is 10.
#' @param runs The number of parallel chains in case of parallel processing. Default is 2.
#' @param cores The number of CPU cores to use for parallel processing. Default is 2.
#' @param verbose If TRUE, print detailed progress information during the fitting process. Default is FALSE.
#' @param ... Additional parameters to be passed to the underlying MCMC fitting functions.
#'
#' @return An object containing the results of the fitted model and MCMC sampling.
#'
#' @examples
#' # Fit a Gaussian multivariate time series model
#' fbms_result <- fbms(
#'  X1 ~ .,
#'  family = "gaussian",
#'  data = data.frame(matrix(rnorm(600), 100)),
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
fbms <- function(formula = NULL, family = "gaussian", data = NULL, transforms = NULL,
                 loglik.pi = gaussian.loglik,
                 loglik.alpha = gaussian.loglik.alpha,
                 P = 10, runs = 10, cores = 1,verbose = FALSE, ...) {
  if (family == "gaussian")
    loglik.pi <- gaussian.loglik
  else if(family == "binomial")
    loglik.pi <- logistic.loglik
  else if(family == "custom")
    loglik.pi <- loglik.pi
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "any")
  X <- model.matrix(formula, data = data)[, -1]
  df <- data.frame(Y, X)
  if (is.null(transforms)) {
    if (cores > 1)
      res <- mjmcmc.parallel(df, loglik.pi, verbose = FALSE, ...)
    else
      res <- mjmcmc(df, loglik.pi, verbose = FALSE, ...)
  } else {
    if (cores > 1)
      res <- gmjmcmc.parallel(runs, cores, data = df, loglik.pi = loglik.pi,
                              loglik.alpha = gaussian.loglik.alpha, transforms = transforms,
                              P = P, ...)
    else
      res <- gmjmcmc(df, loglik.pi, gaussian.loglik.alpha, transforms, 
                     verbose = FALSE, P, ...)
  }
  return(res)
}