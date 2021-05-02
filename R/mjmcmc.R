# Title     : MJMCMC
# Objective : Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-04-27

#' Main algorithm for MJMCMC
#'
#' @param data A matrix containing the data to use in the algorithm,
#' first column should be the dependent variable, second should be the intercept
#' and the rest of the columns should be the independent variables.
#' @param loglik.pi The (log) density to explore
#' @param N The number of iterations to run for
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#'
#' @export mjmcmc
mjmcmc <- function (data, loglik.pi, N, probs, params, sub=F) {
  # Verify that data is well-formed
  data <- check.data(data)
  # Acceptance probability
  accept <- 0

  # Create a population of just the covariates
  S <- gen.covariates(ncol(data)-2)
  complex <- complex.features(S)

  # Initialize first model
  model.cur <- as.logical(rbinom(n = length(S), size = 1, prob = 0.5))
  model.cur <- list(prob=0, model=model.cur, crit=loglik.pre(loglik.pi, model.cur, complex, data, params$loglik), alpha=0)
  best.crit <- model.cur$crit # Set first best criteria value

  # A list of models that have been visited
  models <- vector("list", N)
  # Initialize a vector to contain local opt visited models
  lo.models <- vector("list", 0)
  # If we are running a subsampling strategy, keep a list of best mliks for all models
  if (sub) mliks <- vector("list", (2^(length(S))))
  else mliks <- NULL


  print("MJMCMC begin.")
  progress <- 0
  for (i in 1:N) {
    if (N > 40 && i %% floor(N/40) == 0) progress <- print.progressbar(progress, 40)
    proposal <- mjmcmc.prop(data, loglik.pi, model.cur, S, complex, probs, params)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      cat(paste("\rNew best crit:", best.crit, "\n"))
    }

    # If we did a large jump and visited models to save
    if (!is.null(proposal$models)) {
      lo.models <- c(lo.models, proposal$models)
      # If we are doing subsampling and want to update best mliks
      if (!is.null(mliks)) {
        for (mod in 1:length(proposal$models)) {
          model_idx <- bitsToInt(proposal$models[[mod]]$model)
          # This is a model we have seen before
          if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] < proposal$models[[mod]]$crit) {
            # This is a model which has worse mlik in the previous seen
            mliks[[model_idx]] <- proposal$models[[mod]]$crit
          } else if (is.null(mliks[[model_idx]])) {
            mliks[[model_idx]] <- proposal$models[[mod]]$crit
          }
        }
      }
      proposal$models <- NULL
    }
    if (!is.null(mliks)) {
      model_idx <- bitsToInt(proposal$model)
      # This is a model we have seen before
      if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] < proposal$crit) {
        # This is a model which has worse mlik in the previous seen
        mliks[[model_idx]] <- proposal$crit
      } else if (is.null(mliks[[model_idx]])) {
        mliks[[model_idx]] <- proposal$crit
      }
    }

    if (log(runif(1)) <= proposal$alpha) {
      model.cur <- proposal
      accept <- accept + 1
    }
    # Add the current model to the list of visited models
    models[[i]] <- model.cur
  }
  cat("\nMJMCMC done.\n")
  # Calculate acceptance rate
  accept <- accept / N
  # Return formatted results
  return(list(models=models, populations=S, accept=accept, lo.models=lo.models))
}