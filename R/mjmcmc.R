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
mjmcmc <- function (data, loglik.pi, N, probs, params) {
  # Verify that data is well-formed
  data <- check.data(data)
  # Acceptance probability
  accept <- 0
  # A list of models that have been visited
  models <- vector("list", N)

  # Create a population of just the covariates
  S <- gen.covariates(ncol(data)-2)
  complex <- complex.features(S)

  # Initialize first model of population
  model.cur <- as.logical(rbinom(n = length(S), size = 1, prob = 0.5))
  model.cur <- list(prob=0, model=model.cur, crit=loglik.pre(loglik.pi, model.cur, complex, data, params$loglik), alpha=0)
  best.crit <- model.cur$crit # Set first best criteria value
  # Initialize a vector to contain the models visited in this population
  population.models <- vector("list", N)

  print("MJMCMC begin.")
  progress <- 0
  for (i in 1:N) {
    if (N > 40 && i %% floor(N/40) == 0) progress <- print.progressbar(progress, 40)
    proposal <- mjmcmc.prop(data, loglik.pi, model.cur, S, complex, probs, params)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      cat(paste("\rNew best crit:", best.crit, "\n"))
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
  return(list(models=models, populations=S, accept=accept))
}