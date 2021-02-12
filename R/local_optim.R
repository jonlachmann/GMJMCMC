# Title     : Local optimization
# Objective : Local optimizers for the mode jumps
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing <- function (model, data, loglik.pi, indices, params) {
  # TODO: Is forward and backward local optim depending on each other?

  # TODO: All these should be in params
  # TODO: Further, can any of these be set dynamically based on data?
  t <- 10 # Initial temperature
  t.min <- 0.0001 # Final temperature
  dt <- 3 # Delta temperature
  M <- 12 # Proposals per temperature

  # Calculate current likelihood
  model.lik <- loglik.pi(model, data)
  while (t > t.min) {
    # Make M tries at current temperature
    for (m in 1:M) {
      # Get a modified model as proposal
      proposal <- small.rand(current, indices, type) # TODO: This does not do anything yet
      proposal.lik <- loglik.pi(model, data)

      # Calculate move probability (Bolzmann distribution, see Blum and Roli p. 274)
      alpha <- min(1, exp((model.lik - proposal.lik)/t))
      # Accept move with probability alpha
      if (runif(1) < alpha) {
        model <- proposal
        model.lik <- proposal.lik
      }
    }
    # Update temperature
    t <- t * exp(-dt)
  }
  return(model)
}

greedy.optim <- function () {

}

local.optim <- function (data, loglik.pi, model, features, type) {
  if (type == 1) {
    return(simulated.annealing(model, data, loglik.pi, indices))
  }
  if (type == 2) {
    return(greedy.optim(data, loglik.pi, indices))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}