# Title     : Local optimization
# Objective : Local optimizers for the mode jumps
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing <- function (model, data, loglik.pi, indices, params) {
  # TODO: Can any of these be set dynamically based on data?
  t <- params$t.init # Initial temperature

  # Calculate current likelihood
  model.lik <- loglik.pi(model, data)
  while (t > params$t.min) {
    # Make M tries at current temperature
    for (m in 1:params$M) {
      # Get a modified model as proposal
      proposal <- small.rand(current, indices, type)
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
    t <- t * exp(-params$dt)
  }
  return(model)
}

greedy.optim <- function () {

}

local.optim <- function (model, data, loglik.pi, indices, type, params) {
  if (type == 1) {
    return(simulated.annealing(model, data, loglik.pi, indices, params$sa))
  }
  if (type == 2) {
    return(greedy.optim(model, data, loglik.pi, indices, params$greedy))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}