# Title     : Local optimization
# Objective : Local optimizers for the mode jumps
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing <- function (model, data, loglik.pi, indices, complex, params, loglikparams, kernel=NULL) {
  # Select which kernel to use for the random steps
  if (is.null(kernel)) kernel <- sample.int(n = 6, size = 1, prob = params$kern$probs)

  temp <- params$t.init # Initial temperature

  # Calculate current likelihood
  model.lik <- loglik.pre(loglik.pi, model, complex, data, loglikparams)
  # print(paste("SA Start:", model.lik))
  while (temp > params$t.min) {
    # Make M tries at current temperature
    for (m in 1:params$M) {
      # Get a modified model as proposal and calculate its likelihood
      proposal <- xor(model, gen.proposal(model, params$kern, kernel, indices)$swap)
      proposal.lik <- loglik.pre(loglik.pi, proposal, complex, data, loglikparams)
      # Calculate move probability for negative steps (Bolzmann distribution, see Blum and Roli p. 274)
      if (proposal.lik > model.lik) alpha <- 1
      else alpha <- min(1, exp((proposal.lik - model.lik)/temp))
      # Accept move with probability alpha
      if (runif(1) < alpha) {
        model <- proposal
        model.lik <- proposal.lik
      }
    }
    # Update temperature
    temp <- temp * exp(-params$dt)
  }
  # print(paste("SA Finish:", model.lik))
  return(list(model=model, kern=kernel))
}

greedy.optim <- function (model, data, loglik.pi, indices, complex, params, loglikparams, kernel=NULL) {
  # Select which kernel to use for the random steps
  if (is.null(kernel)) kernel <- sample.int(n = 6, size = 1, prob = params$kern$probs)

  # Calculate current likelihood
  model.lik <- loglik.pre(loglik.pi, model, complex, data, loglikparams)

  # Run the algorithm for the number of steps specified
  for (i in 1:params$steps) {
    # For each step, do the specified number of tries
    proposal.best <- NULL
    proposal.lik.best <- -Inf
    for (j in 1:params$tries) {
      # Get a modified model as proposal and calculate its likelihood
      proposal <- xor(model, gen.proposal(model, params$kern, kernel, indices)$swap)
      proposal.lik <- loglik.pre(loglik.pi, proposal, complex, data, loglikparams)
      if (proposal.lik > proposal.lik.best) {
        proposal.best <- proposal
        proposal.lik.best <- proposal.lik
      }
    }
    # Accept every improvement
    if (proposal.lik > model.lik) {
      model <- proposal
      model.lik <- proposal.lik
    }
  }
  return(list(model=model, kern=kernel))
}

local.optim <- function (model, data, loglik.pi, indices, complex, type, params, kernel=NULL) {
  if (type == 1) {
    return(simulated.annealing(model, data, loglik.pi, indices, complex, params$sa, params$loglik, kernel))
  }
  if (type == 2) {
    return(greedy.optim(model, data, loglik.pi, indices, complex, params$greedy, params$loglik, kernel))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}