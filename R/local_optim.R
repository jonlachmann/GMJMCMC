# Title     : Local optimization
# Objective : Local optimizers for the mode jumps
# Created by: jonlachmann
# Created on: 2021-02-11

simulated.annealing <- function (model, data, loglik.pi, indices, complex, params) {
  # Select which kernel to use for the random steps
  kernel <- sample.int(n = 6, size = 1, prob = params$kern$probs)
  # TODO: Can any of these be set dynamically based on data?
  t <- params$t.init # Initial temperature

  # Calculate current likelihood
  model.lik <- loglik.pre(loglik.pi, model, complex, data)
  while (t > params$t.min) {
    # Make M tries at current temperature
    for (m in 1:params$M) {
      # Get a modified model as proposal and calculate its likelihood
      proposal <- xor(model, gen.proposal(model, params$kern, kernel, indices)$swap)
      proposal.lik <- loglik.pre(loglik.pi, proposal, complex, data)

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

greedy.optim <- function (model, data, loglik.pi, indices, complex, params) {
 return(model)
}

local.optim <- function (model, data, loglik.pi, indices, complex, type, params) {
  if (type == 1) {
    return(simulated.annealing(model, data, loglik.pi, indices, complex, params$sa))
  }
  if (type == 2) {
    return(greedy.optim(model, data, loglik.pi, indices, complex, params$greedy))
  }
  if (type == 3) {
    return("not implemented")
  }
  stop("Invalid local optimizer chosen")
}