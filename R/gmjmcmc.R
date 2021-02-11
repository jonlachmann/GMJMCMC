# Title     : GMJMCMC
# Objective : Genetically Modified Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

# Main algorithm for GMJMCMC
gmjmcmc <- function (data, T, N) {
  # Acceptance probability
  accept <- 0
  # A list of populations that have been visited
  S <- list()
  # A list of models that have been visited, refering to the populations
  models <- list()

  # For every population transition
  for (t in 1:T) {
    model.cur <- model # TODO: FIX ME
    for (i in 1:N) {
      proposal <- mjmcmc.prop(model.cur, S[[t]])
      if (log(runif(1)) <= proposal$alpha) {
        model.cur <- proposal$model
        accept <- accept + 1
      }
      # Add the current model to the list of visited models
      append(models, proposal$model)
    }
    # Generate a new population of features
    S[[t+1]] <- gmjmcmc.transition(S[[t]])
  }
}

# Subalgorithm for generating a proposal and acceptance probability
mjmcmc.prop <- function (model.cur, S, probs) {
  l <- runif(1)
  if (l > probs$large.prob) {
    # TODO: finish this step
    # Select kernels to use for the large jump
    q.l <- sample.int(n = 3, size = 1, prob = probs$largejump.prob) # Select large jump kernel
    q.o <- sample.int(n = 3, size = 1, prob = probs$localopt.prob) # Select optimizer function
    q.r <- sample.int(n = 3, size = 1, prob = probs$random.prob) # Select randomization kernel

    large.jump.ind <- large.jump(q.l) # Get the large jump indices TODO: Implement function


    proposal <- 1
  } else {
    # TODO: generate regular MH proposal
    proposal <- 1
  }
  return(proposal)
}

# Subalgorithm for generating a new population of features
gmjmcmc.transition <- function (S.t) {
  new.feat.count <- 10 # TODO: How to choose this?
  for (i in 1:new.feat.count) {
    proposal <- gen.feature()
    check.collinearity(S.t, proposal)
  }
  return(proposal)
}

