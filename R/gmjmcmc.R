# Title     : GMJMCMC
# Objective : Genetically Modified Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

### Main algorithm for GMJMCMC

# data is the data to use in the algorithm
# loglik.pi is the (log) density to explore
# T is the number of population iterations
# N is the number of iterations per population (total iterations = T*N)
# probs is a list of the various probability vectors to use TODO: specification of this

gmjmcmc <- function (data, loglik.pi, T, N, probs) {
  # Acceptance probability
  accept <- 0
  # A list of populations that have been visited
  S <- list()
  # A list of models that have been visited, refering to the populations
  models <- list()

  # For every population transition
  for (t in 1:T) {
    # Initialize a vector to contain the models visited in this population
    population.models <- list()
    model.cur <- model # TODO: FIX ME
    for (i in 1:N) {
      proposal <- mjmcmc.prop(model.cur, S[[t]], probs)
      if (log(runif(1)) <= proposal$alpha) {
        model.cur <- proposal$model
        accept <- accept + 1
      }
      # Add the current model to the list of visited models
      append(population.models, model.cur)
    }
    # Add the models visited in the current population to the model list
    append(models, population.models)
    # Generate a new population of features for the next iteration
    S[[t+1]] <- gmjmcmc.transition(S[[t]])
  }
}

# Subalgorithm for generating a proposal and acceptance probability
mjmcmc.prop <- function (model.cur, features, probs) {
  l <- runif(1)
  if (l > probs$large.prob) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 3, size = 1, prob = probs$largejump.prob) # Select large jump kernel
    q.o <- sample.int(n = 3, size = 1, prob = probs$localopt.prob) # Select optimizer function
    q.r <- sample.int(n = 3, size = 1, prob = probs$random.prob) # Select randomization kernel

    ### Do large jump and backwards large jump
    large.jump.ind <- large.jump(q.l) # Get the large jump indices to swap TODO: Implement function
    chi.0.star <- xor(model.cur, large.jump.ind) # Swap indices
    chi.k.star <- local.optim(chi.0.star, features, q.o) # Do local optimization
    gamma.star <- small.rand(chi.k.star, q.r) # Randomize around the mode TODO: Implement function
    chi.0 <- xor(gamma.star, large.jump.ind) # Do a backwards large jump
    chi.k <- local.optim(chi.0, features, q.o) # Do backwards local optimization
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # TODO: Implement functions
    # Calculate small.rand probabilties
    prob.gamma_chi.k <- p.small.rand(gamma, chi.k) # Probability of gamma given chi.k
    prob.gamma.star_chi.k.star <- p.small.rand(gamma.star, chi.k.star) # Probability of gamma.star given chi.k.star
  } else {
    ### Regular MH step

    # Small randomization around current model
    gamma.star <- small.rand(cur.mod) # Randomize around the mode
  }
  # Calculate log likelihoods for models
  prob.cur <- loglik.pi(model.cur)
  prob.gamma.star <- loglik.pi(gamma.star)

  ### Calculate acceptance probability
  if (l > probs$large.prob) {
    # Calculate acceptance probability for large jump
    alpha <- min(0, (prob.gamma.star + prob.gamma_chi.k) - (prob.cur + prob.gamma.star_chi.k.star))
  } else {
    # Calculate regular acceptance probability (assuming small rand to be symmetric here)
    alpha <- min(0, (prob.gamma.star - prob.cur))
  }

  ### Format results and return them
  proposal <- list(model=gamma.star, alpha=alpha)
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

