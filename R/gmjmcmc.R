# Title     : GMJMCMC
# Objective : Genetically Modified Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

#' Main algorithm for GMJMCMC
#'
#' @param data The data to use in the algorithm
#' @param loglik.pi The (log) density to explore
#' @param transforms A list of the available nonlinear transformations for feature generation
#' @param T The number of population iterations
#' @param N The number of iterations per population (total iterations = T*N)
#' @param probs A list of the various probability vectors to use TODO: specification of this
#' @param params A list of the various parameters for all the parts of the algorithm TODO: specification of this
#'
#' @export gmjmcmc
gmjmcmc <- function (data, loglik.pi, transforms, T, N, probs, params) {
  # Acceptance probability
  accept <- 0
  # A list of populations that have been visited
  S <- vector("list", T+1)
  # A list of models that have been visited, refering to the populations
  models <- vector("list", T+1)

  # TODO: Initialization of first model
  F.0 <- gen.covariates(ncol(data)-1)
  S[[1]] <- F.0[as.logical(rbinom(n = length(F.0), size = 1, prob = 0.5))]
  model.cur <- as.logical(rbinom(n = length(S[[1]]), size = 1, prob = 0.6))
  marg.probs <- rep(0.5, length(S[[1]]))

  # For every population transition
  for (t in 1:T) {
    # Load a precalculated covariates in data.t
    data.t <- precalc.features(data, S[[t]], transforms)
    # Initialize a vector to contain the models visited in this population
    population.models <- vector("list", N)
    for (i in 1:N) {
      proposal <- mjmcmc.prop(data.t, loglik.pi, model.cur, S[[t]], probs, params)
      if (log(runif(1)) <= proposal$alpha) {
        model.cur <- proposal$model
        accept <- accept + 1
      }
      # Add the current model to the list of visited models
      population.models[[i]] <- model.cur
    }
    # Add the models visited in the current population to the model list
    models[[t]] <- population.models
    # Calculate marginal likelihoods for current features
    marg.probs <- marginal.probs(population.models)
    # Generate a new population of features for the next iteration
    S[[t+1]] <- gmjmcmc.transition(S[[t]], F.0, marg.probs, transforms)
  }
}

#' Subalgorithm for generating a proposal and acceptance probability
#'
#' @param data The data to use in the algorithm
#' @param loglik.pi The the (log) density to explore
#' @param model.cur The current model to make the proposal respective to
#' @param features The features available
#' @param probs A list of the various probability vectors to use TODO: specification of this
#' @param params A list of the various parameters for all the parts of the algorithm TODO: specification of this
#'
mjmcmc.prop <- function (data, loglik.pi, model.cur, features, probs, params) {
  l <- runif(1)
  if (l > probs$large) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 4, size = 1, prob = probs$largejump) # Select large jump kernel
    q.o <- sample.int(n = 2, size = 1, prob = probs$localopt) # Select optimizer function
    q.r <- sample.int(n = 4, size = 1, prob = probs$random) # Select randomization kernel

    # Generate and do large jump
    large.jump <- gen.proposal(model.cur, params$large, q.l) # Get the large jump
    chi.0.star <- xor(model.cur, large.jump$swap) # Swap large jump indices

    # Optimize to find a mode
    chi.k.star <- local.optim(chi.0.star, data, loglik.pi, !large.jump$swap, q.o, params) # Do local optimization

    # Randomize around the mode
    proposal <- gen.proposal(chi.k.star, params$random, q.r, large.jump$swap, prob=T)
    prop.model <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump
    chi.0 <- xor(prop.model, large.jump$swap)

    # Do a backwards local optimization
    chi.k <- local.optim(data, loglik.pi, chi.0, features, q.o)
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # Calculate probaility of current model given chi.k
    swaps <- xor(model.cur, chi.k)

    prob.gamma_chi.k <- model.proposal.1_4.prob(swaps, probs, proposal$S, params$random$min, params$random$max) # Probability of gamma given chi.k
    prob.gamma.star_chi.k.star <- proposal$prob # Probability of gamma.star given chi.k.star
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Small randomization around current model
    proposal <- gen.proposal(model.cur, params$mh, q.g, probs=marg.probs, prob=T)
  }
  # Calculate log likelihoods for models
  prob.cur <- loglik.pi(model.cur, data)
  prob.gamma.star <- loglik.pi(gamma.star, data)

  ### Calculate acceptance probability
  if (l > probs$large) {
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
gmjmcmc.transition <- function (S.t, F.0, marg.probs, transforms, probs) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs/probs$filter, 1)))

  # Generate new features to replace the filtered ones
  feats.replace <- which(!feats.keep)

  for (i in feats.replace) {
    S.t[[i]] <- gen.feature(c(S.t[!feat.keep], F.0), transforms, probs)
  }
  return(S.t)
}

