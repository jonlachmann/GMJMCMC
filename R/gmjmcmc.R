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
  S[[1]] <- F.0[as.logical(rbinom(n = length(F.0), size = 1, prob = 1))] # TODO: How should this be done propely?
  model.cur <- as.logical(rbinom(n = length(S[[1]]), size = 1, prob = 0.6))

  ### Main algorithm loop - Iterate over T different populations
  for (t in 1:T) {
    # TODO: Temporary marginal prob calculator - Look into what Aliaksandr mentioned
    marg.probs <- rep(0.5, length(S[[t]]))

    # Precalculate covariates and put them in data.t
    data.t <- precalc.features(data, S[[t]], transforms)
    # Initialize a vector to contain the models visited in this population
    population.models <- vector("list", N)
    for (i in 1:N) {
      print("Generating proposal")
      proposal <- mjmcmc.prop(data.t, loglik.pi, model.cur, S[[t]], marg.probs, probs, params)
      if (log(runif(1)) <= proposal$alpha) {
        print(paste("Accepted move with alpha ", proposal$alpha))
        model.cur <- proposal$model
        accept <- accept + 1
      }
      # Add the current model to the list of visited models
      population.models[[i]] <- model.cur
    }
    # Set the marginal probabilities for the bare covariates if this is the first run
    if (t == 1) cov.probs <- marginal.probs(population.models)
    # Add the models visited in the current population to the model list
    models[[t]] <- population.models
    # Calculate marginal likelihoods for current features
    marg.probs <- marginal.probs(population.models)
    # Generate a new population of features for the next iteration
    S[[t+1]] <- gmjmcmc.transition(S[[t]], F.0, cov.probs, marg.probs, transforms, probs)
  }
  return(list(models=models, populations=S))
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
mjmcmc.prop <- function (data, loglik.pi, model.cur, features, marg.probs, probs, params) {
  l <- runif(1)
  if (l < probs$large) {
    ### Large jump
    print("Large jump!")

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
    model.prop <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump
    chi.0 <- xor(model.prop, large.jump$swap)

    # Do a backwards local optimization
    chi.k <- local.optim(chi.0, data, loglik.pi, features, q.o, params)
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.min=params$random$min, neigh.max=params$random$max, neigh.size=proposal$S)
    prob.gamma_chi.k <- prob.proposal(model.prop, chi.k, q.r, prop.params, marg.probs) # Get probability of gamma given chi.k
    prob.gamma.star_chi.k.star <- proposal$prob # Probability of gamma.star given chi.k.star
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Small randomization around current model
    model.prop <- xor(gen.proposal(model.cur, params$mh, q.g, probs=marg.probs, prob=T)$swap, model.cur)
  }
  # Calculate log likelihoods for models
  prob.cur <- loglik.pre(loglik.pi, model.cur, data)
  prob.gamma.star <- loglik.pre(loglik.pi, model.prop, data)

  ### Calculate acceptance probability
  if (l < probs$large) {
    # Calculate acceptance probability for large jump
    alpha <- min(0, (prob.gamma.star + prob.gamma_chi.k) - (prob.cur + prob.gamma.star_chi.k.star))
  } else {
    # Calculate regular acceptance probability (assuming small rand to be symmetric here)
    alpha <- min(0, (prob.gamma.star - prob.cur))
  }

  ### Format results and return them
  proposal <- list(model=model.prop, alpha=alpha)
  return(proposal)
}

# Subalgorithm for generating a new population of features
gmjmcmc.transition <- function (S.t, F.0, cov.probs, marg.probs, transforms, probs) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs/probs$filter, 1)))

  # Generate new features to replace the filtered ones
  feats.replace <- which(!feats.keep)

  for (i in feats.replace) {
    if (length(c(F.0, S.t[feats.keep])) != length(c(cov.probs, marg.probs[feats.keep]))) {
      print("Uh-oh!")
    }
    S.t[[i]] <- gen.feature(c(F.0, S.t[feats.keep]), transforms, c(cov.probs, marg.probs[feats.keep]), probs)
  }
  return(S.t)
}

