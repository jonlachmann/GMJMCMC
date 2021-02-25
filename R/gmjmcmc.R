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
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#'
#' @export gmjmcmc
gmjmcmc <- function (data, loglik.pi, transforms, T, N, N.final, probs, params) {
  # Acceptance probability
  accept <- 0
  # A list of populations that have been visited
  S <- vector("list", T)
  # A list of models that have been visited, refering to the populations
  models <- vector("list", T)

  # TODO: Initialization of first model
  F.0 <- gen.covariates(ncol(data)-1)
  S[[1]] <- F.0
  complex <- complex.features(S[[1]])
  model.cur <- as.logical(rbinom(n = length(S[[1]]), size = 1, prob = 0.9))
  model.cur <- list(model=model.cur, crit=loglik.pre(loglik.pi, model.cur, complex, data))

  ### Main algorithm loop - Iterate over T different populations
  for (t in 1:T) {
    # TODO: This is using the frequency, maybe switch to renormalized?
    marg.probs <- rep(0.5, length(S[[t]]))

    # Precalculate covariates and put them in data.t
    data.t <- precalc.features(data, S[[t]], transforms)
    # Initialize a vector to contain the models visited in this population
    population.models <- vector("list", N)

    if (t==T) N <- N.final
    for (i in 1:N) {
      proposal <- mjmcmc.prop(data.t, loglik.pi, model.cur, S[[t]], complex, marg.probs, probs, params)
      if (log(runif(1)) <= proposal$alpha) {
        model.cur <- proposal
        accept <- accept + 1
      }
      # Add the current model to the list of visited models
      population.models[[i]] <- model.cur
    }
    print(paste("Population", t, "done."))
    # Set the marginal probabilities for the bare covariates if this is the first run
    if (t == 1) cov.probs <- marginal.probs(population.models)
    # Add the models visited in the current population to the model list
    models[[t]] <- population.models
    # Calculate marginal likelihoods for current features
    marg.probs <- marginal.probs(population.models)
    # Generate a new population of features for the next iteration (if this is not the last)
    if (t != T) {
      S[[t+1]] <- gmjmcmc.transition(S[[t]], F.0, marg.probs, transforms, probs, params)
      complex <- complex.features(S[[t+1]])
    }
  }
  # Calculate acceptance rate
  accept <- accept / (N*T)
  # Return formatted results
  return(list(models=models, populations=S, accept=accept))
}

#' Subalgorithm for generating a proposal and acceptance probability in (G)MJMCMC
#'
#' @param data The data to use in the algorithm
#' @param loglik.pi The the (log) density to explore
#' @param model.cur The current model to make the proposal respective to
#' @param features The features available
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#'
mjmcmc.prop <- function (data, loglik.pi, model.cur, features, complex, marg.probs, probs, params) {
  l <- runif(1)
  if (l < probs$large) {
    ### Large jump

    ### Select kernels to use for the large jump
    q.l <- sample.int(n = 4, size = 1, prob = probs$large.kern) # Select large jump kernel
    q.o <- sample.int(n = 2, size = 1, prob = probs$localopt) # Select optimizer function
    q.r <- sample.int(n = 4, size = 1, prob = probs$random) # Select randomization kernel

    # Generate and do large jump
    large.jump <- gen.proposal(model.cur$model, params$large, q.l) # Get the large jump
    chi.0.star <- xor(model.cur$model, large.jump$swap) # Swap large jump indices

    # Optimize to find a mode
    chi.k.star <- local.optim(chi.0.star, data, loglik.pi, !large.jump$swap, complex, q.o, params) # Do local optimization

    # Randomize around the mode
    proposal <- gen.proposal(chi.k.star, params$random, q.r, !large.jump$swap, prob=T)
    proposal$model <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump
    chi.0 <- xor(proposal$model, large.jump$swap)

    # Do a backwards local optimization
    chi.k <- local.optim(chi.0, data, loglik.pi, !large.jump$swap, complex, q.o, params)
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.min=params$random$min, neigh.max=params$random$max, neigh.size=proposal$S)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, chi.k, q.r, prop.params) # Get probability of gamma given chi.k
  } else {
    ### Regular MH step
    # Select MH kernel
    q.g <- sample.int(n = 6, size = 1, prob = probs$mh)
    # Generate the proposal
    proposal <- gen.proposal(model.cur$model, params$mh, q.g, prob=T)
    proposal$model <- xor(proposal$swap, model.cur$model)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, model.cur$model, q.g, params$mh)
  }
  # Calculate log likelihoods for the proposed model
  proposal$crit <- loglik.pre(loglik.pi, proposal$model, complex, data)

  # Calculate acceptance probability for proposed model
  proposal$alpha <- min(0, (proposal$crit + model.cur$prob) - (model.cur$crit + proposal$prob))

  ### Format results and return them
  proposal$swap <- NULL; proposal$S <- NULL
  return(proposal)
}

# Subalgorithm for generating a new population of features
gmjmcmc.transition <- function (S.t, F.0, marg.probs, transforms, probs, params) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs/probs$filter, 1)))

  # Generate new features to replace the filtered ones
  feats.replace <- which(!feats.keep)

  for (i in feats.replace) {
    S.t[[i]] <- gen.feature(c(F.0, S.t[feats.keep]), transforms, probs, length(F.0), params$feat)
    feats.keep[i] <- T
  }
  return(S.t)
}

