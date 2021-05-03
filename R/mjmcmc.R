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
mjmcmc <- function (data, loglik.pi, N, probs, params, sub=F) {
  # Verify that data is well-formed
  data <- check.data(data)
  # Acceptance probability
  accept <- 0

  # Create a population of just the covariates
  S <- gen.covariates(ncol(data)-2)
  complex <- complex.features(S)

  # Initialize first model
  model.cur <- as.logical(rbinom(n = length(S), size = 1, prob = 0.5))
  model.cur <- list(prob=0, model=model.cur, crit=loglik.pre(loglik.pi, model.cur, complex, data, params$loglik), alpha=0)
  best.crit <- model.cur$crit # Set first best criteria value

  # A list of models that have been visited
  models <- vector("list", N)
  # Initialize a vector to contain local opt visited models
  lo.models <- vector("list", 0)
  # If we are running a subsampling strategy, keep a list of best mliks for all models
  if (sub) mliks <- vector("list", (2^(length(S))))
  else mliks <- NULL

  cat("\nMJMCMC begin.\n")
  progress <- 0
  for (i in 1:N) {
    if (N > 40 && i %% floor(N/40) == 0) progress <- print.progressbar(progress, 40)
    proposal <- mjmcmc.prop(data, loglik.pi, model.cur, S, complex, probs, params, mliks)
    if (proposal$crit > best.crit) {
      best.crit <- proposal$crit
      cat(paste("\rNew best crit:", best.crit, "\n"))
    }

    # If we did a large jump and visited models to save
    if (!is.null(proposal$models)) {
      lo.models <- c(lo.models, proposal$models)
      # If we are doing subsampling and want to update best mliks
      if (!is.null(mliks)) {
        for (mod in 1:length(proposal$models)) {
          model_idx <- bitsToInt(proposal$models[[mod]]$model)
          # This is a model we have seen before
          if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] < proposal$models[[mod]]$crit) {
            # This is a model which has worse mlik in the previous seen
            print("mlik changed localopt")
            print(proposal$models[[mod]]$crit)
            mliks[[model_idx]] <- proposal$models[[mod]]$crit
          } else if (is.null(mliks[[model_idx]])) {
            mliks[[model_idx]] <- proposal$models[[mod]]$crit
          }
        }
      }
      proposal$models <- NULL
    }
    if (!is.null(mliks)) {
      model_idx <- bitsToInt(proposal$model)
      # This is a model we have seen before
      if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] < proposal$crit) {
        # This is a model which has worse mlik in the previous seen
        mliks[[model_idx]] <- proposal$crit
      } else if (is.null(mliks[[model_idx]])) {
        mliks[[model_idx]] <- proposal$crit
      }
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
  return(list(models=models, populations=S, accept=accept, lo.models=lo.models))
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
mjmcmc.prop <- function (data, loglik.pi, model.cur, features, complex, probs, params, mliks=NULL) {
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
    localopt <- local.optim(chi.0.star, data, loglik.pi, !large.jump$swap, complex, q.o, params) # Do local optimization
    chi.k.star <- localopt$model


    # Randomize around the mode
    proposal <- gen.proposal(chi.k.star, params$random, q.r, !large.jump$swap, prob=T)
    proposal$model <- xor(chi.k.star, proposal$swap)

    # Do a backwards large jump and add in the kernel used in local optim to use the same for backwars local optim.
    chi.0 <- xor(proposal$model, large.jump$swap)

    # Do a backwards local optimization
    localopt2 <- local.optim(chi.0, data, loglik.pi, !large.jump$swap, complex, q.o, params, kern=localopt$kern)
    chi.k <- localopt2$model
    # TODO: We could compare if chi.k is reached by optimising from gamma (model.cur) as "intended"

    ### Calculate acceptance probability
    # Set up the parameters that were used to generate the proposal
    prop.params <- list(neigh.min=params$random$min, neigh.max=params$random$max, neigh.size=proposal$S)

    # Calculate current model probability given proposal
    model.cur$prob <- prob.proposal(proposal$model, chi.k, q.r, prop.params) # Get probability of gamma given chi.k

    # Store models visited during local optimization
    proposal$models <- c(localopt$models, localopt2$models)
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
  proposal$crit <- loglik.pre(loglik.pi, proposal$model, complex, data, params$loglik)

  # TODO: Compare to a list of best mliks for all visited models,
  # TODO: update that list if our estimate is better, otherwise update our estimate.
  # TODO: Save all models visited by local optim, and update the best mliks if we see one during local optim.
  # If we are running with subsampling, check the list for a better mlik
  if (!is.null(mliks)) {
    model_idx <- bitsToInt(proposal$model)
    # This is a model we have not seen before
    if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] > proposal$crit) {
      # This is a model which has better mlik in the previous seen
      proposal$crit <- mliks[[model_idx]]
      print("mlik changed cur")
      print(proposal$crit)
    }
  }

  # Calculate acceptance probability for proposed model
  proposal$alpha <- min(0, (proposal$crit + model.cur$prob) - (model.cur$crit + proposal$prob))

  ### Format results and return them
  proposal$swap <- NULL; proposal$S <- NULL
  return(proposal)
}