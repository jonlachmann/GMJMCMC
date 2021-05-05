# Title     : GMJMCMC
# Objective : Genetically Modified Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

# Allow the package to access Rcpp functions
#' @useDynLib GMJMCMC
#' @importFrom Rcpp sourceCpp
NULL

#' Main algorithm for GMJMCMC
#'
#' @param data A matrix containing the data to use in the algorithm,
#' first column should be the dependent variable, second should be the intercept
#' and the rest of the columns should be the independent variables.
#' @param loglik.pi The (log) density to explore
#' @param loglik.alpha The likelihood function to use for alpha calculation
#' @param transforms A list of the available nonlinear transformations for feature generation
#' @param T The number of population iterations
#' @param N The number of iterations per population (total iterations = (T-1)*N+N.final)
#' @param N.final The number of iterations for the final population (total iterations = (T-1)*N+N.final)
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#'
#' @export gmjmcmc
gmjmcmc <- function (data, loglik.pi, loglik.alpha, transforms, T, N, N.final, probs, params, sub=F) {
  # Verify that data is well-formed
  data <- check.data(data)
  # Acceptance probability
  accept <- 0
  # A list of populations that have been visited
  S <- vector("list", T)
  # A list of models that have been visited, refering to the populations
  models <- vector("list", T)

  # TODO: Initialization of first model
  F.0 <- gen.covariates(ncol(data)-2)
  S[[1]] <- F.0
  complex <- complex.features(S[[1]])

  ### Main algorithm loop - Iterate over T different populations
  for (t in 1:T) {
    # Precalculate covariates and put them in data.t
    if (t != 1) data.t <- precalc.features(data, S[[t]], transforms)
    else data.t <- data
    # Initialize first model of population
    model.cur <- as.logical(rbinom(n = length(S[[t]]), size = 1, prob = 0.5))
    model.cur <- list(prob=0, model=model.cur, crit=loglik.pre(loglik.pi, model.cur, complex, data.t, params$loglik), alpha=0)
    if (t==1) best.crit <- model.cur$crit # Set first best criteria value
    # Initialize a vector to contain the models visited in this population
    population.models <- vector("list", N)
    # Initialize a vector to contain local opt visited models
    population.lo.models <- vector("list", 0)
    # If we are running a subsampling strategy, keep a list of best mliks for all models
    if (sub) mliks <- vector("list", (2^(length(S[[t]]))))
    else mliks <- NULL

    if (t==T) N <- N.final
    cat(paste("Population", t, "begin."))
    progress <- 0
    for (i in 1:N) {
      if (N > 40 && i %% floor(N/40) == 0) progress <- print.progressbar(progress, 40)
      proposal <- mjmcmc.prop(data.t, loglik.pi, model.cur, S[[t]], complex, probs, params, mliks)
      if (proposal$crit > best.crit) {
        best.crit <- proposal$crit
        cat(paste("\rNew best crit:", best.crit, "\n"))
      }

      # If we did a large jump and visited models to save
      if (!is.null(proposal$models)) {
        population.lo.models <- c(population.lo.models, proposal$models)
        # If we are doing subsampling and want to update best mliks
        if (!is.null(mliks)) {
          for (mod in 1:length(proposal$models)) {
            model_idx <- bitsToInt(proposal$models[[mod]]$model)
            # This is a model we have seen before
            if (!is.null(mliks[[model_idx]]) && mliks[[model_idx]] < proposal$models[[mod]]$crit) {
              # This is a model which has worse mlik in the previous seen
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
      population.models[[i]] <- model.cur
    }
    cat(paste("\nPopulation", t, "done.\n"))
    # Add the models visited in the current population to the model list
    models[[t]] <- population.models
    # Calculate marginal likelihoods for current features
    marg.probs <- marginal.probs.renorm(population.models)
    # Print the marginal posterior distribution of the features after MJMCMC
    cat(paste("\rCurrent best crit:", best.crit, "\n"))
    cat("Feature importance:\n")
    print.dist(marg.probs, sapply(S[[t]], print.feature, transforms), probs$filter)
    # Generate a new population of features for the next iteration (if this is not the last)
    if (t != T) {
      S[[t+1]] <- gmjmcmc.transition(S[[t]], F.0, data, loglik.alpha, marg.probs, transforms, probs, params$feat)
      complex <- complex.features(S[[t+1]])
    }
  }
  # Calculate acceptance rate
  accept <- accept / (N*T)
  # Return formatted results
  return(list(models=models, populations=S, accept=accept))
}


#' Subalgorithm for generating a new population of features in GMJMCMC
#'
#' @param S.t The current population of features
#' @param F.0 The initial population of features, i.e. the bare covariates
#' @param data The data used in the model, here we use it to generate alphas for new features
#' @param loglik.alpha The log likelihood function to optimize the alphas for
#' @param marg.probs The marginal inclusion probabilities of the current features
#' @param transforms The nonlinear transformations available
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#'
#' @return The updated population of features, that becomes S.t+1
gmjmcmc.transition <- function (S.t, F.0, data, loglik.alpha, marg.probs, transforms, probs, params) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs/probs$filter, 1)))

  # Always keep original covariates if that setting is on
  if (params$keep.org) feats.keep[1:length(F.0)] <- T

  # Avoid removing too many features
  if (sum(feats.keep) < params$keep.min) {
    feats.add.n <- params$keep.min - sum(feats.keep)
    feats.add <- sample(which(!feats.keep), feats.add.n)
    feats.keep[feats.add] <- T
  }

  # Create a list of which features to replace
  feats.replace <- which(!feats.keep)

  # TODO: Let filtered features become part of new features - tuning parameter
  # TODO: Avoid killing too many features
  # Create a list of inclusion probabilities
  marg.probs.use <- sqrt(c(rep(params$eps, length(F.0)), pmin(pmax(marg.probs[feats.keep], params$eps), (1-params$eps))))

  # Perform the replacements
  for (i in feats.replace) {
    prev.size <- length(S.t)
    print(paste0("Replacing feature ", print.feature(S.t[[i]], transforms)))
    S.t[[i]] <- gen.feature(c(F.0, S.t[feats.keep]), marg.probs.use, data, loglik.alpha, transforms, probs, length(F.0), params)
    if (prev.size > length(S.t)) {
      print("Population shrinking, returning.")
      return(S.t)
    }
    feats.keep[i] <- T
    marg.probs.use <- append(marg.probs.use, mean(marg.probs.use), length(F.0)+i-1)
  }

  # Add additional features if the population is not at max size
  if (length(S.t) < params$pop.max) {
    for (i in (length(S.t)+1):params$pop.max) {
      prev.size <- length(S.t)
      S.t[[i]] <- gen.feature(c(F.0, S.t[feats.keep]), marg.probs.use, data, loglik.alpha, transforms, probs, length(F.0), params)
      if (prev.size == length(S.t)) {
        print("Population not growing, returning.")
        return(S.t)
      }
      marg.probs.use <- c(marg.probs.use, params$eps)
    }
  }
  return(S.t)
}

