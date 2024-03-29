# Title     : GMJMCMC
# Objective : Genetically Modified Mode Jumping MCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

# Allow the package to access Rcpp functions
#' @useDynLib GMJMCMC
#' @importFrom Rcpp sourceCpp
NULL

#' Main algorithm for GMJMCMC
#' TODO: More documentation - borrow from https://github.com/aliaksah/EMJMCMC2016/blob/master/man/EMJMCMC.Rd if applicable.
#'
#' @param data A matrix containing the data to use in the algorithm,
#' first column should be the dependent variable, second should be the intercept
#' and the rest of the columns should be the independent variables.
#' @param loglik.pi The (log) density to explore
#' @param loglik.alpha The likelihood function to use for alpha calculation
#' @param transforms A Character vector including the names of the non-linear functions to be used by the modification 
#' and the projection operator. 
#' @param P The number of generations for GMJMCMC. 
#' The default value is $P = 10$.
#' A larger value like $P = 50$ might be more realistic for more complicated examples where one expects a lot of non-linear structures. 
#' @param N.init The number of iterations per population (total iterations = (T-1)*N.init+N.final)
#' @param N.final The number of iterations for the final population (total iterations = (T-1)*N.init+N.final)
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#' @param sub An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)
#'
#' @export gmjmcmc
gmjmcmc <- function (data, loglik.pi = gaussian.loglik, loglik.alpha = gaussian.loglik.alpha, transforms, P = 10, N.init = 100, N.final = 100, probs = NULL, params = NULL, sub = FALSE) {
  # Verify that the data is well-formed
  data <- check.data(data)

  # Generate default probabilities and parameters if there are none supplied.
  if (is.null(probs)) probs <- gen.probs.gmjmcmc(transforms)
  if (is.null(params)) params <- gen.params.gmjmcmc(data)

  # Extract labels from column names in dataframe
  labels <- get.labels(data)
  # Set the transformations option
  options("gmjmcmc-transformations" = transforms)
  # Acceptance probability per population
  accept <- vector("list", P)
  accept <- lapply(accept, function (x) x <- 0)
  # A list of populations that have been visited
  S <- vector("list", P)
  # A list of models that have been visited, refering to the populations
  models <- vector("list", P)
  lo.models <- vector("list", P)
  # A list of all the marginal probabilities for the features, per population
  marg.probs <- vector("list", P)
  # A list of all the marginal probabilities for the models, per population
  model.probs <- vector("list", P)
  # A list of all the indices of the models which the marginal probabilities for the models refer to, per population
  model.probs.idx <- vector("list", P)
  # A list of all the best marginal model likelihoods, per population
  best.margs <- vector("list", P)

  # Create first population
  F.0 <- gen.covariates(ncol(data) - 2)
  if (is.null(params$feat$prel.filter))
    S[[1]] <- F.0
  else
    S[[1]] <- F.0[params$feat$prel.filter]

  complex <- complex.features(S[[1]])

  ### Main algorithm loop - Iterate over P different populations
  for (p in seq_len(P)) {
    # Set population iteration count
    if (p != P) N <- N.init
    else N <- N.final
    # Precalculate covariates and put them in data.t
    if (length(params$feat$prel.filter) > 0 | p != 1) data.t <- precalc.features(data, S[[p]])
    else data.t <- data
    
    # Initialize first model of population
    model.cur <- as.logical(rbinom(n = length(S[[p]]), size = 1, prob = 0.5))
    model.cur.res <- loglik.pre(loglik.pi, model.cur, complex, data.t, params$loglik)
    model.cur <- list(prob = 0, model = model.cur, coefs = model.cur.res$coefs, crit = model.cur.res$crit, alpha = 0)
    best.crit <- model.cur$crit # Reset first best criteria value

    # Run MJMCMC over the population
    cat(paste("Population", p, "begin."))
    mjmcmc_res <- mjmcmc.loop(data.t, complex, loglik.pi, model.cur, N, probs, params, sub)
    cat(paste("\nPopulation", p, "done.\n"))

    # Add the models visited in the current population to the model list
    models[[p]] <- mjmcmc_res$models
    lo.models[[p]] <- mjmcmc_res$lo.models
    # Store marginal likelihoods for current features
    marg.probs[[p]] <- mjmcmc_res$marg.probs
    # Store marginal likelihoods for the visited models
    model.probs[[p]] <- mjmcmc_res$model.probs
    # Store indices for which the marginal likelihoods for the visited models refer to
    model.probs.idx[[p]] <- mjmcmc_res$model.probs.idx
    # Store best marginal model probability for current population
    best.margs[[p]] <- mjmcmc_res$best.crit
    # Print the marginal posterior distribution of the features after MJMCMC
    cat(paste("\rCurrent best crit:", mjmcmc_res$best.crit, "\n"))
    cat("Feature importance:\n")
    print.dist(marg.probs[[p]], sapply(S[[p]], print.feature, labels = labels, round = 2), probs$filter)
    if (params$rescale.large) prev.large <- params$large
    # Generate a new population of features for the next iteration (if this is not the last)
    if (p != P) {
      S[[p + 1]] <- gmjmcmc.transition(S[[p]], F.0, data, loglik.alpha, marg.probs[[1]], marg.probs[[p]], labels, probs, params$feat)
      complex <- complex.features(S[[p + 1]])
      if (params$rescale.large) params$large <- lapply(prev.large, function(x) x * length(S[[p + 1]]) / length(S[[p]]))
    }
  }
  # Calculate acceptance rate
  accept.tot <- sum(unlist(accept)) / (N.init * (P - 1) + N.final)
  accept <- lapply(accept, function (x) x / N.init)
  accept[[P]] <- accept[[P]] * N.init / N.final
  # Return formatted results
  results <- list(
    models = models,                   # All models per population
    lo.models = lo.models,             # All local optim models per population
    populations = S,                   # All features per population
    marg.probs = marg.probs,           # Marginal feature probabilities per population
    model.probs = model.probs,         # Marginal feature probabilities per population
    model.probs.idx = model.probs.idx, # Marginal feature probabilities per population
    best.margs = best.margs,           # Best marginal model probability per population
    accept = accept,                   # Acceptance rate per population
    accept.tot = accept.tot,           # Overall acceptance rate
    best = max(unlist(best.margs))     # Best marginal model probability throughout the run
  )
  attr(results, "class") <- "gmjmcmc"
  return(results)
}


#' Subalgorithm for generating a new population of features in GMJMCMC
#'
#' @param S.t The current population of features
#' @param F.0 The initial population of features, i.e. the bare covariates
#' @param data The data used in the model, here we use it to generate alphas for new features
#' @param loglik.alpha The log likelihood function to optimize the alphas for
#' @param marg.probs The marginal inclusion probabilities of the current features
#' @param marg.probs.F.0 The marginal inclusion probabilities of the initial population of features
#' @param labels Variable labels for printing
#' @param probs A list of the various probability vectors to use
#' @param params A list of the various parameters for all the parts of the algorithm
#'
#' @return The updated population of features, that becomes S.t+1
gmjmcmc.transition <- function (S.t, F.0, data, loglik.alpha, marg.probs.F.0, marg.probs, labels, probs, params) {
  # Sample which features to keep based on marginal inclusion below probs$filter
  feats.keep <- as.logical(rbinom(n = length(marg.probs), size = 1, prob = pmin(marg.probs / probs$filter, 1)))

  # Always keep original covariates if that setting is on
  if (params$keep.org) {
    if (params$prel.filter > 0) {
      # Do preliminary filtering if turned on
      feats.keep[(seq_along(F.0))[marg.probs.F.0 > params$prel.filter]] <- T
      #removed.count <- sum(marg.probs.F.0 <= params$prel.filter)
      #cat("Preliminary filtering removed",removed.count,"features.")
    } # Keep all if no preliminary filtering
    else feats.keep[seq_along(F.0)] <- T
  }

  # Avoid removing too many features
  if (length(feats.keep) > 0 && mean(feats.keep) < params$keep.min) {
    feats.add.n <- round((params$keep.min - mean(feats.keep)) * length(feats.keep))
    feats.add <- sample(which(!feats.keep), feats.add.n)
    feats.keep[feats.add] <- T
  }

  # Create a list of which features to replace
  feats.replace <- which(!feats.keep)

  # TODO: Let filtered features become part of new features - tuning parameter
  # Create a list of inclusion probabilities
  marg.probs.use <- c(rep(params$eps, length(F.0)), pmin(pmax(marg.probs, params$eps), (1-params$eps)))

  # Perform the replacements
  for (i in feats.replace) {
    prev.size <- length(S.t)
    prev.feat.string <- print.feature(S.t[[i]], labels=labels, round = 2)
    S.t[[i]] <- gen.feature(c(F.0, S.t), marg.probs.use, data, loglik.alpha, probs, length(F.0), params)
    if (prev.size > length(S.t)) {
      cat("Removed feature", prev.feat.string, "\n")
      cat("Population shrinking, returning.\n")
      return(S.t)
    }
    cat("Replaced feature", prev.feat.string, "with", print.feature(S.t[[i]], labels=labels, round = 2), "\n")
    feats.keep[i] <- T
    marg.probs.use[i] <- mean(marg.probs.use)
  }

  # Add additional features if the population is not at max size
  if (length(S.t) < params$pop.max) {
    for (i in (length(S.t)+1):params$pop.max) {
      prev.size <- length(S.t)
      S.t[[i]] <- gen.feature(c(F.0, S.t), marg.probs.use, data, loglik.alpha, probs, length(F.0), params)
      if (prev.size == length(S.t)) {
        cat("Population not growing, returning.\n")
        return(S.t)
      }
      cat("Added feature", print.feature(S.t[[i]], labels=labels, round = 2), "\n")
      marg.probs.use <- c(marg.probs.use, params$eps)
    }
  }
  return(S.t)
}

