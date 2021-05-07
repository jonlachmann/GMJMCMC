# Title     : Results functions
# Objective : Functions to handle the results produced
# Created by: jonlachmann
# Created on: 2021-05-06

#' Merge a list of multiple results from many runs
#' This function will weight the features based on the best mlik in that population
#' and merge the results together, simplifying by merging equivalent features (having high correlation).
#'
#' @param results A list containing multiple results from GMJMCMC.
#' @param transforms A list of the available nonlinear transformations for feature generation used in the runs.
#' @param complex.measure The complex measure to use when finding the simplest equivalent feature,
#' 1=total width, 2=operation count and 3=depth.
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.
#'
#' @export merge.results
merge.results <- function (results, transforms, complex.measure=1, tol=0) {
  res.count <- length(results)

  # Collect all feature populations and save their lengths to be able to map back to the original populations
  populations <- vector("list", res.count)
  pop.lengths <- matrix(NA, res.count, 1)
  for (i in 1:res.count) {
    populations[[i]] <- results[[i]]$populations[[length(results[[i]]$populations)]]
    pop.lengths[i] <- length(populations[[i]])
  }
  features <- unlist(populations, recursive = F)
  feat.count <- length(features)

  # Get complexity for all features
  complex <- complex.features(features)

  # Get the population weigths to be able to weight the features
  pop.weights <- population.weigths(results)

  # Get the max mlik to use for renormalized estimates
  max.crits <- matrix(NA, res.count)
  for (i in 1:res.count) max.crits[i] <- results[[i]]$best
  max.crit <- max(max.crits)

  # Get all the renormalized estimates, renormalizing using the best found mlik
  renorms <- matrix(NA, feat.count, 1)
  start <- 1
  for (i in 1:res.count) {
    renorms[start:(start+pop.lengths[i]-1)] <- pop.weights[i]*marginal.probs.renorm(results[[i]]$models[[length(results[[i]]$models)]], max.crit)
    start <- start+pop.lengths[i]
  }

  ## Detect equivalent features
  # Generate mock data to compare features with
  mock.data <- matrix(runif((feat.count+2)^2, -100, 100), ncol=feat.count+2)
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features(mock.data, features)[,-(1:2)]
  # Calculate the correlation to find equivalent features
  cors <- cor(mock.data.precalc)
  # A map to link equivalent features together,
  # row 1-3 are the simplest equivalent features based on three different complexity measures
  # row 4 is the total weighted density of those features
  feats.map <- matrix(1:feat.count, 4, feat.count, byrow=T)
  for (i in 1:nrow(cors)) {
    equiv.feats <- which(cors[i,] >= (1-tol))
    # Compare equivalent features complexity to find most simple
    equiv.complex <- list(width=complex$width[equiv.feats], oc=complex$oc[equiv.feats], depth=complex$depth[equiv.feats])
    equiv.simplest <- lapply(equiv.complex, which.min)
    #if (length(equiv.feats) > 1) print("Equivalent features:")
    #if (length(equiv.feats) > 1) print(sapply(features[equiv.feats], print.feature))
    feats.map[1:3,equiv.feats] <- c(equiv.feats[equiv.simplest$width], equiv.feats[equiv.simplest$oc], equiv.feats[equiv.simplest$depth])
    feats.map[4,equiv.feats] <- sum(renorms[equiv.feats])
  }
  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- feats.map[complex.measure,unique(feats.map[complex.measure,])]
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4,feats.simplest.ids])]
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4,feats.simplest.ids]
  return(list(feats=feats.simplest, importance=importance))
}

# Function for calculating the weights of different populations based on best mlik (other version not implemented yet).
population.weigths <- function (results, simple=T) {
  pop.count <- length(results)
  modmats <- vector("list", pop.count)
  max.crits <- matrix(NA, pop.count)
  if (simple) {
    for (i in 1:pop.count) max.crits[i] <- results[[i]]$best
    return(exp(max.crits)/sum(exp(max.crits)))
  } else {
    for (i in 1:pop.count) {
      model.size <- length(results[[i]]$models[[1]]$model)
      models <- results[[i]]$models[[length(results[[i]]$models)]]
      modmats[[i]] <- matrix(unlist(models), ncol=model.size+3, byrow=T)
      max.crits[i] <- max(modmats[[i]][,(model.size+2)])
    }
  }
  return(max.crits)
}

#' Function to generate a function string for a model consisting of features
#'
#' @param model A logical vector indicating which features to include
#' @param features The population of features
#' @param link The link function to use, as a string
#'
#' @export model.string
model.string <- function (model, features, link) {
  modelstring <- paste0(sapply(features[model], print.feature, alphas=T), collapse="+")
  modelfun <- set_alphas(modelstring)
  modelfun$formula <- paste0(link, "(", modelfun$formula, ")")
  return(modelfun)
}