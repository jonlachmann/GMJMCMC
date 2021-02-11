# Title     : Features generation for use in GMJMCMC
# Objective : Generate features for use in GMJMCMC at the population transition step
# Created by: jonlachmann
# Created on: 2021-02-10

# Generate a multiplication feature
genMultiplication <- function (features, pool.probs) {
  # generate two features to be multiplied
  feats <- sample.int(n = length(features), size = 2, prob = pool.probs)
  createFeature(0, features[[feats]])
}

# Generate a modification feature
genModification <- function (features, transforms, pool.probs, trans.probs) {
  feat <- sample.int(n = length(features), size = 1, prob = pool.probs)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  createFeature(trans, feat)
}

# Generate a projection feature
# TODO: This is not working according to spec yet
genProjection <- function (features, transforms, pool.probs, trans.probs) {
  feat.count <- sample.int(n = length(features), size = 1) # TODO: Should be a specific distribution?
  feats <- sample.int(n = length(features), size = feat.count, prob = pool.probs)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  # TODO: Generate alphas properly using various methods
  alphas <- rep(1, length(feats)+1)
  createFeature(trans, feats, alphas)
}

# Select a feature to generate and generate it
gen.feature <- function (features, transforms, probs) {
  feat.type <- sample.int(n = 3, size = 1, prob = probs$gen)
  colinear <- T
  while (colinear) {
    if (feat.type == 1) feat <- genMultiplication(features, probs$pool)
    if (feat.type == 2) feat <- genModification(features, transforms, probs$pool, probs$trans)
    if (feat.type == 3) feat <- genProjection(features, transforms, probs$pool, probs$trans)
    # TODO: Check for collinearity etc.
    colinear <- check.collinearity(features, feat)
  }
  return(feat)
}

check.collinearity <- function (features, proposal) {
  # TODO: How can we do this?
  return(T)
}