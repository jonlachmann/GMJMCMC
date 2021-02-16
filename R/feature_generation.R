# Title     : Features generation for use in GMJMCMC
# Objective : Generate features for use in GMJMCMC at the population transition step
# Created by: jonlachmann
# Created on: 2021-02-10

# Generate a multiplication feature
gen.multiplication <- function (features, marg.probs) {
  # generate two features to be multiplied
  feats <- sample.int(n = length(features), size = 2, prob = marg.probs)
  create.feature(0, features[feats])
}

# Generate a modification feature
gen.modification <- function (features, transforms, marg.probs, trans.probs) {
  feat <- sample.int(n = length(features), size = 1, prob = marg.probs)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  create.feature(trans, features[feat])
}

# Generate a projection feature
# TODO: This is not working according to spec yet
gen.projection <- function (features, transforms, marg.probs, trans.probs) {
  feat.count <- sample.int(n = length(features), size = 1) # TODO: Should be a specific distribution?
  feats <- sample.int(n = length(features), size = feat.count, prob = marg.probs)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  # TODO: Generate alphas properly using various methods
  alphas <- rep(1, length(feats)+1)
  create.feature(trans, features[feats], alphas)
}

# Select a feature to generate and generate it
gen.feature <- function (features, transforms, marg.probs, probs) {
  # TODO: Do not generate too advanced features, note max depth and width
  feat.type <- sample.int(n = 3, size = 1, prob = probs$gen)
  colinear <- T
  while (colinear) {
    if (feat.type == 1) feat <- gen.multiplication(features, marg.probs)
    if (feat.type == 2) feat <- gen.modification(features, transforms, marg.probs, probs$trans)
    if (feat.type == 3) feat <- gen.projection(features, transforms, marg.probs, probs$trans)
    # TODO: Check for collinearity etc.
    colinear <- check.collinearity(features, feat)
  }
  print(paste("New feature:", print.feature(feat, transforms)))
  return(feat)
}

check.collinearity <- function (features, proposal) {
  # TODO: How can we do this?
  return(F)
}

# Generate features to represent the covariates, just takes the count needed
gen.covariates <- function (count) {
  features <- list()
  for (i in 1:count) {
    features <- c(features, i)
    class(features[[i]]) <- "feature"
  }
  return(features)
}