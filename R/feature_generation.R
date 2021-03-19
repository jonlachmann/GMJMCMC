# Title     : Features generation for use in GMJMCMC
# Objective : Generate features for use in GMJMCMC at the population transition step
# Created by: jonlachmann
# Created on: 2021-02-10

# Generate a multiplication feature
gen.multiplication <- function (features) {
  # generate two features to be multiplied
  feats <- sample.int(n = length(features), size = 2)
  create.feature(0, features[feats])
}

# Generate a modification feature
gen.modification <- function (features, transforms, trans.probs) {
  feat <- sample.int(n = length(features), size = 1)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  create.feature(trans, features[feat])
}

# Generate a projection feature
# TODO: This is not working according to spec yet
gen.projection <- function (features, transforms, trans.probs) {
  feat.count <- sample.int(n = length(features), size = 1) # TODO: Should be a specific distribution?
  feats <- sample.int(n = length(features), size = feat.count)
  trans <- sample.int(n = length(transforms), size = 1, prob = trans.probs)
  # TODO: Generate alphas properly using various methods
  alphas <- rep(1, length(feats)+1)
  create.feature(trans, features[feats], alphas)
}

# Generate a new features from the initial covariates
gen.new <- function (features, F.0.size) {
  covariate <- sample.int(n = F.0.size, size = 1)
  return(features[[covariate]])
}

# Select a feature to generate and generate it
gen.feature <- function (features, transforms, probs, F.0.size, params) {
  feat.type <- sample.int(n = 4, size = 1, prob = probs$gen)
  colinear <- T
  too.large <- T
  while (colinear || too.large) {
    if (feat.type == 1) feat <- gen.multiplication(features)
    if (feat.type == 2) feat <- gen.modification(features, transforms, probs$trans)
    if (feat.type == 3) feat <- gen.projection(features, transforms, probs$trans)
    if (feat.type == 4) feat <- gen.new(features, F.0.size)
    # Check that the feature is not too wide or deep
    if (depth.feature(feat) <= params$D && width.feature(feat) <= params$L) too.large <- F
    # Check for linear dependence of new the feature
    if (!too.large) {
      if (length(features) == F.0.size) feats <- list()
      else feats <- features[(F.0.size+1):length(features)]
      colinear <- check.collinearity(feat, feats, transforms, F.0.size)
    }
  }
  print(paste("New feature:", print.feature(feat, transforms), "depth:", depth.feature(feat), "width:", width.feature(feat)))
  return(feat)
}

check.collinearity <- function (proposal, features, transforms, F.0.size) {
  # Add the proposal to the feature list for evaluation
  features[[length(features)+1]] <- proposal
  # Generate mock data to test with (avoiding too costly computations)
  mock.data <- matrix(c(runif((F.0.size*2), -1, 1), rep(1,F.0.size*2),
                        runif((F.0.size*2)*(F.0.size), -1, 1)), F.0.size*2, F.0.size+2)
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features(mock.data, features, transforms)
  # Fit a linear model with the mock data precalculated features
  linearmod <- lm(as.data.frame(mock.data.precalc[,-2]))
  # Check if all coefficients were possible to calculate
  if (sum(is.na(linearmod$coefficients)) == 0) return(F)
  else return(T)
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