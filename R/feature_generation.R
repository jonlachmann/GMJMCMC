# Title     : Features generation for use in GMJMCMC (Genetically Modified MJMCMC)
# Objective : Generate features for use in GMJMCMC at the population transition step
# Created by: jonlachmann
# Created on: 2021-02-10

# Generate a multiplication feature
gen.multiplication <- function (features, marg.probs) {
  # Sample two features to be multiplied
  feats <- sample.int(n = length(features), size = 2, prob = marg.probs+0.00001, replace = TRUE)
  create.feature(0, features[feats])
}

# Generate a modification feature
gen.modification <- function (features, marg.probs, trans.probs, trans.priors) {
  feat <- sample.int(n = length(features), size = 1, prob = marg.probs+0.00001)
  trans <- sample.int(n = length(trans.probs), size = 1, prob = trans.probs)
  create.feature(trans, features[feat], trans.priors)
}

# Generate a projection feature
gen.projection <- function (features, marg.probs, trans.probs, max.width, max.size, trans.priors) {
  if (!is.null(max.size)) {
    max.width <- min(max.width, max.size + 1)
  }
  feat.count <- sample.int(n = (min(max.width, (length(features)))-1), size = 1)
  feats <- sample.int(n = length(features), size = feat.count, prob = marg.probs+0.00001)
  trans <- sample.int(n = length(trans.probs), size = 1, prob = trans.probs)
  # TODO: Generate alphas properly using various methods
  alphas <- rep(1, length(feats)+1)
  create.feature(trans, features[feats], trans.priors, alphas)
}

# Generate new features from the initial covariates
gen.new <- function (features, F.0.size) {
  covariate <- sample.int(n = F.0.size, size = 1)
  return(features[[covariate]])
}

# Select a feature to generate and generate it
gen.feature <- function (features, marg.probs, data, loglik.alpha, probs, F.0.size, params, verbose = TRUE) {
  tries <- 0
  feat.ok <- F
  while (!feat.ok && tries < 50) {
    feat.type <- sample.int(n = 4, size = 1, prob = probs$gen)
    if (feat.type == 1) feat <- gen.multiplication(features, marg.probs)
    if (feat.type == 2) feat <- gen.modification(features, marg.probs, probs$trans, probs$trans_priors)
    if (feat.type == 3) feat <- gen.projection(features, marg.probs, probs$trans, params$L, params$max.proj.size, probs$trans_priors)
    if (feat.type == 4) feat <- gen.new(features, F.0.size)
    # Check that the feature is not too wide or deep
    if (!(depth.feature(feat) > params$D || width.feature(feat) > params$L)) {
      # Generate alphas using the strategy chosen
      if (params$alpha != "unit") {
        feat <- gen.alphas(params$alpha, feat, data, loglik.alpha, verbose)
      }
      if (!is.null(feat)) {
        # Check for linear dependence of new the feature
        if (length(features) == F.0.size) feats <- list()
        else feats <- features[(F.0.size + 1):length(features)]
        if (params$check.col && !check.collinearity(feat, feats, F.0.size, data, params$col.check.mock.data))
          feat.ok <- T
        else if (!params$check.col)
          feat.ok <- T
      }
    }
    tries <- tries + 1
    params$eps <- min(params$eps + 0.01, 0.5)
    marg.probs <- pmin(pmax(marg.probs, params$eps), (1 - params$eps))
  }
  if (!feat.ok) return(NULL)
  else return(feat)
}

# Check if there is collinearity present in the current set of features
check.collinearity <- function (proposal, features, F.0.size, data, mock) {
  # Add the proposal to the feature list for evaluation
  features[[length(features) + 1]] <- proposal
  # Generate mock data to test with (avoiding too costly computations)
  n <- F.0.size * 2
  if (mock) {
    mock.data <- list(x = matrix(runif(n * (F.0.size + data$fixed), -100, 100), n, F.0.size + data$fixed),
                      y = matrix(runif(n * (ncol(data$y)), -100, 100), n, ncol(data$y)),
                      fixed = data$fixed)
  } else {
    obs_idx <- seq_len(min(n, nrow(data$x)))
    mock.data <- list(x = data$x[obs_idx, ], y = data$y[obs_idx, ], fixed = data$fixed)
  }
  # Use the mock data to precalc the features
  mock.data.precalc <- precalc.features(mock.data, features)
  # Fit a linear model with the mock data precalculated features
  linearmod <- lm.fit(mock.data.precalc$x, mock.data.precalc$y)
  # Check if all coefficients were possible to calculate
  if (sum(is.na(linearmod$coefficients)) == 0) return(FALSE)
  else return(TRUE)
}

# Generate features to represent the covariates, just takes the count needed
gen.covariates <- function (data) {
  features <- list()
  for (i in seq_len(ncol(data$x) - data$fixed)) {
    features <- c(features, i)
    class(features[[length(features)]]) <- "feature"
  }
  return(features)
}