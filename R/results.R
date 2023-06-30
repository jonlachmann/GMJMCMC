# Title     : Results functions
# Objective : Functions to handle the results produced
# Created by: jonlachmann
# Created on: 2021-05-06

#' Merge a list of multiple results from many runs
#' This function will weight the features based on the best mlik in that population
#' and merge the results together, simplifying by merging equivalent features (having high correlation).
#'
#' @param results A list containing multiple results from GMJMCMC.
#' @param complex.measure The complex measure to use when finding the simplest equivalent feature,
#' 1=total width, 2=operation count and 3=depth.
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.
#' @param data Data to use when comparing features, default is NULL meaning that mock data will be generated,
#' if data is supplied it should be of the same form as is required by gmjmcmc, i.e. with both x, y and an intercept.
#'
#' @export merge.results
merge.results <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL) {
  # Default values
  if (is.null(populations))
    populations <- "last"
  if (is.null(complex.measure))
    complex.measure <- 1
  if (is.null(tol))
    tol <- 0

  res.count <- length(results)

  # Select populations to use
  res.lengths <- vector("list")
  for (i in 1:res.count) res.lengths[[i]] <- length(results[[i]]$populations)
  if (populations == "last") pops.use <- res.lengths
  else if (populations == "all") pops.use <- lapply(res.lengths, function(x) 1:x)
  else if (populations == "best") pops.use <- lapply(1:res.count, function(x) which.max(unlist(results[[x]]$best.marg)))

  # Get the population weigths to be able to weight the features
  pop.weights <- population.weigths(results, pops.use)

  # Collect all features and their renormalized weighted values
  features <- vector("list")
  renorms <- vector("list")
  weight_idx <- 1
  for (i in 1:res.count) {
    results[[i]]$pop.weights <- rep(NA, length(results[[i]]$populations))
    results[[i]]$model.probs <- list()
    for (pop in pops.use[[i]]) {
      features <- append(features, results[[i]]$populations[[pop]])
      renorms <- append(renorms, pop.weights[weight_idx] * results[[i]]$marg.probs[[pop]])
      results[[i]]$pop.weights[pop] <- pop.weights[weight_idx]
      weight_idx <- weight_idx + 1

      model.probs <- marginal.probs.renorm(results[[i]]$models[[pop]], "models")
      results[[i]]$model.probs[[pop]] <- model.probs$probs
      results[[i]]$models[[pop]] <- results[[i]]$models[[pop]][model.probs$idx]
    }
    accept.tot <- results[[i]]$accept.tot
    best <- results[[i]]$best
    results[[i]] <- lapply(results[[i]], function (x) x[pops.use[[i]]])
    results[[i]]$accept.tot <- accept.tot
    results[[i]]$best <- best
  }
  renorms <- unlist(renorms)
  na.feats <- which(is.na(renorms))
  if (length(na.feats) != 0) {
    cat("Underflow occurred,", length(na.feats), "features removed.\n")
    renorms <- renorms[-na.feats]
    features <- features[-na.feats]
  }
  feat.count <- length(features)

  # Get complexity for all features
  complex <- complex.features(features)

  ## Detect equivalent features
  # Generate mock data to compare features with
  if (is.null(data)) mock.data <- matrix(runif((feat.count+2)^2, -100, 100), ncol=feat.count+2)
  else mock.data <- check.data(data)
  mock.data.precalc <- precalc.features(mock.data, features)[,-(1:2)]

  # Calculate the correlation to find equivalent features
  cors <- cor(mock.data.precalc)

  # A map to link equivalent features together,
  # row 1-3 are the simplest equivalent features based on three different complexity measures
  # row 4 is the total weighted density of those features
  feats.map <- matrix(1:feat.count, 4, feat.count, byrow=T)
  for (i in seq_len(nrow(cors))) {
    equiv.feats <- which(cors[i, ] >= (1 - tol))
    # Compare equivalent features complexity to find most simple
    equiv.complex <- list(width=complex$width[equiv.feats], oc=complex$oc[equiv.feats], depth=complex$depth[equiv.feats])
    equiv.simplest <- lapply(equiv.complex, which.min)
    feats.map[1:3,equiv.feats] <- c(equiv.feats[equiv.simplest$width], equiv.feats[equiv.simplest$oc], equiv.feats[equiv.simplest$depth])
    feats.map[4,equiv.feats] <- sum(renorms[equiv.feats])
  }
  # Select the simplest features based on the specified complexity measure and sort them
  feats.simplest.ids <- feats.map[complex.measure,unique(feats.map[complex.measure, ])]
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[1,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]
  merged <- list(features = feats.simplest, marg.probs = importance, counts = counts, results = results)
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}

# Function for calculating the weights of different populations based on best mlik
population.weigths <- function (results, pops.use) {
  max.crits <- vector("list")
  for (i in seq_along(results)) {
    for (pop in pops.use[[i]]) max.crits <- append(max.crits, results[[i]]$best.margs[[pop]])
  }
  max.crits <- unlist(max.crits)
  max.crit <- max(max.crits)
  return(exp(max.crits-max.crit) / sum(exp(max.crits-max.crit)))
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

#' Function to print a quick summary of the results
#'
#' @param results The results to use
#' @param pop The population to print for, defaults to last
#'
#' @export
summary.gmjmcmc <- function (results, pop = "last", tol = 0.0001) {
  if (pop == "last") pop <- length(results$models)
  summary.mjmcmc(list(best = results$best, models = results$models[[pop]], populations = results$populations[[pop]]),tol = tol)
}

#' @export
summary.gmjmcmc_merged <- function (x, tol = 0.0001) {
  best <- max(sapply(x$results, function (y) y$best))
  feats.strings <- sapply(x$features, print)
  summary_internal(best, feats.strings, x$marg.probs, tol = tol)
}

#' Function to print a quick summary of the results
#'
#' @param results The results to use
#'
#' @export
summary.mjmcmc <- function (results, tol = 0.0001) {
  return(summary.mjmcmc_parallel(list(results),tol = tol))
}

#' Function to print a quick summary of the results
#'
#' @param results The results to use
#'
#' @export
summary.mjmcmc_parallel <- function (results , tol = 0.0001) {
  # Get features as strings for printing
  feats.strings <- sapply(results[[1]]$populations, print.feature, round = 2)
  # Get marginal posterior of features
  models <- unlist(lapply(results, function (x) x$models), recursive = FALSE)
  marg.probs <- marginal.probs.renorm(models)$probs
  best <- max(sapply(results, function (x) x$best))
  return(summary_internal(best, feats.strings, marg.probs, tol = tol))
}

summary_internal <- function (best, feats.strings, marg.probs, tol = 0.0001) {
  # Print the final distribution
  keep <- which(marg.probs[1, ] > tol)
  cat("                   Importance | Feature\n")
  print.dist(marg.probs[keep], feats.strings[keep], -1)
  # Print the best marginal likelihood
  cat("\nBest marginal likelihood: ", best, "\n")
  
  
  feats.strings <- feats.strings[keep]
  marg.probs <- marg.probs[1,keep]
  
  ord.marg <- order(marg.probs, decreasing = T) 
  
  
  return(data.frame(feats.strings = feats.strings[ord.marg], marg.probs = marg.probs[ord.marg]))
}

#' Function to plot the results, works both for results from gmjmcmc and
#' merged results from merge.results
#'
#' @param results The results to use
#' @param count The number of features to plot, defaults to all
#' @param pop The population to plot, defaults to last
#'
#' @export
plot.gmjmcmc <- function (results, count="all", pop="last") {
  if (pop == "last") pop <- length(results$populations)
  if (is.null(results$populations)) {
    pops <- results$features
    marg.probs <- results$marg.probs
  } else {
    pops <- results$populations[[pop]]
    marg.probs <- results$marg.probs[[pop]]
  }
  plot.mjmcmc(list(populations = pops, marg.probs = marg.probs), count)
}

#' Function to plot the results, works both for results from gmjmcmc and
#' merged results from merge.results
#'
#' @param results The results to use
#' @param count The number of features to plot, defaults to all
#'
#' @export
plot.mjmcmc <- function (results, count = "all") {
  ## Get features as strings for printing and marginal posteriors
  # If this is a merged results the structure is one way
  if (is.null(results$populations)) {
    feats.strings <- sapply(results$features, print)
    feats.strings <- paste0(feats.strings, ", ", results$count)
    marg.probs <- results$marg.probs
  } # If this is a result that is not merged, it is another way
  else {
    feats.strings <- sapply(results$populations, print)
    marg.probs <- results$marg.probs
  }

  marg.prob.plot(feats.strings, marg.probs, count)
}

marg.prob.plot <- function (feats.strings, marg.probs, count = "all") {
  # Plot the distribution
  feats.strings <- feats.strings[order(marg.probs)]
  marg.probs <- sort(marg.probs)
  tot <- length(marg.probs)
  if (count=="all") count <- tot
  y <- barplot(marg.probs[(tot - count + 1):tot], horiz = T, xlab = "Marginal probability", ylab = "Feature")
  text((max(marg.probs[(tot - count + 1):tot]) / 2), y, feats.strings[(tot - count + 1):tot])
}

#' Plot a mjmcmc_parallel run
#' @export
plot.mjmcmc_parallel <- function (x, count = "all") {
  merged <- merge.mjmcmc_parallel(x)
  marg.prob.plot(merged$features, merged$marg.probs, count)
}

merge.mjmcmc_parallel <- function (x) {
  run.weights <- run.weigths(x)
  marg.probs <- x[[1]]$marg.probs * run.weights[1]
  for (i in seq_along(x[-1])) {
    marg.probs <- marg.probs + x[[i]]$marg.probs * run.weights[i]
  }
  return(structure(
    list(
      features = sapply(x[[1]]$populations, print),
      marg.probs = marg.probs,
      results = x
    ),
    class = "mjmcmc_merged"
  ))
}

run.weigths <- function (results) {
  best.crits <- sapply(results, function (x) x$best.crit)
  max.crit <- max(best.crits)
  return(exp(best.crits - max.crit) / sum(exp(best.crits - max.crit)))
}

#' Plot a gmjmcmc_merged run
#' @export
plot.gmjmcmc_merged <- function (x, count = "all") {
  marg.prob.plot(sapply(x$features, print), x$marg.probs, count = count)
}
