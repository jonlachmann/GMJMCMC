# Title     : Results functions
# Objective : Functions to handle the results produced
# Created by: jonlachmann
# Created on: 2021-05-06

#' Merge a list of multiple results from many runs
#' This function will weight the features based on the best mlik in that population
#' and merge the results together, simplifying by merging equivalent features (having high correlation).
#'
#' @param results A list containing multiple results from GMJMCMC.
#' @param populations Which populations should be merged from the results, can be "all", "last" (default) or "best".
#' @param complex.measure The complex measure to use when finding the simplest equivalent feature,
#' 1=total width, 2=operation count and 3=depth.
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.
#' @param data Data to use when comparing features, default is NULL meaning that mock data will be generated,
#' if data is supplied it should be of the same form as is required by gmjmcmc, i.e. with both x, y and an intercept.
#'
#' @export merge_results
merge_results <- function (results, populations = NULL, complex.measure = NULL, tol = NULL, data = NULL) {
  # Default values
  if (is.null(populations))
    populations <- "last"
  if (is.null(complex.measure))
    complex.measure <- 2
  if (is.null(tol))
    tol <- 0.0000001

  res.count <- length(results)

  # Select populations to use
  res.lengths <- vector("list")
  for (i in 1:res.count) res.lengths[[i]] <- length(results[[i]]$populations)
  if (populations == "last") pops.use <- res.lengths
  else if (populations == "all") pops.use <- lapply(res.lengths, function(x) 1:x)
  else if (populations == "best") pops.use <- lapply(1:res.count, function(x) which.max(unlist(results[[x]]$best.marg)))

  # Get the population weigths to be able to weight the features
  pw <- population.weigths(results, pops.use)
  pop.weights <- pw$weights
  
 
  crit.best <- -Inf
  pop.best <- 1
  thread.best <- 1
  for (i in seq_along(results)) {
    for (pop in 1:(length(results[[1]]$populations))) 
      if(results[[i]]$best.margs[[pop]] > crit.best)
      {
        print(results[[i]]$best.margs[[pop]])
        crit.best <- results[[i]]$best.margs[[pop]]
        pop.best <- pop
        thread.best <- i
      }
  }
  
  
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
  feats.simplest.ids <- unique(feats.map[complex.measure, ])
  feats.simplest.ids <- feats.simplest.ids[order(feats.map[4, feats.simplest.ids])]
  counts <- sapply(feats.simplest.ids, function(x) sum(feats.map[complex.measure,] == x))
  feats.simplest <- features[feats.simplest.ids]
  importance <- feats.map[4, feats.simplest.ids, drop = FALSE]
  merged <- list(features = feats.simplest, marg.probs = importance, counts = counts, results = results, pop.best = pop.best, thread.best = thread.best, crit.best = crit.best, 
                 reported = pw$best, rep.pop = pw$pop.best, rep.thread = pw$thread.best)
  attr(merged, "class") <- "gmjmcmc_merged"
  return(merged)
}

# Function for calculating the weights of different populations based on best mlik
population.weigths <- function (results, pops.use) {
  max.crits <- vector("list")
  max.crit <- -Inf
  pop.best <- 1
  thread.best <- 1
  for (i in seq_along(results)) {
    for (pop in pops.use[[i]]) 
    {
      max.crits <- append(max.crits, results[[i]]$best.margs[[pop]])
      if(results[[i]]$best.margs[[pop]] > max.crit)
      {
        max.crit <- results[[i]]$best.margs[[pop]]
        pop.best <- pop
        thread.best <- i
      }
    }
      
  }
  max.crits <- unlist(max.crits)

  return(list(weights = exp(max.crits-max.crit) / sum(exp(max.crits-max.crit)), best = max.crit, thread.best = thread.best, pop.best = pop.best))
}

#' Function to generate a function string for a model consisting of features
#'
#' @param model A logical vector indicating which features to include
#' @param features The population of features
#' @param link The link function to use, as a string
#' @param round Rounding error for the features in the printed format
#'
#' @export model.string
model.string <- function (model, features, link = "I", round = 2) {
  modelstring <- paste0(sapply(features[model], print.feature, alphas = TRUE, round = round), collapse="+")
  modelfun <- set_alphas(modelstring)
  modelfun$formula <- paste0(link, "(", modelfun$formula, ")")
  return(modelfun)
}

#' Function to print a quick summary of the results
#'
#' @param object The results to use
#' @param pop The population to print for, defaults to last
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param ... Not used.
#'
#' @export
summary.gmjmcmc <- function (object, pop = "last", tol = 0.0001, ...) {
  if (pop == "last") pop <- length(object$models)
  else if (pop == "best") pop <- which.max(unlist(object$best.margs))
  feats.strings <- sapply(object$populations[[pop]], print.feature, round = 2)
  
  summary_internal(
    best = object$best,
    marg.probs = object$marg.probs[[pop]],
    feats.strings = feats.strings,
    best.pop = which.max(unlist(object$best.margs)),
    reported = object$best.margs[[pop]],
    rep.pop = pop,
    tol = tol
  )
}

#' Function to print a quick summary of the results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param ... Not used.
#'
#' @export
summary.gmjmcmc_merged <- function (object, tol = 0.0001, ...) {
  best <- max(sapply(object$results, function (y) y$best))
  feats.strings <- sapply(object$features, print)
  summary_internal(best = object$crit.best, feats.strings, object$marg.probs, 
                   best.pop = object$pop.best, thread.best = object$thread.best,  
                   reported = object$reported, rep.pop = object$rep.pop, rep.thread = object$rep.thread, tol = tol)
}

#' Function to print a quick summary of the results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param ... Not used.
#'
#' @export
summary.mjmcmc <- function (object, tol = 0.0001, ...) {
  return(summary.mjmcmc_parallel(list(object), tol = tol))
}

#' Function to print a quick summary of the results
#'
#' @param object The results to use
#' @param tol The tolerance to use as a threshold when reporting the results.
#' @param ... Not used.
#'
#' @export
summary.mjmcmc_parallel <- function (object, tol = 0.0001, ...) {
  # Get features as strings for printing
  feats.strings <- sapply(object[[1]]$populations, print.feature, round = 2)
  # Get marginal posterior of features
  models <- unlist(lapply(object, function (x) x$models), recursive = FALSE)
  marg.probs <- marginal.probs.renorm(models)$probs
  best <- max(sapply(object, function (x) x$best))
  return(summary_internal(best, feats.strings, marg.probs, tol = tol))
}

summary_internal <- function (best, feats.strings, marg.probs, tol = 0.0001, best.pop = NULL,reported = NULL, rep.pop = NULL, rep.thread = NULL, thread.best = NULL) {
  # Print the final distribution
  keep <- which(marg.probs[1, ] > tol)
  cat("                   Importance | Feature\n")
  print.dist(marg.probs[keep], feats.strings[keep], -1)
  # Print the best log marginal posterior
  if(length(best.pop) > 0){
    
    if(length(thread.best) > 0)
    {
      cat("\nBest   population:", best.pop, " thread:", thread.best,  " log marginal posterior:", best,"\n")
      cat("Report population:", rep.pop," thread:", rep.thread,  " log marginal posterior:", reported,"\n")
      
    }else{
      cat("\nBest   population:", best.pop,  " log marginal posterior:", best,"\n")
      cat("Report population:", rep.pop,  " log marginal posterior:", reported,"\n")
    }
  }else
  {
    cat("\nBest log marginal posterior: ", best,"\n")
  }
  
  cat("\n")
  feats.strings <- feats.strings[keep]
  marg.probs <- marg.probs[1,keep]
  
  ord.marg <- order(marg.probs, decreasing = T) 
  
  
  return(data.frame(feats.strings = feats.strings[ord.marg], marg.probs = marg.probs[ord.marg]))
}

#' Function to get a character respresentation of a list of features 
#' A list or a population of features in a character representation
#'
#' @param x A list of feature objects
#' @param round Rounding precision for parameters of the features
#' @export
string.population <- function(x, round = 2) {
  cbind(sapply(x, print.feature, round = round))
}

#' Function to get a character respresentation of a list of models
#' A list of models in a character representation
#'
#' @param features A list of feature objects on which the models are build
#' @param models A list of model objects
#' @param round Rounding precision for parameters of the features
#' @param link The link function to use, as a string
#' @export
string.population.models <- function(features, models, round = 2, link = "I") {
  cbind(sapply(seq_along(models), FUN = function(x) model.string(features = features, model = (models[[x]]$model), round = round, link = "I")))
}

#' Function to plot the results, works both for results from gmjmcmc and
#' merged results from merge.results
#'
#' @param x The results to use
#' @param count The number of features to plot, defaults to all
#' @param pop The population to plot, defaults to last
#' @param ... Not used.
#'
#' @export
plot.gmjmcmc <- function (x, count = "all", pop = "last", ...) {
  if (pop == "last") pop <- length(x$populations)
  if (is.null(x$populations)) {
    pops <- x$features
    marg.probs <- x$marg.probs
  } else {
    pops <- x$populations[[pop]]
    marg.probs <- x$marg.probs[[pop]]
  }
  plot.mjmcmc(list(populations = pops, marg.probs = marg.probs), count)
}

#' Function to plot the results, works both for results from gmjmcmc and
#' merged results from merge.results
#'
#' @param x The results to use
#' @param count The number of features to plot, defaults to all
#' @param ... Not used.
#'
#' @export
plot.mjmcmc <- function (x, count = "all", ...) {
  ## Get features as strings for printing and marginal posteriors
  # If this is a merged results the structure is one way
  if (is.null(x$populations)) {
    feats.strings <- sapply(x$features, print)
    feats.strings <- paste0(feats.strings, ", ", x$count)
    marg.probs <- x$marg.probs
  } # If this is a result that is not merged, it is another way
  else {
    feats.strings <- sapply(x$populations, print)
    marg.probs <- x$marg.probs
  }

  marg.prob.plot(feats.strings, marg.probs, count)
}

marg.prob.plot <- function (feats.strings, marg.probs, count = "all", ...) {
  # Plot the distribution
  feats.strings <- feats.strings[order(marg.probs)]
  marg.probs <- sort(marg.probs)
  tot <- length(marg.probs)
  if (count=="all") count <- tot
  y <- barplot(marg.probs[(tot - count + 1):tot], horiz = T, xlab = "Marginal probability", ylab = "Feature")
  text((max(marg.probs[(tot - count + 1):tot]) / 2), y, feats.strings[(tot - count + 1):tot])
}

#' Plot a mjmcmc_parallel run
#' @inheritParams plot.mjmcmc
#' @export
plot.mjmcmc_parallel <- function (x, count = "all", ...) {
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
#' @inheritParams plot.gmjmcmc
#' @export
plot.gmjmcmc_merged <- function (x, count = "all", ...) {
  marg.prob.plot(sapply(x$features, print), x$marg.probs, count = count)
}
