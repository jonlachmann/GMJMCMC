# Title     : GMJMCMC Support functions
# Objective : Support functions for GMJMCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

# Function to verify inputs and help the user find if they did anything wrong
verify.inputs <- function (data, loglik.pi, transforms, T, N, N.final, probs, params) {
  # Get information about the data
  n.cov <- ncol(data)-1
  n.obs <- nrow(data)

  # Get information about the transforms
  n.trans <- length(transforms)

  # Verify that the probability list is formatted as it should be
  if (length(probs$large) != 1) error <- c(error, "Too many large jump probabilities (should be 1).")
  else if (probs$large > 1 | probs$large < 0) error <- c(error, "Large jump probability must be in [0,1].")
  if (length(probs$large.kern) != 4) error <- c(error, "There should be 4 large jump kernel probabilities.")
  if (sum(probs$large.kern > 1 | probs$large.kern < 0) != 0) error <- c(error, "Large jump kernel probabilities must be in [0,1].")
}

# Function for calculating marginal inclusion probabilities of features given a list of models
marginal.probs <- function (models) {
  mod.count <- length(models)
  probs <- rep(0, length=length(models[[1]]$model))
  for (i in 1:mod.count) {
    probs <- probs + models[[i]]$model
  }
  probs <- probs / mod.count
  return(probs)
}

# Function for calculating feature importance through renormalized model estimates
marginal.probs.renorm <- function (models, max_mlik=NULL) {
  model.size <- length(models[[1]]$model)
  models.matrix <- matrix(unlist(models), ncol=model.size+3, byrow=T)
  models.matrix <- models.matrix[(!duplicated(models.matrix[,2:(model.size+1)], dim=1, fromLast=T)),]
  if(is.null(max_mlik)) max_mlik <- max(models.matrix[,(model.size+2)])
  crit.sum <- sum(exp(models.matrix[,(model.size+2)]-max_mlik))
  probs <- matrix(NA,1,model.size)
  for (i in 2:(model.size+1)) probs[i-1] <- sum(exp(models.matrix[as.logical(models.matrix[,i]),(model.size+2)]-max_mlik))/crit.sum
  return(probs)
}

# Function for precalculating features for a new feature population
precalc.features <- function (data, features, transforms) {
  precalc <- matrix(NA, nrow(data), length(features)+2)
  precalc[,1:2] <- data[,1:2]
  for (f in 1:length(features)) {
    feature_string <- print.feature(features[[f]], transforms, dataset=T)
    precalc[,(f+2)] <- eval(parse(text=feature_string))
  }
  # Replace any -Inf and Inf values caused by under- or overflow
  precalc <- replace.infinite.data.frame(precalc)
  return(precalc)
}

# TODO: Compare to previous mliks here instead, also add a flag to do that in full likelihood estimation scenarios.
# Function to call the model function
loglik.pre <- function (loglik.pi, model, complex, data, params) {
  # Get the complexity measures for just this model
  complex <- list(width=complex$width[model], depth=complex$depth[model])
  # Call the model estimator with the data and the model, note that we add the intercept to every model
  return(loglik.pi(data[,1], data[,-1], c(T,model), complex, params))
}

#' Summarize results from GMJMCMC
#'
#' @param results The results from GMJMCMC
#' @param populations A list of the populations to include in the summary, defaults to the last one
#'
#' @export summary.gmjresult
summary.gmjresult <- function (results, population="last") {
  if (population=="last") pops <- length(results$models)
  else pops <- population
  feature_strings <- vector("list", length(results$populations[[pops]]))
  for (i in 1:length(feature_strings)) {
    feature_strings[[i]] <- print.feature(results$populations[[pops]][[i]], transforms)
  }
  feature_importance <- marginal.probs.renorm(results$models[[pops]])
  return(list(features=feature_strings, importance=feature_importance))
}

# Function to print all features in a model
print.model <- function (model, features, transforms) {
  # Create a list to store the features in
  model_print <- vector("list", sum(model$model))
  for (i in 1:length(model$model)) {
    if (model$model[i]) model_print[[i]] <- print.feature(features[[i]], transforms)
  }
  return(model_print)
}

# Function to check the data
# Checks that there is an intercept in the data, adds it if missing
# Coerces the data to be of type matrix
check.data <- function (data) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
    print("Data coerced to matrix type")
  }
  if (sum(data[,2] == 1) != nrow(data)) {
    data <- cbind(data[,1],1,data[,-1])
    print("Intercept added to data")
  }
  return(data)
}

# Function to get the dimensions of a dataset, adding an intercept if necessary
data.dims <- function (data) {
  dims <- dim(data)
  if (sum(data[,2] == 1) != nrow(data)) {
    dims[2] <- dims[2] + 1
  }
  return(dims)
}