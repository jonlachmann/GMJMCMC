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
marginal.probs.renorm <- function (models) {
  mod.count <- length(models)
  probs <- rep(0, length=length(models[[1]]$model))
  crit.sum <- 0
  for (i in 1:mod.count) {
    probs <- probs + models[[i]]$model * models[[i]]$crit
    crit.sum <- crit.sum + models[[i]]$crit
  }
  probs <- probs / crit.sum
  return(probs)
}

# Function for precalculating features for a new feature population
precalc.features <- function (data, features, transforms) {
  precalc <- data.frame(matrix(NA, nrow(data), length(features)+1))
  precalc[,1] <- data[,1]
  for (f in 1:length(features)) {
    feature_string <- print.feature(features[[f]], transforms, dataset=T)
    precalc[,(f+1)] <- eval(parse(text=feature_string))
  }
  return(precalc)
}

# Function to call the model function
loglik.pre <- function (loglik.pi, model, data) {
  # Create a formula with only an intercept
  formula <- paste0(colnames(data)[1], " ~ 1 ")
  # Add covariates to formula if we have any
  if (sum(model) != 0) formula <- paste0(formula, "+ ", paste(colnames(data)[c(F,model)], collapse=" + "))
  # Call the model estimator with the data and the formula
  return(loglik.pi(data, model, as.formula(formula)))
}

# Function to summarize results
summary.gmjresult <- function (results, population="last") {
  if (population=="last") pops <- length(results$models)
  feature_strings <- vector("list", length(result$populations[[pops]]))
  for (i in 1:length(feature_strings)) {
    feature_strings[[i]] <- print.feature(result$populations[[pops]][[i]], transforms)
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