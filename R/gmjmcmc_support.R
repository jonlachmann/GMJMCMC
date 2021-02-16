# Title     : GMJMCMC Support functions
# Objective : Support functions for GMJMCMC algorithm
# Created by: jonlachmann
# Created on: 2021-02-11

# Function for calculating marginal probabilities of features given a list of models
marginal.probs <- function (models) {
  mod.count <- length(models)
  probs <- rep(0, length=length(models[[1]]))
  for (i in 1:mod.count) {
    probs <- probs + models[[i]]
  }
  probs <- probs / mod.count
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
  # Call the model estimator with the subset of the data
  return(loglik.pi(data, model, as.formula(formula)))
}