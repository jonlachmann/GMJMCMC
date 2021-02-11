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
