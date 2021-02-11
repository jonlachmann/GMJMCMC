# Title     : MJMCMC Model Generator
# Objective : Generate models for proposals in MJMCMC
# Created by: jonlachmann
# Created on: 2021-02-11

### Generators for proposal models in MJMCMC
### From page 8 in Hubin et. al. 2018
### "Mode Jumping MCMC for Bayesian variable selection in GLMM"

# Random change with random size of the neighborhood (Type 1)
rand.rand.modelgen <- function (model.size, probs, neigh.max, neigh.min) {
  neigh.size <- sample.int(n = neigh.max - neigh.min, size = 1) + neigh.min - 1
  rand.fixed.modelgen(probs, model.size, neigh.size)
}

# Random change with fixed size of the neighborhood (Type 2)
rand.fixed.modelgen <- function (probs, model.size, neigh.size) {
  change <- sample(model.size, size = neigh.size, prob = probs)
  return(which(change))
}

# Swap with random size of the neighborhood (Type 3)
swap.rand.modelgen <- function (model, neigh.max, neigh.min) {
  neigh.size <- sample.int(n = neigh.max - neigh.min, size = 1) + neigh.min - 1
  swap.fixed.modelgen(model, neigh.size)
}

# Swap with fixed size of the neighborhood (Type 4)
swap.fixed.modelgen <- function (model, neigh.size) {
  # Number of currently active features in model
  n.active <- sum(model)
  if (neigh.size > n.active | neigh.size > (length(model) - n.active)) stop("Too many swaps for model")
  # Get active and inactive indices
  active <- which(model)
  inactive <- which(!model)
  # Generate swaps
  active.swap <- sample(active, neigh.size)
  inactive.swap <- sample(inactive, neigh.size)
  return(c(which(active.swap), which(inactive.swap)))
}

# Uniform addition of a covariate (Type 5)
uni.add.modelgen <- function (model) {
  inactive <- which(!model)
  return(sample(inactive, 1))
}

# Uniform deletion of a covariate (Type 6)
uni.del.modelgen <- function () {
  active <- which(model)
  return(sample(active, 1))
}

# Function to generate a small random jump given a current model
small.rand <- function (model, indices, type) {

}