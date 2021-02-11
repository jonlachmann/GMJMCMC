# Title     : MJMCMC Model Generator
# Objective : Generate models for proposals in MJMCMC
# Created by: jonlachmann
# Created on: 2021-02-11

### Generators for proposal models in MJMCMC
### From page 8 in Hubin et. al. 2018
### "Mode Jumping MCMC for Bayesian variable selection in GLMM"

# Random change with random size of the neighborhood (Type 1)
rand.rand.modelgen <- function () {

}

# Random change with fixed size of the neighborhood (Type 2)
rand.fixed.modelgen <- function () {

}

# Swap with random size of the neighborhood (Type 3)
swap.rand.modelgen <- function () {

}

# Swap with fixed size of the neighborhood (Type 4)
swap.fixed.modelgen <- function () {

}

# Uniform addition of a covariate (Type 5)
uni.add.modelgen <- function () {

}

# Uniform deletion of a covariate (Type 6)
uni.del.modelgen <- function () {

}

# Function to generate a small random jump given a current model
small.rand <- function (model, indices, type) {

}