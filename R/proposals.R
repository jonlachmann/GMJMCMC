# Title     : MJMCMC Proposal Generator
# Objective : Generate models for proposals in MJMCMC
# Created by: jonlachmann
# Created on: 2021-02-11

### Generators for proposal models in MJMCMC
### From page 8 in Hubin et. al. 2018
### "Mode Jumping MCMC for Bayesian variable selection in GLMM"

# All the functions are made such that they return the indices that should be swapped in the model
# The idea is then that the function that requested the proposal can just xor the result with the model vector
# We also specify the possibility to include the probability of this exact outcome

# Random change with random size of the neighborhood (Type 1-4)
# With neigh.min=neigh.max, we get a fixed neighborhood size (Type 2 and 4)
# With probs left out, we get a swap instead of a random change (Type 3 and 4)
# Indices tells the sampler which indices it is allowed to sample from
model.proposal.1_4 <- function (model.size, neigh.min, neigh.max, indices, probs=NULL, prob=FALSE) {
  # If no probs, set all to 1 as we are doing a swap
  if (is.null(probs)) probs <- rep(1,model.size)
  # Set neighborhood size, random or fixed
  if (neigh.max == neigh.min) neigh.size <- neigh.min
  else neigh.size <- sample.int(n = neigh.max - neigh.min, size = 1) + neigh.min - 1
  # Select the negihborhood by sampling from the p covariates
  neighborhood <- sample2((1:model.size)[indices], size = neigh.size, prob = probs[indices])

  # Sample which variables to change based on the probs vector
  swaps <- as.logical(rbinom(neigh.size, 1, probs[neighborhood]))
  swaps <- ind.to.log(neighborhood[swaps], model.size)

  if (prob) {
    prob <- model.proposal.1_4.prob(swaps, probs, neigh.size, neigh.min, neigh.max)
    return(list(swap=swaps, S=neigh.size, prob=prob))
  }
  return(list(swap=swaps, S=neigh.size))
}

# Probability for random change with random size of the neighborhood (Type 1)
# By setting neigh.max=neigh.min we get nonrandom neighborhood size (Type 2)
# By setting prob vector to all ones, we get swap instead of random change (Type 3 and 4)
model.proposal.1_4.prob <- function (swaps, probs, neigh.size, neigh.min, neigh.max) {
  p <- length(probs) # Get number of available covariates
  log(prod(probs[swaps]) / (choose(p, neigh.size)*(neigh.max-neigh.min+1)))
}

# Uniform addition and deletion of a covariate (Type 5 and 6)
model.proposal.5_6 <- function (model, addition=TRUE, indices, probs=NULL, prob=FALSE) {
  # If no probs, set all to 1
  if (is.null(probs)) probs <- rep(1,length(model))

  if (addition) change <- which(!model & indices)
  else change <- which(model & indices)

  if (sum(change)==0) swap <- rep(F, length(model)) # Model is full or empty, no change
  else swap <- ind.to.log(change[sample(length(change), 1)], length(model))
  if (prob) {
    prob <- model.proposal.5_6.prob(model, addition)
    return(list(swap=swap, S=1, prob=prob))
  }
  return(list(swap=swap, S=1))
}

# Probability for addition or subtraction of a parameter
model.proposal.5_6.prob <- function (model, addition) {
  p <- length(model)
  modsum <- sum(model)
  if (addition) {
    if (modsum==p) return(log(.Machine$double.eps))
    else return(log(1/(p-modsum)))
  } else {
    if (modsum==0) return(log(.Machine$double.eps))
    else return(log(1/modsum))
  }
}

# Function to generate a proposed model given a current one
gen.proposal <- function (model, params, type, indices=NULL, probs=NULL, prob=FALSE) {
  # If no indices are selected, allow all
  if (is.null(indices)) indices <- rep(T, length(model))
  if (type < 5) {
    # Generate a proposal of type 1, 2, 3 or 4
    if (type == 2 || type == 4) {
      params$neigh.min <- params$neigh.size
      params$neigh.max <- params$neigh.size
    }
    # Generate a proposal of type 3 or 4, i.e. a swap
    if (type > 2) probs <- NULL
    else if (!is.null(probs)) probs[model] <- 1 - probs[model]
    proposal <- model.proposal.1_4(length(model), params$neigh.min, params$neigh.max, indices, probs, prob)
  } else if (type == 5) {
    # Generate a proposal of type 5 (addition of a covariate)
    proposal <- model.proposal.5_6(model, addition=TRUE, indices, probs, prob)
  } else if (type == 6) {
    # Generate a proposal of type 6 (subtraction of a covariate)
    proposal <- model.proposal.5_6(model, addition=FALSE, indices, probs, prob)
  }
  return(proposal)
}

# Calculate the probaility of getting a specified proposal given the current model (i.e. a pdf function)
prob.proposal <- function (proposal, current, type, params, probs=NULL) {
  # Get the difference between the two models
  swaps <- xor(proposal, current)
  if (type < 5) {
    # Prepare parameters for probability calculation
    if (is.null(probs)) probs <- rep(1, length(proposal))
    if (type == 2 || type == 4) {
      params$neigh.min <- params$neigh.size
      params$neigh.max <- params$neigh.size
    }
    prob <- model.proposal.1_4.prob(swaps, probs, params$neigh.size, params$neigh.min, params$neigh.max)
  } else if (type == 5) {
    # Generate a proposal of type 5 (addition of a covariate)
    prob <- model.proposal.5_6.prob(current, addition=TRUE)
  } else if (type == 6) {
    # Generate a proposal of type 6 (subtraction of a covariate)
    prob <- model.proposal.5_6.prob(current, addition=FALSE)
  }
  return(prob)
}
