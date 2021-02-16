# Title     : MJMCMC Model Generator
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
model.proposal.1_4 <- function (model.size, neigh.min, neigh.max, indices=NULL, probs=NULL, prob=F) {
  # If no probs, set all to 1 as we are doing a swap
  if (is.null(probs)) probs <- rep(1,model.size)
  # If no indices are selected, allow all
  if (is.null(indices)) indices <- rep(T,model.size)
  # Set neighborhood size, random or fixed
  if (neigh.max == neigh.min) neigh.size <- neigh.min
  else neigh.size <- sample.int(n = neigh.max - neigh.min, size = 1) + neigh.min - 1
  # Select the negihborhood by sampling from the p covariates
  neighborhood <- sample(model.size, size = neigh.size, prob = probs)

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
model.proposal.1_4.prob <- function (swaps, probs, neigh.size, neigh.max, neigh.min) {
  p <- length(probs) # Get number of available covariates
  print(neigh.size)
  prod(probs[swaps]) / (choose(p, neigh.size)*(neigh.max-neigh.min+1))
}

# Uniform addition and deletion of a covariate (Type 5 and 6)
model.proposal.5_6 <- function (model, addition=T, probs=NULL, prob=F) {
  # If no probs, set all to 1
  if (is.null(probs)) probs <- rep(1,length(model))

  if (addition) change <- which(!model)
  else change <- which(model)

  if (sum(change)==0) stop("No variables to change")

  if (prob) {

  }
  swap <- ind.to.log(sample(change, 1), length(model))
  return(list(swap=swap, S=1))
}

model.proposal.5_6.prob <- function (model, addition) {
  # TODO: Get this finished
  return(0.1)
}

# Function to generate a proposed model given a current one
gen.proposal <- function (model, params, type, indices=NULL, probs=NULL, prob=F) {
  if (type < 5) {
    # Generate a proposal of type 1, 2, 3 or 4
    if (type == 2 || type == 4) {
      params$neigh.min <- params$neigh.size
      params$neigh.max <- params$neigh.size
    }
    proposal <- model.proposal.1_4(length(model), params$neigh.min, params$neigh.max, indices, probs, prob)
  } else if (type == 5) {
    # Generate a proposal of type 5 (addition of a covariate)
    proposal <- model.proposal.5_6(model, addition=T, prob)
  } else if (type == 6) {
    # Generate a proposal of type 6 (subtraction of a covariate)
    proposal <- model.proposal.5_6(model, addition=F, prob)
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
    prob <- model.proposal.5_6.prob(model, addition=T, prob)
  } else if (type == 6) {
    # Generate a proposal of type 6 (subtraction of a covariate)
    prob <- model.proposal.5_6.prob(model, addition=T, prob)
  }
  return(prob)
}

# Function to generate a small random jump given a current model (q.r)
small.rand <- function (model, indices, type, probs=NULL, params, prob=F) {
  if (type == 1 || type == 3) {
    # Load max and min neighborhood sizes from params vector
    neigh.max <- params$small.max
    neigh.min <- params$small.min
  } else if (type == 2 || type == 4) {
    neigh.min <- params$small
    neigh.max <- params$small
  }
  proposal <- model.proposal.1_4(length(model), neigh.min, neigh.max, indices, probs, prob)
  if (prob) return(list(model=xor(proposal$swap, model), prob=proposal$prob)) # Return actual model and probability
  else return(xor(model, proposal)) # Return actual model
}

# Function for generating indices for a large jump given a current model (q.l)
large.jump <- function (model.size, type, probs, params, prob=F) {
  if (type == 1 || type == 3) {
    # Load max and min neighborhood sizes from params vector
    neigh.max <- params$large.max
    neigh.min <- params$large.min
  } else if (type == 2 || type == 4) {
    neigh.min <- params$large
    neigh.max <- params$large
  }
  indices <- model.proposal.1_4(model.size, neigh.min, neigh.max, indices=NULL, probs, prob)
  return(indices) # Return just the indices to be swapped
}

ind.to.log <- function (ind, length) {
  log <- rep(F,length)
  log[ind] <- T
  return(log)
}