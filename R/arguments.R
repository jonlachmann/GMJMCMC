# Title     : Arguments generator
# Objective : Functions intended to help the user generate proper arguments for the algorithm
# Created by: jonlachmann
# Created on: 2021-02-19

#' Generate a probability list for GMJMCMC
#'
#' @param transforms A list of the transformations used (to get the count)
#'
#' @export gen.probs.list
gen.probs.list <- function (transforms) {
  ### Create a probability list for algorithm

  ## Mode jumping algorithm probabilities
  large <- 0.05                         # probability of a large jump
  large.kern <- c(0, 0, 0, 1)           # probability for type of large jump, only allow type 1-4
  localopt.kern <- c(0.5, 0.5)          # probability for each localopt algorithm
  random.kern <- c(0.3, 0.3, 0.2, 0.2)  # probability for random jump kernels
  mh <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1) # probability for regular mh kernels

  ## Feature generation probabilities
  transcount <- length(transforms)
  filter <- 0.6                         # filtration threshold
  gen <- rep(1/4, 4)                    # probability for different feature generation methods
  trans <- c(1/transcount, transcount)  # probability for each different nonlinear transformation

  ## Compile the list and return
  probs <- list(large=large, large.kern=large.kern, localopt.kern=localopt.kern,
                random.kern=random.kern, filter=filter, gen=gen,
                trans=trans, mh=mh)
  return(probs)
}

#' Generate a parameter list for GMJMCMC
#'
#' @export gen.probs.list
gen.params.list <- function () {
  ### Create a list of parameters for the algorithm

  ## Local optimization parameters
  sa_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15),
                  neigh.size=1, neigh.min=1, neigh.max=2)               # Simulated annealing proposal kernel parameters
  sa_params <- list(t.init=10, t.min=0.0001, dt=3, M=12, kern=sa_kern)  # Simulated annealing parameters
  greedy_params <- list()                                               # Greedy algorithm parameters

  ## MJMCMC parameters
  large_params <- list(neigh.size=2, neigh.min=1, neigh.max=2)          # Large jump parameters
  random_params <- list(neigh.size=1, neigh.min=1, neigh.max=2)         # Small random jump parameters
  mh_params <- list(neigh.size=1, neigh.min=1, neigh.max=2)             # Regular MH parameters

  ## GM parameters
  feat_params <- list(D=5, L=15)                                        # Hard limits on feature complexity

  ## Compile the list and return
  params <- list(mh=mh_params, large=large_params, random=random_params,
                 sa=sa_params, greedy=greedy_params, feat=feat_params)
  return(params)
}