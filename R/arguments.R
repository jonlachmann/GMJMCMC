# Title     : Arguments generator
# Objective : Functions intended to help the user generate proper arguments for the algorithm
# Created by: jonlachmann
# Created on: 2021-02-19

#' Generate a probability list for MJMCMC (Mode Jumping MCMC)
#'
#' @return A list of probabilities to be used as input for the mjmcmc function.
#'
#' @examples
#' gen.probs.mjmcmc()
#' 
#' @export gen.probs.mjmcmc
#' 
gen.probs.mjmcmc <- function () {
  ## Mode jumping algorithm probabilities
  large <- 0.05                         # probability of a large jump
  large.kern <- c(0, 0, 0, 1)           # probability for type of large jump, only allow type 1-4
  localopt.kern <- c(0.5, 0.5)          # probability for each localopt algorithm
  random.kern <- c(0.3, 0.3, 0.2, 0.2)  # probability for random jump kernels
  mh <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1) # probability for regular mh kernels

  # Compile the list
  probs <- list(large=large, large.kern=large.kern, localopt.kern=localopt.kern,
                random.kern=random.kern, mh=mh)

  return(probs)
}

#' Generate a probability list for GMJMCMC (Genetically Modified MJMCMC)
#'
#' @param transforms A list of the transformations used (to get the count).
#'
#' @return A list of probabilities to be used as input for the gmjmcmc function.
#'
#' @examples
#' gen.probs.gmjmcmc(c("p0", "exp_dbl"))
#' 
#'
#' @export gen.probs.gmjmcmc
gen.probs.gmjmcmc <- function (transforms) {
  if (!is.character(transforms))
    stop("The argument transforms must be a character vector specifying the transformations.")

  # Get probs for mjmcmc
  probs <- gen.probs.mjmcmc()

  ## Feature generation probabilities
  transcount <- length(transforms)
  filter <- 0.6                             # filtration threshold
  gen <- c(0.40,0.40,0.10,0.10)             # probability for different feature generation methods
  trans <- rep(1 / transcount, transcount)  # probability for each different nonlinear transformation
  trans_priors <- rep(1, transcount)        # Default values assigned to each transformation to be used as "operation count".

  probs$filter <- filter
  probs$gen <- gen
  probs$trans <- trans
  probs$trans_priors <- trans_priors

  return(probs)
}

#' Generate a parameter list for MJMCMC (Mode Jumping MCMC)
#'
#' @param data The dataset that will be used in the algorithm
#'
#' @return A list of parameters to use when running the mjmcmc function.
#'
#' Note that the $loglik item is an empty list, which is passed to the log likelihood function of the model,
#' intended to store parameters that the estimator function should use.
#'
#' @examples
#' gen.params.mjmcmc(matrix(rnorm(600), 100))
#' 
#'
#' @export gen.params.mjmcmc
gen.params.mjmcmc <- function (data) {
  ### Create a list of parameters for the algorithm

  ## Get the dimensions of the data to set parameters based on it
  data.dim <- data.dims(data)
  ncov <- data.dim[2] - 2

  ## Local optimization parameters
  sa_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15),
                  neigh.size=1, neigh.min=1, neigh.max=2)               # Simulated annealing proposal kernel parameters
  sa_params <- list(t.init=10, t.min=0.0001, dt=3, M=12, kern=sa_kern)  # Simulated annealing parameters
  greedy_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15),
                      neigh.size=1, neigh.min=1, neigh.max=2)           # Greedy algorithm proposal kernel parameters
  greedy_params <- list(steps=20, tries=3, kern=greedy_kern)            # Greedy algorithm parameters (60 models default)

  ## MJMCMC parameters
  burn_in <- 100                                                        # TODO

  # Large jump parameters
  large_params <- list(
    neigh.size = min(as.integer(ncov * 0.35),35),
    neigh.min = min(as.integer(ncov * 0.25),25),
    neigh.max = min(as.integer(ncov * 0.45),45)
  )
  random_params <- list(neigh.size = 1, neigh.min = 1, neigh.max = 2)  # Small random jump parameters
  mh_params <- list(neigh.size = 1, neigh.min = 1, neigh.max = 2)      # Regular MH parameters
  ## Compile the list and return
  params <- list(burn_in=burn_in, mh=mh_params, large=large_params, random=random_params,
                 sa=sa_params, greedy=greedy_params, loglik=list())

  return(params)
}

#' Generate a parameter list for GMJMCMC (Genetically Modified MJMCMC)
#'
#' @param data The dataset that will be used in the algorithm
#'
#' @return A list of parameters to use when running the mjmcmc function.
#'
#' @examples
#' gen.params.gmjmcmc(matrix(rnorm(600), 100))
#' 
#'
#' @export gen.params.gmjmcmc
gen.params.gmjmcmc <- function (data) {
  # Get mjmcmc params
  params <- gen.params.mjmcmc(data)

  ncov <- ncol(data) - 2

  feat_params <- list(D = 5, L = 15,                                # Hard limits on feature complexity
                      alpha = 0,                                    # alpha strategy (0 = None, 1,2,3 = strategies as per Hubin et al.) TODO: Fully Bayesian
                      pop.max = min(100,as.integer(ncov * 1.5)),    # Max features population size
                      keep.org = FALSE,                             # Always keep original covariates in every population
                      prel.filter = NULL,                           # Filtration threshold for first population (i.e. filter covariates even if keep.org=TRUE)
                      keep.min = 0.8,                               # Minimum proportion of features to always keep [0,1]
                      eps = 0.05,                                   # Inclusion probability limit for feature generation
                      check.col = TRUE,                             # Whether the colinearity should be checked
                      col.check.mock.data = FALSE,                  # Use mock data when checking for colinearity during feature generation
                      max.proj.size = 15)                           # Maximum projection size
  params$feat <- feat_params
  params$rescale.large <- FALSE
  params$prel.filter <- NULL                                        # Specify which covariates to keep in the first population. See Issue #15.

  return(params)
}
