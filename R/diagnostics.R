# Title     : Diagnostic functions
# Objective : Functions to do diagnostics on a GMJMCMC run
# Created by: jonlachmann
# Created on: 2021-02-24

# TODO: We do not really know if this works as intended...
# TODO: It should be investigated how to properly calculate acf for binary vectors
#' Integrated Auto-Correlation Time
#'
#' @param x A matrix where each row is a model
#'
#' @export gmjmcmc.iact
gmjmcmc.iact <- function(x){
  xlen <- nrow(x)
  tmp <- acf(as.numeric(x),lag.max = xlen,plot = FALSE)$acf
  li <- min(which(tmp<0.05))
  out <- 1 + 2*sum(tmp[1:(li-1)])
  out
}

#' Total explored density
#'
#' @param models A list of models with their criteria values (i.e. marginal likelihoods)
#'
#' @export gmjmcmc.totdens
gmjmcmc.totdens <- function (models) {
  model_size <- length(models[[1]]$model)
  models <- matrix(unlist(models), ncol=3+model_size, byrow=TRUE)
  total_dens <- sum(exp(models[,2+model_size]))
  return(total_dens)
}

gmjmcmc.totdens.plot <- function (result) {
  totdens_pops <- sapply(result$models, gmjmcmc.totdens)
  plot(log(totdens_pops), type="l", xlab="Population", ylab="Total log density")
  return(totdens_pops)
}