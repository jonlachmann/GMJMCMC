# Title     : Diagnostic functions
# Objective : Functions to do diagnostics on a GMJMCMC run
# Created by: jonlachmann
# Created on: 2021-02-24

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