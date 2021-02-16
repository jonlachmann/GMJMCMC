# Title     : General functions
# Objective : General functions used by, but not specific to GMJMCMC
# Created by: jonlachmann
# Created on: 2021-02-16

# Convert a vector of TRUE indices to a logical vector of specified length
ind.to.log <- function (ind, length) {
  log <- rep(F,length)
  log[ind] <- T
  return(log)
}