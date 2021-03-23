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

# Replace non-finite elements of a data frame
replace.infinite.data.frame <- function(df, replacewith=c(.Machine$double.xmin, .Machine$double.xmax)){
  df[df == "-Inf"] <- replacewith[1]
  df[df == "Inf"] <- replacewith[2]
  return(df)
}

# Print a progress bar while iterating over a population
print.progressbar <- function (progress, size=40) {
  cat("\r", "|")
  for (p in 1:size-1) {
    if (progress >= p) cat("=")
    else cat(" ")
  }
  cat("|")
  return(progress+1)
}
