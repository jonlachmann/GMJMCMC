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
replace.infinite.data.frame <- function(df, replacewith = c(.Machine$double.xmin, .Machine$double.xmax)){
  df[df == "-Inf"] <- replacewith[1]
  df[df == "Inf"] <- replacewith[2]
  return(df)
}

# Print a progress bar while iterating over a population
print_progressbar <- function (progress, size=40) {
  cat("\r", "|")
  for (p in 1:size-1) {
    if (progress >= p) cat("=")
    else cat(" ")
  }
  cat("|")
  return(progress+1)
}

# Print a distribution as a horizontal histogram
print_dist <- function(probs, labels, threshold, size=30) {
  threshold <- round((1 - threshold) * size)
  for (i in seq_along(probs)) {
    for (p in 1:size - 1) {
      if (p == threshold) cat("!")
      else if ((1 - probs[i]) * size <= p) cat("#")
      else cat(" ")
    }
    cat("| ")
    cat(labels[i], "\n")
  }
}

# A more intuitive sample function which does not change behaviour when length(x) == 1.
sample2 <- function(x, size, replace = FALSE, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}