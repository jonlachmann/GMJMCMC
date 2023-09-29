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
gmjmcmc.iact <- function (x) {
  xlen <- nrow(x)
  tmp <- acf(as.numeric(x), lag.max = xlen, plot = FALSE)$acf
  li <- min(which(tmp < 0.05))
  out <- 1 + 2 * sum(tmp[1:(li - 1)])
  out
}

# TODO: Inter-thread variance comparison of the marginal log posterior of the best found model.


#' Plot convergence of best/median/mean/other summary log posteriors in time 
#'
#' @param res Object corresponding gmjmcmc output
#' @param FUN The summary statistics to check convergence
#' @param conf which confidence intervals to plot
#' @param burnin how many first populations to skip
#' @param window sliding window for computing the standard deviation
#' @return summary statistics with given confidence intervals
#'
plot.diagn <- function (res, FUN = median, conf = 0.95, burnin = 0, window = 10000) {
  
  if(length(res$thread.best)>0)
    matrix.results <- res$best.log.posteriors
  else
    matrix.results <- as.matrix(unlist(res$best.margs))
  sr <- sapply((1+burnin):dim(matrix.results)[1], FUN = function(x)FUN(matrix.results[x,]))
  sds <- c(0,sapply(2:length(sr), function(x)sd(sr[max(1,x-window):x])))
  
  ub <- sr + qnorm(p = 1-(1-conf)/2)*sds
  lb <- sr - qnorm(p = 1-(1-conf)/2)*sds
  
  plot(y = sr,x = (burnin+1):(dim(matrix.results)[1]), type = "l",col = 1,ylim = c(min(lb), max(ub)), main = "Convergence", xlab = "Population", ylab = "Summary")
  lines(y = ub,x = (burnin+1):(dim(matrix.results)[1]),col = 1,lty = 2)
  lines(y = lb,x = (burnin+1):(dim(matrix.results)[1]),col = 1,lty = 2)
  
  return(list(stat = sr, lower = lb, upper = ub))
}
