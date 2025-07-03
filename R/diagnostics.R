# Title     : Diagnostic functions
# Objective : Functions to do diagnostics on a GMJMCMC run
# Created by: jonlachmann
# Created on: 2021-02-24


#' Plot convergence of best/median/mean/other summary log posteriors in time 
#'
#' @param res Object corresponding gmjmcmc output
#' @param FUN The summary statistics to check convergence
#' @param conf Which confidence intervals to plot
#' @param burnin How many first populations to skip
#' @param window Sliding window for computing the standard deviation
#' @param ylim Limits for the plotting range; if unspecified, min and max of confidence intervals will be used
#' @param ... Additional graphical parameters passed to plot and lines functions, e.g. col, lwd, lty, main, xlab, ylab, ylim
#'
#' @return A list of summary statistics for checking convergence with given confidence intervals
#' 
#' @examples
#' result <- gmjmcmc(y = matrix(rnorm(100), 100),
#'                   x = matrix(rnorm(600), 100), 
#'                   P = 2, 
#'                   transforms =  c("p0", "exp_dbl"))
#' diagnstats <- diagn_plot(result)
#' 
#' @export
diagn_plot <- function(res, FUN = median, conf = 0.95, burnin = 0, window = 5, ylim = NULL, ...) {
  
  args <- list(...)
  args[["..."]] <- NULL  # Remove any "..." element to avoid warning
  
  if(length(res$thread.best) > 0)
    matrix.results <- res$best.log.posteriors
  else
    matrix.results <- as.matrix(unlist(res$best.margs))
  
  sr <- sapply((1 + burnin):dim(matrix.results)[1], function(x) FUN(matrix.results[x, ]))
  sds <- c(0, sapply(2:length(sr), function(x) sd(sr[max(1, x - window):x])))
  
  ub <- sr + qnorm(1 - (1 - conf) / 2) * sds
  lb <- sr - qnorm(1 - (1 - conf) / 2) * sds
  
  if (is.null(ylim) && is.null(args$ylim))
    ylim <- c(min(lb), max(ub))
  else if (is.null(ylim))
    ylim <- args$ylim
  
  main <- if (!is.null(args$main)) args$main else "Convergence"
  xlab <- if (!is.null(args$xlab)) args$xlab else "Population"
  ylab <- if (!is.null(args$ylab)) args$ylab else "Summary"
  
  args$main <- NULL
  args$xlab <- NULL
  args$ylab <- NULL
  args$ylim <- NULL
  
  do.call(plot, c(
    list(
      y = sr,
      x = (burnin + 1):(dim(matrix.results)[1]),
      type = "l",
      col = 1,
      ylim = ylim,
      main = main,
      xlab = xlab,
      ylab = ylab
    ),
    args
  ))
  
  do.call(lines, c(
    list(
      y = ub,
      x = (burnin + 1):(dim(matrix.results)[1]),
      col = 1,
      lty = 2
    ),
    args
  ))
  
  do.call(lines, c(
    list(
      y = lb,
      x = (burnin + 1):(dim(matrix.results)[1]),
      col = 1,
      lty = 2
    ),
    args
  ))
  
  return(list(stat = sr, lower = lb, upper = ub))
}

