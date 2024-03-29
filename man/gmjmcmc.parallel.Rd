% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel.R
\name{gmjmcmc.parallel}
\alias{gmjmcmc.parallel}
\title{Run multiple gmjmcmc runs in parallel returning a list of all results.}
\usage{
gmjmcmc.parallel(
  runs,
  cores = getOption("mc.cores", 2L),
  merge.options = list(populations = "best", complex.measure = 2, tol = 1e-07),
  data,
  loglik.pi = gaussian.loglik,
  loglik.alpha = gaussian.loglik.alpha(),
  transforms,
  ...
)
}
\arguments{
\item{runs}{The number of runs to run}

\item{cores}{The number of cores to run on}

\item{merge.options}{A list of options to pass to the \code{\link[=merge_results]{merge_results()}} function run after the}

\item{data}{A matrix containing the data to use in the algorithm,
first column should be the dependent variable, second should be the intercept
and the rest of the columns should be the independent variables.}

\item{loglik.pi}{The (log) density to explore}

\item{loglik.alpha}{The likelihood function to use for alpha calculation}

\item{transforms}{A Character vector including the names of the non-linear functions to be used by the modification
and the projection operator.}

\item{...}{Further params passed to mjmcmc.}
}
\value{
Results from multiple gmjmcmc runs
}
\description{
Run multiple gmjmcmc runs in parallel returning a list of all results.
}
