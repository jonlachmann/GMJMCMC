% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmjmcmc.R
\name{gmjmcmc}
\alias{gmjmcmc}
\title{Main algorithm for GMJMCMC
TODO: More documentation - borrow from https://github.com/aliaksah/EMJMCMC2016/blob/master/man/EMJMCMC.Rd if applicable.}
\usage{
gmjmcmc(
  data,
  loglik.pi = gaussian.loglik,
  loglik.alpha = gaussian.loglik.alpha,
  transforms,
  P = 10,
  N.init = 100,
  N.final = 100,
  probs = NULL,
  params = NULL,
  sub = FALSE
)
}
\arguments{
\item{data}{A matrix containing the data to use in the algorithm,
first column should be the dependent variable, second should be the intercept
and the rest of the columns should be the independent variables.}

\item{loglik.pi}{The (log) density to explore}

\item{loglik.alpha}{The likelihood function to use for alpha calculation}

\item{transforms}{A Character vector including the names of the non-linear functions to be used by the modification
and the projection operator.}

\item{P}{The number of generations for GMJMCMC.
The default value is $P = 10$.
A larger value like $P = 50$ might be more realistic for more complicated examples where one expects a lot of non-linear structures.}

\item{N.init}{The number of iterations per population (total iterations = (T-1)*N.init+N.final)}

\item{N.final}{The number of iterations for the final population (total iterations = (T-1)*N.init+N.final)}

\item{probs}{A list of the various probability vectors to use}

\item{params}{A list of the various parameters for all the parts of the algorithm}

\item{sub}{An indicator that if the likelihood is inexact and should be improved each model visit (EXPERIMENTAL!)}
}
\description{
Main algorithm for GMJMCMC
TODO: More documentation - borrow from https://github.com/aliaksah/EMJMCMC2016/blob/master/man/EMJMCMC.Rd if applicable.
}
