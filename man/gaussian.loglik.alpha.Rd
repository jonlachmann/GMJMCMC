% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{gaussian.loglik.alpha}
\alias{gaussian.loglik.alpha}
\title{Log likelihood function for gaussian regression for alpha calculation
This function is just the bare likelihood function
Note that it only gives a proportional value and is equivalent to least squares}
\usage{
gaussian.loglik.alpha(a, data, mu_func)
}
\arguments{
\item{a}{A vector of the alphas to be used}

\item{data}{The data to be used for calculation}

\item{mu_func}{The function linking the mean to the covariates,
as a string with the alphas as a\link{i}.}
}
\description{
Log likelihood function for gaussian regression for alpha calculation
This function is just the bare likelihood function
Note that it only gives a proportional value and is equivalent to least squares
}
