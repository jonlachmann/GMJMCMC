\name{abalone}
\alias{abalone}
\docType{data}
\title{
Physical measurements of 4177 abalones, a species of sea snail.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("abalone")}
\format{
  A data frame with 4177 observations on the following 9 variables.
  \describe{
    \item{\code{sex}}{Sex of the abalone, \code{F} is female, \code{M} male, and \code{I} infant, categorical.}
    \item{\code{length}}{Longest shell measurement, continuous}
    \item{\code{diameter}}{Diameter Perpendicular to length, continuous}
    \item{\code{height}}{Height with with meat in shell, continuous.}
    \item{\code{whole_weight}}{Grams whole abalone, continuous.}
    \item{\code{shucked_weight}}{Grams weight of meat, continous.}
    \item{\code{viscera_weight}}{Grams gut weight (after bleeding), continuous.}
    \item{\code{shell_weight}}{Grams after being dried, continuous.}
    \item{\code{rings}}{+1.5 gives the age in years, integer}
  }
}
\details{
See the web page \url{https://archive.ics.uci.edu/ml/datasets/Abalone} for more information about the data set.
}
\source{
Dua, D. and Graff, C. (2019). UCI Machine Learning Repository \url{https://archive.ics.uci.edu/ml/. Irvine, CA: University of California, School of Information and Computer Science.}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(abalone)
## maybe str(abalone) ; plot(abalone) ...
}
\keyword{datasets}
