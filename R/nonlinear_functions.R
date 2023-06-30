# Title     : Nonlinear functions for GMJMCMC
# Objective : Give an example library of nonlinear functions to use for GMJMCMC
# Created by: jonlachmann
# Created on: 2021-03-22

#' Sigmoid function
#'
#' @param x The vector of values
#' @return The sigmoid of x
#'
#' @export sigmoid
sigmoid <- function(x)  1 / (1 - exp(-x))

#' Sine function for degrees
#'
#' @param x The vector of values in degrees
#' @return The sine of x
#'
#' @export sin.rad
sin.rad <- function(x) sin(x/180*pi)

#' Cosine function for degrees
#'
#' @param x The vector of values in degrees
#' @return The cosine of x
#'
#' @export cos.rad
cos.rad <- function(x) cos(x/180*pi)

#' Double exponential function
#'
#' @param x The vector of values
#' @return e^(-abs(x))
#'
#' @export exp.dbl
exp.dbl <- function(x) exp(-abs(x))


#' Square root function
#'
#' @param x The vector of values
#' @return The square root of the absolute value of x
#'
#' @export sqroot
sqroot <- function(x) abs(x)^(1/2)

#' Cube root function
#'
#' @param x The vector of values
#' @return The cube root of x
#'
#' @export troot
troot <- function(x) abs(x)^(1/3)

#' To the 2.3  power function
#'
#' @param x The vector of values
#' @return x^2.3
#'
#' @export troot
to23 <- function(x) abs(x)^(2.3)

#' To the 7/2  power function
#'
#' @param x The vector of values
#' @return x^(7/2)
#'
#' @export troot
to72 <- function(x) abs(x)^(7/2)

#' Gaussian function
#'
#' @param x The vector of values
#' @return e^(-x^2)
#'
#' @export gauss
gauss <- function(x) exp(-x*x)

#' To 2.5 power
#'
#' @param x The vector of values
#' @return x^(2.5)
#'
#' @export to25
to25 <- function(x)abs(x)^(2.5)


#' To 3.5 power
#'
#' @param x The vector of values
#' @return x^(3.5)
#'
#' @export to35
to35 <- function(x)abs(x)^(3.5)

#' p0 polynomial term
#'
#' @param x The vector of values
#' @return log(abs(x) + .Machine$double.eps)
#'
#' @export p0
p0 <- function(x) log(abs(x)+.Machine$double.eps)

#' pm1 polynomial term
#'
#' @param x The vector of values
#' @return sign(x)*(abs(x)+.Machine$double.eps)^(-1)
#'
#' @export pm1
pm1 <- function(x) sign(x)*(abs(x)+.Machine$double.eps)^(-1)

#' pm2 polynomial term
#'
#' @param x The vector of values
#' @return sign(x)*(abs(x)+.Machine$double.eps)^(-2)
#'
#' @export pm2
pm2 <- function(x) sign(x)*(abs(x)+.Machine$double.eps)^(-2)

#' pm05 polynomial term
#'
#' @param x The vector of values
#' @return  (abs(x)+.Machine$double.eps)^(-0.5)
#'
#' @export pm05
pm05 <- function(x) (abs(x)+.Machine$double.eps)^(-0.5)

#' p05 polynomial term
#'
#' @param x The vector of values
#' @return (abs(x)+.Machine$double.eps)^(0.5)
#'
#' @export p05
p05 <- function(x) (abs(x)+.Machine$double.eps)^(0.5)

#' p2 polynomial term
#'
#' @param x The vector of values
#' @return x^(2)
#'
#' @export p2
p2 <- function(x) x^(2)

#' p3 polynomial term
#'
#' @param x The vector of values
#' @return x^(3)
#'
#' @export p3
p3 <- function(x) x^(3)

#' p0p0 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*p0(x)
#'
#' @export p0p0
p0p0 <- function(x) p0(x)*p0(x)

#' p0pm1 polynomial terms
#'
#' @param x The vector of values
#' @return p0(x)*(x+.Machine$double.eps)^(-1)
#'
#' @export p0pm1
p0pm1 <- function(x) p0(x)*(x+.Machine$double.eps)^(-1)

#' p0pm2 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-2)
#'
#' @export p0pm2
p0pm2 <- function(x) p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-2)

#' p0pm05 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*sign(x)*(abs(x)+.Machine$double.eps)^(-0.5)
#'
#' @export p0pm05
p0pm05 <- function(x) p0(x)*(abs(x)+.Machine$double.eps)^(-0.5)

#' p0p05 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*(abs(x)+.Machine$double.eps)^(0.5)
#'
#' @export p0p05
p0p05 <- function(x) p0(x)*(abs(x)+.Machine$double.eps)^(0.5)

#' p0p1 polynomial term
#'
#' @param x The vector of values
#' @return  p0(x)*x
#'
#' @export p0p1
p0p1 <- function(x) p0(x)*x

#' p0p2 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*x^(2)
#'
#' @export p0p2
p0p2 <- function(x) p0(x)*x^(2)

#' p0p3 polynomial term
#'
#' @param x The vector of values
#' @return p0(x)*x^(3)
#'
#' @export p0p3
p0p3 <- function(x) p0(x)*x^(3)


#' ReLu function
#'
#' @param x The vector of values
#' @return max(x,0)
#'
#' @export relu
relu <- function(x) max(x,0)

#' negative ReLu function
#'
#' @param x The vector of values
#' @return max(x,0)
#'
#' @export nrelu
nrelu <- function(x) max(x,0)


#' negative ReLu function
#'
#' @param x The vector of values
#' @return max(-x,0)
#'
#' @export nrelu
nrelu <- function(x) max(-x,0)

#' GELU function
#'
#' @param x The vector of values
#' @return x*pnorm(x)
#'
#' @export gelu
gelu <- function(x)x *pnorm(x)


#' Negative GELU function
#'
#' @param x The vector of values
#' @return -x*pnorm(-x)
#'
#' @export ngelu
ngelu <- function(x) -x*pnorm(-x)

#' erf function
#'
#' @param x The vector of values
#' @return 2 * pnorm(x * sqrt(2)) - 1
#'
#' @export erf
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


#' heavy side function
#'
#' @param x The vector of values
#' @return as.integer(x>0)
#'
#' @export hs
hs <- function(x) as.integer(x>0)


#' negative heavy side function
#'
#' @param x The vector of values
#' @return as.integer(x<0)
#'
#' @export nhs
nhs <- function(x) as.integer(x<0)

#' not x
#'
#' @param x The vector of binary values
#' @return 1-x
#'
#' @export not
not <- function(x) (1-x)

