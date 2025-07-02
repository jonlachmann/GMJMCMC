#######################################################
#
# Example 2 (Section 4.1):
#
# Simulated data without any nonlinearities, only using fbms function
#
# This is the valid version for the JSS Paper
#
#######################################################


library(mvtnorm)
library(FBMS)



n <- 100  # sample size
p <- 20   # number of covariates
p.vec <- 1:p


k <- 5    #size of the data generating model


correct.model <- 1:k
beta.k <- (1:5)/5   # Coefficents of the correct submodel

beta <- c(rep(0, p))
beta[correct.model] <- beta.k

set.seed(123)

x <- rmvnorm(n, rep(0, p))
y <- x %*% beta    + rnorm(n)
X <- as.matrix(x)

y<-scale(y)
X<-scale(X)/sqrt(n)


df <- as.data.frame(cbind(y, X))
colnames(df) <- c("Y", paste0("X", seq_len(ncol(df) - 1)))

correct.model
beta.k


to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")

set.seed(1)
  result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, P = 40)
 
  summary(result)
 
set.seed(2)
  result2 <- result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                            N = 1000, P = 40)
summary(result2)




########################################################
#
#  Model which includes no non-linear effects
#
#


set.seed(1)
  result.lin <- fbms(data = df, N = 5000)
  summary(result.lin)

  set.seed(2)
  result.lin <- fbms(data = df, N = 1000)
  summary(result.lin)
  
  plot(result.lin)
 
  


#The default value of  N = 100 by accident yields the correct results 
set.seed(3)

  result.lindef <- fbms(data = df)

plot(result.lindef)
summary(result.lindef)
