#######################################################
#
# Example 4 (Section 4.2):
#
# Simulated data with interactions
#
# This is the valid version for the JSS Paper
#
#######################################################

library(mvtnorm)
library(FBMS)
use.fbms <- TRUE  
stronger.singal <- FALSE

n <- 100*ifelse(stronger.singal,10,1)  # sample size
p <- 20   # number of covariates

# Model:  
# X1: Pure Main effect
# X2 : X3: Pure interaction effect
# X4 * X5: Main effects plus interaction effect


set.seed(1003)

x = rmvnorm(n, rep(0, p))
X <- as.matrix(x)
X <- scale(X)/sqrt(n)

y <- (1.2 * x[,1] + 1.5 * x[,2]* x[,3] - x[,4] + 1.1*x[,5] - 1.3 * x[,4]*x[,5])+ rnorm(n)
y<-scale(y)

df <- data.frame(y = y, X)


transforms <- c("")
params <- gen.params.gmjmcmc(ncol(df) - 1)
#params$mlpost$var = "unknown" #this will set the variance to unknwon
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,0,0,1)            #Include interactions and mutations

####################################################
#
# single thread analysis (two different runs)
#
####################################################

set.seed(123)
if (use.fbms) {
  result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, beta_prior = list(type = "EB-local"),
                  probs = probs, params = params, P=40)
} else {
  result <- gmjmcmc(x = df[, -1], y = df[, 1], mlpost_params = list(family = "gaussian", beta_prior = list(type = "EB-local")), transforms = transforms, params = params, probs = probs, P=40)
}
summary(result)


set.seed(123)
if (use.fbms) {
  result2 <- fbms(data = df, method = "gmjmcmc", transforms = transforms, beta_prior = list(type = "EB-local"),
                 probs = probs, params = params, P=40)
} else {
  result2 <- gmjmcmc(x = df[, -1], y = df[, 1], mlpost_params = list(family = "gaussian", beta_prior = list(type = "EB-local")), transforms = transforms, N = 1000, N.final = 5000,
                     probs = probs, params = params,  P = 40)
}
summary(result2, tol = 0.01)


####################################################
#
# multiple thread analysis
#
####################################################


set.seed(123)

if (use.fbms) {
  result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms,beta_prior = list(type = "EB-local"),
                 runs = 40, cores = 10,
                 probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 40, cores = 10, x = df[, -1], y = df[, 1],
                            transforms = transforms, probs = probs, params = params, P=25)
}

summary(result_parallel, tol = 0.01)



# Using longer more iterations of MJMCMC chains does not change results substantially
set.seed(123)

if (use.fbms) {
  result_parallel2 <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms,beta_prior = list(type = "EB-local"),
                 runs = 40, cores = 10, N=1000, N.final=2000,
                 probs = probs, params = params, P=25)
} else {
  result_parallel2 =  gmjmcmc.parallel(runs = 40, cores = 10, x = df[, -1], y = df[, 1],
                                transforms = transforms,mlpost_params = list(family = "gaussian", beta_prior = list(type = "EB-local")), probs = probs, params = params, P=25, 
                                N=1000, N.final=2000)
}
summary(result_parallel2, tol = 0.01)




########################################################
#
#  Model which includes no interactions effects
#
#


set.seed(123)


if (use.fbms) {
  result.lin <- fbms(data = df, beta_prior = list(type = "EB-local"), N = 5000)
} else {
  result.lin <- mjmcmc(x = df[, -1], y = df[, 1],mlpost_params = list(family = "gaussian", beta_prior = list(type = "EB-local")), N = 5000)
}

plot(result.lin)
summary(result.lin)



