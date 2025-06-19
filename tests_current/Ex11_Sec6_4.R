#######################################################
#
# Example 11 (Section 6.4):
#
# Subsampling
#
# Heart Disease Health Indicators Dataset”
#
# This is the valid version for the JSS Paper
#
#######################################################

library(tictoc)
library(FBMS)
#library(devtools)
#devtools::install_github("jonlachmann/irls.sgd", force=T, build_vignettes=F)
library(irls.sgd)

use.fbms <- T

df <- read.csv2(file = "/Users/aliaksandrhome/GMJMCMC/tests/heart_disease_health_indicators_BRFSS2015.csv",sep = ",",dec = ".")

summary(df)
dim(df)

#number of observations in the data

n <- dim(df)[1] 

#number of covariates

p <- dim(df)[2] - 1   

params <- gen.params.gmjmcmc(ncol(df) - 1)
transforms <- c("sigmoid","pm1","p0","p05","p2","p3")
probs <- gen.probs.gmjmcmc(transforms)

logistic.posterior.bic.irlssgd <- function (y, x, model, complex, mlpost_params) 
{
  if (!is.null(mlpost_params$crit)) {
    mod <- glm.sgd(x[,model], y, binomial(), sgd.ctrl = list(start=mlpost_params$coefs, subs=mlpost_params$subs, maxit=10, alpha=0.00008, decay=0.99, histfreq=10))
    mod$deviance <- get_deviance(mod$coefficients, x[,model], y, binomial())
    mod$rank <- length(mod$coefficients)
  } else {
    mod <- irls.sgd(as.matrix(x[,model]), y, binomial(),
                  irls.control=list(subs=mlpost_params$subs, maxit=20, tol=1e-7, cooling = c(1,0.9,0.75), expl = c(3,1.5,1)),
                  sgd.control=list(subs=mlpost_params$subs, maxit=250, alpha=0.001, decay=0.99, histfreq=10))
  }
  
  # logarithm of marginal likelihood
  mloglik <- -mod$deviance / 2 - 0.5 * log(length(y)) * (mod$rank - 1)
    
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  crit <- mloglik + lp

  if (!is.null(mlpost_params$crit) && mlpost_params$crit > crit) {
    return(list(crit = mlpost_params$crit, coefs = mlpost_params$coefs))
  }
  
  return(list(crit = crit, coefs = mod$coefficients))
}



############################
#
# Testing runtime
#
############################

set.seed(100001)
tic()
if (use.fbms) {
  result1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "custom", loglik.pi = logistic.posterior.bic.irlssgd, method = "gmjmcmc", 
                   model_prior = list(r = 0.5, subs = 0.01), transforms = transforms,
                   params = params, P = 2, sub  = T)
} else { 
  result1 <- gmjmcmc(x = df[, -1], y = df[, 1], loglik.pi = logistic.posterior.bic.irlssgd,
                  mlpost_params = list(r = 0.5, subs = 0.01), transforms = transforms,
                  params = params, P = 2, sub  = T) 
}
time1 <- toc() 

set.seed(100002)
tic()
# regular analysis

if (use.fbms) {
  result2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "binomial", method = "gmjmcmc", 
               model_prior = list(r = 0.5, subs = 0.01), 
               beta_prior = list(type = "Jeffreys-BIC"),
               transforms = transforms, params = params, P = 2, sub  = T)
} else { 
  result2 <- gmjmcmc(x = df[, -1], y = df[, 1],
               mlpost_params = list(r = 0.5, family = "binomial", 
               beta_prior = list(type = "Jeffreys-BIC")), 
               params = params, P = 2)
}

time2 <- toc() 

c(time1, time2)

############################
#
# More serious analysis
#
############################


# with subsampling

set.seed(100003)

tic()
if (use.fbms) {
  result_parallel_1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "custom", loglik.pi = logistic.posterior.bic.irlssgd, 
              method = "gmjmcmc.parallel", runs = 10, cores = 10,
              model_prior = list(r = 0.5, subs = 0.01), transforms = transforms,
              params = params, P = 2, sub  = T)
} else { 
  result_parallel_1 <- gmjmcmc.parallel(runs = 10, cores = 10, x = df[, -1], y = df[, 1],
              loglik.pi = logistic.posterior.bic.irlssgd,
              mlpost_params = list(r = 0.5, subs = 0.01), 
              transforms = transforms, params = params,  P = 2, sub = T)
}
time3 <- toc() 

summary(result_parallel_1)

# without subsampling


set.seed(100004)

tic()
if (use.fbms) {
  result_parallel_2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "binomial", 
                            method = "gmjmcmc.parallel", runs = 10, cores = 10,
                            model_prior = list(r = 0.5),
                            beta_prior = list(type = "Jeffreys-BIC"), 
                            transforms = transforms,
                            params = params, P = 2, sub  = T)
} else { 
  result_parallel_2 <- gmjmcmc.parallel(runs = 10,cores = 10, x = df[, -1], y = df[, 1],
                            mlpost_params = list(r = 0.5, family = "binomial", 
                            beta_prior = list(type = "Jeffreys-BIC")), 
                            transforms = transforms, params = params,  P = 2)
}
time4 <- toc() 

summary(result_parallel_2)



############################
#
# Probably too ambitious analysis
#
############################

# with subsampling

set.seed(100005)
tic()
if (use.fbms) {
  result_parallel_long_1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "custom", loglik.pi = logistic.posterior.bic.irlssgd, 
                            method = "gmjmcmc.parallel", runs = 40, cores = 40,
                            model_prior = list(r = 0.5, subs = 0.01), transforms = transforms,
                            params = params, P = 10, sub  = T)
} else { 
  result_parallel_long_1 <- gmjmcmc.parallel(formula = HeartDiseaseorAttack ~ 1 + ., runs = 40, cores = 40, x = df[, -1], y = df[, 1],
                                        loglik.pi = logistic.posterior.bic.irlssgd,
                                        mlpost_params = list(r = 0.5, subs = 0.01), 
                                        transforms = transforms, params = params,  P = 10, sub = T)
}

time5 <- toc() 
summary(result2)



# regular analysis


set.seed(100006)

tic()
if (use.fbms) {
  result_parallel_long_2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, family = "binomial", 
                            method = "gmjmcmc.parallel", runs = 40, cores = 40,
                            model_prior = list(r = 0.5),
                            beta_prior = list(type = "Jeffreys-BIC"), 
                            transforms = transforms,
                            params = params, P = 10, sub  = T)
} else { 
  result_parallel_long_2 <- gmjmcmc.parallel(runs = 40,cores = 40, x = df[, -1], y = df[, 1],
                                        mlpost_params = list(r = 0.5, family = "binomial", 
                                                             beta_prior = list(type = "Jeffreys-BIC")), 
                                        transforms = transforms, params = params,  P = 10)
}
time6 <- toc() 


summary(result_parallel_long_2)


############################################################################

C = cor(df, use = "everything",
        method = "spearman")

corrplot::corrplot(C)

apply((abs(C - diag(diag(C)))), 2, max)

save.image("Ex11_Results.RData")

