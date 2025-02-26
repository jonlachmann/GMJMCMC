#######################################################
#
# Example 8 (Section 5.2):
#
# Logic regression with a different model prior
#
#  DATA - simulated 
#
#
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
use.fbms = FALSE  

set.seed(1)
X2 <- as.data.frame(array(data = rbinom(n = 50*2000,size = 1,prob = runif(n = 50*2000,0,1)),dim = c(2000,50)))
Y2 <- rnorm(n = 2000,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*(X2$V37),sd = 1)
df <- data.frame(Y2,X2)
summary(df)

str(df)

#number of observations in the data

n = dim(df)[1] 



# remove variable sex because gmjmcmc cannot handle factor variables
df.training <- df[1:1000,]
df.test <- df[1001:n,]
df.test$Mean <- (1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*(X2$V37))[1001:n]

#FBMS unlike EMJMCMC package does not explicitly have GMJMCMC for logic regression, but we can easily run it without an 
# "or" operator as "and" and "not" allow to automatically handle "or" through de Morgan law.
transforms <- c("not")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) #No projections allowed
probs$filter <- 0.6

params <- gen.params.gmjmcmc(df.training)
params$loglik$m <- 50
params$loglik$n <- n
params$loglik$r <- 1
params$feat$pop.max <- 31
params$large$neigh.max <- 21
params$large$neigh.size <- 10
##############################################
#
#  New function for model prior
#
#

estimate.logic.lm = function(y, x, model, complex, params)
{
  
  if (length(params) == 0) 
    params <- list(r = 1/dim(x)[1]) 
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  sj <- complex$width
 
  lp <- sum(log(factorial(sj))) - sum(sj*log(params$m) + (2*sj-2)*log(2))
  
  #print(lp)
  
  mloglik <- -(mod$aic + (log(length(y))-2) * (mod$rank))/2 
  
  logpost <- mloglik + lp + params$n #just add a constant to all models
  
  if(logpost==-Inf)
    logpost = -10000
  
  return(list(crit = logpost + lp, coefs = mod$coefficients))
}



#############################################################################
#
#   Logic regression training
#
#############################################################################

set.seed(5001)


#result <- gmjmcmc(df.training, logistic.loglik, logistic.loglik.alpha, transforms 
#                  = transforms, probs = probs, params = params)

if (use.fbms) {
  result <- fbms(data = df.training, family = "custom", loglik.pi = estimate.logic.lm,N.init = 500,N.final = 500, P = 25,
                 method = "gmjmcmc", transforms = transforms, 
                 probs = probs, params = params)
} else {
  #  result <- gmjmcmc(df.training, transforms = transforms, probs = probs)
  
  result <- gmjmcmc(df.training, loglik.pi = estimate.logic.lm,N.init = 500,N.final = 500, , P = 25,
                    transforms = transforms,params = params, probs = probs)
  
}
summary(result)
mpm <- get.mpm.model(result,y = df.training$Y2,x = df.training[,-1],family = "custom", loglik.pi = estimate.logic.lm,params = params$loglik)
mbest <- get.best.model(result)


pred <- predict(result, x =  df.test[,-1], link = function(x)(x))  
pred_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_best - df.test$Y2)^2))
sqrt(mean((pred_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_best - df.test$Mean)^2))
sqrt(mean((pred_mpm - df.test$Mean)^2))



plot(pred$aggr$mean, df.test$Y2)
points(pred$aggr$mean,df.test$Mean,col = 2)
points(pred_best,df.test$Mean,col = 3)
points(pred_mpm,df.test$Mean,col = 4)




#############################################################################
#
#   Parallel version just 16 chains on 8 cores
#
#############################################################################


set.seed(5002)

if (use.fbms) {
  result_parallel <- fbms(data = df.training, family = "custom", loglik.pi = estimate.logic.lm,N.init = 500,N.final = 500,
                          method = "gmjmcmc.parallel", runs = 16, cores = 8,
                          transforms = transforms, probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 16, cores = 8, data = df.training, 
                                      loglik.pi = estimate.logic.lm,N.init = 500,N.final = 500,
                                      transforms = transforms, probs = probs, params = params, P=25)
}
summary(result_parallel)
mpm <- get.mpm.model(result_parallel,y = df.training$Y2,x = df.training[,-1],family = "custom", loglik.pi = estimate.logic.lm,params = params$loglik)
mbest <- get.best.model(result_parallel)


pred_parallel <- predict(result_parallel, x =  df.test[,-1], link = function(x)(x))  
pred_par_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_par_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred_parallel$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_par_best - df.test$Y2)^2))
sqrt(mean((pred_par_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred_parallel$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_par_best - df.test$Mean)^2))
sqrt(mean((pred_par_mpm - df.test$Mean)^2))



plot(pred_parallel$aggr$mean, df.test$Y2)
points(pred_parallel$aggr$mean,df.test$Mean,col = 2)
points(pred_par_best,df.test$Mean,col = 3)
points(pred_par_mpm,df.test$Mean,col = 4)

