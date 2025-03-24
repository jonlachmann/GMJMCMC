#######################################################
#
# Example 8 (Section 5.2):
#
# Prediction using non-linear Projections with different model prior
#
#  DATA - abalone data set
#
#  Data set is available at https://www.kaggle.com/datasets/rodolfomendes/abalone-dataset
#
#  For convenience we provide the file abalone.csv which contains already the names of the variables
#
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
use.fbms = FALSE  

data("abalone")

df = abalone
df$Sex_F_vs_I = as.numeric(df$Sex == "F")
df$Sex_M_vs_I = as.numeric(df$Sex == "M")
df$Sex = as.factor(df$Sex)
df$Rings = as.numeric(df$Rings)

summary(df)


#number of observations in the data

n = dim(df)[1] 

# Create Training and Test dataset

# Sam Waugh (1995) "Extending and benchmarking Cascade-Correlation", PhD
# thesis, Computer Science Department, University of Tasmania.

#-- Test set performance (final 1044 examples, first 3133 used for training):
#  24.86% Cascade-Correlation (no hidden nodes)
#  26.25% Cascade-Correlation (5 hidden nodes)
#  21.5%  C4.5
#  0.0%  Linear Discriminate Analysis
#  3.57% k=5 Nearest Neighbour
#  (Problem encoded as a classification task)


# remove variable sex because gmjmcmc cannot handle factor variables
df.training = df[1:3133,-2]
df.test = df[3134:n,-2]


summary(df.training)


pred.RMSE = rep(0,5)   # Collect the results of prediction RMSE from the five different methods



transforms <- c("sigmoid")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(0,0,1,1) #Only projections!

params <- gen.params.gmjmcmc(df.training)
params$mlpost$r = 0.2

##############################################
#
#  New function for model prior
#
#

gaussian.loglik.gap.open.prior = function(y, x, model, complex, params){
  {
    if (length(params) == 0) 
      params <- list(r = 1/dim(x)[1]) #list(r = 0.5)
    suppressWarnings({
      mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
    })

    # logarithm of marginal likelihood
    mloglik <- -mod$deviance/2 +  log(length(y)) * (mod$rank - 1) 
              
   # logarithm of model prior
#    lp <-  log(params$r) * (sum(complex$depth) + 0.1 * sum(complex$oc))
    lp <-  log(params$r) * (0.9*sum(complex$depth) + 0.1 * sum(complex$width))
    return(list(crit = mloglik + lp, coefs = mod$coefficients))
  }

}




#############################################################################
#
#   Using method 0 for alpha (simply set to 1, default)
#
#############################################################################

set.seed(5001)


#result <- gmjmcmc(df.training, logistic.loglik, logistic.loglik.alpha, transforms 
#                  = transforms, probs = probs, params = params)

if (use.fbms) {
  result <- fbms(data = df.training, family = "custom", loglik.pi = gaussian.loglik.gap.open.prior,
                 method = "gmjmcmc", transforms = transforms, 
                 probs = probs, params = params)
} else {
#  result <- gmjmcmc(df.training, transforms = transforms, probs = probs)

  result <- gmjmcmc(df.training, loglik.pi = gaussian.loglik.gap.open.prior, 
                        transforms = transforms, probs = probs)
  
  }
summary(result)

pred = predict(result, x =  df.test[,-1], link = function(x)(x))  

pred.RMSE[1] = sqrt(mean((pred$aggr$mean - df.test$Rings)^2))

plot(pred$aggr$mean, df.test$Rings)





#############################################################################
#
#   Parallel version
#
#############################################################################


set.seed(5002)

if (use.fbms) {
  result_parallel <- fbms(data = df.training, family = "custom", loglik.pi = gaussian.loglik.gap.open.prior, 
                          loglik.alpha = gaussian.loglik.alpha,
                          method = "gmjmcmc.parallel", runs = 40, cores = 40,
                          transforms = transforms, probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 40, cores = 40, data = df.training, 
                                   loglik.pi = gaussian.loglik.gap.open.prior, 
                                   loglik.alpha = gaussian.loglik.alpha,
                                    transforms = transforms, probs = probs, params = params, P=25)
}
summary(result_parallel)



pred_parallel = predict(result_parallel, x =  df.test[,-1], link = function(x)(x))  

pred.RMSE[2] = sqrt(mean((pred_parallel$aggr$mean - df.test$Rings)^2))

plot(pred_parallel$aggr$mean, df.test$Rings)
abline(0,1)



#############################################################################
#
#   Using method 3 to estimate alpha
#
#############################################################################

params$feat$alpha = 3


set.seed(5003)


if (use.fbms) {
  result.a3 <- fbms(data = df.training, family = "custom", loglik.pi = gaussian.loglik.gap.open.prior, 
                    method = "gmjmcmc", transforms = transforms, 
                 probs = probs, params = params)
} else {
  result.a3 <- gmjmcmc(df.training, transforms = transforms, probs = probs, params = params)
}
summary(result.a3)



pred.a3 = predict(result.a3, x =  df.test[,-1], link = function(x)(x))  

pred.RMSE[3] = sqrt(mean((pred.a3$aggr$mean - df.test$Rings)^2))

plot(pred.a3$aggr$mean, df.test$Rings)





#############################################################################
#
#   Parallel version  params$feat$alpha = 3
#
#############################################################################

params$feat$alpha = 3

set.seed(5004)

if (use.fbms) {
  result_parallel.a3 <- fbms(data = df.training, family = "custom", loglik.pi = gaussian.loglik.gap.open.prior, 
                             loglik.alpha = gaussian.loglik.alpha,
                             method = "gmjmcmc.parallel", runs = 40, cores = 40,
                          transforms = transforms, probs = probs, params = params, P=25)
} else {
  result_parallel.a3 =  gmjmcmc.parallel(runs = 40, cores = 40, data = df.training, 
                                    loglik.pi = gaussian.loglik.gap.open.prior, 
                                    loglik.alpha = gaussian.loglik.alpha,
                                    transforms = transforms, probs = probs, params = params, P=25)
}
summary(result_parallel.a3)



pred_parallel.a3 = predict(result_parallel.a3, x =  df.test[,-1], link = function(x)(x))  

pred.RMSE[4] = sqrt(mean((pred_parallel.a3$aggr$mean - df.test$Rings)^2))

plot(pred_parallel.a3$aggr$mean, df.test$Rings)
abline(0,1)







#############################################################################
#
#   Parallel version with fractional polynomials
#
#############################################################################

transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(0,1,0,1) #Only modifications!

set.seed(50005)

if (use.fbms) {
  result.fp <- fbms(data = df.training, method = "gmjmcmc.parallel", runs = 40, cores = 40,
                             transforms = transforms, probs = probs, params = params, P=25)
} else {
  result.fp <- gmjmcmc.parallel(runs = 40, cores = 40, data = df.training, 
                              loglik.pi =gaussian.loglik,loglik.alpha = gaussian.loglik.alpha, 
                              transforms = transforms, probs = probs, params = params, P=25)
}
summary(result.fp)



pred_fp = predict(result.fp, x =  df.test[,-1], link = function(x)(x))  

pred.RMSE[5] = sqrt(mean((pred_fp$aggr$mean - df.test$Rings)^2))

plot(pred_fp$aggr$mean, df.test$Rings)



pred.RMSE
