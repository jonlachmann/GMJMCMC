#######################################################
#
# Example 9 (Section 6.1):
#
# Logistic Regression
#
# This is the valid version for the JSS Paper
#
#######################################################

library(devtools)
devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
use.fbms = FALSE  

#setwd("/home/florian/FBMS/")

df = read.csv2(file = "/Users/aliaksandrhome/EMJMCMC/supplementaries/BGNLM/spam/spam.data",sep = " ",dec = ".")[,c(58,1:57)]

summary(df)


#number of observations in the data

n = dim(df)[1] 

#number of covariates

p = dim(df)[2] - 1   


colnames(df) =  c("y", paste0("x",1:p))


to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,1,1) 

params <- gen.params.gmjmcmc(df)
params$feat$check.col <- F
####################################################
#
# single thread analysis
#
####################################################


set.seed(6001)

df$y <- df$x1 + df$x10*df$x11 + df$x5 + 1

df$y <- as.numeric(df$y > mean(df$y))

params$loglik = list(r = 1/dim(df)[1], family = "binomial", betaprior = g.prior(100), laplace = F)

result2 <- fbms(data = df, method = "gmjmcmc", family = "custom",
               transforms = transforms, probs = probs, loglik.pi = glm.logpost.bas, params = params,P=3)

# Perform analysis with logistic.loglik
if (use.fbms) {
  result <- fbms(data = df, method = "gmjmcmc", family = "binomial",
                 transforms = transforms, probs = probs, params = params)
} else {
  result <- gmjmcmc(df, logistic.loglik, transforms 
                  = transforms, probs = probs, params = params)
}
# Default tuning parameters for logistic regression:
#
# params$loglik$r = 1/n

summary(result)


# IMPORTANT: specify correct link function for predict

pred <- predict(result, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred$aggr$mean)==df$y)

bm <- get.best.model(result = result)
preds <-  predict(object = bm, df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)

mpm <- get.mpm.model(result = result,family = "binomial",y = df$y,x=df[,-1])
preds <-  predict(mpm, df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)



plot(pred$aggr$mean)
points(pred$aggr$quantiles[1,], col = 2)
points(pred$aggr$quantiles[3,], col = 3)


head(cbind(pred$aggr$mean, pred$aggr$quantiles[1,],pred$aggr$quantiles[3,]))


####################################################
#
# multiple thread analysis
#
####################################################

set.seed(6002)

if (use.fbms) {
  result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", family = "binomial",
                          runs = 40, cores = 40, transforms = transforms, 
                          probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 40, cores = 40, data = df, 
                                      loglik.pi = logistic.loglik, transforms = transforms, 
                                      probs = probs, params = params, P=25)
}
summary(result_parallel)

# IMPORTANT: specify correct link function for predict

pred_parallel = predict(result_parallel, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred_parallel$aggr$mean)==df$y)

preds <-  predict(get.best.model(result_parallel), df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)


preds <-  predict(get.mpm.model(result = result_parallel,family = "binomial",y = df$y,x=df[,-1]), df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)

