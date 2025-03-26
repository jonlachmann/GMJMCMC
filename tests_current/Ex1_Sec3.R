#################################################
#
# Example 1:
#
# Kepler Example with the most recent database update
#
# This is the valid version for the JSS paper
#
##################################################

#install.packages("FBMS")
#install.packages("devtools")
#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)
library(FBMS)

data <- read.csv("https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/master/comma_separated/open_exoplanet_catalogue.txt")
data <- na.omit(data[,c("semimajoraxis","mass","radius","period","eccentricity","hoststar_mass","hoststar_radius","hoststar_metallicity","hoststar_temperature","binaryflag")]) 
summary(data)


te.ind <- 540:939
df.train = data[-te.ind,]
df.test = data[te.ind,]
n <- dim(df.train)[1]

to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")

# Logical to decide whether to perform analysis with fbms function
# If FALSE then gmjmcmc or gmjmcmc.parallel function is used
use.fbms = FALSE  

####################################################
#
# single thread analysis (default values, Section 3.1)
#
####################################################

params <- gen.params.gmjmcmc(df.train)
params$loglik$var <- "unknown"

params <- gen.params.gmjmcmc(df.train)
params$loglik$prior_beta <- "Jeffreys-BIC"
params$loglik$var <- 1
params$loglik$family <- "gaussian"
params$loglik$r <- 1/n
result <- gmjmcmc.parallel(data = df.train,loglik.pi = fbms.mlik.master, transforms = transforms, params = params,cores = 10,runs = 40, P = 25)
summary(result)


params <- gen.params.gmjmcmc(df.train)
params$loglik$prior_beta <- "g-prior"
params$loglik$g <- n
params$loglik$family <- "gaussian"
params$loglik$r <- 1/n
result <- gmjmcmc.parallel(data = df.train,loglik.pi = fbms.mlik.master, transforms = transforms, params = params,cores = 10,runs = 40, P = 25)
summary(result)


params <- gen.params.gmjmcmc(df.train)
params$loglik$prior_beta <- "EB-local"
params$loglik$family <- "gaussian"
params$loglik$r <- 1/n
result <- gmjmcmc.parallel(data = df.train,loglik.pi = fbms.mlik.master, transforms = transforms, params = params,cores = 10,runs = 40, P = 25)
summary(result)

params <- gen.params.gmjmcmc(df.train)
params$loglik$prior_beta <- "hyper-g"
params$loglik$family <- "gaussian"
params$loglik$a <- 3
params$loglik$r <- 1/n
result <- gmjmcmc.parallel(data = df.train,loglik.pi = fbms.mlik.master, transforms = transforms, params = params,cores = 10,runs = 40, P = 25)
summary(result)


if (use.fbms) {
 result.default <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, method = "gmjmcmc", transforms = transforms, params = params)
} else {
 result.default <- gmjmcmc(df.train, transforms = transforms, params = params)
}
summary(result.default,labels = F)


preds <- predict(result.default, df.test[,-1], link = function(x) x)
sqrt(mean((preds$aggr$mean - df.test$semimajoraxis)^2))

#new additional ways to predict using MPM and best model
get.best.model(result = result.default)
preds <- predict(get.best.model(result.default), df.test[,-1])
sqrt(mean((preds - df.test$semimajoraxis)^2))

get.mpm.model(result = result.default,y = df.test$semimajoraxis,x=df.test[,-1])
preds <- predict(get.mpm.model(result.default,y = df.test$semimajoraxis,x=df.test[,-1]), df.test[,-1])
sqrt(mean((preds - df.test$semimajoraxis)^2))

####################################################
#
# single thread analysis (more iterations, Section 3.2)
#
####################################################


set.seed(123)

if (use.fbms) {
 result.P50 <- fbms(data = df.train, method = "gmjmcmc", transforms = transforms,
                    P=50, N.init=1000, N.final=1000, params = params)
} else {
 result.P50 <- gmjmcmc(df.train,  transforms = transforms,
                       P=50, N.init=1000, N.final=1000, params = params)
}
summary(result.P50, labels = names(df.train)[-1])

####################################################
#
# multiple thread analysis (Section 3.3)
#
####################################################

set.seed(123)
if (use.fbms) {
 result_parallel <- fbms(data = df.train, method = "gmjmcmc.parallel", transforms = transforms,
                         runs = 40, cores = 10, P=25,params = params)
} else {
 result_parallel <- gmjmcmc.parallel(runs = 40, cores = 10, data = df.train, loglik.pi = gaussian.loglik, 
                                     transforms = transforms, P=25,params = params)
}
summary(result_parallel, tol = 0.01)


####### fixed variance
params$loglik$var <- 1
set.seed(124)
if (use.fbms) {
  result_parallel_unitphi <- fbms(data = df.train, method = "gmjmcmc.parallel", transforms = transforms,
                                  runs = 40, cores = 10, P=25,params = params)
} else {
  result_parallel_unitphi <- gmjmcmc.parallel(runs = 40, cores = 10, data = df.train, loglik.pi = gaussian.loglik, 
                                              transforms = transforms, P=25,params = params)
}
summary(result_parallel_unitphi, tol = 0.01)


#g prior with g = n is perfect 
gaussian.loglik.g <- function (y, x, model, complex, params)
{
  
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  # Calculate R-squared
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  Rsquare <- 1 - (RSS / TSS)
  
  # logarithm of marginal likelihood
  mloglik <- 0.5*(log(1.0 + params$g) * (dim(x)[1] - mod$rank)  - log(1.0 + params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1))*(mod$rank!=1)
  
  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(params, complex)
  
  return(list(crit = mloglik + lp, coefs = mod$coefficients))
}


#default for N.final = N.init
params$loglik$g = n
set.seed(1)
if (use.fbms) {
  result_parallel_g <- fbms(data = df.train,family = "custom", method = "gmjmcmc.parallel",loglik.pi = gaussian.loglik.g, transforms = transforms,
                          runs = 40, cores = 10, P=25,params = params)
} else {
  result_parallel_g <- gmjmcmc.parallel(runs = 40, cores = 10, data = df.train, loglik.pi = gaussian.loglik.g, 
                                      transforms = transforms, P=25,params = params)
}
summary(result_parallel_g, tol = 0.01)
####################################################
#
# Inspection of Results (Section 3.4)
#
####################################################

######################
# summary

summary(result.default)
summary(result.default, labels = names(df.train)[-1])

summary(result.P50)

summary(result_parallel)
summary(result_parallel, tol = 0.01,labels = names(df.train)[-1])


######################
# plot

pdf.train("result.pdf.train") 
plot(result.default)
dev.off()

plot(result.default)



pdf.train("result.P50.pdf.train") 
plot(result.P50)
dev.off()

plot(result.P50)



pdf.train("result_parallel.pdf.train") 
plot(result_parallel)
dev.off()

plot(result_parallel)
plot(result_parallel, 12)

pdf.train("result_parallel_unitphi.pdf.train") 
plot(result_parallel_unitphi)
dev.off()

plot(result_parallel_unitphi)
plot(result_parallel_unitphi, 12)


######################
# Prediction

#preds <-  predict(result, df.test[,-1], link = function(x) x)  
preds <-  predict(result.default, df.test[,-1])

pdf.train("prediction.pdf.train") 
plot(preds$aggr$mean, df.test$semimajoraxis)
dev.off()

plot(preds$aggr$mean, df.test$semimajoraxis)

rmse.default <- sqrt(mean((preds$aggr$mean - df.test$semimajoraxis)^2))

###############################


#preds.P50 = predict(result.P50, df.test[,-1], link = function(x) x)  
preds.P50 = predict(result.P50, df.test[,-1])  

pdf.train("prediction.P50.pdf.train") 
plot(preds.P50$aggr$mean, df.test$semimajoraxis)
dev.off()

plot(preds.P50$aggr$mean, df.test$semimajoraxis)

rmse.P50 <-  sqrt(mean((preds.P50$aggr$mean - df.test$semimajoraxis)^2))


###############################


preds.multi <- predict(result_parallel , df.test[,-1], link = function(x) x)

pdf.train("pred_parallel.pdf.train") 
plot(preds.multi$aggr$mean, df.test$semimajoraxis)
dev.off()

rmse.parallel <- sqrt(mean((preds.multi$aggr$mean - df.test$semimajoraxis)^2))


###############################


preds_unitphi <- predict(result_parallel_unitphi , df.test[,-1], link = function(x) x)

pdf.train("pred_parallel.pdf.train") 
plot(preds_unitphi$aggr$mean, df.test$semimajoraxis)
dev.off()

rmse_unitphi <- sqrt(mean((preds_unitphi$aggr$mean - df.test$semimajoraxis)^2))


###############################


preds_g <- predict(result_parallel_g , df.test[,-1], link = function(x) x)

pdf.train("pred_parallel.pdf.train") 
plot(preds_g$aggr$mean, df.test$semimajoraxis)
dev.off()

rmse_g <- sqrt(mean((preds_g$aggr$mean - df.test$semimajoraxis)^2))

c(rmse.default, rmse.P50, rmse.parallel,rmse_unitphi,rmse_g)



#let us test all priors from BAS ,see prior in ?bas.lm 

library(tictoc)
#just testing all priors I now added, time, etc.
for(prior in c("g-prior",
               "hyper-g",
               "hyper-g-laplace",
               "hyper-g-n",
               "AIC",
               "BIC",
               "ZS-null",
               "ZS-full",
               "EB-local",
               "EB-global",
               "JZS"))
{
  print(paste0("testing ",prior))
  params$loglik <- list(r =  1/dim(df.train)[1], betaprior = prior,alpha = max(dim(df.train)[1],(dim(df.train)[2])^2))
  
  
  #ours are stil a bit faster than the BAS ones, but BAS are relatively fine too
  
  tic()
  result.default <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, method = "gmjmcmc.parallel",cores = 10, runs = 10, transforms = transforms, loglik.pi = lm.logpost.bas, params = params, P = 50)
  time.res = toc()
  preds <- predict(result.default, df.test[,-1], link = function(x) x)
  print(summary(result.default))
  print(sqrt(mean((preds$aggr$mean - df.test$semimajoraxis)^2)))
  print(time.res)
}


#default for N.final = N.init
params <- gen.params.gmjmcmc(df.train)
params$loglik$g <- dim(df.train)[1]
tic()
result.default <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, method = "gmjmcmc.parallel",cores = 10, runs = 10, transforms = transforms, loglik.pi = gaussian.loglik.g, params = params, P = 50)
time.res = toc()
preds <- predict(result.default, df.test[,-1], link = function(x) x)
print(summary(result.default))
print(sqrt(mean((preds$aggr$mean - df.test$semimajoraxis)^2)))
print(time.res)



#testing a bit BAS based stuff vs our implementation, g prior
lm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 10),params = list(r =  1/dim(df.train)[1], betaprior = "g-prior",alpha = min(dim(df.train)[1],(dim(df.train)[2])^2)))
gaussian.loglik.g(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 10),params = list(r =  1/dim(df.train)[1],g = min(dim(df.train)[1],(dim(df.train)[2])^2)))
#perfect agreement 

library(tictoc)
tic()
mean(sapply(1:100000,function(i)lm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 1),params = list(r =  1/dim(df.train)[1], betaprior = "g-prior",alpha = min(dim(df.train)[1],(dim(df.train)[2])^2)))$crit))
toc()

tic()
mean(sapply(1:100000,function(i)gaussian.loglik.g(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 1),params = list(r =  1/dim(df.train)[1],g = min(dim(df.train)[1],(dim(df.train)[2])^2)))$crit))
toc()

tic()
mean(sapply(1:100000,function(i)gaussian.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 1),params = list(r =  1/dim(df.train)[1], var = 1))$crit))
toc()

#BAS version is in fact quicker even than Jeffreys prior based implementation! 


#testing a bit BAS based stuff vs our implementation, Jeffreys prior aka BIC() in BAS
lm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], betaprior = "BIC",alpha = (dim(df.train)[1])))$crit -
  lm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], betaprior = "BIC",alpha = (dim(df.train)[1])))$crit

#var of 1
gaussian.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1],var = 1))$crit - 
  gaussian.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1],var = 1))$crit


#var unknown
gaussian.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1],var = "unknown"))$crit - 
  gaussian.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1],var = "unknown"))$crit

-BIC(lm(semimajoraxis~.,df.train))/2 + BIC(lm(semimajoraxis~.,df.train[,-c(2,3,4,5)]))/2



for(prior in c("CH", "Hyper-g", "Uniform", "Jeffreys", "Beta-prime", "Benchmark", "TruncGamma", "ZS adapted", "Robust", "Hyper-g/n", "Intrinsic"))
{
  print(prior)
  print(gaussian_tcch_log_likelihood(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], prior_beta = prior))$crit-
          gaussian_tcch_log_likelihood(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], prior_beta = prior))$crit)
  
}

lm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T,T),complex = list(oc = 10),params = list(r =  1/dim(df.train)[1], betaprior = "g-prior",alpha = min(dim(df.train)[1],(dim(df.train)[2])^2)))





#let us quickly test the Beroulli responses

df.train$semimajoraxis = as.numeric(df.train$semimajoraxis>mean(df.train$semimajoraxis))


glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = Jeffreys(),laplace = 1))$crit -
  glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = Jeffreys(),laplace = 1))$crit

# laplace or not does not matter, does not fully correspond to ours
glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = bic.prior((dim(df.train)[1])),laplace = 1))$crit -
  glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = bic.prior((dim(df.train)[1])),laplace = 1))$crit


logistic.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1]))$crit - 
  logistic.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,F,F,F,F,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1]))$crit


-BIC(glm(semimajoraxis~.,df.train,family = "binomial"))/2 + BIC(glm(semimajoraxis~.,df.train[,-c(2,3,4,5)],family = "binomial"))/2


tic()
mean(sapply(1:10000,function(i)glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = Jeffreys(),laplace = 1))$crit))
toc()

tic()
mean(sapply(1:10000,function(i)glm.logpost.bas(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1], family = "binomial", betaprior = bic.prior((dim(df.train)[1])),laplace = 1))$crit))
toc()

tic()
mean(sapply(1:10000,function(i)logistic.loglik(y = df.train$semimajoraxis,x = cbind(1,df.train[,-1]),model = c(T,T,T,T,T,T,T,T,T),complex = list(oc = 0),params = list(r =  1/dim(df.train)[1]))$crit))
toc()


