#################################################
#
# Example 1:
#
# Kepler Example from the JAIR paper
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
                         runs = 40, cores = 10, P=30,params = params)
} else {
 result_parallel <- gmjmcmc.parallel(runs = 40, cores = 10, data = df.train, loglik.pi = gaussian.loglik, 
                                     transforms = transforms, P=30,params = params)
}
summary(result_parallel, tol = 0.01)
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

c(rmse.default, rmse.P50, rmse.parallel)


