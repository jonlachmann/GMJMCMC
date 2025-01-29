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
library(devtools)
devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)
library(FBMS)

setwd("/home/florian/FBMS/")

X <- read.csv("exa1.csv")
df <- as.data.frame(cbind(MajorAxis = X[,5], X[,-5]))


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


set.seed(123)

if (use.fbms) {
 result.default <- fbms(formula = MajorAxis ~ 1 + . , data = df, method = "gmjmcmc", transforms = transforms)
} else {
 result.default <- gmjmcmc(df, transforms = transforms)
}


####################################################
#
# single thread analysis (more iterations, Section 3.2)
#
####################################################


set.seed(123)

if (use.fbms) {
 result.P50 <- fbms(data = df, method = "gmjmcmc", transforms = transforms,
                    P=50, N.init=1000, N.final=5000)
} else {
 result.P50 <- gmjmcmc(df,  transforms = transforms,
                       P=50, N.init=1000, N.final=5000)
}

####################################################
#
# multiple thread analysis (Section 3.3)
#
####################################################

set.seed(124)

# Actual parallel analysis works currently only under Linux or Mac
# result_mm =  gmjmcmc.parallel(runs = 4, cores = 4,df, gaussian.loglik, gaussian.loglik.alpha, transforms)

if (use.fbms) {
 result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms,
                         runs = 40, cores = 10, P=25)
} else {
 result_parallel <- gmjmcmc.parallel(runs = 40, cores = 10,data = df, loglik.pi = gaussian.loglik, 
                                     transforms = transforms, P=25)
}

####################################################
#
# Inspection of Results (Section 3.4)
#
####################################################

######################
# summary

summary(result.default)
summary(result.default, labels = names(df))

summary(result.P50)

summary(result_parallel)
summary(result_parallel, tol = 0.01)


######################
# plot

pdf("result.pdf") 
plot(result.default)
dev.off()

plot(result.default)



pdf("result.P50.pdf") 
plot(result.P50)
dev.off()

plot(result.P50)



pdf("result_parallel.pdf") 
plot(result_parallel)
dev.off()

plot(result_parallel)
plot(result_parallel, 12)


######################
# Prediction

#preds <-  predict(result, df[,-1], link = function(x) x)  
preds <-  predict(result.default, df[,-1])  

pdf("prediction.pdf") 
plot(preds$aggr$mean, df$MajorAxis)
dev.off()

plot(preds$aggr$mean, df$MajorAxis)

rmse.default <- sqrt(mean((preds$aggr$mean - df$MajorAxis)^2))

###############################


#preds.P50 = predict(result.P50, df[,-1], link = function(x) x)  
preds.P50 = predict(result.P50, df[,-1])  

pdf("prediction.P50.pdf") 
plot(preds.P50$aggr$mean, df$MajorAxis)
dev.off()

plot(preds.P50$aggr$mean, df$MajorAxis)

rmse.P50 <-  sqrt(mean((preds.P50$aggr$mean - df$MajorAxis)^2))


###############################


preds.multi <- predict(result_parallel , df[,-1], link = function(x) x)  

pdf("pred_parallel.pdf") 
plot(preds.multi$aggr$mean, df$MajorAxis)
dev.off()

rmse.parallel <- sqrt(mean((preds.multi$aggr$mean - df$MajorAxis)^2))



c(rmse.default, rmse.P50, rmse.parallel)