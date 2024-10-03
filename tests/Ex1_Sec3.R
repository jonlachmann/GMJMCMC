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
library(FBMS)

#setwd("/home/florian/FBMS/")

#X <- read.csv("exa1.csv")
#X = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/kepler%20and%20mass/exa1.csv")
#df <- as.data.frame(cbind(MajorAxis = X[,5], X[,-5]))
data(exoplanet)
df <- as.data.frame(cbind(MajorAxis = exoplanet[,5], exoplanet[,-5]))


to3 <- function(x) x^3

transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")


####################################################
#
# single thread analysis (default values, Section 3.1)
#
####################################################


set.seed(123)

result.default <- gmjmcmc(df, transforms=transforms)

params = gen.params.gmjmcmc(df)
params$feat$pop.max = 50
result <- gmjmcmc(df, gaussian.loglik, gaussian.loglik.alpha, transforms,
                  P=200,N.init=1000,N.final=1000)



####################################################
#
# single thread analysis (more iterations, Section 3.2)
#
####################################################


set.seed(123)

result.P50 <- gmjmcmc(df, gaussian.loglik, transforms  = transforms,
                  P=50, N.init=1000, N.final=5000)


####################################################
#
# multiple thread analysis (Section 3.3)
#
####################################################

set.seed(124)

# Actual parallel analysis works currently only under Linux or Mac
# result_mm =  gmjmcmc.parallel(runs = 4, cores = 4,df, gaussian.loglik, gaussian.loglik.alpha, transforms)


result.parallel <- gmjmcmc.parallel(runs = 40, cores = 10,
                                    data = df, loglik.pi = gaussian.loglik, 
                                    transforms = transforms, P=50)


####################################################
#
# Inspection of Results (Section 3.4)
#
####################################################

######################
# summary

summary(result)

summary(result.default, labels = names(df)[-1])

summary(result.P50)

summary(result.parallel)

summary(result.parallel, tol = 0.01)


######################
# plot

pdf("result.pdf") 
plot(result.default)
dev.off()

plot(result)



pdf("result.P50.pdf") 
plot(result.P50)
dev.off()

plot(result.P50)



pdf("result_parallel.pdf") 
plot(result.parallel)
dev.off()

plot(result.parallel)
plot(result.parallel, 12)


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
plot(preds.P50$aggr$mean, df$MajorAxis)


pdf("prediction.P50.pdf") 
plot(preds.P50$aggr$mean, df$MajorAxis)
dev.off()


rmse.P50 <-  sqrt(mean((preds.P50$aggr$mean - df$MajorAxis)^2))


###############################


preds.multi <- predict(result.parallel , df[,-1], link = function(x) x)  

pdf("pred_parallel.pdf") 
plot(preds.multi$aggr$mean, df$MajorAxis)
dev.off()

rmse.parallel <- sqrt(mean((preds.multi$aggr$mean - df$MajorAxis)^2))



c(rmse.default, rmse.P50, rmse.parallel)
