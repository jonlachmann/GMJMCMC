# Title     : A file to keep a running example of the algorithm, not really anything to keep around when it is done
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-12

library(roxygen2)
roxygenize()
library(devtools)
build()

test()

install.packages("../GMJMCMC_1.0.tar.gz")
library(GMJMCMC)

sigmoid <- function(x)exp(-x)
sini <- function(x)sin(x/180*pi)
expi <- function(x)exp(-abs(x))
logi <- function(x)log(abs(x)+.Machine$double.eps)
troot <- function(x)abs(x)^(1/3)
to25 <- function(x)abs(x)^(2.5)
to35 <- function(x)abs(x)^(3.5)

transforms <- c("to25","expi","logi","to35","troot","sigmoid")

probs <- gen.probs.list(transforms)
params <- gen.params.list()

data("breastcancer")

bc <- breastcancer[,c(ncol(breastcancer),2:(ncol(breastcancer)-1))]

result <- gmjmcmc(bc, logistic.loglik, transforms, 30, 1000, 5000, probs, params)

result$accept

summ <- summary.gmjresult(result)

importance <- data.frame(c(summ$importance))
names(summ$importance) <- summ$features
barplot(summ$importance, las=2)

summ$features

library(profvis)
profvis({result <- gmjmcmc(bc, loglik.test, transforms, 30, 20, 50, probs, params)})

x1 <- rnorm(100, 1, 2)
x2 <- rnorm(100, 3, 2)
x3 <- rnorm(100, 5,2)
x4 <- rnorm(100, 1, 200)
x5 <- rnorm(100)
x6 <- rnorm(100)
x7 <- rnorm(100)
x8 <- rnorm(100)
y <- 0.5*sini(x1)+2*logi(x2+x8)+7*troot(x3)-5*x5*x6 -15
y2 <- rbinom(100, 1, 1/(1+exp(-y)))

testdata <- data.frame(cbind(y2,x1,x2,x3,x4,x5,x6,x7,x8))


loglik.tester <- function (data, model, formula, complex) {
  summod <- sum(which(model)) - sum(complex$width)*5
  return(summod)
}

result <- gmjmcmc(testdata, logistic.loglik, transforms, 50, 200, 500, probs, params)

summ <- summary.gmjresult(result)

importance <- data.frame(c(summ$importance))
names(summ$importance) <- summ$features
barplot(summ$importance, las=2)

greedy_kern_test <- list(probs=c(1,0,0,0,0,0),
                      neigh.size=1, neigh.min=1, neigh.max=2)           # Greedy algorithm proposal kernel parameters
greedy_params_test <- list(steps=30, kern=greedy_kern_test)

greedy.optim(c(F,F,F,F,F,F,F,F,F,F), NULL, loglik.tester, c(T,T,T,T,T,T,T,T,T,T), NULL, greedy_params_test)

sa_kern_test <- list(probs=c(1,0,0,0,0,0), neigh.size=1, neigh.min=1, neigh.max=2)
sa_params_test <- list(t.init=100, t.min=0.00001, dt=3, M=12, kern=sa_kern_test)

simulated.annealing(c(F,F,F,F,F,F,F,F,F,F), NULL, loglik.tester, c(T,T,T,T,T,T,T,T,T,T), NULL, sa_params_test)


for(i in 1:10) sales2[,i] <- as.numeric(sales2[,i])

result <- gmjmcmc(sales2, loglik.test, transforms, 100, 200, 500, probs, params)

result$accept

print.model(result$models[[30]][[1]], result$populations[[30]], transforms)

summ <- summary.gmjresult(result)
barplot(summ$importance)

importance <- as.data.frame(summ$importance)
rownames(importance) <- summ$features

importance


### Visualization experiments below

resmat <- matrix(unlist(result$models), ncol=12, byrow=T)

hist(resmat[,12], breaks=100)

# Give each variable its own angle
var.angles <- seq(0, 1*pi, length.out=10+2)[2:11]
x.positions <- cos(var.angles)
y.positions <- sin(var.angles)

# Calculate aggregated angle based measures for each model
mods <- matrix(as.numeric(unlist(resmat[2031:2530,1:10])), ncol=10)
x.pos <- rowSums(mods*x.positions)
y.pos <- rowSums(mods*y.positions)

dff <- as.data.frame(cbind(x.pos, y.pos))

plot(dff)
library(ggplot2)

# Create a plot
ggplot(dff, aes(x=x.pos, y=y.pos)) +
  stat_density2d(aes(fill=..level..), geom = "polygon", colour="white", show.legend=F)

library(gnlm)

attach(sales2)

gnlr(y=price, mu=~cos(b0+b1*size+b2*room_count), pmu=rnorm(3, 0, 0.0001), pshape=1)
