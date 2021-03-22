# Title     : A file to keep a running example of the algorithm, not really anything to keep around when it is done
# Objective : A scratchpad for me to try out different things while developing
# Created by: jonlachmann
# Created on: 2021-02-12

library(roxygen2)
roxygenize()
library(devtools)
build(vignettes=FALSE)
build()
test()

install.packages("../GMJMCMC_1.0.tar.gz")
detach("package:GMJMCMC", unload=TRUE)
library(GMJMCMC)

{
  files.sources <- list.files("R", full.names=T)
  sapply(files.sources, source)

  sigmoid <- function(x)exp(-x)
  sini <- function(x)sin(x/180*pi)
  expi <- function(x)exp(-abs(x))
  logi <- function(x)log(abs(x)+.Machine$double.eps)
  troot <- function(x)abs(x)^(1/3)
  to25 <- function(x)abs(x)^(2.5)
  to35 <- function(x)abs(x)^(3.5)

  transforms <- c("sini", "to25","expi","logi","to35","troot","sigmoid")

  probs <- gen.probs.list(transforms)
  params <- gen.params.list()
}

library(GenSA)
sourceCpp("src/set_alphas.cpp")
params$feat$alpha <- 3

glm.logistic.loglik <- function (y, x, model, complex) {
  r <- 20/223
  suppressWarnings({mod <- glm.fit(as.matrix(x[,model]), y, family=binomial())})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$width)))/2
  return(ret)
}

data("breastcancer")

bc <- breastcancer[,c(ncol(breastcancer),2:(ncol(breastcancer)-1))]

system.time(result2 <- gmjmcmc(testdata, logistic.loglik, logistic.loglik.alpha, transforms, 10, 100, 100, probs, params))
system.time(result2 <- gmjmcmc(testdata, glm.logistic.loglik, transforms, 2, 100, 100, probs, params))

result2$accept

summ <- summary.gmjresult(result)

importance <- data.frame(c(summ$importance))
names(summ$importance) <- summ$features
barplot(summ$importance, las=2)

summ$features

library(profvis)
profvis({result <- gmjmcmc(bc, loglik.test, transforms, 30, 20, 50, probs, params)})

{
  nobs <- 1000
  x1 <- rnorm(nobs, 1, 2)
  x2 <- rnorm(nobs, 3, 2)
  x3 <- rnorm(nobs, 5,2)
  x4 <- rnorm(nobs, 1, 4)
  x5 <- rnorm(nobs)
  x6 <- rnorm(nobs)
  x7 <- rnorm(nobs)
  x8 <- rnorm(nobs)
  y <- 5*sini(x1)+2*logi(x2+x8)+3*troot(x3)-3*x5*x6 -15
  y2 <- rbinom(nobs, 1, 1/(1+exp(-y)))

  testdata <- matrix(cbind(y2,1,x1,x2,x3,x4,x5,x6,x7,x8), nrow=nobs)
}

mattt <- matrix(c(2,2,1,1,1,1), 2, 3)
fett <- list(3, mattt)
class(fett) <- "feature"
print(fett, transforms, dataset=T, alphas=T)


system.time(result3 <- gmjmcmc(testdata, logistic.loglik, logistic.loglik.alpha, transforms, 50, 500, 1000, probs, params))

summ <- summary.gmjresult(result3)

importance <- data.frame(c(summ$importance))
names(summ$importance) <- summ$features
par(mar=c(14,3,3,3))
barplot(summ$importance, las=2)

totdens <- gmjmcmc.totdens.plot(result2)

result$accept

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

# Alpha generation stuff
library(gnlm)
library(Rcpp)


feat <- alpha_3(featt, transforms, "g", loglik)
print(feat, transforms)


featt <- create.feature(1, list(feature2, feature1),c(1,2,3))

feature2 <- result$populations[[9]][[1]]

print(featt, transforms)

feat2 <- update.alphas(featt, c(10.1,9.1,8.1,7.1,6.1,5.1,4.1,3.1,2.1,1.1))
print(feat2, transforms)

# Create an environment that gnlr can work with
testenv <- attach(testdata)

# Create the formula for gnlr
mufcn <- model.function(c(F,F,F,F,T,F,F,F), result$populations[[8]], transforms, "g")

mufcn2 <- mufcn


for (i in 1:20) {
  sarestmp <- GenSA(rep(0,mufcn$count), loglik, rep(-i*10,mufcn$count), rep(i*10,mufcn$count), control=list(max.call=1e4), mufcn$formula)
  sares2[i,] <- c(sarestmp$par, sarestmp$value)
}


loglik(rep(2,5),mufcn)

# Run gnlr to estimate alpha parameters
nlm <- gnlr(y=y2, mu=formul, distribution = "binomial", envir=testenv, pmu=rep(0.3,5))

install.packages("nloptr")
library(nloptr)

gnlr(y=y2, mu=formul, distribution = "binomial", envir=testenv, pmu=rep(0.3,5))

matt <- matrix(NA,200,5)
liks <- matrix(NA, 200, 1)
for(i in 1:200) {
  matt[i,] <- gnlr(y=y2, mu=formul, distribution = "binomial", envir=testenv, pmu=rep(i/50,5), steptol=1e-10)$coefficients
  liks[i] <- gnlr(y=y2, mu=formul, distribution = "binomial", envir=testenv, pmu=rep(i/50,5), steptol=1e-10)$maxlik
}

plot(sort(round(liks, digits=5)), type="l")

plot(sort(liks), type="l", ylim=c(55,70))
lines(sort(samat[,6]), col="red")

min(samat[,6])
min(liks)

cor(matt)

for (i in 2:20) {
  matt[i,] <- matt[i,] / matt[1,]
}

for (i in 2:5) matt[,i] <- matt[,i] / matt[,1]

plot(matt[,1], type="l", col="red")
lines(matt[,2], col="blue")
lines(matt[,3], col="yellow")
lines(matt[,4], col="green")
lines(matt[,5])

samat <- matrix(NA, 20, 6)

for (i in 1:20) {
  satmp <- GenSA(rep(as.numeric(i),5), fcnn, rep(-100,5), rep(100,5))
  samat[i,] <- c(satmp$par, satmp$value)
}


