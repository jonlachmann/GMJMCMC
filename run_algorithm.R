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

  transforms <- c("sini", "to25","expi","logi","to35","troot","sigmoid")
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

  testdata <- matrix(cbind(y,1,x1,x2,x3,x4,x5,x6,x7,x8), nrow=nobs)
}


xx1 <- 5*sini(x1)
xx2 <- 2*logi(x2+x8)
xx3 <- 3*troot(x3)
xx4 <- 3*x5*x6

testdata2 <- cbind(xx1, xx2, xx3, xx4)
glmmodel <- glm.fit(testdata2[1:100,], y[1:100])


mattt <- matrix(c(2,2,1,1,1,1), 2, 3)
fett <- list(3, mattt)
class(fett) <- "feature"
print(fett, transforms, dataset=T, alphas=T)

probs$filter <- 0.8
params$feat$pop.max <- 15

probs <- gen.probs.list(transforms)
params <- gen.params.list(testdata)

system.time(result3 <- gmjmcmc(testdata[1:100,], gaussian.loglik, gaussian.loglik.alpha, transforms, 50, 500, 1000, probs, params))

par(mar=c(3,3,3,3))
gmjmcmc.totdens.plot(result3)

pop <- unlist(print.model(list(model=rep(T,15)), result3$populations[[48]], transforms))
pop[matt[23837]]
marginal.probs.renorm(result3$models[[45]])

matt <- matrix(unlist(result3$models[2:50]), ncol=18, byrow=T)

mat <- matrix(unlist(result3$models[[45]]), ncol=18, byrow=T)

(exp(-10) + exp(-8)) / (exp(-9) + exp(-8))
(exp(-20) + exp(-16)) / (exp(-18) + exp(-16))


log(exp(-10) + exp(-8))

logspace_add <- function(logx,logy) {
    pmax(logx,logy) + log1p(exp(-abs(logx - logy)))
}

install.packages("Rmpfr")
library(Rmpfr)

m5000 <- mpfr(-5001, 128)

exp(m5000)

testx <- rep(c(5,4),500)
testy <- rep(1,1000)

glmod1 <- glm.fit(testx, testy)
glmod <- glm.fit(testx[1:500], testy[1:500])





logspace_add(-4300,-4000)




exp(-10+(-10/-9))



summ <- summary.gmjresult(result3,50)

importance <- data.frame(c(summ$importance))
names(summ$importance) <- summ$features
par(mar=c(14,3,3,3))
barplot(summ$importance, las=2)

totdens <- gmjmcmc.totdens.plot(result2)

result3$accept

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

fcnnn <- function (p)
{
    sum((y - p[1]+p[2]*x1+p[3]*x2+p[4]*x3)^2)
}

library(GenSA)
gensa <- GenSA(rep(0,5), fcnnn, rep(-15,5), rep(15,5))

lm(y ~ x1+x2+x3-1)

gnlr(y=y, mu=formula(~x1+x2+x3), distribution = "normal", pmu=rep(0,4), pshape = rep(0,1))


planets <- read.csv("../DBRM/exa1.csv")
planets <- planets[,-1]

result_planets2 <- gmjmcmc(planets, gaussian.loglik, gaussian.loglik.alpha, transforms, 20, 5000, 10000, probs, params)

gmjmcmc.totdens.plot(result_planets2)

unlist(print.model(list(model=rep(T,15)), result_planets2$populations[[14]], transforms))
unlist(print.model(list(model=rep(T,15)), result_planets2$populations[[15]], transforms))
marginal.probs.renorm(result_planets2$models[[14]])



pop <- unlist(print.model(list(model=rep(T,15)), result_planets$populations[[25]], transforms))
pop[matt[23837]]


matt <- matrix(unlist(result3$models[2:50]), ncol=18, byrow=T)

library(RCurl)

#prepare data
simx <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-x1.txt"),sep = ",")
simy <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-y1.txt"))
data.example <- cbind(simy,simx)
names(data.example)[1]="Y"

params$loglik$g <- 47
params$random$neigh.size <- 2
params$random$neigh.min <- 2
params$mh$neigh.size <- 2
params$mh$neigh.min <- 2
params$large$neigh.size <- 4
params$large$neigh.min <- 4
params$large$neigh.max <- 4


result_crime <- gmjmcmc(data.example, linear.g.prior.loglik, gaussian.loglik.alpha, transforms, 1, 200, 3276, probs, params)
renorm <- marginal.probs.renorm(result_crime$models[[1]])
names(renorm) <- paste0("y",1:15)

modds <- matrix(unlist(result_crime$models), ncol=18, byrow=T)

mcmc <- marginal.probs(result_crime$models[[1]])
names(mcmc) <- paste0("y",1:15)
sort(mcmc)
sort(renorm)

simxx <- cbind(1,simx)

logliks <- matrix(NA,32768, 17)
for (i in 1:32768) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  loglik <- linear.g.prior.loglik(as.matrix(simy), as.matrix(simxx), modelvector, NULL, list(g=47))
  logliks[i,] <- c(loglik, modelvector)
}
unimodds <- modds[(!duplicated(modds[,2:16], dim=1)),]
renormtruth <- matrix(NA,1,15)
for (i in 3:17) renormtruth[i] <- sum(exp(logliks[as.logical(logliks[,i]),1]))/sum(exp(logliks[,1]))
renormmj <- matrix(NA,1,15)
for (i in 2:16) renormmj[i] <- sum(exp(unimodds[as.logical(unimodds[,i]),17]))/sum(exp(unimodds[,17]))

sum(exp(unimodds[,17]))
sum(exp(logliks[,1]))

sort(renormmj)
sort(renormtruth)



# Hassan
x1<-rnorm(20,0,1)
x2<-abs(rnorm(20,2,4))
e<-rnorm(20,0,1)
y<-20+3*x1+4*sqrt(x2)+10*x1^3+e

result <- gmjmcmc(cbind(y,x1,x2), loglik.test,NULL, transforms, 30, 20, 50, probs, params)