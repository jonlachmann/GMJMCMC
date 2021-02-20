# Title     : A file to keep a running example of the algorithm, not really anything to keep around when it is done
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-12

library(roxygen2)
roxygenize()
library(devtools)
build()

install.packages("../GMJMCMC_1.0.tar.gz")
library(GMJMCMC)

# install.packages("RSQLite")
library(RSQLite)

con <- dbConnect(drv=RSQLite::SQLite(), dbname= "../scraper/dalen.db")
tables <- dbListTables(con)
sales <- dbGetQuery(conn=con, statement=paste0("SELECT * FROM sales"))

sales2 <- sales[,3:8]
sales2$x4 <- sales[,4]
sales2$x5 <- sales[,5]
sales2$x6 <- sales[,6]
sales2$x7 <- sales[,7]
sales2$x8 <- sales[,8]

covmat <- cor(sales2)

sigmoid<-function(x)exp(-x)
logg <- function (x) log(abs(x)+.Machine$double.eps)

transforms <- c("logg", "logg")

# Generate lists of parameters and probabilities
probs <- gen.probs.list(transforms)
params <- gen.params.list()

loglik.test <- function (data, model, formula, complex) {
  linmod <- lm(formula = formula, data=data)
  ret <- -AIC(linmod)
  if (ret > 0) {
    print(formula)
    print(ret)
  }
  return(ret)
}

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
