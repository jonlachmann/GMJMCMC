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

# Create probability list for algorithm
large <- 0.05 # probability of a large jump
large.kern <- c(0, 0, 0, 1) # probability for type of large jump, only allow type 1-4
localopt.kern <- c(0.5, 0.5) # probability for each localopt algorithm
random.kern <- c(0.3, 0.3, 0.2, 0.2) # probability for random jump kernels
mh <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1) # probability for regular mh kernels

filter <- 0.6 # filtration threshold
gen <- c(1/4, 1/4, 1/4, 1/4) # probability for different feature generation methods
trans <- c(0.5, 0.5) # probability for each different nonlinear transformation

probs <- list(large=large, large.kern=large.kern, localopt.kern=localopt.kern,
              random.kern=random.kern, filter=filter, gen=gen,
              trans=trans)

# Create the list of parameters
sa_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15), neigh.size=1, neigh.min=1, neigh.max=2) # Simulated annealing proposal kernel parameters
sa_params <- list(t.init=10, t.min=0.0001, dt=3, M=12, kern=sa_kern) # Simulated annealing parameters
greedy_params <- list() # Greedy algorithm parameters
large_params <- list(neigh.size=2, neigh.min=1, neigh.max=2) # Large jump parameters
random_params <- list(neigh.size=1, neigh.min=1, neigh.max=2) # Small random jump parameters
mh_params <- list(neigh.size=1, neigh.min=1, neigh.max=2) # Regular MH parameters
jump_params <- list(large=0.5, large.max=0.6, large.min=0.4, small=0.2, small.max=0.3, small.min=0.1) # Neighborhood size parameters
feat_params <- list(D=5, L=15)

params <- list(mh=mh_params, large=large_params, random=random_params,
               sa=sa_params, greedy=greedy_params, feat=feat_params)

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

logg <- function (x) log(abs(x)+.Machine$double.eps)

transforms <- c("logg", "logg")

loglik.test <- function (data, model, formula) {
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
