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
largejump <- c(0, 0, 0, 1) # probability for type of large jump, only allow type 1-4
localopt <- c(0.5, 0.5) # probability for each localopt algorithm
random <- c(0.3, 0.3, 0.2, 0.2) # probability for random jump kernels
mh <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1) # probability for regular mh kernels

filter <- 0.3 # filtration threshold
gen <- c(1/3, 1/3, 1/3) # probability for different feature generation methods
trans <- c(0.5, 0.5) # probability for each different nonlinear transformation

probs <- list(large=large, largejump=largejump, localopt=localopt,
              random=random, filter=filter, gen=gen,
              trans=trans)

# Create the list of parameters
sa_kern <- list(probs=c(0.1, 0.05, 0.2, 0.3, 0.2, 0.15), neigh.size=1, neigh.min=1, neigh.max=2) # Simulated annealing proposal kernel parameters
sa_params <- list(t.init=10, t.min=0.0001, dt=3, M=12, kern=sa_kern) # Simulated annealing parameters
greedy_params <- list() # Greedy algorithm parameters
large_params <- list(neigh.size=2, neigh.min=1, neigh.max=2) # Large jump parameters
random_params <- list(neigh.size=1, neigh.min=1, neigh.max=2) # Small random jump parameters
mh_params <- list(neigh.size=1, neigh.min=1, neigh.max=2) # Regular MH parameters
jump_params <- list(large=0.5, large.max=0.6, large.min=0.4, small=0.2, small.max=0.3, small.min=0.1) # Neighborhood size parameters

params <- list(mh=mh_params, sa=sa_params, greedy=greedy_params, large=large_params, random=random_params)

# install.packages("RSQLite")
library(RSQLite)

con <- dbConnect(drv=RSQLite::SQLite(), dbname= "../scraper/dalen.db")
tables <- dbListTables(con)
sales <- dbGetQuery(conn=con, statement=paste0("SELECT * FROM sales"))

sales2 <- sales[,3:8]
trans <- c("sin", "cos")

loglik.test <- function (data) {
  linmod <- lm(data)
  AIC(linmod)
}

options(warn=2)
gmjmcmc(sales2, loglik.test, trans, 1, 100, probs, params)





