library(RCurl)

# Download data
crime_x <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-x1.txt"),sep = ",")
crime_y <- read.table(text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-y1.txt"))
# Add intercept
crime_x <- cbind(1, crime_x)

full_model_count <- 32768

# Calculate the full model set
crime_full <- vector("list", full_model_count)
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:15]))
  loglik <- linear.g.prior.loglik(as.matrix(crime_y), as.matrix(crime_x), modelvector, NULL, list(g=47))
  crime_full[[i]] <- list(prob=NA, model=modelvector[-1], crit=loglik, alpha=NA)
}

# Calculate the renormalized values for the full model space (ground truth)
full_renorm <- marginal.probs.renorm(crime_full)
names(full_renorm) <- paste0("y",1:15)

# Set up data and parameters for MJMCMC
colnames(crime_y) <- "Y"
crime_data <- cbind(crime_y, crime_x)
crime_probs <- gen.probs.list(NA)
crime_params <- gen.params.list(crime_data)
crime_params$loglik$g <- 47

# Generate models using MJMCMC
crime_result <- gmjmcmc(crime_data, linear.g.prior.loglik, NA, transforms, 1, 60000, 60000, crime_probs, crime_params)

# Calculate the renormalized values for the MJMCMC result
mjmcmc_renorm <- marginal.probs.renorm(crime_result$models[[1]])

rmse <- function (full, sim, iters) {
  sim_renorm <- marginal.probs.renorm(sim$models[[1]][1:iters])
  names(sim_renorm) <- paste0("y",1:length(sim_renorm))
  sim_renorm <- sim_renorm[order(full)]
  full <- sort(full)
  rmse <- sqrt((sim_renorm - full)^2)
  return(rmse*100)
}

# Calculate the rmse for 1000-60000 iterations as the MJMCMC baseline
rmse_converge <- matrix(NA, 60, 15)
for(i in 1:60) {
  rmse_converge[i,] <- rmse(full_renorm, crime_result, i*1000)
}

# Plot the convergence
plot(rowMeans(rmse_converge), type="l")

### Download simulated logistic data as per example 2
logistic_x <- read.csv(header=F, text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20(Example%203)/sim3-X.txt"))
logistic_y <- read.csv(header=F, text=getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20(Example%203)/sim3-Y.txt"))
# Modify data
logistic_x$V2<-(logistic_x$V10+logistic_x$V14)*logistic_x$V9
logistic_x$V5<-(logistic_x$V11+logistic_x$V15)*logistic_x$V12

# Add intercept
logistic_x <- cbind(1, logistic_x)

logistic.loglik.aic <- function (y, x, model, complex, params) {
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial())})
  ret <- -(mod$deviance/2) - mod$rank
  return(ret)
}

full_model_count <- 2^20

# Calculate the full model set
logistic_full2 <- vector("list", full_model_count)
progress <- 0
for (i in 1:full_model_count) {
  modelvector <- as.logical(c(T,intToBits(i)[1:20]))
  loglik <- logistic.loglik.aic(as.matrix(logistic_y), as.matrix(logistic_x), modelvector, NULL, list(g=47))
  logistic_full2[[i]] <- list(prob=NA, model=modelvector[-1], crit=loglik, alpha=NA)
  if (i %% floor(full_model_count/40) == 0) progress <- print.progressbar(progress, 40)
}

save(logistic_full2, file="logistic_full2")

# Calculate the renormalized values for the full model space (ground truth)
logistic_full.renorm <- marginal.probs.renorm(logistic_full2)
names(logistic_full.renorm) <- paste0("y",1:20)

# Set up data and parameters for MJMCMC
colnames(logistic_y) <- "Y"
logistic_data <- cbind(logistic_y, logistic_x)
logistic_probs <- gen.probs.list(NA)
logistic_params <- gen.params.list(logistic_data)

# Generate models using MJMCMC
logistic_result <- gmjmcmc(logistic_data, logistic.loglik.aic, NA, transforms, 1, 60000, 60000, logistic_probs, logistic_params)

logistic_rmse_converge <- matrix(NA, 60, 20)
for(i in 1:60) {
  logistic_rmse_converge[i,] <- rmse(logistic_full.renorm, logistic_result, i*1000)
}

# Plot the convergence
plot(rowMeans(logistic_rmse_converge), type="l")

# NAIVE SUBSAMPLING EXPERIMENTS BELOW
naive.linear.g.prior.loglik <- function (y, x, model, complex, params) {
  subsample <- sample.int(nrow(x), 20)
  out <- lm.fit(as.matrix(x[subsample,model]),y[subsample])
  rsquared <- 1-sum(var(out$residuals))/sum(var(y[subsample]))
  p <- out$rank
  n <- 20
  logmarglik <- 0.5*(log(1+params$g)*(n-p) - log(1+params$g*(1-rsquared))*(n-1))*(p!=1)
  return(logmarglik)
}

crime_result_naive <- gmjmcmc(crime_data, naive.linear.g.prior.loglik, NA, transforms, 1, 60000, 60000, crime_probs, crime_params)

rmse_converge_naive <- matrix(NA, 60, 15)
for(i in 1:60) {
  rmse_converge_naive[i,] <- rmse(full_renorm, crime_result_naive, i*1000)
}

plot(rowMeans(rmse_converge_naive), type="l")

test_renorm <- marginal.probs.renorm(crime_result_naive$models[[1]])
names(test_renorm) <- paste0("y",1:15)
barplot(rbind(test_renorm, full_renorm), beside=T)