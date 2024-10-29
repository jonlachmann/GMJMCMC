#######################################################
#
# Example 8 (Section 5.2): Sanger data again
#
# High dimensional analysis without nonlinearities
#
# Now using g prior for coefficients
#
# This is the valid version for the JSS Paper
#
#######################################################

library(FBMS)
use.fbms = FALSE  

setwd("/home/florian/FBMS/")
load("Sangerdata.Rdata")

df = as.data.frame(cbind(as.numeric(data[24266,-1]),
                         t(as.matrix(data[-24266,-1]))
))

names(df) = c("y",paste0("x",1:47292))

# Candidates for the first MJMCMC round based on marginal p values
p.vec = unlist(mclapply(2:47293, function(x)cor.test(df[,1],df[,x])$p.value))
ids = sort(order(p.vec)[1:50])          

transforms = c("")
params = gen.params.gmjmcmc(df[,ids])
params$feat$check.col <- F
params$feat$pop.max = 60
params$feat$prel.filter <- ids
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,0,1)


####################################################
#
# Here begin the changes to use Zellers g-prior 
#
####################################################


params$loglik$g <- dim(df)[1]   # Using sample size for g in g-prior

#this will be added to the package
log.prior <- function(params,complex){
  
  pl <-  log(params$r) * (sum(complex$oc))
  return(pl)
}

gaussian.loglik.g <- function (y, x, model, complex, params) 
{
  
  suppressWarnings({
     mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
   # Calculate R-squared
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  Rsquare <- 1 - (RSS / TSS)
  
  # logarithm of marginal likelihood
  mloglik <- 0.5*(log(1.0 + params$g) * (dim(x)[1] - mod$rank)  - log(1.0 + params$g * (1.0 - Rsquare)) * (dim(x)[1]  - 1))*(mod$rank!=1)
  
  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  
  return(list(crit = mloglik + lp, coefs = mod$coefficients))
}


#result <- mjmcmc(params = params, data = df, gaussian.loglik.g, N = 5000)

set.seed(66)


if (use.fbms) {
  result1 <- fbms(data = df, family = "custom", loglik.pi = gaussian.loglik.g, method = "gmjmcmc", 
                  transforms = transforms, probs = probs, params = params, P=25)
} else {
  result1 <-  gmjmcmc(data = df, loglik.pi = gaussian.loglik.g, transforms = transforms, 
                     probs = probs, params = params, P=25)
}
summary(result1)


plot(result1,17)

#Correlation analysis

S = summary(result1)
names.S = S$feats.strings

X.best = df[,names.S]

cor(X.best)
corrplot::corrplot(cor(X.best))
hist(cor(X.best))


# Correlation with results from Example 3




