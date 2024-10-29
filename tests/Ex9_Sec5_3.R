#######################################################
#
# Example 9 (Section 5.3): Zambia data set from the cAIC4 package
#
# Linear Mixed Model with Fractional Polynomials (one might consider other non-linear features)
#
# Added dummy variables for Second and Third visit
#
# This is the valid version for the JSS Paper
#
#######################################################


library(FBMS)
#this will be added to the package

log.prior <- function(params,complex){
  
  pl =  -log(params$r) * (sum(complex$oc))
  return(pl)
}

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

#options(repos=c( inlabruorg = "https://inlabru-org.r-universe.dev", INLA = "https://inla.r-inla-download.org/R/testing", CRAN = "https://cran.rstudio.com") )
#install.packages("fmesher") 

library(INLA)

#install.packages("cAIC4") 
#library(cAIC4)

data(Zambia, package = "cAIC4")

df = Zambia[1:5]


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) # Modifications and interactions!

params <- gen.params.gmjmcmc(df)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$loglik$r = 1/dim(df)[1]

#specify indices for a random effect
params$loglik$dr = Zambia$dr # district ids for repeated measurements 
params$loglik$INLA.num.threads = 10 # Number of threads used by INLA

params$greedy$steps = 2
params$greedy$tries = 1
params$sa$t.min = 0.1
params$sa$dt = 10



#estimator function
mixed.model.loglik.inla <- function (y, x, model, complex, params) 
{
  
  if(sum(model)>1)
  {
    data1 = data.frame(y, as.matrix(x[,model]), params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~",paste0(names(data1)[3:(dim(data1)[2]-1)],collapse = "+"),"+ f(params.dr,model = \"iid\")"))
  } else
  {
    data1 = data.frame(y, params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~","1 + f(params.dr,model = \"iid\")"))
  }
 
  #to make sure inla is not stuck
  inla.setOption(inla.timeout=30)
  inla.setOption(num.threads=params$INLA.num.threads) 
  
  
  mod<-NULL
  #importance with error handling for unstable libraries that one does not trust 100%
  tryCatch({
    mod <- inla(family = "gaussian",silent = 1L,safe = F, data = data1,formula = formula1)
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  

  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  
  if(length(mod)<3||length(mod$mlik[1])==0) {
     return(list(crit = -10000 + lp, coefs = rep(0,dim(data1)[2]-1)))
  } else {
     mloglik <- mod$mlik[1]
     return(list(crit = mloglik + lp, coefs = mod$summary.fixed$mode))
    
  }
}


set.seed(03052024)

result <- gmjmcmc(data = df, loglik.pi = mixed.model.loglik.inla, transforms = transforms, probs = probs, params = params, P = 25)

plot(result)
summary(result)



set.seed(22052024)


params$loglik$INLA.num.threads = 1 # Number of threads used by INLA set to 1
result2 <- gmjmcmc.parallel(runs = 20, cores = 20, data = df, loglik.pi = mixed.model.loglik.inla, transforms = transforms, probs = probs, params = params, P = 25)

plot(result2)
summary(result2)


