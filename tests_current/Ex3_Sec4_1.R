#######################################################
#
# Example 7 (Section 5.1): Sanger data again
#
# High dimensional analysis without nonlinearities
#
# Now using g prior for coefficients
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
library(xtable)
use.fbms <- TRUE  
run.parallel <- T

data(SangerData2)
df = SangerData2
#Rename columns
colnames(df) = c("y",paste0("x",1:(ncol(df)-1)))

#Only linear terms/mutations
transforms = c("")
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,0,1)


# Candidates for the first MJMCMC round based on correlation with response
c.vec = unlist(mclapply(2:ncol(df), function(x)abs(cor(df[,1],df[,x]))))
ids = sort(order(c.vec,decreasing=TRUE)[1:50])
params = gen.params.gmjmcmc(ncol(df) - 1)
params$feat$prel.filter <- ids

params$feat$check.col <- T
params$feat$pop.max <- 50

####################################################
#
# Here we shall use Zellners g-prior with g = max(n,p^2)
#
####################################################

##Parallel runs
if(run.parallel)
{
  set.seed(123)
  if (use.fbms) {
    result_parallel1=fbms(data=df,loglik.pi=gaussian.loglik.g,transforms=transforms,
                          probs=probs,params=params,
                          method="gmjmcmc.parallel",
                          P=50,N.init=1000,N.final=1000,runs=10,cores=10)
  }else {
    start = Sys.time()
    result_parallel1=gmjmcmc.parallel(x = df[, -1], y = df[, 1], loglik.pi=gaussian.loglik.g,transforms=transforms,
                                      probs=probs,params=params,
                                      P=50,N.init=1000,N.final=1000,runs=10,cores=10)
    end = Sys.time()
    print(end-start)
  }
  save(result_parallel1,file="Ex3_parallel1.RData")
  #load("Ex3_parallel1.RData")

  set.seed(1234)
  if (use.fbms) {
    result_parallel2=fbms(data=df,loglik.pi=gaussian.loglik.g,transforms=transforms,
                          probs=probs,params=params,
                          method="gmjmcmc.parallel",
                          P=50,N.init=1000,N.final=1000,runs=10,cores=10)
  } else {
    result_parallel2=gmjmcmc.parallel(x = df[, -1], y = df[, 1], loglik.pi=gaussian.loglik.g,transforms=transforms,
                                     probs=probs,params=params,
                                     P=50,N.init=1000,N.final=1000,runs=10,cores=10)
  }
  save(result_parallel2,file="Ex3_parallel2.RData")
  #load("Ex3_parallel2.RData")

  set.seed(123456)
  if (use.fbms) {
    result_parallel3=fbms(data=df,loglik.pi=gaussian.loglik.g,transforms=transforms,
                          probs=probs,params=params,
                          method="gmjmcmc.parallel",
                          P=50,N.init=1000,N.final=1000,runs=10,cores=10)
  } else {
    result_parallel3=gmjmcmc.parallel(x = df[, -1], y = df[, 1], loglik.pi=gaussian.loglik.g,transforms=transforms,
                                     probs=probs,params=params,
                                     P=50,N.init=1000,N.final=1000,runs=10,cores=10)
  }
  save(result_parallel3,file="Ex3_parallel3.RData")
  #load("Ex3_parallel2.RData")


  ## Combine results from three runs
  res1 = summary(result_parallel1,tol=0.01)
  res1$marg.probs = round(res1$marg.probs,3)
  res2 = summary(result_parallel2,tol=0.01)
  res2$marg.probs = round(res2$marg.probs,3)
  res3 = summary(result_parallel3,tol=0.01)
  res3$marg.probs = round(res3$marg.probs,3)
  names.best = unique(c(res1$feats.strings,res2$feats.strings,res3$feats.strings))
  m = max(nrow(res1),nrow(res2),nrow(res3))
  while(nrow(res1)<m)
    res1 = rbind(res1,c("",""))
  while(nrow(res2)<m)
    res2 = rbind(res2,c("",""))
  while(nrow(res3)<m)
    res3 = rbind(res3,c("",""))
  names(res1) = c("feats","prob")
  names(res2) = c("feats","prob")
  names(res3) = c("feats","prob")
  foo = cbind(res1[1:m,],res2[1:m,],res3[1:m,])
  show(print(xtable(foo),include.rownames=FALSE))

  ind = order(as.numeric(substring(names.best,first=2)))
  names.best = names.best[ind]
  X.best = df[,names.best]
  pdf("crossplot_Sanger_best.pdf")
  corrplot::corrplot(cor(X.best))
  dev.off()

}
