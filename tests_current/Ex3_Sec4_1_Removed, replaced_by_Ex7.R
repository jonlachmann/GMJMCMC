#######################################################
#
# Example 3: Sanger data (Section 4.1)
#
# High dimensional analysis without nonlinearities
#
# This is the valid version for the JSS Paper
#
#######################################################

# Logical to decide whether to perform analysis with fbms function
# If FALSE then gmjmcmc or gmjmcmc.parallel function is used
use.fbms = FALSE

library(FBMS)

data(SangerData2)
df = SangerData2
#Rename columns
colnames(df) = c("y",paste0("x",1:(ncol(df)-1)))

# Candidates for the first MJMCMC round based on correlation with response
c.vec = unlist(mclapply(2:ncol(df), function(x)abs(cor(df[,1],df[,x]))))
ids = sort(order(c.vec,decreasing=TRUE)[1:50])


####################################################
#
# single thread analysis (four different runs)
#
# Comparison of gmjmcmc.parallel with one thread and gmjmcmc
#
####################################################

params = gen.params.gmjmcmc(df)
params$feat$check.col <- F
params$feat$pop.max = 60
params$prel.select <- ids

transforms = c("")
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,0,1)
probs$filter=0.8
params$loglik$var = "unknown"
#params$loglik$r = 1
set.seed(123)

start.time=Sys.time()
if (use.fbms) {
  result1 <- fbms(data = df, method = "gmjmcmc", transforms = transforms,
                  probs = probs, params = params, P=25)
} else {
  result1 =  gmjmcmc(data = df, transforms = transforms,
                     probs = probs, params = params, P=25)
}
show(Sys.time()-start.time)
summary(result1)


################################
set.seed(124)   #Same analysis using a different seed

if (use.fbms) {
  result2 <- fbms(data = df, method = "gmjmcmc", transforms = transforms,
                  probs = probs, params = params, P=25)
} else {
  result2 =  gmjmcmc(data = df, transforms = transforms,
                     probs = probs, params = params, P=25)
}

## Combine results from two runs
res1 = summary(result1,tol=0.05)
res1$marg.probs = round(res1$marg.probs,3)
res2 = summary(result2,tol=0.05)
res2$marg.probs = round(res2$marg.probs,3)
names.best = unique(c(res1$feats.strings,res2$feats.strings))
m = max(nrow(res1),nrow(res2))
while(nrow(res1)<m)
  res1 = rbind(res1,c("",""))
while(nrow(res2)<m)
  res2 = rbind(res2,c("",""))
cbind(res1[1:m,],res2[1:m,])

##
ind = order(as.numeric(substring(names.best,first=2)))
names.best = names.best[ind]
X.best = df[,names.best]
pdf("crossplot_Sanger_best.pdf")
corrplot::corrplot(cor(X.best))
dev.off()

cor(X.best)
min(cor(X.best))


################################
#
# Comparing results
#




################################

#Same analysis but using slightly different initial population
# Candidates for the first MJMCMC round based on marginal p values
ids3 = ids

transforms = c("")
params = gen.params.gmjmcmc(df[,ids3])
params$feat$check.col <- F
params$feat$pop.max = 60
params$prel.select <- ids3
params$loglik$var <- "unknown"
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,0,1)


set.seed(123)

if (use.fbms) {
  result3 <- fbms(data = df, method = "gmjmcmc", transforms = transforms,
                  probs = probs, params = params, P=25)
} else {
  result3 =  gmjmcmc(data = df, transforms = transforms,
                     probs = probs, params = params, P=25)
}


# And again for the sake of comparison
summary(result3,tol = 0.01)
summary(result1,tol = 0.01)
summary(result2,tol = 0.01)









####################################################
#
# multiple thread analysis
#
####################################################

set.seed(123)

if (use.fbms) {
  result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", runs = 10, cores = 10,
                                      transforms = transforms, probs = probs, params = params,
                                      P=25, N.init=500, N.final=500)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 40, cores = 40,data = df,
                                      transforms = transforms, probs = probs, params = params,
                                      P=25, N.init=500, N.final=500)
}
save(result_parallel,file="Ex3_parallel.RData")
load("Ex3_parallel.RData")

plot(result_parallel)
summary(result_parallel,tol = 0.01)

S = summary(result_parallel,tol=0.05)
names.best = S$feats.strings

X.best = df[,names.best]

cor(X.best)
min(cor(X.best))
corrplot::corrplot(cor(X.best))
hist(cor(X.best))


######################################



# repeat same analysis with different seed
set.seed(1234)

if (use.fbms) {
  result_parallel2 <- fbms(data = df, method = "gmjmcmc.parallel", runs = 40, cores = 10,
                          transforms = transforms, probs = probs, params = params,
                          P=25, N.init=500, N.final=500)
} else {
  result_parallel2 =  gmjmcmc.parallel(runs = 40, cores = 10,data = df,
                         transforms = transforms, probs = probs, params = params,
                         P=25, N.init=500, N.final=500)
}
save(result_parallel2,file="Ex3_parallel2.RData")
load("Ex3_parallel2.RData")
plot(result_parallel2)
summary(result_parallel2)

## Combine results from two runs
res1 = summary(result_parallel,tol=0.05)
res1$marg.probs = round(res1$marg.probs,3)
res2 = summary(result_parallel2,tol=0.05)
res2$marg.probs = round(res2$marg.probs,3)
names.best = unique(c(res1$feats.strings,res2$feats.strings))
m = max(nrow(res1),nrow(res2))
while(nrow(res1)<m)
  res1 = rbind(res1,c("",""))
while(nrow(res2)<m)
  res2 = rbind(res2,c("",""))
cbind(res1[1:m,],res2[1:m,])

##
ind = order(as.numeric(substring(names.best,first=2)))
names.best = names.best[ind]
X.best = df[,names.best]
pdf("crossplot_Sanger_best_par.pdf")
corrplot::corrplot(cor(X.best))
dev.off()


S2 = summary(result_parallel2,tol=0.05)
names.best2 = S2$feats.strings

X.best = df[,names.best2]

cor(X.best)
min(cor(X.best))
corrplot::corrplot(cor(X.best))
hist(cor(X.best))

Cor.check = cor(X.best)
diag(Cor.check) = 0
max(Cor.check)
which.max(Cor.check)

# Comparison of the two parallel runs
sum(is.element(names.best,names.best2))
sum(is.element(names.best2,names.best))

length(intersect(names.best,names.best2))

cbind(sort(names.best),sort(names.best2))

##Combine results
feats = unique(c(S[,1],S2[,1]))
n = length(feats)
Scomb = data.frame(feats.strings=feats,marg.probs1=rep(0,n),marg.probs2=rep(0,n))
Scomb$marg.probs1[which(S$feats.strings %in%Scomb$feats.strings)] = S$marg.probs
Scomb$marg.probs2[which(S2$feats.strings %in% Scomb$feats.strings)] = S2$marg.probs
