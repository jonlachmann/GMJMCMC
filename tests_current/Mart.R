
#Loadrequiredlibraries
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(lavaan)
#library(FBMS)
library(psych)
#-------------------------------
#Step1:LoadandRenameData
#-------------------------------
#Setfilepath
file_path<-"DataAnalysis202501.xlsx"

#Loaddataset
data<-read_excel("/Users/aliaksandrhome/Rprojects/Mary SEM/DataAnalysis 202501.xlsx",sheet="Data")

data$Ps1<-8-data$Ps1
data$Ps3<-8-data$Ps3
data$Ps5<-8-data$Ps5
#Loadandpreprocessdata
cleaned_data<-data%>%select(Pi1,Pi2,Pi3,Ps1,PS2,Ps3,Ps4,Ps5,Ps6,Ps7,
                            Bo1,Bo2,Bo3,Bo4,Bo5,Bo6,Bo7,Bo8,
                            Vo1,Vo2,Vo3,Vo4,Vo5,Vo6,
                            Vo7,Vo8,Vo9,
                            Wp1,Wp2,Wp3,Wp4,Wp5)

#Computecompositescores
data_scores<-cleaned_data%>%mutate(
  PerceivedImpact=rowMeans(select(.,Pi1,Pi2,Pi3),na.rm=TRUE),
  PsychologicalSafety=rowMeans(select(.,Ps1,PS2,Ps3,Ps4,Ps5,Ps6,Ps7),na.rm=TRUE),
  Burnout=rowMeans(select(.,Bo1,Bo2,Bo3,Bo4,Bo5,Bo6,Bo7,Bo8),na.rm=TRUE),
  Voice=rowMeans(select(.,Vo1,Vo2,Vo3,Vo4,Vo5,Vo6),na.rm=TRUE),
  Silence=rowMeans(select(.,Vo7,Vo8,Vo9),na.rm=TRUE),
  WorkPerformance=rowMeans(select(.,Wp1,Wp2,Wp3,Wp4,Wp5),na.rm=TRUE)
)

#Bartlett’sTestandKMOMeasuretoassesssamplingadequacy
bt<-bartlett.test(cleaned_data)
print(bt)
kmo_result<-KMO(cleaned_data)
print(kmo_result)#Shouldbe>0.6forgoodfactorability

#ExploratoryFactorAnalysis(EFA)withPrincipalAxisFactoringandVarimaxRotation
efa_result<-fa(cleaned_data,nfactors=6,rotate="varimax",fm="pa")
print(efa_result$loadings,cutoff=0.4)#Displayloadingsabove0.4
print(efa_result$Vaccounted)

print(sum(efa_result$Vaccounted["Proportion Var",]))

barplot(efa_result$Vaccounted["Proportion Var",],
        main="Variance Explainedby Each Factor",
        ylab="Proportionof Variance",
        names.arg=paste("F",1:6,sep=""),
        col="skyblue")

#Interpretation
#-Highfactorloadings(>0.4)indicatestrongassociationwiththefactor.
#-Bartlett’stestshouldbesignificant(p<0.05)indicatingthatthecorrelationmatrixisnotanidentitymatrix.
#-KMOvalues>0.7indicatesamplingadequacyforfactoranalysis.
#-Cronbach’salpha>0.7indicatesgoodinternalconsistency.
#-PCAshouldexplainagoodportionofthevariance(>60%).

#Definesubsetsforeachfactor
factor_groups<-list(
  PerceivedImpact=cleaned_data%>%select(Pi1,Pi2,Pi3),
  PsychologicalSafety=cleaned_data%>%select(Ps1,PS2,Ps3,Ps4,Ps5,Ps6,Ps7),
  Burnout=cleaned_data%>%select(Bo1,Bo2,Bo3,Bo4,Bo5,Bo6,Bo7,Bo8),
  Voice=cleaned_data%>%select(Vo1,Vo2,Vo3,Vo4,Vo5,Vo6),
  Silence=cleaned_data%>%select(Vo7,Vo8,Vo9),
  WorkPerformance=cleaned_data%>%select(Wp1,Wp2,Wp3,Wp4,Wp5)
)

#FunctiontocomputeKMOandBartlett'stest
compute_tests<-function(data_subset,factor_name){
  cat("\n----------------------\n")
  cat("Factor:",factor_name,"\n")
  
  #Bartlett'sTest
  bt<-bartlett.test(data_subset)
  #KMOMeasure
  kmo_result<-KMO(data_subset)
  
  al<-alpha(data_subset)
  
  return(c(bt$statistic,bt$parameter,bt$p.value,kmo_result$MSA,al$total$raw_alpha))
}

#Applyfunctiontoeachfactor
stats.factors<-t(sapply(names(factor_groups),function(f)compute_tests(factor_groups[[f]],f)))
colnames(stats.factors)<-c("b-stat","b-df","b-pval","KMO","alpha")
print(stats.factors)


model<-'
#MeasurementModel:Definelatentvariables
PerceivedImpact=~Pi1+Pi2+Pi3
PsychologicalSafety=~Ps1+PS2+Ps3+Ps4+Ps5+Ps6+Ps7
Burnout=~Bo1+Bo2+Bo3+Bo4+Bo5+Bo6+Bo7+Bo8
Voice=~Vo1+Vo2+Vo3+Vo4+Vo5+Vo6
Silence=~Vo7+Vo8+Vo9
WorkPerformance=~Wp1+Wp2+Wp3+Wp4+Wp5


#StructuralModel:Relationshipsamongfactors
WorkPerformance~Voice+Silence
Burnout~Voice+Silence
Voice~PsychologicalSafety+PerceivedImpact+PsychologicalSafety*PerceivedImpact
Silence~PsychologicalSafety+PerceivedImpact+PsychologicalSafety*PerceivedImpact
'
fit<-sem(model,data=cleaned_data,estimator="MLR")

sf=summary(fit,standardized=TRUE,fit.measures=TRUE)

BIC(fit)

names(cleaned_data)<-c("Pi1","Pi2","Pi3","Ps1","Ps2","Ps3","Ps4","Ps5","Ps6","Ps7","Bo1","Bo2","Bo3","Bo4","Bo5","Bo6","Bo7","Bo8","Vo1","Vo2","Vo3","Vo4","Vo5","Vo6","Si1","Si2","Si3","Wp1","Wp2","Wp3","Wp4","Wp5")

cleaned_data <- data.frame(cleaned_data)


estimator.sem.fbms.p<- function (y, x, model, complex, params) 
{
  
  fparam<-params$strings[which(model==1)]
  
  if(length(fparam) == 1)
  {
    return(list(crit = -10000, coefs = c(0,rep(1,length(fparam)))))
  }
  
  
  data <- data.frame(x[,which(model==1)[-1]])
  names(data) <- fparam[-1]
  
  
  PI <- fparam[which(stringi::stri_startswith_fixed(fparam,"Pi"))]
  PS <- fparam[which(stringi::stri_startswith_fixed(fparam,"Ps"))]
  BU <- fparam[which(stringi::stri_startswith_fixed(fparam,"Bo"))]
  VO <- fparam[which(stringi::stri_startswith_fixed(fparam,"Vo"))]
  SI <- fparam[which(stringi::stri_startswith_fixed(fparam,"Si"))]
  WP <- fparam[which(stringi::stri_startswith_fixed(fparam,"Wp"))]
  
  if(length(PI)==0|length(BU)==0|length(PS)==0|length(VO)==0|length(SI)==0|length(WP)==0)
  {
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  }
  
  PI <- paste0("PerceivedImpact=~", paste0(PI,collapse = "+"))
  PS <- paste0("PsychologicalSafety=~", paste0(PS,collapse = "+"))
  BU <- paste0("Burnout=~", paste0(BU,collapse = "+"))
  VO <- paste0("Voice=~", paste0(VO,collapse = "+"))
  SI <- paste0("Silence=~", paste0(SI,collapse = "+"))
  WP <- paste0("WorkPerformance=~", paste0(WP,collapse = "+"))
  
  
  SEM <- '
  WorkPerformance~Voice+Silence 
  Burnout~Voice+Silence 
  Voice~PsychologicalSafety+PerceivedImpact
  Silence~PsychologicalSafety+PerceivedImpact 
  '
  
  model<- paste0(c(PI,PS,BU,VO,SI,WP,SEM,collapse = "\n"))
  
  
  out<-sem(model,data=data,estimator="WLSMV")
  
  if(!inspect(out,"converged"))
  {
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  }
  meas <- fitMeasures(out, c("cfi", "tli", "rmsea", "srmr","pvalue"))
  
  
  if(meas[5] > 0.05)
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  
  logmarglik <- sum(meas[1:2]) - sum(meas[3:4]) 
  
  
  return(list(crit = logmarglik, coefs = c(0,rep(1,length(fparam))),model = out))
  
  
}


estimator.sem.fbms<- function (y, x, model, complex, params) 
{
  
  fparam<-params$strings[which(model==1)]
  
  if(length(fparam) == 1)
  {
    return(list(crit = -10000, coefs = c(0,rep(1,length(fparam)))))
  }
  
  
  data <- data.frame(x[,which(model==1)[-1]])
  names(data) <- fparam[-1]
  
  
  PI <- fparam[which(stringi::stri_startswith_fixed(fparam,"Pi"))]
  PS <- fparam[which(stringi::stri_startswith_fixed(fparam,"Ps"))]
  BU <- fparam[which(stringi::stri_startswith_fixed(fparam,"Bo"))]
  VO <- fparam[which(stringi::stri_startswith_fixed(fparam,"Vo"))]
  SI <- fparam[which(stringi::stri_startswith_fixed(fparam,"Si"))]
  WP <- fparam[which(stringi::stri_startswith_fixed(fparam,"Wp"))]
  
  if(length(PI)==0|length(BU)==0|length(PS)==0|length(VO)==0|length(SI)==0|length(WP)==0)
  {
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  }
  
  PI <- paste0("PerceivedImpact=~", paste0(PI,collapse = "+"))
  PS <- paste0("PsychologicalSafety=~", paste0(PS,collapse = "+"))
  BU <- paste0("Burnout=~", paste0(BU,collapse = "+"))
  VO <- paste0("Voice=~", paste0(VO,collapse = "+"))
  SI <- paste0("Silence=~", paste0(SI,collapse = "+"))
  WP <- paste0("WorkPerformance=~", paste0(WP,collapse = "+"))
  
  
  SEM <- '
  WorkPerformance~Voice+Silence 
  Burnout~Voice+Silence 
  Voice~PsychologicalSafety+PerceivedImpact
  Silence~PsychologicalSafety+PerceivedImpact 
  '
  
  model<- paste0(c(PI,PS,BU,VO,SI,WP,SEM,collapse = "\n"))
  
  
  out<-sem(model,data=data,estimator="WLSMV")
  
  if(!inspect(out,"converged"))
  {
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  }
  meas <- fitMeasures(out, c("cfi", "tli", "rmsea", "srmr","pvalue"))
  
  
  if(meas[5] > 0.05)
    return(list(crit = -10000 + rnorm(1), coefs = c(0,rep(1,length(fparam)))))
  
  logmarglik <- sum(meas[1:2]) - sum(meas[3:4]) 
  
  
  return(list(crit = logmarglik, coefs = c(0,rep(1,length(fparam)))))
  
  
}

cleaned_data$Y = 1
params <- gen.params.mjmcmc(data = cleaned_data)
params$loglik$strings = c("1","Pi1","Pi2","Pi3","Ps1","Ps2","Ps3","Ps4","Ps5","Ps6","Ps7","Bo1","Bo2","Bo3","Bo4","Bo5","Bo6","Bo7","Bo8","Vo1","Vo2","Vo3","Vo4","Vo5","Vo6","Si1","Si2","Si3","Wp1","Wp2","Wp3","Wp4","Wp5")

probs <- gen.probs.mjmcmc()
probs$large <- 0.005

params$mh$neigh.size <- 3
params$mh$neigh.min <- 1
params$mh$neigh.max <- 5

params$sa$kern$neigh.size <- 2
params$sa$kern$neigh.max <- 3

params$greedy$neigh.size <- 2
params$greedy$neigh.max <- 3
res <- FBMS::fbms(Y~1+.,data = cleaned_data,family = "custom",loglik.pi = estimator.sem.fbms,method = "mjmcmc",params = params,probs = probs, verbose = T, N = 5000)


set.seed(2)
res3 <- FBMS::fbms(Y~1+.,data = cleaned_data,family = "custom",loglik.pi = estimator.sem.fbms,method = "mjmcmc.parallel",params = params,probs = probs,runs = 64, cores = 8, verbose = F, N = 200)

save(res3,file = "res3")

sr3 <- summary(res3,labels=names(cleaned_data)[-33])
sr2 <- summary(res2,labels=names(cleaned_data)[-33])

sr3$feats.strings
sr2$feats.strings
plot(res,count = 20)

sum(sr3$feats.strings%in%sr2$feats.strings)/length(sr3$feats.strings)

sum(sr2$feats.strings%in%sr3$feats.strings)/length(sr2$feats.strings)

res4<-list(res)
bestmod2 = list()
bestcrit = rep(0,1)
res4$best.crit
i = 1
for(res in res4){
  bestmod2[[i]] <- which.max(sapply(res$models,function(x)x$crit))
  bestcrit[i] = res$best.crit
  i = i + 1
}

bestcrit

bc <- which.max(bestcrit)
bestmod = bestmod2[[bc]]

out<-estimator.sem.fbms.p(y = cleaned_data$Y,x = as.matrix(cleaned_data[,c(33,1:32)]),model = c(T,res4[[bc]]$models[[bestmod]]$model),params = params$loglik)

out$model

out$crit


params$loglik$strings[c(T,res4[[bc]]$models[[bestmod]]$model)]

summary(out$model, standardized = TRUE, fit.measures = TRUE)


fitMeasures(out$model, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

BIC(out)

fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

BIC(fit)







estimator.sem<-function (formula, data) 
{
  fmla.proc<-as.character(formula)[2:3]
  fobserved<-fmla.proc[1]
  fmla.proc[2]<-stringi::stri_replace_all(str=fmla.proc[2],
                                          fixed=" ",replacement="")
  fmla.proc[2]<-stringi::stri_replace_all(str=fmla.proc[2],
                                          fixed="\n",replacement="")
  fparam<-stringi::stri_split_fixed(str=fmla.proc[2],pattern="+",
                                    omit_empty=FALSE)[[1]]
  
  PI <- paste0("PerceivedImpact=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Pi"))],collapse = "+"))
  PS <- paste0("PsychologicalSafety=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Ps"))],collapse = "+"))
  BU <- paste0("Burnout=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Bo"))],collapse = "+"))
  VO <- paste0("Voice=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Vo"))],collapse = "+"))
  SI <- paste0("Silence=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Si"))],collapse = "+"))
  WP <- paste0("WorkPerformance=~", paste0(fparam[which(stringi::stri_startswith_fixed(fparam,"Wp"))],collapse = "+"))
  
  SEM <- '
  WorkPerformance~Voice+Silence 
  Burnout~Voice+Silence 
  Voice~PsychologicalSafety+PerceivedImpact+PsychologicalSafety*PerceivedImpact  
  Silence~PsychologicalSafety+PerceivedImpact+PsychologicalSafety*PerceivedImpact  
  '
  
  model<- paste0(c(PI,PS,BU,VO,SI,WP,SEM,collapse = "\n"))
  
  out<-sem(model,data=data,estimator="MLR")
  
  logmarglik<--stats::BIC(out)
  
  return(list(mlik=logmarglik,waic=stats::AIC(out),dic=stats::BIC(out),
              summary.fixed=list(mean=c(0,rep(1,length(fparam))))))
}


estimator.sem(formula = formula,data = cleaned_data)

EMJMCMC::estimate.logic.lm

EMJMCMC::pinferunemjmcmc()

fitMeasures(fit,c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
modificationindices(fit,sort=TRUE,minimum.value=10)