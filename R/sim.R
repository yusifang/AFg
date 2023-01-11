
### generating reference dist for simulation ###

rm(list=ls())
setwd("/home/yuf31/sparse_signals/results")
n=10^3
B=10^5
core_num=60
library(parallel)
library(Rcpp)
library(TFisher)

#-----HC-------------#
cppFunction("NumericVector HCtmpC(NumericVector y,int d0, int d){
            NumericVector out(d0);
            NumericVector out2(d0);
            
            for (int i = 0; i < d0; i++) {
            double dev_temp=i+1;
            double dev= dev_temp/d-y[i];
            double up=sqrt(d)*dev;
            double down=sqrt(y[i]*(1-y[i]));
            out[i]= up/down;
            out2[i]=dev_temp/d;
            }
            return out;
            }")


HC=function(x){
  d=length(x)
  y=sort(x)
  d0=ceiling(d/2)
  res=max(HCtmpC(y,d0,d))
  return(res)
}


#---------AFz------------#
AFz.stat = function(p.values) { 
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  v.stat.standard = rep(NA, num.study)
  for (i in 1:num.study) {
    v.stat = -sum(sort.p.log[1:i])
    w = c(rep(1, i), i/((i+1):(num.study)))
    v.stat.standard[i] = (v.stat - sum(w))/sqrt(sum(w^2))
  }
  v.AFz = max(abs(v.stat.standard))
  return(list(AFz.stat = v.AFz))
}


#----------TFishersoftomni-------#
TFsoft.omni=function(p.values){
  TAU1 = c(0.01, 0.05, 0.5, 1)
  return(stat.soft.omni(p.values, TAU1=TAU1)$omni)
}
#---------TF omni---------------#
TF.omni=function(p.values){
  TAU1 = c(0.01, 0.05, 0.5, 1)
  TAU2 = c(0.1, 0.2, 0.5, 1)
  return(stat.tfisher.omni(p.values, TAU1=TAU1, TAU2=TAU2)$omni)
}
#---------AFp------------#
AFp.stat = function(p.values) {
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  v.stat.pval = rep(NA, num.study)
  for (i in 1:num.study) {
    v.stat = -2 * sum(sort.p.log[1:i])
    v.stat.pval[i] = pchisq(v.stat, 2*i, lower.tail = FALSE)
  }
  v.AFp = min(v.stat.pval)
  return(list(AFp.stat = v.AFp))
}

#--------minP------------#

minP.stat = function(p.values) {
  v.stat = min(p.values)
  return(list(minP.stat = v.stat))
}

#-----------------------AFg----------------------------------#

geomSeries <- function(n,q) {
  max=n^q
  base=(1+1/log(n))
  
  res=ceiling(c(log(n)*base^(0:log(max/(log(n))^2, base)),n))
  res=unique(c(1,res,n))
  return(res)
}

## Primary Reference


### N times (col for study, row for times)


null <-  do.call(rbind, mclapply(1:B, function(x) runif(n), mc.cores = core_num))
sort.null <-  t(log(apply(null, 1, sort)))

AFg.primary_refdistN <- mclapply(geomSeries(n,1), function(x)  apply(as.matrix(sort.null[,1:x]), 1, function(t) -2 * sum(t)), mc.cores = core_num)
names(AFg.primary_refdistN)=geomSeries(n,1)
save(AFg.primary_refdistN,  file = paste0("AFg_primaryrefdistN",n,"final.RData"))



## Secondary Reference


### N times


AFg.stat = function(p.values,n) {
  vec=geomSeries(n,1)
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  position.p = rank(p.values, ties.method = "first")
  v.stat.pval = rep(NA, length(vec))
  for (i in 1:length(vec)) {
    
    v.stat = -2 * sum(sort.p.log[1:vec[i]])
    v.stat.pval[i] = mean(AFg.primary_refdistN[[i]] > v.stat)
  }
  v.AFp = min(v.stat.pval)
  
  return(v.AFp)
}


Fisher=function(p.values){
  stat=sum(-2*log(p.values))
  return(stat)
}

null = do.call(cbind, mclapply(1:B,function(x) runif(n), mc.cores = core_num))
lib_ref=list()
lib_ref[[1]]= unlist(mclapply(1:B, function(x) AFg.stat(null[,x],n),mc.cores=core_num))

lib_ref[[2]]=unlist(mclapply(1:B, function(x) minP.stat(null[,x]),mc.cores=core_num))

lib_ref[[3]]=unlist(mclapply(1:B, function(x) AFp.stat(null[,x]),mc.cores=core_num))
lib_ref[[4]]=unlist(mclapply(1:B, function(x) TF.omni(null[,x]),mc.cores=core_num))
lib_ref[[5]]=unlist(mclapply(1:B, function(x) TFsoft.omni(null[,x]),mc.cores=core_num))
lib_ref[[6]]=unlist(mclapply(1:B, function(x) AFz.stat(null[,x]),mc.cores=core_num))

lib_ref[[7]]=unlist(mclapply(1:B, function(x) HC(null[,x]),mc.cores=core_num))



rm(null)
rm(AFg.primary_refdistN)

names(lib_ref)=c("AFg","minP","AFp","TF.omni","TFsoft.omni","AFz","HC")
save(lib_ref,file="lib_reffinal.Rdata")
cut=0.05
cutoff=list()
cutoff=lapply(1:5, function(x) quantile(lib_ref[[x]],probs=cut,na.rm=T))
cutoff[[6]]=quantile(lib_ref[[6]],probs = (1-cut),na.rm = T)

cutoff[[7]]=quantile(lib_ref[[7]],probs = (1-cut),na.rm = T)


save(cutoff,file=paste0("cutoff_cut=",cut,"final.Rdata"))




rm(list=ls())
setwd("/home/yuf31/sparse_signals/results")
n=10^3
B=10^5
core_num=60
cut=0.05
library(parallel)
library(Rcpp)
library(TFisher)
library(mvtnorm)
load(paste0("AFg_primaryrefdistN",n,"final.RData"))

load(paste0("cutoff_cut=",cut,"final.Rdata"))
names(cutoff)=c("AFg","minP","AFp","TF.omni","TFsoft.omni","AFz","HC")
#-----HC-------------#
cppFunction("NumericVector HCtmpC(NumericVector y,int d0, int d){
            NumericVector out(d0);
            NumericVector out2(d0);
            
            for (int i = 0; i < d0; i++) {
            double dev_temp=i+1;
            double dev= dev_temp/d-y[i];
            double up=sqrt(d)*dev;
            double down=sqrt(y[i]*(1-y[i]));
            out[i]= up/down;
            out2[i]=dev_temp/d;
            }
            return out;
            }")


HC=function(x){
  d=length(x)
  y=sort(x)
  d0=ceiling(d/2)
  res=max(HCtmpC(y,d0,d))
  return(res)
}


#---------AFz------------#
AFz.stat = function(p.values) { 
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  v.stat.standard = rep(NA, num.study)
  for (i in 1:num.study) {
    v.stat = -sum(sort.p.log[1:i])
    w = c(rep(1, i), i/((i+1):(num.study)))
    v.stat.standard[i] = (v.stat - sum(w))/sqrt(sum(w^2))
  }
  v.AFz = max(abs(v.stat.standard))
  return(list(AFz.stat = v.AFz))
}


#----------TFishersoftomni-------#
TFsoft.omni=function(p.values){
  TAU1 = c(0.01, 0.05, 0.5, 1)
  return(stat.soft.omni(p.values, TAU1=TAU1)$omni)
}
#---------TF omni---------------#
TF.omni=function(p.values){
  TAU1 = c(0.01, 0.05, 0.5, 1)
  TAU2 = c(0.1, 0.2, 0.5, 1)
  return(stat.tfisher.omni(p.values, TAU1=TAU1, TAU2=TAU2)$omni)
}
#---------AFp------------#
AFp.stat = function(p.values) {
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  v.stat.pval = rep(NA, num.study)
  for (i in 1:num.study) {
    v.stat = -2 * sum(sort.p.log[1:i])
    v.stat.pval[i] = pchisq(v.stat, 2*i, lower.tail = FALSE)
  }
  v.AFp = min(v.stat.pval)
  return(list(AFp.stat = v.AFp))
}

#--------minP------------#

minP.stat = function(p.values) {
  v.stat = min(p.values)
  return(list(minP.stat = v.stat))
}

#-----------------------AFg----------------------------------#

geomSeries <- function(n,q) {
  max=n^q
  base=(1+1/log(n))
  
  res=ceiling(c(log(n)*base^(0:log(max/(log(n))^2, base)),n))
  res=unique(c(1,res,n))
  return(res)
}

AFg.stat = function(p.values,n) {
  vec=geomSeries(n,1)
  num.study = length(p.values)
  sort.p = sort(p.values)
  sort.p.log = log(sort.p)
  position.p = rank(p.values, ties.method = "first")
  v.stat.pval = rep(NA, length(vec))
  for (i in 1:length(vec)) {
    
    v.stat = -2 * sum(sort.p.log[1:vec[i]])
    v.stat.pval[i] = mean(AFg.primary_refdistN[[i]] > v.stat)
  }
  v.AFp = min(v.stat.pval)
  
  return(v.AFp)
}



Fisher=function(p.values){
  stat=sum(-2*log(p.values))
  return(stat)
}

#-----alt---gen-----------------------------------------------#
detection_boundary2=function(beta,n,shift){
  S=round(n^(1-beta))
  r=0
  if(beta>0.75){r=(1-sqrt(1-beta))^2}
  if(beta<=0.75){r=beta-0.5}
  return(sqrt(2*(r+shift)*log(n)))
}


alt_gen_sub3=function(beta,n,shift,B,sigma){
  S=round(n^(1-beta))
  mu0=detection_boundary2(beta,n,shift)
  mu=rep(mu0,S)+rnorm(S,0,sigma)
 
  y=do.call(rbind,lapply(1:B, function(x) return(rmvnorm(1,mean=mu,sigma=diag(S)))))
  
  return(2*(1-pnorm(abs(y))))
}
alt_gen3=function(beta,n,shift,B,sigma){
  p_alt=alt_gen_sub3(beta,n,shift,B,sigma)
  S=round(n^(1-beta))
  n_temp=n-S
  p_ref=matrix(runif(B*n_temp),ncol = n_temp)
  
  return(cbind(p_alt,p_ref))
}




#---------simulation2-----------------------#
shift_vec=c(0.2,0.05,0.1)
sigma_vec=c(0.05,0.1,0.2)
beta_vec=seq(0.55,0.95,0.05)
B_alt=10^4
R=30
methods=c("AFg","minP","AFp","TF.omni","TFsoft.omni","AFz","HC")
cut=0.05
for(i1 in 1:length(shift_vec)){
  shift=shift_vec[[i1]]
  for(i2 in 1:length(sigma_vec)){
    
    
    sigma=sigma_vec[[i2]]
    power=list()
    
    for(s in 1:R){
      print(s)
      set.seed(s)
      null= mclapply(beta_vec,function(x) alt_gen3(x,n,shift =shift,B_alt,sigma=sigma), mc.cores = core_num)
      
      res_temp=matrix(,nrow = length(beta_vec),ncol = length(methods))
      colnames(res_temp)=methods
      for(j in 1:length(beta_vec)){
        
        res=list()
        res[[1]]= mean(unlist(mclapply(1:B_alt,function(id) (AFg.stat(null[[j]][id,],n)),mc.cores = core_num))<=cutoff[[1]])
        
        
        
        res[[2]]= mean(unlist(mclapply(1:B_alt,function(id) (minP.stat(null[[j]][id,])),mc.cores = core_num))<=cutoff[[2]])
        
        res[[3]]= mean(unlist(mclapply(1:B_alt,function(id) (AFp.stat(null[[j]][id,])),mc.cores = core_num))<=cutoff[[3]])
        
        res[[4]]= mean(unlist(mclapply(1:B_alt,function(id) (TF.omni(null[[j]][id,])),mc.cores = core_num))<=cutoff[[4]])
        
        res[[5]]= mean(unlist(mclapply(1:B_alt,function(id) (TFsoft.omni(null[[j]][id,])),mc.cores = core_num))<=cutoff[[5]])
        
        res[[6]]= mean(unlist(mclapply(1:B_alt,function(id) (AFz.stat(null[[j]][id,])),mc.cores = core_num))>=cutoff[[6]])
        
        res[[7]]= mean(unlist(mclapply(1:B_alt,function(id) (HC(null[[j]][id,])),mc.cores = core_num))>=cutoff[[7]])
        
        res=unlist(res)
        
        names(res)=c("AFg","minP","AFp","TF.omni","TFsoft.omni","AFz","HC")
        res_temp[j,]=res
      }
      power[[s]]=res_temp
    }
    
    save(power,file = paste0("random_HCsetup_power_shift=",shift,"sigma=",sigma,"_n=",n,"_finalv2.Rdata"))
  }
  
}
