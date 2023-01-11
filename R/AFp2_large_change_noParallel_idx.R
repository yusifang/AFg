library(Rfast)
library(parallel)
library(pracma)



mtop=function(sort.p,m=10){
  sort.p.log=-2*log(sort.p)
  return(Rfast::rowsums(sort.p.log[,1:m]))
}
geomSeries <- function(n,q) {
  max=n^q
  base=(1+1/log(n))
  
  res=ceiling(c(log(n)*base^(0:log(max/(log(n))^2, base)),n))
  res=unique(c(1,res,n))
  return(res)
}

AFg=function(sort.p,f_list){
  m_vec=geomSeries(dim(sort.p)[2],1)
  sort.p.log=-2*log(sort.p)
  
  mtx=do.call(cbind,lapply(1:length(m_vec),function(x) f_list[[x]](rowsums(sort.p.log[,1:m_vec[x]]))))
  mtx=cbind(mtx,-log(pbeta(sort.p.log[,1],shape1 = 1,shape2 = dim(sort.p)[2],lower.tail = T))
            ,-log(pchisq(rowsums(sort.p.log),df=2*dim(sort.p)[2],lower.tail = F)))
  return(rowMaxs(mtx,value = T))
  
  
}


#################CE mixed gradient#######################3
CE.AFg.mixed.grad=function(K=50,theta,N1=1000,f_list,idx) {
  grad=1/((1:K)^idx+1)
  
  
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  #apply(mix.prop,1,mean)
  mix.prop=transpose(mix.prop)
  #apply(mix.prop,2,mean)
  
  # set.seed(1)
  x=matrix(Rnorm(K*N1,0,theta),nrow = N1)^{(mix.prop)}+matrix(Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  p=2*pnorm(abs(x),mean=0,sd=1,lower.tail = FALSE)
  AFg_input=Rfast::rowSort(p)
  
  stat.val=AFg(sort.p = AFg_input,f_list)
  
  
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  
  log_g_mix <- rowsums(log(transpose(transpose(dnorm(x,0,theta))*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  
  
  w=exp(log_g_norm-log_g_mix)
  return(list(stat.val=stat.val,weight=w,x=x))
}
mixed.par.update.grad=function(object,ro,K,idx) {
  
  
  
  #stat.val=unlist(lapply(object,"[[",1))
  stat.val=object[[1]]
  gamma=quantile(stat.val,1-ro)
  #weight=unlist(lapply(object,"[[",2))
  weight=object[[2]]
  #x.mat=transpose(do.call(rbind,lapply(object, "[[",3)))
  x.mat=transpose(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)
  
   opt=optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t)+(1-grad)*dnorm(x,0,1))))),
              interval = c(0,5),maximum = TRUE)
   par=opt$maximum
  
  # ft=function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,
  #                                                      function(x) sum(log(1/((1:K)+1)*dnorm(x,0,t)+(1-1/((1:K)+1))*dnorm(x,0,1)))))
  # opt=optim(2,ft,gr = function(t) pracma::grad(ft, t),
  #           control=list(fnscale=-1),method = "L-BFGS-B",lower=1,upper = 6)
  # 
  # #par=opt$maximum
  # par=opt$par
  # 
  
  
  return(list(gamma=gamma,par=par))
}



CE.mixed.grad=function(K,ro,N=10^5,q,f_list,idx) {
  theta=1
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.AFg.mixed.grad(K=K,theta=theta,N1=N,f_list,idx)
    
    par=mixed.par.update.grad(object = obj,ro=ro,K=K,idx)
    gamma=par$gamma
    theta=par$par
    message(paste0(gamma,", ",theta))
  }
  message("stop")
  obj.final=CE.AFg.mixed.grad(K=K,theta=theta,N1=N,f_list,idx)
  
  stat.val.final=obj.final[[1]]
  
  weight.final=obj.final[[2]]
  
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  return(list(p.val=p.val,iter=t,theta=theta))
}

libs_AFg.CE.mixGrad=function(q.val.set=exp(seq(log(3),log(22),length.out = 200)),K=50,ro=0.01,N=10^5,f_list,idx){
  temp.lib=mclapply(q.val.set, function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=x,f_list,idx),mc.cores = 30)
  temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}

AFg.pvals.CE.mixGrad=function(q.val.set=exp(seq(log(100),log(170),length.out = 200)),K=50,ro=0.01,N=10^5,M,f_list) {
  nc=length(q.val.set)/M
  pvals_all=c()
  pi=1
  for (i in 1:M) {
    pvals_i=mclapply(q.val.set[pi:(nc+pi-1)],function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=x,f_list,idx), mc.cores = nc)
    pvals_i=unlist(lapply(pvals_i,"[[",1))
    pvals_all=c(pvals_all,pvals_i)
    pi=nc+pi
  }
  return(pvals_all)
}