library(Rfast)
library(parallel)
library(pracma)



mtop=function(sort.p,m=10){
  sort.p.log=-2*log(sort.p)
  return(Rfast::rowsums(sort.p.log[,1:m]))
}
#################CE mixed gradient#######################3
CE.mtop.mixed.grad=function(K=50,theta,N1=1000,m,idx) {
  grad=1/((1:K)^idx+1)
  
  
  mix.prop=matrix(rbinom(K*N1,size=1,prob=grad),nrow=K)
  
  mix.prop=transpose(mix.prop)
  
  theta_idx_vec=1/sqrt((log(1:K+2)))
  theta_vec=theta^theta_idx_vec
  
  
  x=transpose(matrix(rnorm(K*N1,0,theta_vec),nrow=K))^{(mix.prop)}+matrix(Rnorm(K*N1),nrow = N1)^{1-(mix.prop)}-1
  
  
  
  p=2*pnorm(abs(x),mean=0,sd=1,lower.tail = FALSE)
  mtop_input=Rfast::rowSort(p)
  
  stat.val=mtop(sort.p = mtop_input,m)
  
  
  log_g_norm <- rowsums(log(dnorm(x))) # numerator part
  
  
  
  
  log_g_mix <- rowsums(log(transpose(dnorm(transpose(x),0,theta_vec)*grad)+transpose(transpose(dnorm(x))*(1-grad))))
  
  
  w=exp(log_g_norm-log_g_mix)
  return(list(stat.val=stat.val,weight=w,x=x))
}
mixed.par.update.grad=function(object,ro,K,idx) {
  
  
  
 
  stat.val=object[[1]]
  gamma=quantile(stat.val,1-ro)
  weight=object[[2]]
 
  x.mat=transpose(object[[3]])
  B=length(stat.val)
  grad=1/((1:K)^idx+1)
  theta_idx_vec=1/sqrt((log(1:K+2)))
 
  opt=optimize(function(t) 1/B*sum(weight[stat.val>=gamma]*apply(x.mat[,stat.val>=gamma],2,function(x) sum(log(grad*dnorm(x,0,t^theta_idx_vec)+(1-grad)*dnorm(x,0,1))))),
               interval = c(0,1000),maximum = TRUE)
  par=opt$maximum
  
  return(list(gamma=gamma,par=par))
}



CE.mixed.grad=function(K,ro,N=10^5,q,m,idx) {
  theta=1
  gamma=-Inf
  t=0
  while (gamma<q) {
    t=t+1
    obj=CE.mtop.mixed.grad(K=K,theta=theta,N1=N,m,idx)
   
    par=mixed.par.update.grad(object = obj,ro=ro,K=K,idx)
    gamma=par$gamma
    theta=par$par
    message(paste0(gamma,", ",theta,",","q=",q))
  }
  message("stop")
  obj.final=CE.mtop.mixed.grad(K=K,theta=theta,N1=N,m,idx)
  
  stat.val.final=obj.final[[1]]
  
  weight.final=obj.final[[2]]
  
  p.val=1/N*sum(weight.final[stat.val.final>=q])
  return(list(p.val=p.val,iter=t,theta=theta))
}

libs_mtop.CE.mixGrad=function(q.val.set=exp(seq(log(3),log(22),length.out = 200)),K=50,ro=0.01,N=10^5,m,idx){
  temp.lib=mclapply(q.val.set, function(x) CE.mixed.grad(K=K,ro=ro,N=N,q=x,m,idx),mc.cores = 40)
  temp.p=sapply(temp.lib,"[[",1)
  return(temp.p)
}


