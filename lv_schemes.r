### Schemes used in Aphid example

## Load helper functions and required library

#install.packages("deSolve") ## Uncomment to install "deSolve" package if required 
library(deSolve)
source("general_helpers.r")
source("lv_helpers.r")

## Load in data and tuning matrix

sigtuneLVBench = matrix(scan("sigtuneLVBench.dat"),ncol=3)

dataLV = scan("dataLV.dat")

### RB^- iter standard (C)PMMH
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusCPMMH=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLV,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,
	sige=1,rho,F=diag(c(1,1)),Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  loglikecurr=-1000000 #log-likelihood of current value (to accept first candidate)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else{
    can = rmvn(curr,sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    loglikecan=0
    solarray=lv_lna_filter2(data=x,inter=inter,deltat=deltat,param=exp(can), a=x0[1,],B=diag(c(0.0,0.0)), sige=sige, F=F)$solarray
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j], solarray[,,j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)
    u = runif(1)
    
    if (log(u) < laprob)
    { 
      curr = can
      loglikecurr=loglikecan 
      bmarray=bmarrayprop
      uvec=uprop
      count=count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/(iters-1))
  return(mat[2:iters,])
}



### RB^- iter simplified MALA (C)PMMH
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#eps is the step size for MALA
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusMALA_CPMMH=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLV,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,sige=1,rho,F=diag(c(1,1)),
	eps=1,Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  loglikecurr=-1000000 #log-likelihood of current value (to accept first candidate)
  gradcurr=c(0,0,0) #gradient at current value (first step is rwm step)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else{
    can = rmvn(curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    loglikecan=0
    MEGASOLUTION=lv_lna_MEGASOL_filter3(data=x,inter=inter,deltat=deltat,param=exp(can),m = a, sd = b, a=x0[1,],B=diag(c(0.0,0.0)), sige=sige,F=F)
    gradprop=MEGASOLUTION$grad
    solarray=MEGASOLUTION$solarray
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j], solarray[,,j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)+lgpdf(curr,can+0.5*eps*sigma%*%gradprop,eps*sigma)-lgpdf(can,curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    u = runif(1)
    if (log(u) < laprob)
    { 
      curr = can
      loglikecurr=loglikecan 
      gradcurr=gradprop  
      bmarray=bmarrayprop
      uvec=uprop
      count=count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/(iters-1))
  return(mat[2:iters,])
}



### RB^- iter delayed acceptance (C)PMMH (with RWM or simplified MALA)
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#eps is the step size for MALA (if used)
#MALA is an indicator variable indicating whether MALA is used
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusdaMALA_CPMMH=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLV,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,sige=1,rho,F=diag(c(1,1)),
	eps=1, MALA=FALSE,Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  lna_loglikecurr=-1000000 #approximate log-likelihood of current value (to accept first candidate)
  loglikecurr=-10000000 #log-likelihood of current value (to accept first candidate)
  gradcurr=c(0,0,0) #gradient at current value (first step is rwm step)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count1=0
  count2=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else if(MALA==TRUE){
    can = rmvn(curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    else{
      can = rmvn(curr,sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    MEGASOLUTION=lv_lna_MEGASOL_filter3(data=x,inter=inter,deltat=deltat,param=exp(can),m = a, sd = b, a=x0[1,],B=diag(c(0.0,0.0)), sige=sige,F=F, MALA=MALA)
    gradprop=MEGASOLUTION$grad
    solarray=MEGASOLUTION$solarray
    lna_loglikecan=MEGASOLUTION$llike
    if(MALA==TRUE){
      laprob1 = lna_loglikecan+logprior(can,a,b)-lna_loglikecurr-logprior(curr,a,b)+lgpdf(curr,can+0.5*eps*sigma%*%gradprop,eps*sigma)-lgpdf(can,curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    else{
      laprob1 = lna_loglikecan+logprior(can,a,b)-lna_loglikecurr-logprior(curr,a,b)
    }
    u1 = runif(1)
    if (log(u1) < laprob1){
     loglikecan=0
     count1=count1+1 #track no. Stage 1 acceptances
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j], solarray[,,j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob2 = loglikecan-loglikecurr+lna_loglikecurr-lna_loglikecan
    u2 = runif(1)
    if (log(u2) < laprob2)
    { 
      curr = can
      loglikecurr=loglikecan 
      lna_loglikecurr=lna_loglikecan
      gradcurr=gradprop  
      bmarray=bmarrayprop
      uvec=uprop
      count2=count2+1  #track no. Stage 2 acceptances
    }
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(paste("Stage 1 acceptance rate: ", count1/(iters-1)))
  print(paste("Stage 2 acceptance rate: ", count2/count1))
  return(mat[2:iters,])
}

##############################################################################################


### RB^- part standard (C)PMMH
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusCPMMH2=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLVnoise,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,sige=1,rho,F=diag(c(1,1)),
	Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  loglikecurr=-1000000 #log-likelihood of current value (to accept first candidate)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else{
    can = rmvn(curr,sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    loglikecan=0
    #No filter as we are resolving at each particle, with no MALA or DA
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge2(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)
    u = runif(1)
    if (log(u) < laprob)
    { 
      curr = can
      loglikecurr=loglikecan 
      bmarray=bmarrayprop
      uvec=uprop
      count=count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/(iters-1))
  return(mat[2:iters,])
}

### RB^- part simplified MALA (C)PMMH
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#eps is the step size for MALA
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusMALA_CPMMH2=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLV,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,sige=1,rho,F=diag(c(1,1)),
	eps=1,Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  loglikecurr=-1000000 #log-likelihood of current value (to accept first candidate)
  gradcurr=c(0,0,0) #gradient at current value (first step is rwm step)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else{
    can = rmvn(curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    loglikecan=0
    MEGASOLUTION=lv_lna_MEGASOL_filter(data=x,inter=inter,deltat=deltat,param=exp(can),m = a, sd = b, a=x0[1,],B=diag(c(0.0,0.0)), sige=sige,F=F)
    gradprop=MEGASOLUTION$grad
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge2(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)+lgpdf(curr,can+0.5*eps*sigma%*%gradprop,eps*sigma)-lgpdf(can,curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    u = runif(1)
    if (log(u) < laprob)
    { 
      curr = can
      loglikecurr=loglikecan 
      gradcurr=gradprop  
      bmarray=bmarrayprop
      uvec=uprop
      count=count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/(iters-1))
  return(mat[2:iters,])
}


### RB^- part delayed acceptance (C)PMMH (with RWM or simplified MALA)
#
#Arguments:
#iters is number of iterations
#sigma is the innovation covariance matrix
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#x is the data
#a and b are the parameter prior mean and variance
#init is the initialisation of the Markov chain
#deltat is the intermediate time step between observations
#inter is the total time between observations
#afun and bfun are the drift and diffusion functions
#sige is the observation noise
#rho is the correlation parameter for CPMMH (rho = 0 corresponds to PMMH)
#F is the observation matrix (set for complete observation of both species)
#eps is the step size for MALA (if used)
#MALA is an indicator variable indicating whether MALA is used
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
lv_resminusdaMALA_CPMMH2=function(iters,sigma=sigtuneLVBench,N,x0,x=dataLV,a=0,b=10,
	init=c(0.5,0.0025,0.3),deltat,inter=1,afun=alphaLV,bfun=betaLV,sige=1,rho,F=diag(c(1,1)),
	eps=1, MALA=FALSE,Ntuning=0)
{
  n=dim(x)[1] #no. obs
  d=dim(x0)[2] #total no. components - observed and unobserved
  n2=(inter/deltat) #no. BM increments in each interval 
  bmarray=array(rnorm(N*d*n2*(n-1),0,sqrt(deltat)),dim=c(N,d,n2,n-1))
  uvec=runif(n-1)
  p=length(init)  #no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  curr=log(init) #current param value (log-scale)
  lna_loglikecurr=-1000000 #approximate log-likelihood of current value (to accept first candidate)
  loglikecurr=-10000000 #log-likelihood of current value (to accept first candidate)
  gradcurr=c(0,0,0) #gradient at current value (first step is rwm step)
  mat[1,]=c(curr,0) #params on log-scale and candidate log-like
  count1=0
  count2=0
  for (i in 2:iters) {
    xi=x0
    if(Ntuning==1){
    can=curr #only need for looking at variance of estimator
    }
    else if(MALA==TRUE){
    can = rmvn(curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    else{
      can = rmvn(curr,sigma)
    }
    bmarrayprop=rho*bmarray+sqrt( 1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
    uprop=pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
    MEGASOLUTION=lv_lna_MEGASOL_filter(data=x,inter=inter,deltat=deltat,param=exp(can),m = a, sd = b, a=x0[1,],B=diag(c(0.0,0.0)), sige=sige,F=F, MALA=MALA)
    gradprop=MEGASOLUTION$grad
    lna_loglikecan=MEGASOLUTION$llike
    if(MALA==TRUE){
      laprob1 = lna_loglikecan+logprior(can,a,b)-lna_loglikecurr-logprior(curr,a,b)+lgpdf(curr,can+0.5*eps*sigma%*%gradprop,eps*sigma)-lgpdf(can,curr+0.5*eps*sigma%*%gradcurr,eps*sigma)
    }
    else{
      laprob1 = lna_loglikecan+logprior(can,a,b)-lna_loglikecurr-logprior(curr,a,b)
    }
    u1 = runif(1)
    if (log(u1) < laprob1){
     loglikecan=0
     count1=count1+1 #track no. Stage 1 acceptances
    for(j in 1:(n-1))
    {
     wrstuff=lv_wrresminusbridge2(N,xi,x[j+1,],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     }
    }
    laprob2 = loglikecan-loglikecurr+lna_loglikecurr-lna_loglikecan
    u2 = runif(1)
    if (log(u2) < laprob2)
    { 
      curr = can
      loglikecurr=loglikecan 
      lna_loglikecurr=lna_loglikecan
      gradcurr=gradprop  
      bmarray=bmarrayprop
      uvec=uprop
      count2=count2+1  #track no. Stage 2 acceptances
    }
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(paste("Stage 1 acceptance rate: ", count1/(iters-1)))
  print(paste("Stage 2 acceptance rate: ", count2/count1))
  return(mat[2:iters,])
}
