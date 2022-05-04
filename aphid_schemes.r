### Schemes used in Aphid example

## Load helper functions and required library

#install.packages("deSolve") ## Uncomment to install "deSolve" package if required 
library(deSolve)
source("general_helpers.r")
source("aphid_helpers.r")

## Load in data and tuning matrix

sigtuneAphidBench = matrix(scan("sigtuneAphidBench.dat"),ncol=2)

dataAphid = scan("dataAphid.dat")

### Myopic simulation
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
#F is the observation matrix (set so that C_t is unobserved)
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
aphid_myopicCPMMH=function(iters,sigma=sigtuneAphidBench,N,x0,
	x=dataAphid,a=0,b=10,init=c(1.75,0.001),deltat,inter=1,afun=alphaAphid,
	bfun=betaAphid,sige=1,rho,F=matrix(c(1,0),ncol=1,nrow=2), Ntuning=0)
{
  n=length(x) #no. obs
  d=2 #total no. components - observed and unobserved
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
    for(j in 1:(n-1))
    {
 wrstuff=aphid_wr_myopic(N,xi,x[j+1],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j])
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



### RBpart
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
#F is the observation matrix (set so that C_t is unobserved)
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
aphid_resCPMMHnoiseMV=function(iters,sigma=sigtuneAphidBench,N,x0,x=dataAphid,
	a=0,b=10,init=c(1.75,0.001),deltat,inter=1,afun=alphaAphid,bfun=betaAphid,
	sige=1,rho,F=matrix(c(1,0),ncol=1,nrow=2), Ntuning=0)
{
  n=length(x) #no. obs
  d=2 #total no. components - observed and unobserved
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
    for(j in 1:(n-1))
    {
     wrstuff=aphid_wrresbridgenoiseMVInc(N,xi,x[j+1],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j])
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

### RBiter
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
#F is the observation matrix (set so that C_t is unobserved)
#Ntuning is an indicator variable indicating whether the parameter proposals should be fixed to tune the number of bridges
aphid_resCPMMHnoiseMV2=function(iters,sigma,N,x0,x=dataAphid,
	a=0,b=10,init=c(1.75,0.001),deltat,inter=1,afun=alphaAphid,bfun=betaAphid,
	sige=1,rho,F=matrix(c(1,0),ncol=1,nrow=2), Ntuning=0)
{
  n=length(x) #no. obs
  d=2 #total no. components - observed and unobserved
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
    solarray=aphid_lna_filter2(data=x,inter=inter,deltat=deltat,param=exp(can), a=x0[1,],B=diag(c(0.0,0.0)), sige=sige, F=matrix(c(1,0),ncol=1,nrow=2))$solarray
    for(j in 1:(n-1))
    {
     wrstuff=aphid_wrresbridgenoiseMVInc2(N,xi,x[j+1],F,deltat,inter,afun,bfun,exp(can),sige,bmarrayprop[,,,j],uprop[j],solarray[,,j])
     if(wrstuff[[5]]==1)
     {
      loglikecan=-1e8
      j=n
     }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
     #print(x0)
     #print(log(wrstuff[[2]]))
     }
    }
    laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)
    u = runif(1)
    #print(c(loglikecan,loglikecurr))
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
