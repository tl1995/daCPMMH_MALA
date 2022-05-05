### Helper functions used in Aphid example

### The deSolve package is required

## Drift and diffusion

#Drift (alpha) at x, given parameters theta
alphaAphid=function(x,theta=c(1.75,0.001))
{
 c(theta[1]*x[1]-theta[2]*x[1]*x[2],theta[1]*x[1])
}

#Diffusion (beta) at x, given parameters theta
betaAphid=function(x,theta=c(1.75,0.001))
{
 cov=theta[1]*x[1]
 var1=theta[1]*x[1]+theta[2]*x[1]*x[2]
 var2=theta[1]*x[1]
 matrix(c(var1,cov,cov,var2),ncol=2,byrow=T)
}


###################################################

## Helper functions for different bridge constructs

###################################################

### Used for all schemes

##Evaluate density under Euler-Maruyama approximation
#
#Arguments:
#x is a bridge
#yT is the observation being bridged towards
#T is the time from the start of the bridge to the observation
#afun and bfun are the drift and diffusion functions
#deltat is the intermediate time step
#theta is a vector of the proposed model parameters
#sige is the observation noise
#F is the observation matrix (set so that C_t is unobserved)
EMjointdensitynoiseMV = function(x, yT, T=1, afun=alphaAphid, bfun=betaAphid, deltat, theta=c(1.75,0.001),sige=1, F=matrix(c(1,0),ncol=1,nrow=2))
{
  n = T/deltat
  ll = 0
  for(i in 1:n){
    ll = ll+lgpdf(x[i+1,], x[i,] + afun(x[i,],theta)*deltat, bfun(x[i,],theta)*deltat)
  }
  #Observation equation for state dependent noise regime
  ll=ll+sum(dnorm(yT,t(F)%*%x[n+1,],sige*sqrt(x[n+1,1]),log=T))
  
  return(ll)
}

##############################################################################################

### Myopic forward simulation


##Generate forward simulation from one time point to the next
#
#Arguments:
#x0 is the initial condition of the particle
#deltat is the intermediate time step
#T is the length of time to generate the bridge for
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#bmmat is a d*n matrix of N(0,deltat) quantities (see function for d and n)
aphid_eulersim = function(x0, deltat, T = 1, afun = alphaAphid, bfun = betaAphid, theta = c(1.75,0.001),bmmat)
{
  
  n = T/deltat
  d=length(x0)
  x=matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0

  for(i in 1:n)
  {
    x[i+1,] = x[i,] + afun(x[i,],theta)*deltat + t(chol(bfun(x[i,],theta)*deltat))%*%bmmat[,i]
    for(j in 1:d){
     if(x[i+1,j]<0.0001){
      x[i+1,j]=0.0001
     }
    }
  }
  return(x)
}

##One step of the Particle Filter - forward simulation and weighted resampling
#
#Arguments:
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#yT is the observation (used to compute weights)
#F is the observation matrix (set so that C_t is unobserved)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmarray is an N*d*n array of N(0,deltat) quantities (see function for d and n)
#uni is a U(0,1) quantity to be used in the resampling step
aphid_wr_myopic=function(N,x0,yT=27,
	F=matrix(c(1,0),ncol=1,nrow=2),deltat,T=1,afun=alphaAphid,bfun=betaAphid,
	theta=c(1.75,0.001),sige=1,bmarray,uni)
{
  rejFlag=0
  #bmarray is N*d*n array of N(0,deltat) quantities
  d=dim(x0)[2] #cols = no. components 
  n = T/deltat
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  #propose and calculate weight
  if(N==1){
    bridge=aphid_eulersim(x0,deltat,T,afun,bfun,theta,bmarray)
    wts=exp(sum(dnorm(yT,t(F)%*%bridge[n+1,],sige*sqrt(bridge[n+1,1]),log=T)))
    mat[1,]=bridge[n+1,]
    indices=1
  }else{
   for(i in 1:N)
   {
    bridge=aphid_eulersim(x0[i,],deltat,T,afun,bfun,theta,bmarray[i,,])
    wts[i]=exp(sum(dnorm(yT,t(F)%*%bridge[n+1,],sige*sqrt(bridge[n+1,1]),log=T)))
    mat[i,]=bridge[n+1,]
   }
   if(sum(wts)<=0)
   {
    rejFlag=1; indices=1:N
   }else{ 
   sorted=esort(mat) #Euclidean sorting
   #print(sorted)
   mat=mat[sorted,]; wts=wts[sorted]
   #systematic resampling
   indices=sysresamp2(wts,N,uni)
   #print(indices) 
   }
 } 
 return(list(mat[indices,],mean(wts),mat,ess(wts),rejFlag))
}

##############################################################################################

### RBpart

## Standard (C)PMMH

#Model function for eta
aphidEtaFunc <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11 #eta_{N,t}
    dy12 <- c1*y11          #eta_{C,t}

    return(list(c(dy11, dy12)))
  })
}


##Solver for eta
#
#Arguments:
#xi is the initial condition
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
aphid_eta_solve=function(xi,inter,deltat,param=c(1.75,0.001))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2])
  odepars = c(c1 = param[1], c2 = param[2])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, aphidEtaFunc, odepars)
  return(odeOut)
}


##Generate bridge and evaluate density under proposal
#
#Arguments:
#x0 is the initial condition of the particle
#yT is the observation to bridge towards
#F is the observation matrix (set so that C_t is unobserved)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmmat is a d*n matrix of N(0,deltat) quantities (see function for d and n)
#odeOut is output from aphid_eta_solve
aphid_rbridgenoiseMVInc = function(x0, yT, F=matrix(c(1,0),ncol=1,nrow=2), deltat, T = 1, afun = alphaAphid, bfun = betaAphid, theta = c(1.75,0.001),sige=1.0,bmmat, odeOut)
{
  
  n = T/deltat
  d=length(x0)
  ll = 0
  x=matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0
  R=matrix(0,ncol=d,nrow=n+1)
  #ODE setup from solution argument
  eta=matrix(odeOut[1:(n+1),c("y11","y12")],ncol=2)

  for(i in 1:n)
  {
    #initial setup of time increments and quantities that don't require ODE solution
    tm = (i-1)*deltat
    delT=T-tm
    ai=afun(x[i,],theta); bi=bfun(x[i,],theta)
    #Inverse variance matrix for state dependent noise regime
    inv=solve(t(F)%*%bi%*%F*delT+eta[n+1,1]*sige^2)
    
    #Put it all together
    mu=ai-((eta[i+1,]-eta[i,])/deltat)+bi%*%F%*%inv%*%(yT-t(F)%*%(eta[n+1,]+(R[i,]+(ai-((eta[i+1,]-eta[i,])/deltat))*delT)))
    
    psi=bi-bi%*%F%*%inv%*%t(F)%*%bi*deltat
    
    R[i+1,] = R[i,] + mu*deltat + t(chol(psi))%*%bmmat[,i]
    x[i+1,] = eta[i+1,]+R[i+1,]
    for(j in 1:d){
     if(x[i+1,j]<0.0001){ #Check that bridge never gets too close to 0
      x[i+1,j]=0.0001
     }
    }
    ll = ll+lgpdf(R[i+1,], R[i,] + mu*deltat, psi*deltat)
  }
  return(list(x=x, ll=ll)) 
}


##One step of the Particle Filter - bridge propagation and weighted resampling
#
#Arguments:
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#yT is the observation to bridge towards
#F is the observation matrix (set so that C_t is unobserved)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmarray is an N*d*n array of N(0,deltat) quantities (see function for d and n)
#uni is a U(0,1) quantity to be used in the resampling step
aphid_wrresbridgenoiseMVInc=function(N,x0,yT,
	F=matrix(c(1,0),ncol=1,nrow=2),deltat,T=1,afun=alphaAphid,bfun=betaAphid,
	theta=c(1.75,0.001),sige=1,bmarray,uni)
{
  rejFlag=0
  
  d=dim(x0)[2] #cols = no. components 
  n = T/deltat
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  #propose and calculate weight
  if(N==1){
    #Resolving eta at particle
    newEta=aphid_eta_solve(x0,T,deltat,theta)[,c("y11","y12")]

    temp=aphid_rbridgenoiseMVInc(x0,yT,F,deltat,T,afun,bfun,theta,sige,bmarray, newEta)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts=exp(EMjointdensitynoiseMV(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[1,]=bridge[n+1,]
    indices=1
  }else{
   for(i in 1:N)
   {
    #Resolving eta at each particle
    newEta=aphid_eta_solve(x0[i,],T,deltat,theta)[,c("y11","y12")]

    temp=aphid_rbridgenoiseMVInc(x0[i,],yT,F,deltat,T,afun,bfun,theta,sige,bmarray[i,,], newEta)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts[i]=exp(EMjointdensitynoiseMV(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[i,]=bridge[n+1,]
   }
   if(sum(wts)<=0) #Rare case where all particles have negligible weight, resampling breaks
   {
    rejFlag=1; indices=1:N
   }else{ 
   	sorted=esort(mat) #Euclidean sorting
   	mat=mat[sorted,]; wts=wts[sorted]
   	#systematic resampling
   	indices=sysresamp2(wts,N,uni) 
   }
 } 
 return(list(mat[indices,],mean(wts),mat,ess(wts),rejFlag))
}

##############################################################################################

### RBiter

## Standard (C)PMMH

#Model function for eta and V
aphidFunc <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11 #eta_{N,t}
    dy12 <- c1*y11          #eta_{C,t}
    
    dy21 <- 2*(c1-c2*y12)*y21 - 2*c2*y11*y2C + c1*y11 + c2*y11*y12 #V_11
    dy2C <- c1*(y11+y21) - c2*y11*y22 + y2C*(c1-c2*y12)            #V_12 = V_21
    dy22 <- c1*(y11+2*y2C)                                         #V_22
    
    return(list(c(dy11, dy12, dy21, dy2C, dy22)))
  })
}


#LNA solver for eta and V
#
#Arguments:
#xi is the initial condition for eta
#B0 is the initial condition for V
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
aphid_lna_solve_correct=function(xi,B0,inter,deltat,param=c(1.75,0.001))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = B0[1,1], y2C = B0[1,2], y22 = B0[2,2])
  odepars = c(c1 = param[1], c2 = param[2])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, aphidFunc, odepars)
  return(odeOut)
}

#Forward filter for the LNA - stores solutions in an array
#
#Arguments:
#data is all observations
#inter is the total time between two observations
#deltat is the intermediate time step between two observations
#param is a vector of model parameters
#a is the initial condition for eta
#B is the inital condition for V
#sige is the observation noise
#F is the observation matrix (set so that C_t is unobserved)
aphid_lna_filter2=function(data,inter,deltat,param=c(1.75,0.001), a=c(5,5),B=diag(c(1.0,1.0)), sige=1, F=matrix(c(1,0),ncol=1,nrow=2))
{
  n=length(data) #number of observations
  n2=inter/deltat+1 #number of timepoints between observations (including endpoints)
  solarray=array(0,dim=c(n2,5,(n-1)))
  colnames(solarray)=c("y11","y12","y21","y2C","y22")
  #update posterior
  #State dependent noise regime
  a = a + B%*%F%*%solve(t(F)%*%B%*%F+a[1]*sige^2)%*%(data[1]-t(F)%*%a)
  B = B - B%*%F%*%solve(t(F)%*%B%*%F+a[1]*sige^2)%*%t(F)%*%B

  #loop
  for(i in 2:n)
  {
    #update mu, sigma, sensitivities and G, and store solution
    sensOut=aphid_lna_solve_correct(a,B,inter,deltat,param)
    solarray[,,(i-1)]=sensOut[,c("y11","y12","y21","y2C","y22")]
    mu=sensOut[n2,c("y11","y12")]
    sigma=matrix(sensOut[n2,c("y21", "y2C", "y2C", "y22")], ncol = 2)

    #update posterior
    #State dependent noise regime
    a = mu + sigma%*%F%*%solve(t(F)%*%sigma%*%F+t(F)%*%mu*sige^2)%*%(data[i]-t(F)%*%mu)
    B = sigma - sigma%*%F%*%solve(t(F)%*%sigma%*%F+t(F)%*%mu*sige^2)%*%t(F)%*%sigma

  }
  return(list(solarray = solarray))
}

##Generate bridge and evaluate density under proposal
#
#Arguments:
#x0 is the initial condition of the particle
#yT is the observation to bridge towards
#F is the observation matrix (set so that C_t is unobserved)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmmat is a d*n matrix of N(0,deltat) quantities (see function for d and n)
#odeOut is output from aphid_lna_filter2
aphid_rbridgenoiseMVInc2 = function(x0, yT, F=matrix(c(1,0),ncol=1,nrow=2), deltat, T = 1, afun = alphaAphid, bfun = betaAphid, theta = c(1.75,0.001),sige=1.0,bmmat, odeOut)
{
  
  n = T/deltat
  d=length(x0)
  ll = 0
  x=matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0
  R=matrix(0,ncol=d,nrow=n+1)
  #ODE setup from solution argument
  eta=matrix(odeOut[1:(n+1),c("y11","y12")],ncol=2)
  R[1,]=x[1,]-eta[1,]

  for(i in 1:n)
  {
    #initial setup of time increments and quantities that don't require ODE solution
    tm = (i-1)*deltat
    delT=T-tm
    ai=afun(x[i,],theta); bi=bfun(x[i,],theta)
    #Inverse variance matrix for state dependent noise regime
    inv=solve(t(F)%*%bi%*%F*delT+eta[n+1,1]*sige^2)
    

    #Put it all together
    mu=ai-((eta[i+1,]-eta[i,])/deltat)+bi%*%F%*%inv%*%(yT-t(F)%*%(eta[n+1,]+(R[i,]+(ai-((eta[i+1,]-eta[i,])/deltat))*delT)))
    
    psi=bi-bi%*%F%*%inv%*%t(F)%*%bi*deltat
    
    R[i+1,] = R[i,] + mu*deltat + t(chol(psi))%*%bmmat[,i]
    x[i+1,] = eta[i+1,]+R[i+1,]
    for(j in 1:d){
     if(x[i+1,j]<0.0001){ #Check that bridge never gets too close to 0
      x[i+1,j]=0.0001
     }
    }
    ll = ll+lgpdf(R[i+1,], R[i,] + mu*deltat, psi*deltat)
  }
  return(list(x=x, ll=ll))
}

##One step of the Particle Filter - bridge propagation and weighted resampling
#
#Arguments:
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#yT is the observation to bridge towards
#F is the observation matrix (set so that C_t is unobserved)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmarray is an N*d*n array of N(0,deltat) quantities (see function for d and n)
#uni is a U(0,1) quantity to be used in the resampling step
#odeOut is output from aphid_lna_filter2
aphid_wrresbridgenoiseMVInc2=function(N,x0,yT,F=matrix(c(1,0),ncol=1,nrow=2),deltat,T=1,
	afun=alphaAphid,bfun=betaAphid,theta=c(1.75,0.001),sige=1,bmarray,uni, odeOut)
{
  rejFlag=0
  
  d=dim(x0)[2] #cols = no. components 
  n = T/deltat
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  #propose and calculate weight
  if(N==1){
    temp=aphid_rbridgenoiseMVInc2(x0,yT,F,deltat,T,afun,bfun,theta,sige,bmarray, odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts=exp(EMjointdensitynoiseMV(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[1,]=bridge[n+1,]
    indices=1
  }else{
   for(i in 1:N)
   {
    temp=aphid_rbridgenoiseMVInc2(x0[i,],yT,F,deltat,T,afun,bfun,theta,sige,bmarray[i,,], odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts[i]=exp(EMjointdensitynoiseMV(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[i,]=bridge[n+1,]
   }
   if(sum(wts)<=0) #Rare case where all particles have negligible weight, resampling breaks
   {
    rejFlag=1; indices=1:N
   }else{ 
   sorted=esort(mat) #Euclidean sorting
   
   mat=mat[sorted,]; wts=wts[sorted]
   #systematic resampling
   indices=sysresamp2(wts,N,uni)
   }
 } 
 return(list(mat[indices,],mean(wts),mat,ess(wts),rejFlag))
}

