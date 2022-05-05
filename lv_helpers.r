### Helper functions used in Lotka-Volterra example

### The deSolve package is required

## Drift and diffusion

#Drift (alpha) at x, given parameters theta
alphaLV=function(x,theta=c(0.5,0.0025,0.3))
{
 c(theta[1]*x[1]-theta[2]*x[1]*x[2],theta[2]*x[1]*x[2]-theta[3]*x[2])
}

#Diffusion (beta) at x, given parameters theta
betaLV=function(x,theta=c(0.5,0.0025,0.3))
{
 cov=-theta[2]*x[1]*x[2]
 var1=theta[1]*x[1]+theta[2]*x[1]*x[2]
 var2=theta[2]*x[1]*x[2]+theta[3]*x[2]
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
#F is the observation matrix (set for complete observation of both species)
EMjointdensity = function(x, yT, T=1, afun=alphaLV, bfun=betaLV, deltat, theta=c(0.5,0.0025,0.3),sige=1, F=diag(c(1,1)))
{
  n = T/deltat
  ll = 0
  for(i in 1:n){
    ll = ll+lgpdf(x[i+1,], x[i,] + afun(x[i,],theta)*deltat, bfun(x[i,],theta)*deltat)
  }
  ll=ll+sum(dnorm(yT,t(F)%*%x[n+1,],sige,log=T))
  return(ll)
}

##############################################################################################

### RB^- iter

## Standard (C)PMMH

#Model function for eta, V and G
lvFuncG <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11 #eta_1
    dy12 <- (c2*y11-c3)*y12 #eta_2
    
    dy21 <- 2*(c1-c2*y12)*y21 - 2*c2*y11*y2C + c1*y11 + c2*y11*y12          #V_11
    dy2C <- (c1 + c2*(y11-y12) - c3)*y2C + c2*(y12*y21 - y11*y22 - y11*y12) #V_12 == V_21
    dy22 <- 2*(c2*y11-c3)*y22 + 2*c2*y12*y2C + c3*y12 + c2*y11*y12          #V_22
    
    dG1 <- (c1-c2*y12)*G1-c2*y11*G2   
    dG2 <- c2*y12*G1 + (c2*y11-c3)*G2 
    dG3 <- (c1-c2*y12)*G3-c2*y11*G4   
    dG4 <- c2*y12*G3 + (c2*y11-c3)*G4

    return(list(c(dy11, dy12, dy21, dy2C, dy22, dG1, dG2, dG3, dG4)))
  })
}

#LNA solver for eta, V and G
#
#Arguments:
#xi is the initial condition for eta
#B0 is the initial condition for V
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
lv_lna_solve2=function(xi,B0,inter,deltat,param=c(0.5,0.0025,0.3))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = B0[1,1], y2C = B0[1,2], y22 = B0[2,2],
              G1 = 1, G2 = 0, G3 = 0, G4 = 1)
  odepars = c(c1 = param[1], c2 = param[2], c3 = param[3])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, lvFuncG, odepars)
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
#F is the observation matrix (set for complete observation of both species)
lv_lna_filter2=function(data,inter,deltat,param=c(0.5,0.0025,0.3), a=c(100,100),B=diag(c(0.0,0.0)), sige=1.0, F=diag(c(1,1)))
{
  n=dim(data)[1] #number of observations
  n2=inter/deltat+1 #number of timepoints between observations (including endpoints)
  solarray=array(0,dim=c(n2,9,(n-1)))
  colnames(solarray)=c("y11","y12","y21","y2C","y22","G1","G2","G3","G4")
  #update posterior
  a = a + B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%(data[1,]-t(F)%*%a)
  B = B - B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%B
  #loop
  for(i in 2:n)
  {
    #update mu, sigma, and store solution
    sensOut=lv_lna_solve2(a,B,inter,deltat,param)
    solarray[,,(i-1)]=sensOut[,c("y11","y12","y21","y2C","y22","G1","G2","G3","G4")]
    mu=sensOut[n2,c("y11","y12")]
    sigma=matrix(sensOut[n2,c("y21", "y2C", "y2C", "y22")], ncol = 2)

    #update posterior
    a = mu + sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%(data[i,]-t(F)%*%mu)
    B = sigma - sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%sigma
  }
  return(list(solarray = solarray))
}

##Generate bridge and evaluate density under proposal
#
#Arguments:
#x0 is the initial condition of the particle
#yT is the observation to bridge towards
#F is the observation matrix (set for complete observation of both species)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmmat is a d*n matrix of N(0,deltat) quantities (see function for d and n)
#odeOut is output from lv_lna_filter2
lv_rminusbridge = function(x0, yT, F=diag(c(1,1)), deltat, T = 1, afun = alphaLV, bfun = betaLV, theta = c(0.5,0.0025,0.3),sige=1.0,bmmat, odeOut)
{
  
  n = T/deltat
  d=length(x0)
  ll = 0
  x=matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0
  R=matrix(0,ncol=d,nrow=n+1)
  #ODE setup from solution argument
  eta=matrix(odeOut[1:(n+1),c("y11","y12")],ncol=2)
  rmean=matrix(0,ncol=d,nrow=n+1)
  rmean[1,]=x[1,]-eta[1,]
  R[1,]=0

  V = array(as.vector(t(odeOut[1:(n+1),c("y21", "y2C", "y2C", "y22")])), dim = c(2,2,n+1))
  G = array(as.vector(t(odeOut[1:(n+1),c("G1", "G2", "G3", "G4")])), dim=c(2,2,n+1))
  
  rmean[n+1,]=G[,,n+1]%*%rmean[1,]+V[,,n+1]%*%F%*%solve(t(F)%*%V[,,n+1]%*%F+diag(rep(sige*sige,2)))%*%(yT-t(F)%*%eta[n+1,]-t(F)%*%G[,,(n+1)]%*%rmean[1,])
  for(i in 1:n)
  {
    #initial setup of time increments and quantities that don't require ODE solution
    tm = (i-1)*deltat
    delT=T-tm
    ai=afun(x[i,],theta); bi=bfun(x[i,],theta)
    inv=solve(t(F)%*%bi%*%F*delT+diag(rep(sige*sige,d)))

    #rmean at next time point
    rmean[i+1,]=G[,,i+1]%*%rmean[1,]+V[,,i+1]%*%t(solve(G[,,i+1]))%*%t(G[,,n+1])%*%F%*%solve(t(F)%*%V[,,n+1]%*%F+diag(rep(sige*sige,2)))%*%(yT-t(F)%*%eta[n+1,]-t(F)%*%G[,,(n+1)]%*%rmean[1,])

    #Put it all together
    mu=ai-((eta[i+1,]-eta[i,])/deltat)-((rmean[i+1,]-rmean[i,])/deltat)+bi%*%F%*%inv%*%(yT-t(F)%*%(eta[n+1,]+rmean[n+1,]+(R[i,]+(ai-((eta[i+1,]-eta[i,])/deltat)-((rmean[i+1,]-rmean[i,])/deltat))*delT)))
    psi=bi-bi%*%F%*%inv%*%t(F)%*%bi*deltat
    R[i+1,] = R[i,] + mu*deltat + t(chol(psi))%*%bmmat[,i]
    x[i+1,] = eta[i+1,]+rmean[i+1,]+R[i+1,]
    for(j in 1:d){
     if(x[i+1,j]<0.0001){
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
#F is the observation matrix (set for complete observation of both species)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmarray is an N*d*n array of N(0,deltat) quantities (see function for d and n)
#uni is a U(0,1) quantity to be used in the resampling step
#odeOut is output from lv_lna_filter2
lv_wrresminusbridge=function(N,x0,yT,F=diag(c(1,1)),deltat,T=1,afun=alphaLV,bfun=betaLV,theta=c(0.5,0.0025,0.3),sige=1.0,bmarray,uni, odeOut)
{
  rejFlag=0
  d=dim(x0)[2] #cols = no. components 
  n = T/deltat
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  #propose and calculate weight
  if(N==1){
    temp=lv_rminusbridge(x0,yT,F,deltat,T,afun,bfun,theta,sige,bmarray, odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts=exp(EMjointdensity(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[1,]=bridge[n+1,]
    indices=1
  }else{
   for(i in 1:N)
   {
    temp=lv_rminusbridge(x0[i,],yT,F,deltat,T,afun,bfun,theta,sige,bmarray[i,,], odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts[i]=exp(EMjointdensity(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
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


## additional helpers for MALA and delayed acceptance

#Model function for eta, V, G, and drift Sensitivities
lvFuncGS <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11 #eta_1
    dy12 <- (c2*y11-c3)*y12 #eta_2
    
    dy21 <- 2*(c1-c2*y12)*y21 - 2*c2*y11*y2C + c1*y11 + c2*y11*y12          #V_11
    dy2C <- (c1 + c2*(y11-y12) - c3)*y2C + c2*(y12*y21 - y11*y22 - y11*y12) #V_12 == V_11
    dy22 <- 2*(c2*y11-c3)*y22 + 2*c2*y12*y2C + c3*y12 + c2*y11*y12          #V_22
        
    dG1 <- (c1-c2*y12)*G1-c2*y11*G2
    dG2 <- c2*y12*G1 + (c2*y11-c3)*G2
    dG3 <- (c1-c2*y12)*G3-c2*y11*G4
    dG4 <- c2*y12*G3 + (c2*y11-c3)*G4

    dS11_1 <- (c1-c2*y12)*S11_1 - c2*y11*S12_1 + y11 #Sensitivities for eta w.r.t. c1
    dS12_1 <- c2*y12*S11_1 + (c2*y11-c3)*S12_1
    
    dS11_2 <- (c1-c2*y12)*S11_2 - c2*y11*S12_2 - y12*y11 #Sensitivities for eta w.r.t. c2
    dS12_2 <- c2*y12*S11_2 + (c2*y11-c3)*S12_2 + y11*y12
        
    dS11_3 <- (c1-c2*y12)*S11_3 - c2*y11*S12_3       #Sensitivities for eta w.r.t. c3
    dS12_3 <- c2*y12*S11_3 + (c2*y11-c3)*S12_3 - y12

    return(list(c(dy11, dy12, dy21, dy2C, dy22, dG1, dG2, dG3, dG4,
		  dS11_1, dS12_1, dS11_2, dS12_2, dS11_3, dS12_3)))
  })
}

#LNA solver for eta, V, G, and drift sensitivities
#
#Arguments:
#xi is the initial condition for eta
#B0 is the initial condition for V
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
lv_lna_solve_MALA2=function(xi,B0,inter,deltat,param=c(0.5,0.0025,0.3))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = B0[1,1], y2C = B0[1,2], y22 = B0[2,2],
              G1 = 1, G2 = 0, G3 = 0, G4 = 1,
               S11_1 = 0, S12_1 = 0, S11_2 = 0, S12_2 = 0, S11_3 = 0, S12_3 = 0)
  odepars = c(c1 = param[1], c2 = param[2], c3 = param[3])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, lvFuncGS, odepars)
  return(odeOut)
}

#Forward filter for the LNA - stores gradient and likelihood under LNA, and an array of ODE solutions
#
#Arguments:
#data is all observations
#inter is the total time between two observations
#deltat is the intermediate time step between two observations
#param is a vector of model parameters
#m and sd are the prior mean of the (log) params, used for gradient calculation
#a is the initial condition for eta
#B is the inital condition for V
#sige is the observation noise
#F is the observation matrix (set for complete observation of both species)
#MALA is a boolean indicating whether the gradient should be calculated
lv_lna_MEGASOL_filter3=function(data,inter,deltat,param=c(0.5,0.0025,0.3),m = 0, sd = 10, a=c(100,100),B=diag(c(0.0,0.0)), sige=1.0, F=diag(c(1,1)), MALA=TRUE)
{
  n=dim(data)[1] #number of observations
  n2=inter/deltat+1 #number of timepoints between observations (including endpoints)

  solarray=array(0,dim=c(n2,9,(n-1)))
  colnames(solarray)=c("y11","y12","y21","y2C","y22","G1","G2","G3","G4")

  if(MALA==TRUE){
    gradthet= -(log(param)-m)/sd^2 #prior
  }
  else{
    gradthet=0
  }
  llike = 0
  #update marginal likelihood
  llike = llike + lgpdf(data[1,], t(F)%*%a, t(F)%*%B%*%F + diag(rep(sige*sige,2)))
  #update posterior
  a = a + B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%(data[1,]-t(F)%*%a)
  B = B - B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%B
  #loop
  for(i in 2:n)
  {
    #update mu, sigma, and G (and sensitivities if MALA)
    if(MALA==TRUE){
      sensOut=lv_lna_solve_MALA2(a,B,inter,deltat,param)
    }
    else{
      sensOut=lv_lna_solve2(a,B,inter,deltat,param)
    }
    solarray[,,(i-1)]=sensOut[,c("y11","y12","y21","y2C","y22","G1","G2","G3","G4")]
    mu=sensOut[n2,c("y11","y12")]
    sigma=matrix(sensOut[n2,c("y21", "y2C", "y2C", "y22")], ncol = 2)
    sigmaInv=solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2))) #Only used in gradient contribution so may be partially observed

    #update marginal likelihood
    llike = llike + lgpdf(data[i,], t(F)%*%mu, t(F)%*%sigma%*%F + diag(rep(sige*sige,2)))

    if(MALA==TRUE){
      #add gradient contribution
      geta1=c(sensOut[n2,"S11_1"],sensOut[n2,"S12_1"])
      geta2=c(sensOut[n2,"S11_2"],sensOut[n2,"S12_2"])
      geta3=c(sensOut[n2,"S11_3"],sensOut[n2,"S12_3"])
      C=sigmaInv%*%(data[i,]-t(F)%*%mu)
      gradthet=gradthet+c(t(C)%*%geta1,t(C)%*%geta2,t(C)%*%geta3)
    }

    #update posterior
    a = mu + sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%(data[i,]-t(F)%*%mu)
    B = sigma - sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%sigma
  }
  grad.out = c(gradthet[1]*param[1],gradthet[2]*param[2], gradthet[3]*param[3])
  return(list(grad = grad.out, llike = llike, solarray = solarray))
}

##############################################################################################

### RB^- part

## Standard (C)PMMH

#LNA solver for eta, V, and G that initialises at the particle
#
#Arguments:
#xi is the initial condition for eta
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
lv_lna_solve3=function(xi,inter,deltat,param=c(0.5,0.0025,0.3))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = 0, y2C = 0, y22 = 0,
              G1 = 1, G2 = 0, G3 = 0, G4 = 1)
  odepars = c(c1 = param[1], c2 = param[2], c3 = param[3])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, lvFuncG, odepars)
  return(odeOut)
}

##One step of the Particle Filter - bridge propagation and weighted resampling
#
#Arguments:
#N is the number of bridges
#x0 is an N*d matrix of initial conditions for the bridges (see function for d)
#yT is the observation to bridge towards
#F is the observation matrix (set for complete observation of both species)
#deltat is the intermediate time step
#T is the time from the initial condition to the observation
#afun and bfun are the drift and diffusion functions
#theta is a vector of the proposed model parameters
#sige is the observation noise
#bmarray is an N*d*n array of N(0,deltat) quantities (see function for d and n)
#uni is a U(0,1) quantity to be used in the resampling step
lv_wrresminusbridge2=function(N,x0,yT,F=diag(c(1,1)),deltat,T=1,afun=alphaLV,bfun=betaLV,theta=c(0.5,0.0025,0.3),sige=1.0,bmarray,uni)
{
  rejFlag=0
  d=dim(x0)[2] #cols = no. components 
  n = T/deltat
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  #propose and calculate weight
  if(N==1){
    #Resolving LNA ODEs at particle
    odeOut=lv_lna_solve3(x0,T,deltat,theta)

    temp=lv_rminusbridge(x0,yT,F,deltat,T,afun,bfun,theta,sige,bmarray, odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts=exp(EMjointdensity(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
    mat[1,]=bridge[n+1,]
    indices=1
  }else{
   for(i in 1:N)
   {
    #Resolving LNA ODEs at particle
    odeOut=lv_lna_solve3(x0[i,],T,deltat,theta)

    temp=lv_rminusbridge(x0[i,],yT,F,deltat,T,afun,bfun,theta,sige,bmarray[i,,], odeOut)
    #temp$x is the proposed bridge, temp$ll is its likelihood under this proposal
    bridge=temp$x
    wts[i]=exp(EMjointdensity(bridge,yT,T,afun,bfun,deltat,theta,sige,F)-temp$ll)
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


## additional helpers for MALA and delayed acceptance

#Model function for eta, V, and drift Sensitivities - used for MALA
lvFuncS <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11 #eta_1
    dy12 <- (c2*y11-c3)*y12 #eta_2
    
    dy21 <- 2*(c1-c2*y12)*y21 - 2*c2*y11*y2C + c1*y11 + c2*y11*y12          #V_11
    dy2C <- (c1 + c2*(y11-y12) - c3)*y2C + c2*(y12*y21 - y11*y22 - y11*y12) #V_12 == V_21
    dy22 <- 2*(c2*y11-c3)*y22 + 2*c2*y12*y2C + c3*y12 + c2*y11*y12          #V_22
    
    dS11_1 <- (c1-c2*y12)*S11_1 - c2*y11*S12_1 + y11     #Sensitivities for eta w.r.t. c1
    dS12_1 <- c2*y12*S11_1 + (c2*y11-c3)*S12_1
    
    dS11_2 <- (c1-c2*y12)*S11_2 - c2*y11*S12_2 - y12*y11 #Sensitivities for eta w.r.t. c2
    dS12_2 <- c2*y12*S11_2 + (c2*y11-c3)*S12_2 + y11*y12
        
    dS11_3 <- (c1-c2*y12)*S11_3 - c2*y11*S12_3           #Sensitivities for eta w.r.t. c3
    dS12_3 <- c2*y12*S11_3 + (c2*y11-c3)*S12_3 - y12

    return(list(c(dy11, dy12, dy21, dy2C, dy22,
		  dS11_1, dS12_1, dS11_2, dS12_2, dS11_3, dS12_3)))
  })
}

#LNA solver for eta, V, and drift sensitivities - used for MALA
#
#Arguments:
#xi is the initial condition for eta
#B0 is the initial condition for V
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
lv_lna_solve_MALA=function(xi,B0,inter,deltat,param=c(0.5,0.0025,0.3))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = B0[1,1], y2C = B0[1,2], y22 = B0[2,2],
               S11_1 = 0, S12_1 = 0, S11_2 = 0, S12_2 = 0, S11_3 = 0, S12_3 = 0)
  odepars = c(c1 = param[1], c2 = param[2], c3 = param[3])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, lvFuncS, odepars)
  return(odeOut)
}


#Model function for eta and V - used for delayed acceptance without MALA
lvFunc <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dy11 <- (c1-c2*y12)*y11
    dy12 <- (c2*y11-c3)*y12
    
    dy21 <- 2*(c1-c2*y12)*y21 - 2*c2*y11*y2C + c1*y11 + c2*y11*y12
    dy2C <- (c1 + c2*(y11-y12) - c3)*y2C + c2*(y12*y21 - y11*y22 - y11*y12)
    dy22 <- 2*(c2*y11-c3)*y22 + 2*c2*y12*y2C + c3*y12 + c2*y11*y12
    
    return(list(c(dy11, dy12, dy21, dy2C, dy22)))
  })
}

#LNA solver for eta and V - used for delayed acceptance without MALA
#
#Arguments:
#xi is the initial condition for eta
#B0 is the initial condition for V
#inter is the time to be integrated to
#deltat is the size of the intermediate time steps
#param is a vector of model parameters
lv_lna_solve=function(xi,B0,inter,deltat,param=c(0.5,0.0025,0.3))
{
  odeinit <- c(y11 = xi[1], y12 = xi[2], y21 = B0[1,1], y2C = B0[1,2], y22 = B0[2,2])
  odepars = c(c1 = param[1], c2 = param[2], c3 = param[3])
  times=seq(0,inter,by=deltat)
  odeOut = ode(odeinit, times, lvFunc, odepars)
  return(odeOut)
}

#Forward filter for the LNA - stores gradient and likelihood under LNA, and an array of ODE solutions
#
#Arguments:
#data is all observations
#inter is the total time between two observations
#deltat is the intermediate time step between two observations
#param is a vector of model parameters
#m and sd are the prior mean of the (log) params, used for gradient calculation
#a is the initial condition for eta
#B is the inital condition for V
#sige is the observation noise
#F is the observation matrix (set for complete observation of both species)
#MALA is a boolean indicating whether the gradient should be calculated
lv_lna_MEGASOL_filter=function(data,inter,deltat,param=c(0.5,0.0025,0.3),m = 0, sd = 10, a=c(100,100),B=diag(c(0.0,0.0)), sige=1.0, F=diag(c(1,1)), MALA=TRUE)
{
  n=dim(data)[1] #number of observations
  n2=inter/deltat+1 #number of timepoints between observations (including endpoints)

  if(MALA==TRUE){
    gradthet= -(log(param)-m)/sd^2 #prior
  }
  else{
    gradthet=0
  }
  llike = 0
  #update marginal likelihood
  llike = llike + lgpdf(data[1,], t(F)%*%a, t(F)%*%B%*%F + diag(rep(sige*sige,2)))
  #update posterior
  a = a + B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%(data[1,]-t(F)%*%a)
  B = B - B%*%F%*%solve(t(F)%*%B%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%B
  #loop
  for(i in 2:n)
  {
    #update mu, sigma, and G (and sensitivities if MALA)
    if(MALA==TRUE){
      sensOut=lv_lna_solve_MALA(a,B,inter,deltat,param)
    }
    else{
      sensOut=lv_lna_solve(a,B,inter,deltat,param)
    }
    #solarray[,,(i-1)]=sensOut[,c("y11","y12","y21","y2C","y22","G1","G2","G3","G4")]
    mu=sensOut[n2,c("y11","y12")]
    sigma=matrix(sensOut[n2,c("y21", "y2C", "y2C", "y22")], ncol = 2)
    sigmaInv=solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2))) #Only used in gradient contribution so partially observed

    #update marginal likelihood
    llike = llike + lgpdf(data[i,], t(F)%*%mu, t(F)%*%sigma%*%F + diag(rep(sige*sige,2)))

    if(MALA==TRUE){
      #add gradient contribution
      geta1=c(sensOut[n2,"S11_1"],sensOut[n2,"S12_1"])
      geta2=c(sensOut[n2,"S11_2"],sensOut[n2,"S12_2"])
      geta3=c(sensOut[n2,"S11_3"],sensOut[n2,"S12_3"])
      C=sigmaInv%*%(data[i,]-t(F)%*%mu)
      gradthet=gradthet+c(t(C)%*%geta1,t(C)%*%geta2,t(C)%*%geta3)
      #print(gradthet)
    }

    #update posterior
    a = mu + sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%(data[i,]-t(F)%*%mu)
    B = sigma - sigma%*%F%*%solve(t(F)%*%sigma%*%F+diag(rep(sige*sige,2)))%*%t(F)%*%sigma
  }
  grad.out = c(gradthet[1]*param[1],gradthet[2]*param[2], gradthet[3]*param[3])
  return(list(grad = grad.out, llike = llike))
}
