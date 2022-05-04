### Helper functions used in both applications


#Generates a single innovation from a multivariate random normal distribution, given a mean vector m and covariance matrix V (used for parameter proposal)
rmvn=function(m,V)
{
  p=length(m)
  z=rnorm(p)
  return(m+t(chol(V))%*%z)
}

#Evaluates the (log) prior density given a set of log-parameters, param, prior mean a and prior standard deviation b
logprior=function(param,a=0,b=10)
{
 #assumes log-params follow a N(a,b^2) a priori
 return(sum(dnorm(param,a,b,log=T)))
}

#Evaluates the log of a multivariate Gaussian (mmean, mvar) pdf at x
lgpdf=function(x,mmean,mvar) 
{
 d=length(mmean)
 -0.5*log(det(mvar)) -0.5*t(x-mmean)%*%solve(mvar)%*%(x-mmean) -0.5*d*log(2*pi)
}

#Euclidean distance between two vectors (vec is the difference between two vectors)
euclid=function(vec)
{
	return(sqrt(sum(vec^2)))
}

#returns indices of Euclidean sorted (smallest to largest) vectors, where each vector is a row of xmat
esort=function(xmat)
{
 N=dim(xmat)[1]
 indices=rep(0,N)
 iset=1:N
 indices[1]=which.min(xmat[,1])
 for(j in 2:N)
 {
  xstar=xmat[indices[j-1],]
  iset=iset[iset!=indices[(j-1)]]
  dist=rep(0,length(iset))
  for(i in 1:length(iset))
  {
   dist[i]=euclid(xstar-xmat[iset[i],])
  }
  ind=which.min(dist)
  indices[j]=iset[ind]
 }
 indices
}

# Systematic weighted resampling:
# wts are the weights
# N is the number of times to resample
# uni is a U(0,1) random variate
sysresamp2=function(wts,N,uni)
{
 vec=rep(0,N)
 wsum=sum(wts)
 k=1
 u=uni/N
 wsumtarg=u
 wsumcurr=wts[k]/wsum
 delta=1/N
 for(i in 1:N)
 {
  while (wsumcurr<wsumtarg)
  {
   k=k+1
   wsumcurr=wsumcurr+wts[k]/wsum
  }   
  vec[i]=k 
  wsumtarg=wsumtarg+delta
 }
 vec
}

#Approximate effective sample size of weighted resamples (given weights)
ess=function(wts)
{
(sum(wts)^2)/sum(wts^2) #ess
}

