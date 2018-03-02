#### provides functions to compute premalignant-to-cancer sojourn time distribution

## for testing require 3-stage model parameters (e.g. Meza et al 2010) 
source('~/parameters_males_right_colon.R')

require(statmod)

## for fast Gauss-Legendre integration
gauss = gauss.quad(16)
# normalize legauss to [0.1]
gauss$nodes = 0.5*(gauss$nodes+1); gauss$weights = 0.5*gauss$weights

## 3-stage model hazard function 
h3  = function(a,mu0,mu1,X,alpha,p,q) {
  mu0X = mu0*X
  return(mu0X*(1-((q-p)/(q*exp(-p*a)-p*exp(-q*a)))^(mu1/alpha)))  
}

## 3-stage model survival function (simple integration)
S3 = function(a,mu0,mu1,X,alpha,p,q) {
  mu0X = mu0*X
  t = seq(0,a,.1)
  dum = sum(h3(t,mu0,mu1,X,alpha,p,q))*.1
  return(exp(-dum))  
}

## unconditional survival function for premaligant clones onset f_A(u)
S.adenoma = function(t,mu0,mu1,X) {
  mu0X = mu0*X
  return(exp(mu0X*((1-exp(-mu1*t))/mu1 - t)))
}

## unconditional density function for premalignant clones
f.adenoma = function(t,mu0,mu1,X) {
  mu0X = mu0*X
  return(mu0X*(1-exp(-mu1*t))*S.adenoma(t,mu0,mu1,X))
}

## 1-stage survival function for premalignant clone (p-clone)
S1 = function(t,alpha,p,q) {
  return(1+p*q*(exp(-p*t)-exp(-q*t))/(q*exp(-p*t)-p*exp(-q*t))/alpha)
}

## corresponding 1-stage density function
f1 = function(t,alpha,p,q) {
  w = q*exp(-p*t)-p*exp(-q*t)
  return(q*((-p*exp(-p*t)+q*exp(-q*t))/w + p*q*((exp(-p*t)-exp(-q*t))/w)^2))
}

## Jeon 3.10
phiY0 = function(y,t,a,alphaI,pI,qI) {
  au1 = a-t
  w=alphaI*(y-1)
  return(1+(pI*(qI-w)*exp(-pI*au1)-qI*(pI-w)*exp(-qI*au1))/
           ((qI-w)*exp(-pI*au1)-(pI-w)*exp(-qI*au1))/alphaI)
}

## 1-stage density function for non-extinction of p-clone and no cancer
p0 = function(x,alpha,p,q) {
  xi = (exp(-p*x)-exp(-q*x))/((q+alpha)*exp(-p*x)-(p+alpha)*exp(-q*x))
  out = xi*(alpha+p)*(alpha+q)*(q*exp(-p*x)-p*exp(-q*x))/(q*(alpha+p)*exp(-p*x)-p*(alpha+q)*exp(-q*x))
  return(out)
  # return(gamI*(exp(-p*x)-1)/(exp(-p*x) - gamI))
}

## construction of p-clone survival function (surviving a malignant transformation)
S.ad0 = function(a,mu0, mu1,X,alphaI,pI,qI) {
  U0 = a*gauss$nodes
  dum1 = 0; j=1
  for(u0 in U0) {
    x = (a-u0)*gauss$nodes
    dum2 = -mu1*(a-u0)*sum((1-p0(x,alphaI,pI,qI))*gauss$weights)
#    dum2 = -mu1*pinfI*(a-u0)
    dum1 = dum1+(1-exp(dum2))*gauss$weights[j];j=j+1
  }
  return(exp(-mu0*X*(a*dum1)))
}

## example
a =200; t = matrix(seq(a,0,-2),ncol=1)
y = apply(t,1,function(x) {S.ad0(a-x,mu0, mu1,X,alphaI,pI,qI)})

plot(a-t,y,type='l')
lines(a-t,y,col=2)

## density function by numerical derivative  
f.ad0 = function(a,mu0, mu1,X,alphaI,pI,qI) {
  if(a > .5) {
    a2 = a+.5; a1 = a-.5
    return(S.ad0(a1,mu0, mu1,X,alphaI,pI,qI)-S.ad0(a2,mu0, mu1,X,alphaI,pI,qI))
  } else {
    a2 = a+.5; a1 = 0
    return((1-S.ad0(a2,mu0, mu1,X,alphaI,pI,qI))*.5)
    
  }
}

## example for f.ad0 as define above
t  = matrix(seq(1,200,1),ncol=1)
aux = apply(t,1,function(x) {f.ad0(x,mu0, mu1,X,alphaI,pI,qI)}) 
plot(t,aux,type='l')
lines(t,aux,col=2,lwd=2)

## 1-stage density function conditioned on time < a (age malignancy occurs)
f1.cond = function(t,a,alpha,p,q) {
  # if(max(t)>a) {warning('max(t)>a'); exit}
  out = numeric(length(t)); cond = (t <= a); t = t[cond]
  w = q*exp(-p*t)-p*exp(-q*t)
  F1 = q*((-p*exp(-p*t)+q*exp(-q*t))/w + p*q*((exp(-p*t)-exp(-q*t))^2)/w/w)
  out[cond] = pinf*F1/(1-S1(a,alpha,p,q))
  return(out)
}

## normalization
gnorm = function(a,mu0,mu1,X,alpha,p,q) {
  x = matrix(a*gauss$nodes,ncol=1)
  dum = apply(x,1,function(x) {f1(a-x,alpha,p,q)*f.ad0(x,mu0,mu1,X,alpha,p,q)})
  return(a*sum(dum*gauss$weights))
}

## full convolution
g.adenoma = function(t,a,mu0,mu1,X,alpha,p,q) {
  if(max(t)>a) {warning('max(t)>a'); break}
  dum1 = f1(a-t,alpha,p,q)*f.ad0(t,mu0,mu1,X,alpha,p,q)
  return(dum1/gnorm(a,mu0,mu1,X,alpha,p,q))
}

## apply/vectorize convolution to given time vector u1 < a (=max(u1))
f.adenoma = function(u1) {
  U1 = matrix(u1,ncol=1)
  return(apply(U1,1,function(x) {g.adenoma(x,max(u1),mu0,mu1,X,alphaI,pI,qI)}))
}

## example: cancer at age 50 
u1 = seq(0,44,1); tmp = f.adenoma(u1)
dwellT = 50 - u1 
plot(dwellT,tmp,type='l')
lines(u1,tmp,col=4)

#### code below is for testing the impact of a stochastic clonal expansion of malignant cells 
#### convolve with malignant exansion. Tricky - requires high precision in integration
# fullConv = function(a) {
#   source('./Simulations/CRC_parameters.R')
#   U1 = seq(0,a,1)
#   dens = numeric(length(U1)); j=1
#   for(u1 in U1) {
#     u2  = matrix((a-u1)*gauss$nodes+u1,ncol=1)
#     dum = apply(u2,1,function(x) {f1(a-x,alphaM,pM,qM)*g.adenoma(u1,x,mu0,mu1,X,alphaI,pI,qI)})
#     dens[j] = (a-u1)*sum(dum*gauss$weights); dens[is.na(dens)]=0
#     j=j+1
#   }
#   return(dens)        
# }














