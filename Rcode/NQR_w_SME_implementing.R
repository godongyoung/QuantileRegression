rm(list = ls())

setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')
# library(mvtnorm)
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)
library(Brq)
library(bayesQR)
library(quantreg)

########
set.seed(20210317)

# Function--------------------------------------------------------------------------------
smooth.y=function(knots,g.tau,xout,version=1){
  if(version==1){
    mspline=spline(x = knots,y = g.tau,xout = xout)
    y.est=mspline$y
  }
  if(version==2){
    msmooth.spline=smooth.spline(x = knots,y = g.tau,cv = NA,lambda = lambda.t)
    mspline=predict(msmooth.spline,xout)
    y.est=mspline$y
  }
  if(version==3){
    g.tau=as.numeric(g.tau)
    fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
    y.est=predict(fit.ns, data.frame(knots=xout))
  }
  return(y.est)
}

# Make data--------------------------------------------------------------------------------
p0=0.25
n=1000
inp.version=3
multiply_c=1
N=30

Mu_x=5
sigma2_xx=1
sigma2_22=1

# x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sigma2_xx)
x1i=runif(n=n,min=0,max=2*Mu_x)
X=cbind(1,x1i)
y=2+sin(x1i)+rnorm(n,0,0.1)
X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)

#generate w1,w2
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

W2=X[,2]+delta2


# Fitting --------------------------------------------------------------------------------
plot(X[,2],y)
for(p0 in c(0.1,0.25,0.5,0.75,0.9)){
  NQR_SME_res=NQR_w_SME(y,W2,p0)
  
  g.est=colMeans(NQR_SME_res$g_trace,na.rm = T)
  y.est=smooth.y(NQR_SME_res$Knots,g.est,X_range)
  points(X_range,y.est,type = 'l')    
}
naive_NQR_res=NQR(y,cbind(1,W2),p0)
g.est=colMeans(naive_NQR_res$g_trace,na.rm = T)
y.est=smooth.y(naive_NQR_res$Knots,g.est,X_range)
points(X_range,y.est,type = 'l')    

NQR_SME_res=NQR_w_SME(y,W2,0.25,is.plot = T)
plot(colMeans(NQR_SME_res$X_trace),y)
plot(W2,y)

