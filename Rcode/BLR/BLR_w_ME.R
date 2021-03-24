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
library(MASS)



# Define True parameter--------------------------------------------------------------------------------
n=1000
multiply_c=2
beta=c(3,5)
alpha=c(4,3)

Mu_x=5
sigma2=1
sigma2_xx=1
sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
sim_idx=1
nmax=500
is.plot=F

for(sim_idx in 100:nmax){
  # Make data--------------------------------------------------------------------------------
  set.seed(sim_idx)
  x1i=runif(n=n,min=0,max=2*Mu_x)
  x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
  X=cbind(1,x1i)
  y=X%*%beta+rnorm(n,0,1)
  X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
  
  #generate W1,W2
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  if(is.plot){
    plot(X[,2],y)
  }
  
  BLR_res=BLR_w_MME(y,W1,W2)
  
  X.est=colMeans(BLR_res$X_trace)
  beta.est=colMeans(BLR_res$beta_trace)
  alpha.est=colMeans(BLR_res$alpha_trace)
  if(is.plot){
    plot(X.est,X[,2]);abline(0,1)
    plot(X.est,W1);abline(alpha.est)
    plot(X.est,y);abline(beta.est)
  }
  save(BLR_res, file=sprintf('../debugging/BLR_%s_v2.RData',sim_idx))
}
