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

# for(sim_idx in 1:nmax){
for(sim_idx in c(1,  4,  7,  8, 10, 11, 14, 16) ){
# for(sim_idx in c( 1,   5,   7,  11,  17,  18,  19,  21,  24,  25,  26,  32,  33,  34,  35,  36,  37,  39,  40,  42,  44,  46,  47,  48,  50,  54,  61,  65,  66,  73,  74,  75,  76,  77,  78,  79,  81,  84,  86,  89,  90,  92,  94,  95,  96, 100)){
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
  
  # If scale --------------
  X=scale(X,center = T,scale = F)
  W1=scale(W1,center = T,scale = F)
  W2=scale(W2,center = T,scale = F)
  # If scale --------------
  
  if(is.plot){
    plot(X[,2],y)
  }
  
  BLR_res=BLR_w_MME(y,W1,W2,multiply_c = 20)
  
  X.est=colMeans(BLR_res$X_trace)
  beta.est=colMeans(BLR_res$beta_trace)
  alpha.est=colMeans(BLR_res$alpha_trace)
  if(is.plot){
    plot(X.est,X[,2]);abline(0,1)
    plot(X.est,W1);abline(alpha.est)
    plot(X.est,y);abline(beta.est)
  }
  save(BLR_res, file=sprintf('../debugging/BLR_%s_scale.RData',sim_idx))
}

