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
# Define True parameter--------------------------------------------------------------------------------
n=1000
alpha=c(4,3)

Mu_x=5
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
p0=0.25
sim_idx=1
nmax=500
for(sim_idx in 1:nmax){
  set.seed(sim_idx)
  # Make data--------------------------------------------------------------------------------
  
  # x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sigma2_xx)
  x1i=runif(n=n,min=0,max=2*Mu_x)
  x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
  X=cbind(1,x1i)
  X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
  y=2+sin(x1i)+rnorm(n,0,0.1)
  
  #generate W1,W2
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  NQR_res=NQR_w_MME(y,W1,W2,p0,inp.min = 0,inp.max = 2*Mu_x,inp.version = 1,multiply_c = 3)
  save(NQR_res, file=sprintf('../debugging/NQR_old_%s_%s.RData',p0,sim_idx)) #old is W1 version!
  
}
# 
# 
# tau.i=NQR_res$Knots
# g.est=colMeans(NQR_res$g_trace,na.rm = T)
# plot(X[,2],y)
# points(tau.i,g.est,type='l')
# 
# #Result check-----------------------------------------------------------------------------------------------
# 
# idx=2
# ts.plot(X_trace[,idx]);abline(h=X[idx,2])
# 
# g.est=colMeans(g_trace,na.rm = T)
# plot(X[,2],y)
# points(tau.i,g.est,type='l')
# 
# X.est=colMeans(X_trace,na.rm = T)
# alpha.est=colMeans(alpha_trace,na.rm = T)
# plot(X.est,y)
# points(tau.i,g.est,type='l')
# plot(X.est,X[,2]);abline(0,1);abline(v=tau.i)
# plot(X.est,W1);abline(alpha.est)
# plot(X.est,W2);abline(0,1)
# 
# mean(sigma2_11_trace,na.rm = T)
# mean(sigma2_22_trace,na.rm = T)
# mean(sigma2_xx_trace,na.rm = T)

