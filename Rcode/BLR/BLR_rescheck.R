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
sim_idx=500
nmax=500
is.plot=F

# Simulation check --------------------------------------------------------------------------------
tic()
beta_save = matrix(NA,ncol=2,nrow=nmax)
alpha_save = matrix(NA,ncol=2,nrow=nmax)
X_save=matrix(NA,ncol=n,nrow=nmax)
mux_save=rep(NA,nmax)
sigma2_11_save=rep(NA,nmax)
sigma2_22_save=rep(NA,nmax)
sigma2_xx_save=rep(NA,nmax)
sigma2_save=rep(NA,nmax)

for(sim_idx in 1:nmax){
  load(file=sprintf('../debugging/BLR_%s.RData',sim_idx))
  beta.est=colMeans(BLR_res$beta_trace)
  alpha.est=colMeans(BLR_res$alpha_trace)
  X.est=colMeans(BLR_res$X_trace)
  mux.est=mean(BLR_res$mux_trace)
  sigma2.est=mean(BLR_res$sigma2_trace)
  sigma2_11.est=mean(BLR_res$sigma2_11_trace)
  sigma2_22.est=mean(BLR_res$sigma2_22_trace)
  sigma2_xx.est=mean(BLR_res$sigma2_xx_trace)
  
  alpha_save[sim_idx,]=alpha.est
  beta_save[sim_idx,]=beta.est
  X_save[sim_idx,]=X.est
  mux_save[sim_idx]=mux.est
  sigma2_save[sim_idx]=sigma2.est
  sigma2_11_save[sim_idx]=sigma2_11.est
  sigma2_22_save[sim_idx]=sigma2_22.est
  sigma2_xx_save[sim_idx]=sigma2_xx.est
  if(is.plot){
    par(mfrow=c(2,2))
    plot(X.est,X[,2]);abline(0,1)
    plot(X.est,W1);abline(alpha.est)
    plot(X.est,y);abline(beta.est)
    par(mfrow=c(1,1))
  }
}
toc()

# Make data--------------------------------------------------------------------------------
set.seed(sim_idx)
n=1e4
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

# Check result --------------------------------------------------------------------------------
par(mfrow=c(3,2))
hist(beta_save[,1],nclass=100);abline(v=beta[1],col=2,lwd=3)
hist(beta_save[,2],nclass=100);abline(v=beta[2],col=2,lwd=3)

hist(alpha_save[,1],nclass=100);abline(v=alpha[1],col=2,lwd=3)
hist(alpha_save[,2],nclass=100);abline(v=alpha[2],col=2,lwd=3)

plot(X[,2],y,main='X vs Y, with beta.est');abline(colMedians(beta_save,na.rm=T))
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMedians(alpha_save,na.rm=T))

hist(sigma2_save,nclass=100);abline(v=sigma2,col=2,lwd=3)
hist(mux_save,nclass=100);abline(v=Mu_x,col=2,lwd=3)
hist(sigma2_11_save,nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(sigma2_22_save,nclass=100);abline(v=sigma2_22,col=2,lwd=3)
hist(sigma2_xx_save,nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
par(mfrow=c(1,1))

# Debugging weired result --------------------------------------------------------------------------------
condition=c( 1,   5,   7,  11,  17,  18,  19,  21,  24,  25,  26,  32,  33,  34,  35,  36,  37,  39,  40,  42,  44,  46,  47,  48,  50,  54,  61,  65,  66,  73,  74,  75,  76,  77,  78,  79,  81,  84,  86,  89,  90,  92,  94,  95,  96, 100)
condition=which(beta_save[,1]>10)
sim_idx = condition[2] # 5,7 are weired!
# sim_idx = 5

n=1000
# load(file=sprintf('../debugging/BLR_%s_v2.RData',sim_idx))
load(file=sprintf('../debugging/BLR_%s.RData',sim_idx))
beta.est=colMeans(BLR_res$beta_trace)
alpha.est=colMeans(BLR_res$alpha_trace)
X.est=colMeans(BLR_res$X_trace)
mux.est=mean(BLR_res$mux_trace)
sigma2.est=mean(BLR_res$sigma2_trace)
sigma2_11.est=mean(BLR_res$sigma2_11_trace)
sigma2_22.est=mean(BLR_res$sigma2_22_trace)
sigma2_xx.est=mean(BLR_res$sigma2_xx_trace)

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


par(mfrow=c(4,2))
plot(X.est,X[,2],main='X.est Vs X with y=x line');abline(0,1)
plot(X.est,W1,main='X.est Vs W1 with alpha.est');abline(alpha.est)
plot(X.est,W2,main='X.est Vs W2 with y=x line');abline(0,1)
plot(X.est,y,main='X.est Vs Y with beta.est');abline(beta.est)
hist(BLR_res$sigma2_trace,nclass=100)
hist(BLR_res$sigma2_11_trace,nclass=100)
hist(BLR_res$sigma2_22_trace,nclass=100)
hist(BLR_res$sigma2_xx_trace,nclass=100)
par(mfrow=c(1,1))

par(mfrow=c(4,2))
ts.plot(BLR_res$beta_trace[,1])
ts.plot(BLR_res$alpha_trace[,1])
ts.plot(BLR_res$X_trace[,1])
ts.plot(BLR_res$sigma2_trace)
ts.plot(BLR_res$sigma2_11_trace)
ts.plot(BLR_res$sigma2_22_trace)
ts.plot(BLR_res$sigma2_xx_trace)
par(mfrow=c(1,1))
sigma2.est
sigma2_11.est
sigma2_22.est
sigma2_xx.est
