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
alpha=c(4,3)

Mu_x=5
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
sim_idx=4
p0=0.25
nmax=500
is.plot=F




# Simulation check --------------------------------------------------------------------------------
tic()

alpha_save = matrix(NA,ncol=2,nrow=nmax)
g_save=matrix(NA,ncol=30,nrow=nmax)
X_save=matrix(NA,ncol=n,nrow=nmax)
mux_save=rep(NA,nmax)
sigma2_11_save=rep(NA,nmax)
sigma2_22_save=rep(NA,nmax)
sigma2_xx_save=rep(NA,nmax)
accept_g_save = rep(NA,nmax)

load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
tau.i=NQR_res$Knots
tau.i=seq(from = 0,to = 2*Mu_x,length.out = 30)



for(sim_idx in 1:nmax){
  load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
  alpha.est=colMeans(NQR_res$alpha_trace)
  g.est=colMeans(NQR_res$g_trace)
  X.est=colMeans(NQR_res$X_trace)
  mux.est=mean(NQR_res$mux_trace)
  sigma2_11.est=mean(NQR_res$sigma2_11_trace)
  sigma2_22.est=mean(NQR_res$sigma2_22_trace)
  sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
  
  alpha_save[sim_idx,]=alpha.est
  g_save[sim_idx,]=g.est
  X_save[sim_idx,]=X.est
  mux_save[sim_idx]=mux.est
  sigma2_11_save[sim_idx]=sigma2_11.est
  sigma2_22_save[sim_idx]=sigma2_22.est
  sigma2_xx_save[sim_idx]=sigma2_xx.est
  accept_g_save[sim_idx]=NQR_res$g_accept_ratio
  
  if(is.plot){
    # Make data--------------------------------------------------------------------------------
    set.seed(sim_idx)
    
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
    
    par(mfrow=c(2,2))
    plot(X.est,X[,2]);abline(0,1)
    plot(X.est,W1);abline(alpha.est)
    plot(X.est,y);points(tau.i,g.est,type='l')
    par(mfrow=c(1,1))
    
    
    # mean(NQR_res$sigma2_11_trace)
    # mean(NQR_res$sigma2_22_trace)
    # mean(NQR_res$sigma2_xx_trace)
    # 
    # ts.plot(NQR_res$sigma2_11_trace)
    # ts.plot(NQR_res$sigma2_22_trace)
    # ts.plot(NQR_res$sigma2_xx_trace)
    # 
  }
}
toc()


# Make larger data for Ground Truth--------------------------------------------------------------------------------
set.seed(sim_idx)
n=1e4
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


# Check result --------------------------------------------------------------------------------
nconverge_idx=which(accept_g_save<0.1)


par(mfrow=c(4,2))


plot(X[,2],y,main='X vs Y, with g.est');points(tau.i,colMedians(g_save,na.rm = T),type='l',col=2,lwd=3)
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMedians(alpha_save,na.rm=T),col=2,lwd=3)

hist(alpha_save[-nconverge_idx,1],nclass=100);abline(v=alpha[1],col=2,lwd=3)
hist(alpha_save[-nconverge_idx,2],nclass=100);abline(v=alpha[2],col=2,lwd=3)

hist(mux_save[-nconverge_idx],nclass=100);abline(v=Mu_x,col=2,lwd=3)
hist(sigma2_11_save[-nconverge_idx],nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(sigma2_22_save[-nconverge_idx],nclass=100);abline(v=sigma2_22,col=2,lwd=3)
hist(sigma2_xx_save[-nconverge_idx],nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
par(mfrow=c(1,1))


# mean(mux_save,na.rm = T)
# median(sigma2_11_save,na.rm = T)
# hist(sigma2_11_save,nclass=100)
# mean(sigma2_22_save,na.rm = T)
# mean(sigma2_xx_save,na.rm = T)

# Debugging weired result --------------------------------------------------------------------------------
{sim_idx=11
sim_idx=which(sigma2_xx_save>4)[2]

load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
alpha.est=colMeans(NQR_res$alpha_trace)
g.est=colMeans(NQR_res$g_trace)
X.est=colMeans(NQR_res$X_trace)
mux.est=mean(NQR_res$mux_trace)
sigma2_11.est=mean(NQR_res$sigma2_11_trace)
sigma2_22.est=mean(NQR_res$sigma2_22_trace)
sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)


set.seed(sim_idx)
n=1000
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

par(mfrow=c(4,2))
plot(X.est,X[,2],main='X.est Vs X with y=x line');abline(0,1)
plot(X.est,W1,main='X.est Vs W1 with alpha.est');abline(alpha.est)
plot(X.est,W2,main='X.est Vs W2 with y=x line');abline(0,1)
plot(X.est,y,main='X.est Vs Y with beta.est');points(tau.i,g.est,type='l',col=2,lwd=3)

hist(NQR_res$sigma2_11_trace,nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(NQR_res$sigma2_22_trace,nclass=100);abline(v=sigma2_22,col=2,lwd=3)
hist(NQR_res$sigma2_xx_trace,nclass=100);abline(v=sigma2_xx,col=2,lwd=3)

ts.plot(NQR_res$sigma2_22_trace)
par(mfrow=c(1,1))
}
NQR_res$g_accept_ratio
NQR_res$l_accept_ratio
NQR_res$x_accept_ratio
sigma2_11.est
sigma2_22.est
sigma2_xx.est

par(mfrow=c(3,2))
ts.plot(NQR_res$g_trace[,1])
ts.plot(NQR_res$alpha_trace[,1])
ts.plot(NQR_res$X_trace[,1])
ts.plot(NQR_res$sigma2_11_trace)
ts.plot(NQR_res$sigma2_22_trace)
ts.plot(NQR_res$sigma2_xx_trace)
par(mfrow=c(1,1))

accept_g_save[which(sigma2_xx_save>4)]

median(accept_g_save)
