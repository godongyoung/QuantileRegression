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

# Function --------------------------------------------------------------------------------
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
N=30

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
is.cal_quantile=T
if.short=F


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
y.quantile_save = rep(NA,nmax)
y.est.quantile_save = rep(NA,nmax)
bias.quantile_save = rep(NA,nmax)

y.est.quantile_save2 = rep(NA,nmax)

load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
tau.i=NQR_res$Knots
tau.i=seq(from = 0,to = 2*Mu_x,length.out = 30)



for(sim_idx in 1:nmax){
  # Make data--------------------------------------------------------------------------------
  set.seed(sim_idx)
  x1i=runif(n=n,min=0,max=2*Mu_x)
  x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
  X=cbind(1,x1i)
  X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
  y=2+sin(x1i)+rnorm(n,0,0.1)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  # Load MCMC result--------------------------------------------------------------------------------
  load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
  if(!if.short){
    alpha.est=colMeans(NQR_res$alpha_trace)
    g.est=colMeans(NQR_res$g_trace)
    X.est=colMeans(NQR_res$X_trace)
    mux.est=mean(NQR_res$mux_trace)
    sigma2_11.est=mean(NQR_res$sigma2_11_trace)
    sigma2_22.est=mean(NQR_res$sigma2_22_trace)
    sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
    Knots = NQR_res$Knots
  }
  if(if.short){
    alpha.est = NQR_res_short$alpha.est
    g.est = NQR_res_short$g.est
    X.est = NQR_res_short$X.est
    mux.est = NQR_res_short$mux.est
    sigma2_11.est = NQR_res_short$sigma2_11.est
    sigma2_22.est = NQR_res_short$sigma2_22.est
    sigma2_xx.est = NQR_res_short$sigma2_xx.est
    Knots = NQR_res_short$Knots
  }
  
  alpha_save[sim_idx,]=alpha.est
  g_save[sim_idx,]=g.est
  X_save[sim_idx,]=X.est
  mux_save[sim_idx]=mux.est
  sigma2_11_save[sim_idx]=sigma2_11.est
  sigma2_22_save[sim_idx]=sigma2_22.est
  sigma2_xx_save[sim_idx]=sigma2_xx.est
  accept_g_save[sim_idx]=NQR_res$g_accept_ratio
  
  if(is.cal_quantile){
    ordered_idx=order(X[,2])
    X=X[ordered_idx,]
    y=y[ordered_idx]
    
    y.est=smooth.y(Knots, g.est, X[,2], version=3)
    y.est.quantile_save[sim_idx] = quantile(y.est,p0)
    y.est.quantile_save2[sim_idx] = mean(y.est) #quantile(y.est,0.5)
    y.quantile_save[sim_idx] = quantile(y,p0)
    bias.quantile_save[sim_idx] = quantile(y-y.est,p0)
    # plot(X[,2],y)
    # points(X[,2],y.est,type='l')    
  }

  if(is.plot){
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

# plot(X[,2],y)
# points(X[,2],y.est,type='l')
# boxplot(y,y.est);abline(h=mean(y.quantile_save[-nconverge_idx]))
# boxplot(y.est.quantile_save[-nconverge_idx]);abline(h=mean(y.quantile_save[-nconverge_idx]))
# boxplot(y.est.quantile_save2[-nconverge_idx]);abline(h=mean(y.quantile_save[-nconverge_idx]))
# hist(y.est.quantile_save[-nconverge_idx],nclass=100);abline(v=mean(y.quantile_save[-nconverge_idx]))
# boxplot(bias.quantile_save[-nconverge_idx]);abline(h=0)
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
hist(accept_g_save[-nconverge_idx])


par(mfrow=c(4,2))
plot(X[,2],y,main='X vs Y, with g.est');points(tau.i,colMedians(g_save[-nconverge_idx,],na.rm = T),type='l',col=2,lwd=3)
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMedians(alpha_save[-nconverge_idx,],na.rm=T),col=2,lwd=3)

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
sim_idx=254
sim_idx=3
condition = which(abs(alpha_save[,2])>10)
condition = which(abs(alpha_save[,2])<10)

for(idx in 1:length(condition)){
  sim_idx=condition[idx]
  # sim_idx=idx
  if(sim_idx %in% nconverge_idx){next}
  
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
  
  # par(mfrow=c(4,2))
  # plot(X.est,X[,2],main='X.est Vs X with y=x line');abline(0,1)
  # plot(X.est,W1,main='X.est Vs W1 with alpha.est');abline(alpha.est)
  # plot(X.est,W2,main='X.est Vs W2 with y=x line');abline(0,1)
  # plot(X.est,y,main='X.est Vs Y with beta.est');points(tau.i,g.est,type='l',col=2,lwd=3)
  # 
  # hist(NQR_res$sigma2_11_trace,nclass=100);abline(v=sigma2_11,col=2,lwd=3)
  # hist(NQR_res$sigma2_22_trace,nclass=100);abline(v=sigma2_22,col=2,lwd=3)
  # hist(NQR_res$sigma2_xx_trace,nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
  # 
  # ts.plot(NQR_res$sigma2_22_trace)
  # par(mfrow=c(1,1))
  
  NQR_res$g_accept_ratio
  NQR_res$l_accept_ratio
  NQR_res$x_accept_ratio
  sigma2_11.est
  sigma2_22.est
  sigma2_xx.est
  
  par(mfrow=c(3,2))
  ts.plot(NQR_res$g_trace[,1],main=sim_idx)
  ts.plot(NQR_res$alpha_trace[,1])
  ts.plot(NQR_res$X_trace[,1])
  ts.plot(NQR_res$sigma2_11_trace)
  ts.plot(NQR_res$sigma2_22_trace)
  ts.plot(NQR_res$sigma2_xx_trace)
  par(mfrow=c(1,1))
  
  accept_g_save[condition]
  
  median(accept_g_save,na.rm = T)

}
# sim_idx=condition[1]

