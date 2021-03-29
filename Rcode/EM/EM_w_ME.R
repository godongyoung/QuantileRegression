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

# Sim start ----------------------------------------------------------------------------------------



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

rst=matrix(NA,nmax,9) 
colnames(rst)=c('beta0','beta1','alpha0','alpha1','Mux','sigma2_xx','sigma2','sigma2_11','sigma2_22')

for(sim_idx in 1:nmax){
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
  
  EM_res=EM_w_MME(y,W1,W2)

  rst[sim_idx,]=EM_res
}

colMeans(rst)

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

# If scale --------------
X=scale(X,center = T,scale = F)
W1=scale(W1,center = T,scale = F)
W2=scale(W2,center = T,scale = F)
# If scale --------------

if(is.plot){
  plot(X[,2],y)
}

# Check result --------------------------------------------------------------------------------
par(mfrow=c(3,2))
hist(rst[,'beta0'],nclass=100);abline(v=beta[1],col=2,lwd=3)
hist(rst[,'beta1'],nclass=100);abline(v=beta[2],col=2,lwd=3)

hist(rst[,'alpha0'],nclass=100);abline(v=alpha[1],col=2,lwd=3)
hist(rst[,'alpha1'],nclass=100);abline(v=alpha[2],col=2,lwd=3)

plot(X[,2],y,main='X vs Y, with beta.est');abline(colMeans(rst[,c('beta0','beta1')]))
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMeans(rst[,c('alpha0','alpha1')]))

hist(rst[,'sigma2'],nclass=100);abline(v=sigma2,col=2,lwd=3)
hist(rst[,'Mux'],nclass=100);abline(v=Mu_x,col=2,lwd=3)
hist(rst[,'sigma2_11'],nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(rst[,'sigma2_xx'],nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
hist(rst[,'sigma2_22'],nclass=100);abline(v=sigma2_22,col=2,lwd=3)
par(mfrow=c(1,1))