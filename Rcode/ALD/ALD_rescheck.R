rm(list = ls())
# tryCatch( # set wording directory to current source file folder
#   {setwd(getSrcDirectory()[1])}
#   ,error=function(e){setwd(dirname(rstudioapi::getActiveDocumentContext()$path))}
# )
# setwd('C:/Users/admin/Documents/dongyoung/Rcode')
# setwd('D:/dy_result/QuantileRegression/Rcode')
setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_w_ME.R')
source('fn_wo_ME.R')
# library(mvtnorm)
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)

# Function --------------------------------------------------------------------------------
gen_error=function(type){
  if (type=='N'){
    ei=rnorm(n,0,1)
  }
  if (type=='t'){
    ei=rt(n,df = 3)
  }
  if (type=='locshit'){
    ei=rnorm(n,0,1)*(1+x1i)
  }  
  if (type=='MixN'){
    if(p0==0.1){to_subtract=1.43}
    if(p0==0.25){to_subtract=0.72}
    if(p0==0.5){to_subtract=-0.02}
    if(p0==0.75){to_subtract=-0.77}
    if(p0==0.9){to_subtract=-1.555}
    
    g_indx=rbinom(n = n,size = 1,prob = c(0.1))
    n1=sum(g_indx==0)
    n2=sum(g_indx==1)
    samp1=rnorm(n=n1,mean=to_subtract,sd=1)
    samp2=rnorm(n=n2,mean=(to_subtract+1),sd=5)
    samp=c(samp1,samp2)
    ei=samp
  }
  if (type=='BayesQR'){
    ei=rnorm(n=n, mean=0, sd=locshit_c*x1i)
  }
  return(ei)
}

# Calculate True beta based on input types & p0
gen_true.beta=function(type){
  if (type=='N'){
    beta.true=beta+c(qnorm(p0,0,1),0)
  }
  if (type=='t'){
    beta.true=beta+c(qt(p0,df = 3),0)
  }
  if (type=='locshit'){
    beta.true=c(beta[1]+qnorm(p0,0,1),beta[2]+qnorm(p0,0,1))
  }  
  if (type=='MixN'){
    beta.true=beta
  }
  if(type=='BayesQR'){
    beta.true=c(beta[1]+0*qnorm(p0,0,1),beta[2]+locshit_c*qnorm(p0,0,1))
  }
  return(beta.true)
}

# Define true parameter--------------------------------------------------------------------------------
type='BayesQR'
p0=0.25
n=1000
Mu_x=10
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1
beta=c(3,5)
alpha=c(4,5)
if(type=='BayesQR'){
  locshit_c=0.6
}




# Simulation start --------------------------------------------------------------------------------
sim_idx=266
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
sigma_save=rep(NA,nmax)

for(sim_idx in 1:nmax){
  load(file=sprintf('../debugging/ALD_%s_%s_%s_W1.RData',type,p0,sim_idx))
  X.est=colMeans(ALD_res$X_trace)
  beta.est=colMeans(ALD_res$beta_trace)
  alpha.est=colMeans(ALD_res$alpha_trace)
  X.est=colMeans(ALD_res$X_trace)
  mux.est=mean(ALD_res$mux_trace)
  sigma.est=mean(ALD_res$sigma_trace)
  sigma2_11.est=mean(ALD_res$sigma2_11_trace)
  sigma2_22.est=mean(ALD_res$sigma2_22_trace)
  sigma2_xx.est=mean(ALD_res$sigma2_xx_trace)
  
  alpha_save[sim_idx,]=alpha.est
  beta_save[sim_idx,]=beta.est
  X_save[sim_idx,]=X.est
  mux_save[sim_idx]=mux.est
  sigma_save[sim_idx]=sigma.est
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

# Make larger data for Ground Truth--------------------------------------------------------------------------------
set.seed(sim_idx)

# x1i=runif(n=n,min=0,max=2*Mu_x)
n=1e4
x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
X=cbind(1,x1i)

#generate W1,W2
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

W1=X%*%alpha+delta1
W2=X[,2]+delta2

ei=gen_error(type = type)
y=(X%*%beta)+ei
beta.true=gen_true.beta(type)

# Check result --------------------------------------------------------------------------------
par(mfrow=c(3,2))
hist(beta_save[,1],nclass=100);abline(v=beta.true[1],col=2,lwd=3)
hist(beta_save[,2],nclass=100);abline(v=beta.true[2],col=2,lwd=3)

hist(alpha_save[,1],nclass=100);abline(v=alpha[1],col=2,lwd=3)
hist(alpha_save[,2],nclass=100);abline(v=alpha[2],col=2,lwd=3)

plot(X[,2],y,main='X vs Y, with beta.est');abline(colMedians(beta_save,na.rm=T),col=2,lwd=3)
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMedians(alpha_save,na.rm=T),col=2,lwd=3)


hist(sigma_save,nclass=100);abline(v=sqrt(sigma2),col=2,lwd=3)
hist(mux_save,nclass=100);abline(v=Mu_x,col=2,lwd=3)
hist(sigma2_11_save,nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(sigma2_22_save,nclass=100);abline(v=sigma2_22,col=2,lwd=3)
hist(sigma2_xx_save,nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
par(mfrow=c(1,1))

# Debugging weired result --------------------------------------------------------------------------------
sim_idx=1
condition = which(abs(beta_save[,2])>10)
sim_idx=condition[3] # 26 & 29 are weired

load(file=sprintf('../debugging/ALD_%s_%s_%s_W1.RData',type,p0,sim_idx))
X.est=colMeans(ALD_res$X_trace)
beta.est=colMeans(ALD_res$beta_trace)
alpha.est=colMeans(ALD_res$alpha_trace)
X.est=colMeans(ALD_res$X_trace)
mux.est=mean(ALD_res$mux_trace)
sigma.est=mean(ALD_res$sigma_trace)
sigma2_11.est=mean(ALD_res$sigma2_11_trace)
sigma2_22.est=mean(ALD_res$sigma2_22_trace)
sigma2_xx.est=mean(ALD_res$sigma2_xx_trace)

n=1000
set.seed(sim_idx)
x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
X=cbind(1,x1i)
#generate W1,W2
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))
W1=X%*%alpha+delta1
W2=X[,2]+delta2
ei=gen_error(type = type)
y=(X%*%beta)+ei
beta.true=gen_true.beta(type = type)

par(mfrow=c(4,2))
plot(X.est,X[,2],main='X.est Vs X with y=x line');abline(0,1)
plot(X.est,W1,main='X.est Vs W1 with alpha.est');abline(alpha.est)
plot(X.est,W2,main='X.est Vs W2 with y=x line');abline(0,1)
plot(X.est,y,main='X.est Vs Y with beta.est');abline(beta.est)
hist(ALD_res$sigma_trace,nclass=100)
hist(ALD_res$sigma2_11_trace,nclass=100)
hist(ALD_res$sigma2_22_trace,nclass=100)
hist(ALD_res$sigma2_xx_trace,nclass=100)
par(mfrow=c(1,1))


par(mfrow=c(4,2))
ts.plot(ALD_res$beta_trace[,1])
ts.plot(ALD_res$alpha_trace[,1])
ts.plot(ALD_res$X_trace[,1])
ts.plot(ALD_res$sigma_trace)
ts.plot(ALD_res$sigma2_11_trace)
ts.plot(ALD_res$sigma2_22_trace)
ts.plot(ALD_res$sigma2_xx_trace)
par(mfrow=c(1,1))


sigma.est
sigma2_11.est
sigma2_22.est
sigma2_xx.est
# 
# mean(mux_save,na.rm = T)
# mean(sigma2_save,na.rm = T)
# hist(sigma2_save,nclass=100)
# median(sigma2_11_save,na.rm = T)
# hist(sigma2_11_save,nclass=100)
# mean(sigma2_22_save,na.rm = T)
# mean(sigma2_xx_save,na.rm = T)
