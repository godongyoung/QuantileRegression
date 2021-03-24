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
sim_idx=1
nmax=500
is.plot=F

for(sim_idx in 1:nmax){
  set.seed(sim_idx)
  
  
  # Make data --------------------------------------------------------------------------------
  # x1i=runif(n=n,min=0,max=2*Mu_x)
  x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
  X=cbind(1,x1i)
  
  #generate W1,W2
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  # Generate error based on input types & p0 ---------------------------------------
  ei=gen_error(type = type)
  y=(X%*%beta)+ei
  
  beta.true=gen_true.beta(type = type)
  if(is.plot){
    plot(X[,2],y,main='Groud Truth')
    abline(beta.true)
  }
  
  ALD_res=ALD_w_MME(y,W1,W2,p0)
  
  X.est=colMeans(ALD_res$X_trace)
  beta.est=colMeans(ALD_res$beta_trace)
  alpha.est=colMeans(ALD_res$alpha_trace)
  
  if(is.plot){
    plot(X.est,X[,2]);abline(0,1)
    plot(X.est,W1);abline(alpha.est)
    plot(X.est,y);abline(beta.est)
  }
  save(ALD_res, file=sprintf('../debugging/ALD_%s_%s_%s.RData',type,p0,sim_idx))
}

