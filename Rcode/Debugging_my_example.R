rm(list = ls())
setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)
library(Brq)
library(bayesQR)
library(SuppDists)
library(quantreg)

p0_list=c(0.1,0.5,0.9)
res_list=list()
for(p0 in p0_list){
  nmax=100
  beta_save1=matrix(NA,ncol=2,nrow=nmax)
  beta_save2=matrix(NA,ncol=2,nrow=nmax)
  beta_save3=matrix(NA,ncol=2,nrow=nmax)
  beta_save4=matrix(NA,ncol=2,nrow=nmax)
  
  sigma_save1=rep(NA,nmax)
  sigma_save2=rep(NA,nmax)
  sigma_save3=rep(NA,nmax)
  tic()
  for(ii in 1:nmax){
    n=100
    ui=rnorm(n,0,1)
    zi=rexp(n,1)
    A=(1-2*p0)/(p0*(1-p0))
    B=2/(p0*(1-p0))
    beta=c(4,1)
    x1i=seq(0,20,length.out = n)
    X=cbind(1,x1i)
    sigma=3
    y=X%*%beta+sigma*A*zi+sigma*sqrt(B)*sqrt(zi)*ui
    
    GAL_res=GAL_wo_ME(y,X,p0)
    beta_save1[ii,]=colMeans(GAL_res[['beta_trace']])
    sigma_save1[ii]=mean(GAL_res[['sigma_trace']])
    
    ALD_res=ALD_wo_ME(y,X,p0)
    beta_save2[ii,]=colMeans(ALD_res[['beta_trace']])
    sigma_save2[ii]=mean(ALD_res[['sigma_trace']])
    
    BQR_res=mBayesQR(y,X,p0,F)
    beta_save3[ii,]=colMeans(BQR_res[['beta_trace']])
    sigma_save3[ii]=mean(BQR_res[['sigma_trace']])
    
    res=rq(y~X[,-1],p0)
    beta_save4[ii,]=res$coefficients
    
  }
  toc()
  
  # colMeans(beta_save1)
  # colMeans(beta_save2)
  # colMeans(beta_save3)
  # 
  # colVars(beta_save1)
  # colVars(beta_save2)
  # colVars(beta_save3)
  # 
  # sigma
  # mean(sigma_save1)
  # mean(sigma_save2)
  # mean(sigma_save3)
  # 
  # var(sigma_save1)
  # var(sigma_save2)
  # var(sigma_save3)
  # 
  # hist(GAL_res[['sigma_trace']],nclass=100)
  # hist(ALD_res[['sigma_trace']],nclass=100)
  # ts.plot(ALD_res[['sigma_trace']])
  # ts.plot(GAL_res[['sigma_trace']])
  # ts.plot(GAL_res[['gamma_trace']])
  # 
  res_list[[as.character(p0)]]=list()
  true_beta=beta
  df=data.frame(matrix(NA,nrow=4,ncol=((dim(X)[2]*3+2))),row.names = c('GAL','ALD','BayesQR','QR'))
  df[1,]=c((colMeans(beta_save1)-true_beta)/1, NA, (colVars(beta_save1)), NA, (colMeans(beta_save1)-true_beta)^2+ colVars(beta_save1))
  df[2,]=c((colMeans(beta_save2)-true_beta)/1, NA, (colVars(beta_save2)), NA, (colMeans(beta_save2)-true_beta)^2+ colVars(beta_save2))
  df[3,]=c((colMeans(beta_save3)-true_beta)/1, NA, (colVars(beta_save3)), NA, (colMeans(beta_save3)-true_beta)^2+ colVars(beta_save3))
  df[4,]=c((colMeans(beta_save4)-true_beta)/1, NA, (colVars(beta_save4)), NA, (colMeans(beta_save4)-true_beta)^2+ colVars(beta_save4))
  colnames(df)=c(paste0('bias',1:dim(X)[2]),NA,paste0('var',1:dim(X)[2]),NA,paste0('mse',1:dim(X)[2]))
  res_list[[as.character(p0)]][['beta']]=df
  
  
  true_sigma=sigma
  df=data.frame(matrix(NA,nrow=4,ncol=((1*3+2))),row.names = c('GAL','ALD','BayesQR','QR'))
  df[1,]=c((mean(sigma_save1)-true_sigma)/1, NA, (var(sigma_save1)), NA, (mean(sigma_save1)-true_sigma)^2+ var(sigma_save1))
  df[2,]=c((mean(sigma_save2)-true_sigma)/1, NA, (var(sigma_save2)), NA, (mean(sigma_save2)-true_sigma)^2+ var(sigma_save2))
  df[3,]=c((mean(sigma_save3)-true_sigma)/1, NA, (var(sigma_save3)), NA, (mean(sigma_save3)-true_sigma)^2+ var(sigma_save3))
  colnames(df)=c(paste0('bias',1),NA,paste0('var',1),NA,paste0('mse',1))
  res_list[[as.character(p0)]][['sigma']]=df
}

# save.image(file='../debugging/Debugging_w_my_example_data.RData')