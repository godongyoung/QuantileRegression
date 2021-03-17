# The result says that it is not reliable... need fixing but stopped implenenting 0316
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




# Set Defaults --------------------------------------------------------------------------------
nmax=1000
p0_list=c(0.1,0.5,0.9)
sim_idx=1

GAL_list=list()
naive_list=list()
for(p0 in p0_list){
  GAL_list[[as.character(p0)]]=matrix(NA,ncol=2,nrow=nmax)
  naive_list[[as.character(p0)]]=matrix(NA,ncol=2,nrow=nmax)
}

# Iteration starts --------------------------------------------------------------------------------
par(mfrow=c(1,3))
for(sim_idx in 1:nmax){
  set.seed(sim_idx)
  # Make Data ------------------------------------------------
  n=500
  
  X=rnorm(n,4,1)
  sigma_delta=1/4
  ui=rnorm(n,0,sigma_delta)
  wi=X+ui
  W1=cbind(1,wi)
  beta=c(0,2)
  y=cbind(1,X)%*%beta + (cbind(1,X)%*%c(0,0.5))*rnorm(n,0,1)
  
  alpha=c(0,1)
  
  rm(X)
  
  for(p0 in p0_list){
    
    true.beta = beta + c(0,0.5*qnorm(p0))
    res_list=GAL_w_SME(y,W1,p0)

    if(sim_idx<101){
      save(res_list,file=sprintf('../debugging/ME_single_nalpha_%s_%s.RData',p0,sim_idx))  
    }
    
    
    GAL_list[[as.character(p0)]][sim_idx,] = colMeans(res_list[['beta_trace']]) - true.beta
    rq_res=rq(y~wi,p0)
    naive_list[[as.character(p0)]][sim_idx, ] = rq_res$coefficients - true.beta
    
    if(sim_idx%%5==1){
      p0=as.character(p0)
      data=cbind(GAL_list[[p0]],naive_list[[p0]])
      colnames(data)=c('GAL_b0','GAL_b1','naive_b0','naive_b1')
      boxplot(data,main=sprintf('%s result at %s simulations',p0,sim_idx));abline(h=0)
      }
  }
}
par(mfrow=c(1,1))

# GAL_list=list()
# naive_list=list()
# for(p0 in p0_list){
#   true.beta = beta + c(0,0.5*qnorm(p0))
#   beta_save=matrix(NA,ncol=2,nrow=nmax)
#   naive_beta_save=matrix(NA,ncol=2,nrow=nmax)
#   for(sim_idx in 1:nmax){
#     load(file=sprintf('../debugging/ME_single_nalpha_%s_%s.RData',p0,sim_idx))
#     beta_save[sim_idx,]=colMeans(beta_trace) - true.beta
#     
#     res=rq(y~wi,p0)
#     naive_beta_save[sim_idx,]=res$coefficients - true.beta
#   }
#   GAL_list[[as.character(p0)]]=beta_save
#   naive_list[[as.character(p0)]]=naive_beta_save
# }
# 
# p0='0.9'
# data=cbind(GAL_list[[p0]],naive_list[[p0]])
# colnames(data)=c('GAL_b0','GAL_b1','naive_b0','naive_b1')
# boxplot(data,main=p0);abline(h=0)
# 
# # set.seed(sim_idx)
# # n=200
# # X=rnorm(n,4,1)
# # sigma_delta=1/4
# # ui=rnorm(n,0,sigma_delta)
# # wi=X+ui
# # W1=cbind(1,wi)
# # beta=c(0,2)
# # y=cbind(1,X)%*%beta + (cbind(1,X)%*%c(0,0.5))*rnorm(n,0,1)
# # plot(X,y)
# # for(p0 in p0_list){
# #   load(file=sprintf('../debugging/ME_single_nalpha_%s_%s.RData',p0,sim_idx))
# #   abline(colMeans(beta_trace))
# #   Sys.sleep(1)
# # }
# # ts.plot(beta_trace[,2])
# # colMeans(beta_trace)
# # true.beta = beta + c(0,0.5*qnorm(p0))
# # p0=0.1
