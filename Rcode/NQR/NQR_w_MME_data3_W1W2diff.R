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


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
start.idx = as.numeric(args[1])
end.idx = as.numeric(args[2])
print(sprintf("start:%s / end:%s",start.idx,end.idx))


###########################################################################################################
################################## Simulated Data #########################################################
###########################################################################################################

# Define True parameter--------------------------------------------------------------------------------
if.short = T
if.NQR_wo_ME = F
inp.N.Knots = 30
inp.mul = 10
n=1000
alpha=c(4,3)
X_demean = T

sigma2_11=1
sigma2_22=1

n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)

make_data = function(X,W1.v2=F,W2.v2=F,version=1){
  set.seed(sim_idx)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  if(W1.v2){
    if(version==1){
      tmp.seed=8
      set.seed(tmp.seed)
      a = runif(1,1/15,1/10)
      b = runif(1,1/6,1/3)
      d = runif(1,-5,5)
      c = runif(1,3,6)
      W1 = (a*X[,2]**3 + b*X[,2]**2 + d*X[,2] -c)*1/4 + delta1
      }
    if(version==2){
      tmp.seed=5
      set.seed(tmp.seed) #1, 9
      a = runif(1,min = 3,max = 5)
      b = -1*runif(1,1/10,1/5)
      c = runif(1,-3,3)
      W1 = b*(X[,2]-a)**2+c*X[,2]+delta1
    }
    if(version==3){
      W1 = 1/5*X[,1] + 9/5*X[,2] + 1/5*(X[,2])**2 + delta1
    }
  }
  if(W2.v2){
    if(version==1){
      tmp.seed=8
      set.seed(tmp.seed) #
      a = runif(1,min = 9.25,max = 9.35)
      b = runif(1,9.5,10.5)
      c = -runif(1,20.5,21.5)
      W2 = a*log(X[,2]+b) + c+ delta2
      }
    
    if(version==2){
      tmp.seed=8
      set.seed(tmp.seed) #1, 8
      a = runif(1,min = 5,max = 15)
      b = runif(1,5,25)
      tmp = a*log(X[,2]+b) + delta2
      c = min(tmp)+ifelse(tmp.seed==1,6,8)
      W2 = tmp-c
    }
    if(version==3){
      W2 = -21 + 9.3*log(X[,2]+10) + delta2
    }
  }

  return(list('W1'=W1, 'W2'=W2))
}


# version=2
# if(version==1){tmp.seed=9}
# if(version==2){tmp.seed=1}
# 
# set.seed(tmp.seed) #1, 9
# a = runif(1,min = 3,max = 5)
# b = -1*runif(1,1/10,1/5)
# c = runif(1,-3,3)
# W1 = b*(X[,2]-a)**2+c*X[,2]+delta1
# cat(a**2*b,-2*a*b+c,b)
# 
# 
# if(version==1){tmp.seed=1}
# if(version==2){tmp.seed=8}
# 
# set.seed(tmp.seed) #1, 8
# a = runif(1,min = 5,max = 15)
# b = runif(1,5,25)
# tmp = a*log(X[,2]+b) + delta2
# c = min(tmp)+ifelse(tmp.seed==1,6,8)
# W2 = tmp-c
# cat(-c, a,b)

# par(mfrow=c(1,1))
# set.seed(sim_idx)
# x1 = runif(n,0,1)
# inp.sd = 0.5
# y = x1*sin(2.5*pi*x1) + rnorm(n,0,inp.sd)
# X=cbind(1,x1*Xmul-X_shit)
# Xrange = seq(-5, 5,length.out = 100)
# x1range = seq(0, 1,length.out = 100)
# plot(X[,2],y)
# for(p0.tmp in p0_list){
#   y.p0 = x1range*sin(2.5*pi*x1range) + qnorm(p0.tmp,0,inp.sd)
#   points(Xrange,y.p0,col=2,lwd=2,type='l')
# }
# W_list = make_data(X,W1.V2,W2.V2)
# plot(X[,2],W_list$W2);abline(0,1)
# plot(X[,2],W_list$W1);abline(lm(X[,2]~W_list$W1)$coeff)
# 
# # for deviating linear relationship------
# tmp.seed=1
# set.seed(5) #1, 9
# a = runif(1,min = 3,max = 5)
# b = -1*runif(1,1/10,1/5)
# c = runif(1,-3,3)
# tmp = b*(X[,2]-a)**2+c*X[,2]+delta1
# cat(b,-2*a*b+c,a**2*b)
# plot(X[,2],tmp);abline(lm(tmp~X[,2])$coeff)
# tmp.seed=tmp.seed+1
# 
# tmp.seed=1
# set.seed(5)  #5
# a = runif(1,min = 3,max = 5)
# b = runif(1,1/10,1/5)
# c = runif(1,-3,3)
# tmp = b*(X[,2]-a)**2+c*X[,2]+delta1
# cat(b,-2*a*b+c,a**2*b)
# plot(X[,2],tmp);abline(lm(tmp~X[,2])$coeff)
# tmp.seed=tmp.seed+1
# 
# # for deviating y = x relationship------
# tmp.seed=1
# set.seed(tmp.seed) #1, 8
# a = runif(1,min = 5,max = 15)
# b = runif(1,5,25)
# tmp = a*log(X[,2]+b) + delta2
# c = min(tmp)+ifelse(tmp.seed==1,6,8)
# cat(a,b,c)
# plot(X[,2],tmp-c);abline(0,1)
# tmp.seed=tmp.seed+1


# for making nonlinear relationship
# plot(X[,2],W1)
# plot(X[,2],(1/5*X[,1]+9/5*X[,2]+1/5*(X[,2])**2+delta1))
# y.tmp = (1/5*X[,1]+9/5*X[,2]+1/5*(X[,2])**2+delta1)
# lm1 = lm(y.tmp~X[,2])
# abline(lm1$coefficients)
# 
# plot(X[,2],W2)
# lm(X[,2] ~ I(X[,2]^2) + I(X[,2]^3))
# plot(X[,2],0.08*X[,1] + -0.0025*X[,2]**2 + 0.056*X[,2]**3+delta2);abline(0,1)
# plot(X[,2],-21+9.3*log(X[,2]+10)+delta2);abline(0,1)

make_X_shit = function(Mu_x, Xmul){
  if(X_demean){
    X_shit = Mu_x*Xmul
  }
  else{X_shit = 0}
  
  return(X_shit)
}

sim_idx = 15
p0 = 0.1
is.t=''
data.type=3
W1.V2=T;W2.V2=F
inp.version = 2
# Loop start #############################################################################################
# plot(X[,2],W_list$W2);abline(0,1)
# plot(X[,2],W_list$W1);abline(lm(W_list$W1~X[,2])$coeff)
for(sim_idx in start.idx:end.idx){
  for(p0 in p0_list){
    for(inp.version in c(1)){ 
      # if(W1.W2.case==1){W1.V2=T;W2.V2=F}
      # if(W1.W2.case==2){W1.V2=F;W2.V2=T}
      # if(W1.W2.case==3){W1.V2=T;W2.V2=T}
      
      ## Data Gen #############################################################################################
      Xmul = 10
      X_shit = make_X_shit(0.5,Xmul)
      if(data.type==3){
        is.t=''
        # # An Introduction to Kernel and Nearest-Neighbor Nonparametric Regression #############################################################################################
        set.seed(sim_idx)
        x1 = runif(n,0,1)
        inp.sd = 0.5
        y = x1*sin(2.5*pi*x1) + rnorm(n,0,inp.sd)
        X=cbind(1,x1*Xmul-X_shit)
        Xrange = seq(-5, 5,length.out = 100)
        x1range = seq(0, 1,length.out = 100)
        plot(X[,2],y)
        for(p0.tmp in p0_list){
          y.p0 = x1range*sin(2.5*pi*x1range) + qnorm(p0.tmp,0,inp.sd)
          points(Xrange,y.p0,col=2,lwd=2,type='l')
        }
        W_list = make_data(X,W1.V2,W2.V2,version = inp.version)
      }
      
      tryCatch(
        {
          ## Modeling for each type of data #############################################################################################
          f_name = sprintf('../debugging/NQR_data%s_W1%sW2%s_ver%s_%swME_%s_%s.RData',data.type,W1.V2,W2.V2,inp.version,is.t,p0,sim_idx)
          if(!file.exists(file=f_name)){
            NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
            NQR_res$X_trace = colMeans(NQR_res$X_trace)
            NQR_res$mux_trace = mean(NQR_res$mux_trace)
            NQR_res$alpha_trace = colMeans(NQR_res$alpha_trace)
            save(NQR_res, file=f_name)
          }
          
          # f_name = sprintf('../debugging/NQR_data%s_W1%s_W2%s_%swoME_%s_%s.RData',data.type,W1.V2,W2.V2,is.t,p0,sim_idx)
          # if(!file.exists(file=f_name)){
          #   NQR_wo_ME_res=NQR(y,X,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
          #   save(NQR_wo_ME_res, file=f_name)
          # }
          
          if(W2.V2==T){
            f_name = sprintf('../debugging/NQR_data%s_W1%sW2%s_ver%s_%sW2_%s_%s.RData',data.type,W1.V2,W2.V2,inp.version,is.t,p0,sim_idx) 
            if(!file.exists(file=f_name)){
              NQR_W2_ME_res=NQR(y,cbind(1,W_list$W2),p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
              save(NQR_W2_ME_res, file=f_name)
              }
          }
          # if(file.exists(file=f_name)){
          #   load(f_name)
          #   if(NQR_W2_ME_res$g_accept_ratio<0.1){
          #     NQR_W2_ME_res=NQR(y,cbind(1,W_list$W2),p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = 30,N.Knots = inp.N.Knots)
          #     save(NQR_W2_ME_res, file=f_name)
          #   }
          # }
        },
        error = function(e) cat('Error \n'))      
    }
  }
}
