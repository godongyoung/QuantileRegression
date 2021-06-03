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

# sigma2_11=1
# sigma2_22=1
# sigma2_11=0.8**2
# sigma2_22=1.2**2
sigma2_11=1
sigma2_22=2

n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)

make_data = function(X,W1.v2=F,W2.v2=F){
  set.seed(sim_idx)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  if(W1.v2){
    W1 = 1/5*X[,1] + 9/5*X[,2] + 1/5*(X[,2])**2 + delta1
  }
  if(W2.v2){
    W2 = -21 + 9.3*log(X[,2]+10) + delta2
  }
  
  return(list('W1'=W1, 'W2'=W2))
}

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
p0 = 0.9
is.t=''
data.type=1
if.W1.scaled=F
# Loop start #############################################################################################

for(sim_idx in start.idx:end.idx){
  if(sim_idx>100){next}
  for(p0 in p0_list){
    # for(data.type in c(1,2,3)){
    for(data.type in c(3)){
      
      ## Data Gen #############################################################################################
      if(data.type==1){
        is.t='t_'
        
        # ESL data1 #############################################################################################
        set.seed(sim_idx)
        inp.sd = 1
        Xmul = 10
        X_shit = make_X_shit(0.5,Xmul)
        x1 = runif(n,0,1)
        y = sin(12*(x1+0.2))/(x1+0.2) + rt(n,df=3)
        X=cbind(1,x1*Xmul-X_shit)
        Xrange = seq(0-X_shit, 1*Xmul-X_shit,length.out = 100)
        x1range = seq(0, 1,length.out = 100)
        plot(X[,2],y)
        for(p0.tmp in p0_list){
          y.p0 = sin(12*(x1range+0.2))/(x1range+0.2) + qt(p0.tmp,df=3)
          points(Xrange,y.p0,col=2,lwd=2,type='l')
        }
        W_list = make_data(X)
      }

      if(data.type==2){
        is.t=''
        # data2 #############################################################################################
        # Fast Nonparametric Quantile Regression With Arbitrary Smoothing Methods
        # Heteroscedasticity
        set.seed(sim_idx)
        x1 = (seq(1,n)-1)/n
        inp.sd = 0.07
        y = sin(10*x1) + (x1+0.25)/(0.1)*rnorm(n,0,sd=inp.sd)
        X_shit = make_X_shit(0.5,Xmul)
        X=cbind(1,x1*Xmul-X_shit)
        Xrange = seq(0-X_shit, 1*Xmul-X_shit,length.out = 100)
        x1range = seq(0, 1,length.out = 100)
        
        plot(X[,2],y)
        for(tmp.p0 in p0_list){
          y.p0 = sin(10*x1range) + (x1range+0.25)/(0.1)*qnorm(tmp.p0,0,inp.sd)
          points(Xrange,y.p0,type='l',col=2,lwd=3)
        }
        W_list = make_data(X)
      }
      if(data.type==3){
        is.t=''
        Xmul = 10
        X_shit = make_X_shit(0.5,Xmul)
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
        W_list = make_data(X)
      }
      
      tryCatch(
        {
          ## Modeling for each type of data #############################################################################################
          
          # Our model --------------------------
          {if(if.W1.scaled){
            f_name = sprintf('../debugging/NQR_data%s_%swME_W1scaled_%s_%s.RData',data.type,is.t,p0,sim_idx)
            W_list$W1 = (scale(W_list$W1)*sd(W_list$W2))
          }
          else{
            f_name = sprintf('../debugging/NQR_data%s_%swME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          }}
          if(sigma2_22!=sigma2_11){
            f_name = sprintf('../debugging/NQR_data%s_%swME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
          }
          if(!file.exists(file=f_name)){
            NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
            NQR_res$X_trace = colMeans(NQR_res$X_trace)
            NQR_res$mux_trace = mean(NQR_res$mux_trace)
            NQR_res$alpha_trace = colMeans(NQR_res$alpha_trace)
            save(NQR_res, file=f_name)
          }
          
          # # woME model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%swoME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(!file.exists(file=f_name)){
          #   NQR_wo_ME_res=NQR(y,X,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
          #   save(NQR_wo_ME_res, file=f_name)
          # }
          # 
          # naive model --------------------------
          f_name = sprintf('../debugging/NQR_data%s_%sW2_%s_%s.RData',data.type,is.t,p0,sim_idx)
          if(sigma2_22!=sigma2_11){
            f_name = sprintf('../debugging/NQR_data%s_%sW2_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
          }
          if(!file.exists(file=f_name)){
            NQR_W2_ME_res=NQR(y,cbind(1,W_list$W2),p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
            save(NQR_W2_ME_res, file=f_name)
          }
          if(file.exists(file=f_name)){
            load(f_name)
            if(NQR_W2_ME_res$g_accept_ratio<0.1){
              NQR_W2_ME_res=NQR(y,cbind(1,W_list$W2),p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = 30,N.Knots = inp.N.Knots)
              save(NQR_W2_ME_res, file=f_name)
            }
          }

          # # SME model --------------------------
          f_name = sprintf('../debugging/NQR_data%s_%swSME_cvar_%s_%s.RData',data.type,is.t,p0,sim_idx)
          if(!file.exists(file=f_name)){
            
            NQR_SME_res = list()
            NQR_SME_res$g_accept_ratio = 0
            # Repeat Fitting until it have higher convergence ratio
            while(NQR_SME_res$g_accept_ratio>0.1){
              NQR_SME_res = NQR_w_SME(y = y,W1 = W_list$W1,W2 = W_list$W2,p0 = p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots,Knots.direct = NA,inp.sigma.g = 0.05)
            }
          }
          NQR_SME_res$X_trace = colMeans(NQR_SME_res$X_trace)
          NQR_SME_res$mux_trace = mean(NQR_SME_res$mux_trace)
          save(NQR_SME_res, file=f_name)
          
          # # MME_cvar model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%swME_cvar_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(!file.exists(file=f_name)){
          #   NQR_cvar_res = NQR_w_MME_cvar(y = y,W1 = W_list$W1,W2 = W_list$W2,p0 = p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots,Knots.direct = NA)
          #   NQR_cvar_res$X_trace = colMeans(NQR_cvar_res$X_trace)
          #   NQR_cvar_res$mux_trace = mean(NQR_cvar_res$mux_trace)
          #   NQR_cvar_res$alpha_trace = colMeans(NQR_cvar_res$alpha_trace)
          #   save(NQR_cvar_res, file=f_name)
          # }
        },
        error = function(e) cat('Error \n'))      
    }
  }
}
