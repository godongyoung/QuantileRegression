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
# start.idx = 1
# end.idx = 500
print(sprintf("start:%s / end:%s",start.idx,end.idx))


########
# Define True parameter--------------------------------------------------------------------------------
if.short = T
if.NQR_wo_ME = T
inp.N.Knots = 20
n=1000
alpha=c(4,3)
inp.mul = 6

Mu_x=3*pi/2
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1

Xmul = 5
Xmul_list = c(2,3,5)
# Simulation start --------------------------------------------------------------------------------
sim_idx=1
nmax=500
p0=0.1

for(inp.N.Knots in c(30)){
  for(sim_idx in start.idx:end.idx){
    for(Xmul in Xmul_list){
      tryCatch(
        {
          set.seed(sim_idx)
          # Make data--------------------------------------------------------------------------------
          
          # x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sqrt(sigma2_xx))
          x1i=runif(n=n,min=0,max=2*Mu_x)
          # x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
          
          inp.sd=1
          y=2+sin(x1i)+rnorm(n,0,inp.sd)    
          
          
          X=cbind(1,x1i*Xmul)
          
          #generate W1,W2
          delta1=rnorm(n,0,sd=sqrt(sigma2_11))
          delta2=rnorm(n,0,sd=sqrt(sigma2_22))
          
          W1=X%*%alpha+delta1
          W2=X[,2]+delta2
          
          for(p0 in c(0.1,0.25,0.5,0.75,0.9)){
            
            # fit NQR w MME--------------------------------------------------------------------------------
            NQR_res=NQR_w_MME(y,W1,W2,p0,inp.min = 0,inp.max = 2*Mu_x*Xmul,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
            
            if(if.short){
              NQR_res_short = list() 
              NQR_res_short[['g.est']] = colMeans(NQR_res$g_trace)
              NQR_res_short[['lambda.est']] = mean(NQR_res$lambda_trace)
              NQR_res_short[['g_accept_ratio']] =NQR_res$g_accept_ratio
              NQR_res_short[['l_accept_ratio']] = NQR_res$l_accept_ratio
              NQR_res_short[['X.est']] = colMeans(NQR_res$X_trace)
              NQR_res_short[['x_accept_ratio']] = NQR_res$x_accept_ratio
              NQR_res_short[['alpha.est']] = colMeans(NQR_res$alpha_trace)
              NQR_res_short[['mux.est']] = mean(NQR_res$mux_trace)
              NQR_res_short[['sigma2_22.est']] = mean(NQR_res$sigma2_22_trace)
              NQR_res_short[['sigma2_11.est']] = mean(NQR_res$sigma2_11_trace)
              NQR_res_short[['sigma2_xx.est']] = mean(NQR_res$sigma2_xx_trace)
              NQR_res_short[['Knots']] = NQR_res$Knots
              NQR_res_short[['inp.version']]=NQR_res$inp.version
              save(NQR_res_short, file=sprintf('../debugging/NQR_short_%s_%s_sd%s_NKnots%s_mul%s_Xrange%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul,Xmul))
              
              if(sim_idx%%20==1){
                save(NQR_res, file=sprintf('../debugging/NQR_%s_%s_Xrange%s.RData',p0,sim_idx,Xmul))
              }
            }
            else{
              save(NQR_res, file=sprintf('../debugging/NQR_%s_%s_Xrange%s.RData',p0,sim_idx,Xmul))
            }
            
            # fit NQR wo MME--------------------------------------------------------------------------------
            if(if.NQR_wo_ME){
              NQR_wo_ME_res=NQR(y,X,p0,inp.min = 0,inp.max = 2*Mu_x,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
              NQR_res_woME_short = list() 
              NQR_res_woME_short[['g.est']] = colMeans(NQR_wo_ME_res$g_trace)
              NQR_res_woME_short[['lambda.est']] = mean(NQR_wo_ME_res$lambda_trace)
              NQR_res_woME_short[['g_accept_ratio']] =NQR_wo_ME_res$g_accept_ratio
              NQR_res_woME_short[['l_accept_ratio']] = NQR_wo_ME_res$l_accept_ratio
              NQR_res_woME_short[['Knots']] = NQR_wo_ME_res$Knots
              NQR_res_woME_short[['inp.version']]=NQR_wo_ME_res$inp.version
              save(NQR_res_woME_short, file=sprintf('../debugging/NQR_woME_short_%s_%s_sd%s_NKnots%s_mul%s_Xrange%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul,Xmul))
            }
            
          }
        },
        error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))      
    }
  }  
}

g.est = colMeans(NQR_res$g_trace)
Knots = NQR_res$Knots

plot(X[,2],y)
points(Knots,g.est,type='l',col=2,lwd=3)
# 
# 
# tau.i=NQR_res$Knots
# g.est=colMeans(NQR_res$g_trace,na.rm = T)
# plot(X[,2],y)
# points(tau.i,g.est,type='l')
# 
# #Result check-----------------------------------------------------------------------------------------------
# 
# idx=2
# ts.plot(X_trace[,idx]);abline(h=X[idx,2])
# 
# g.est=colMeans(g_trace,na.rm = T)
# plot(X[,2],y)
# points(tau.i,g.est,type='l')
# 
# X.est=colMeans(X_trace,na.rm = T)
# alpha.est=colMeans(alpha_trace,na.rm = T)
# plot(X.est,y)
# points(tau.i,g.est,type='l')
# plot(X.est,X[,2]);abline(0,1);abline(v=tau.i)
# plot(X.est,W1);abline(alpha.est)
# plot(X.est,W2);abline(0,1)
# 
# mean(sigma2_11_trace,na.rm = T)
# mean(sigma2_22_trace,na.rm = T)
# mean(sigma2_xx_trace,na.rm = T)
