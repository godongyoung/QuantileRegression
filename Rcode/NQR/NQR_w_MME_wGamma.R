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


########
# Define True parameter--------------------------------------------------------------------------------
if.short = T
if.NQR_wo_ME = F # we already did this in NQR_w_MME.R !
inp.N.Knots = 30
inp.mul = 10
n=1000
alpha=c(4,3)
gamma=c(6,5)

# Mu_x=5
Mu_x=3*pi/2
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
sim_idx=1
nmax=500

for(inp.N.Knots in c(30)){
  for(sim_idx in start.idx:end.idx){
    tryCatch(
      {
        set.seed(sim_idx)
        # Make data--------------------------------------------------------------------------------
        
        # x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sqrt(sigma2_xx))
        x1i=runif(n=n,min=0,max=2*Mu_x)
        # x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
        X=cbind(1,x1i)
        inp.sd=1
        y=2+sin(x1i)+rnorm(n,0,inp.sd)    
        
        #generate W1,W2
        delta1=rnorm(n,0,sd=sqrt(sigma2_11))
        delta2=rnorm(n,0,sd=sqrt(sigma2_22))
        
        W1=X%*%alpha+delta1
        W2=X%*%gamma+delta2
        
        for(p0 in c(0.1,0.25,0.5,0.75,0.9)){
          # fit NQR w MME_Gamma--------------------------------------------------------------------------------
          NQR_Gamma_res=NQR_w_MME_Gamma(y,W1,W2,p0,inp.min = 0,inp.max = 2*Mu_x,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
          
          if(if.short){
            NQR_Gamma_res_short = list() 
            NQR_Gamma_res_short[['g.est']] = colMeans(NQR_Gamma_res$g_trace)
            NQR_Gamma_res_short[['lambda.est']] = mean(NQR_Gamma_res$lambda_trace)
            NQR_Gamma_res_short[['g_accept_ratio']] =NQR_Gamma_res$g_accept_ratio
            NQR_Gamma_res_short[['l_accept_ratio']] = NQR_Gamma_res$l_accept_ratio
            NQR_Gamma_res_short[['X.est']] = colMeans(NQR_Gamma_res$X_trace)
            NQR_Gamma_res_short[['x_accept_ratio']] = NQR_Gamma_res$x_accept_ratio
            NQR_Gamma_res_short[['alpha.est']] = colMeans(NQR_Gamma_res$alpha_trace)
            NQR_Gamma_res_short[['gamma.est']] = colMeans(NQR_Gamma_res$gamma_trace)
            NQR_Gamma_res_short[['mux.est']] = mean(NQR_Gamma_res$mux_trace)
            NQR_Gamma_res_short[['sigma2_22.est']] = mean(NQR_Gamma_res$sigma2_22_trace)
            NQR_Gamma_res_short[['sigma2_11.est']] = mean(NQR_Gamma_res$sigma2_11_trace)
            NQR_Gamma_res_short[['sigma2_xx.est']] = mean(NQR_Gamma_res$sigma2_xx_trace)
            NQR_Gamma_res_short[['Knots']] = NQR_Gamma_res$Knots
            NQR_Gamma_res_short[['inp.version']]=NQR_Gamma_res$inp.version
            save(NQR_Gamma_res_short, file=sprintf('../debugging/NQR_wGamma_short_%s_%s_sd%s_NKnots%s_mul%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul)) #old is W1 version!
            if(sim_idx%%20==1){
              save(NQR_Gamma_res, file=sprintf('../debugging/NQR_wGamma_%s_%s.RData',p0,sim_idx)) #old is W1 version! 
            }
          }
          else{
            save(NQR_Gamma_res, file=sprintf('../debugging/NQR_wGamma_%s_%s.RData',p0,sim_idx)) #old is W1 version!
          }
          
          # fit NQR wo MME--------------------------------------------------------------------------------
          if(if.NQR_wo_ME){
            NQR_wo_ME_res=NQR(y,X,p0,inp.min = 0,inp.max = 2*Mu_x,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
            NQR_Gamma_res_woME_short = list() 
            NQR_Gamma_res_woME_short[['g.est']] = colMeans(NQR_wo_ME_res$g_trace)
            NQR_Gamma_res_woME_short[['lambda.est']] = mean(NQR_wo_ME_res$lambda_trace)
            NQR_Gamma_res_woME_short[['g_accept_ratio']] =NQR_wo_ME_res$g_accept_ratio
            NQR_Gamma_res_woME_short[['l_accept_ratio']] = NQR_wo_ME_res$l_accept_ratio
            NQR_Gamma_res_woME_short[['Knots']] = NQR_wo_ME_res$Knots
            NQR_Gamma_res_woME_short[['inp.version']]=NQR_wo_ME_res$inp.version
            save(NQR_Gamma_res_woME_short, file=sprintf('../debugging/NQR_woME_short_%s_%s_sd%s_NKnots%s_mul%s_Xunif.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul))
            
          }
        }
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))

  }  
}


# 
# #Result check-----------------------------------------------------------------------------------------------
# 
# ts.plot(NQR_Gamma_res$gamma_trace[,1])
# ts.plot(NQR_Gamma_res$gamma_trace[,2])
# ts.plot(NQR_Gamma_res$alpha_trace[,1])
# ts.plot(NQR_Gamma_res$alpha_trace[,2])
# ts.plot(NQR_Gamma_res$g_trace[,10])
# ts.plot(NQR_Gamma_res$X_trace[,10]);abline(h=X[10,2])
# colMeans(NQR_Gamma_res$gamma_trace);gamma
# colMeans(NQR_Gamma_res$alpha_trace);alpha
# mean(NQR_Gamma_res$sigma2_11_trace)
# mean(NQR_Gamma_res$sigma2_22_trace)
# mean(NQR_Gamma_res$sigma2_xx_trace)
# plot(X[,2],colMeans(NQR_Gamma_res$X_trace));abline(0,1)
# g.est=colMeans(NQR_Gamma_res$g_trace)
# plot(X[,2],y)
# lines(NQR_Gamma_res$Knots,g.est)
# y.true=2+sin(NQR_Gamma_res$Knots)+qnorm(p0,0,inp.sd)
# lines(NQR_Gamma_res$Knots,y.true,col=2)
