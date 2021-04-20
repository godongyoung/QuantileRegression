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

make_data = function(X){
  set.seed(sim_idx)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  return(list('W1'=W1, 'W2'=W2))
}

make_X_shit = function(Mu_x, Xmul){
  if(X_demean){
    X_shit = Mu_x*Xmul
  }
  else{X_shit = 0}
  
  return(X_shit)
}

sim_idx = 1

# Loop start #############################################################################################

for(sim_idx in start.idx:end.idx){
  for(p0 in p0_list){

    
    # ESL data1 #############################################################################################
    
    set.seed(sim_idx)
    
    inp.sd = 1
    Xmul = 10
    x1 = runif(n,0,1)
    y = sin(12*(x1+0.2))/(x1+0.2) + rnorm(n,0,inp.sd)
    X_shit = make_X_shit(0.5,Xmul)
    X=cbind(1,x1*Xmul-X_shit)
    Xrange = seq(0-X_shit, 1*Xmul-X_shit,length.out = 100)
    x1range = seq(0, 1,length.out = 100)
    plot(X[,2],y)
    for(p0 in p0_list){
      y.p0 = sin(12*(x1range+0.2))/(x1range+0.2) + qnorm(p0,0,inp.sd)
      points(Xrange,y.p0,col=2,lwd=2,type='l')
    }
    
    
    W_list = make_data(X)
    NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
    
    # save result -----------------------------------------------------------------------------------------
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
    save(NQR_res_short, file=sprintf('../debugging/NQR_data1_short_%s_%s_sd%s_NKnots%s_mul%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul))
    if(sim_idx%%20==1){
      save(NQR_res, file=sprintf('../debugging/NQR_data1_%s_%s.RData',p0,sim_idx))
    }
    
    # ESL data2 #############################################################################################
    # X = runif(n,0,1)
    # y = sin(4*X) + rnorm(n,0,inp.sd)
    # plot(X,y)
    # for(p0 in p0_list){
    #   y.p0 = sin(4*Xrange) + qnorm(p0,0,inp.sd)
    #   points(Xrange,y.p0,col=2,lwd=2,type='l')
    # }
    
    # Design Adaptive Nonparametric Regression #############################################################################################
    
    x1 = runif(n,-4,2.5)
    # X = cbind(rnorm(n/2,-1,1), rnorm(n/2,1.75,0.25));hist(X,nclass=100)
    y = sin(2.5*x1) + 0.4*rnorm(n,0,inp.sd)
    x1range = seq(-4,2.5,length.out = 100)
    Xrange = seq(-5,5,length.out = 100)
    X=cbind(1,((x1+4)/6.5)*Xmul-X_shit)
    plot(X[,2],y)
    for(p0 in p0_list){
      y.p0 = sin(2.5*x1range) + 0.4*qnorm(p0,0,inp.sd)
      points(Xrange,y.p0,col=2,lwd=2,type='l')
    }
    
    W_list = make_data(X)
    NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
    
    # save result -----------------------------------------------------------------------------------------
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
    save(NQR_res_short, file=sprintf('../debugging/NQR_data2_short_%s_%s_sd%s_NKnots%s_mul%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul))
    if(sim_idx%%20==1){
      save(NQR_res, file=sprintf('../debugging/NQR_data2_%s_%s.RData',p0,sim_idx)) #old is W1 version! 
    }
    
    # An Introduction to Kernel and Nearest-Neighbor Nonparametric Regression #############################################################################################
    X = runif(n,0,1)
    inp.sd = 0.5
    y = X*sin(2.5*pi*X) + rnorm(n,0,inp.sd)
    Xrange = seq(0,1,length.out = 100)
    # plot(X,y)
    # for(p0 in p0_list){
    #   y.p0 = Xrange*sin(2.5*pi*Xrange) + qnorm(p0,0,inp.sd)
    #   points(Xrange,y.p0,col=2,lwd=2,type='l')
    # }
    X = cbind(1,X)
    W_list = make_data(X)
    NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
    
    # save result -----------------------------------------------------------------------------------------
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
    save(NQR_res_short, file=sprintf('../debugging/NQR_data3_short_%s_%s_sd%s_NKnots%s_mul%s.RData',p0,sim_idx,inp.sd,inp.N.Knots,inp.mul))
    if(sim_idx%%20==1){
      save(NQR_res, file=sprintf('../debugging/NQR_data3_%s_%s.RData',p0,sim_idx)) #old is W1 version! 
    }
    
    # # DPpackage: Bayesian Semi- and Nonparametric Modeling in R ------------------------------------------------------------------------------------------
    # set.seed(0)
    # n <- 500
    # x1i <- runif(n)
    # X = cbind(1,x1i)
    # y1 <- x1i + rnorm(n, 0, sqrt(0.01))
    # y2 <- x1i^4 + rnorm(n, 0, sqrt(0.04))
    # u <- runif(n)
    # prob <- exp(-2 * x1i)
    # y <- ifelse(u < prob, y1, y2)
    # plot(X[,2],y)
    
  }
}

# ###########################################################################################################
# ################################## Real Data #############################################################
# ###########################################################################################################
# 
# 
# # ESL-----------------------------------------------------------------------------------------
# library(ElemStatLearn)
# 
# data("bone") #bone mineral
# head(bone);dim(bone)
# plot(bone$age,bone$spnbmd)
# 
# # Bayesian nonparametric quantile regression ------------------------------------------------------------------------------------------
# # I will follow the procedure they described. 
# # They firstly fit the ncs on mcycle data (didn't specified the knots)
# # And then based on the fitted value of ncs at 30 equally spaced on X, 
# # they generated 1000 replicated with gaussian error on each spaces
# data("mcycle")
# y=mcycle$accel
# x1=mcycle$times
# n=length(y)
# 
# # x1=x1+rnorm(n,0,0.1)
# # order_idx=order(x1)
# # x1=x1[order_idx]
# # y=y[order_idx]
# 
# N=10
# Basis.Knots = seq(min(x1),max(x1),length.out = N)
# plot.Knots = seq(min(x1),max(x1),length.out = 1000)
# plot(x1,y)
# 
# data = data.frame(y,x1)
# colnames(data) = c('y','x1')
# model = lm(y~
#              ns(x1,knots=Basis.Knots[-c(1,length(Basis.Knots))],
#                 # Boundary.knots = Basis.Knots[c(1,length(Basis.Knots))]
#              ),
#            data = data )
# y.est = predict(model,newdata = list(x1=plot.Knots))
# points(plot.Knots,y.est,type='l')
# 
# # start generating data
# N=30
# Knots = seq(min(x1),max(x1),length.out = N)
# y.est.smooth = predict(model,newdata = list(x1=Knots))
# y.generate = sapply(X = y.est.smooth,FUN = function(x) x+rnorm(1000,0,sd = 20))
# y.generate = as.numeric(y.generate)
# x.generate = rep(Knots,each = 1000)
# plot(x.generate,y.generate)
# 
# y.est.smooth.true = predict(model,newdata = list(x1=plot.Knots)) + qnorm(0.95,0,20)
# points(plot.Knots,y.est.smooth.true,type='l',col=4,lwd=3)
# 
# 
# 
# # {
# #   # Chekcing the fitted result of NQR paper -----------------------------------------------------------------
# #   # X.generate = cbind(1,x.generate)
# #   # NQR_res = NQR(y.generate,X.generate,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 1,multiply_c = 6)
# #   # save(NQR_res, file=sprintf('./NQR_v1_fitted_%s.RData',0.95))
# #   
# #   # NQR_res = NQR(y.generate,X.generate,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 3,multiply_c = 6)
# #   # save(NQR_res, file=sprintf('./NQR_v3_fitted_%s.RData',0.95))
# #   
# #   load(file=sprintf('./NQR_v1_fitted_%s.RData',0.95))
# #   # ts.plot(NQR_res$g_trace[,10])
# #   g.est = colmeans(NQR_res$g_trace[3000:5000,])
# #   points(NQR_res$Knots,g.est,type='l',col=2,lwd=3)
# #   
# #   
# #   
# #   # Chekcing the fitted result of Our NQR w MME -----------------------------------------------------------------
# #   n=length(y.generate)
# #   sigma2_11 = 1
# #   sigma2_22 = 1
# #   delta1=rnorm(n,0,sd=sqrt(sigma2_11))
# #   delta2=rnorm(n,0,sd=sqrt(sigma2_22))
# #   
# #   alpha=c(4,3)
# #   W1=X.generate%*%alpha+delta1
# #   W2=X.generate[,1]+delta2
# #   
# #   # NQR_w_MME_res = NQR_w_MME(y.generate,W1,W2,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 1,multiply_c = 6)
# #   # save(NQR_w_MME_res, file=sprintf('./NQRwMME_v1_fitted_%s.RData',0.95))
# #   load(file=sprintf('./NQRwMME_v1_fitted_%s.RData',0.95))
# #   
# #   # ts.plot(NQR_w_MME_res$g_trace[,10])
# #   g.est = colmeans(NQR_w_MME_res$g_trace[3000:5000,])
# #   points(NQR_w_MME_res$Knots,g.est,type='o',col=3,lwd=3)
# #   
# #   plot(x.generate,colMeans(NQR_w_MME_res$X_trace));abline(0,1)
# # }
# 
# 
# 
# 
# 
# # # An Introduction to Kernel and Nearest-Neighbor Nonparametric Regression ------------------------------------------------------------------------------------------
# # chicago = read.table('./DATA/T67.1.txt')[,-c(1,2,3)]
# # colnames(chicago) = c('Zip Code', 'Racial comp', 'Fire', 'Theft' ,'Age', 'Voluntary market activity', 'Involuntary market activity', 'Income')
# # head(chicago)
# # X = chicago$Theft
# # y = chicago$'Voluntary market activity'
# # 
# # adj.const = 10
# # X = rep(X,adj.const)
# # y = rep(y,adj.const)
# # n = length(y)
# # X = X + rnorm(n,0,4)
# # y = y + rnorm(n)
# # X = cbind(1,X)
# # plot(X[,2],y)