rm(list = ls())

setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')
source('./income_analysis/income_preprocess.R')



#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
out.const = as.numeric(args[1])
W2.const = as.numeric(args[2])
print(sprintf("out.const:%s / W2.const:%s",out.const,W2.const))


out.const = 0.001
W2.const = 1000

data = load_data( out.const = out.const,W2.const = W2.const )
W1 = data$W1
W2 = data$W2
y = data$y
log_W1 = data$log_W1
log_W2 = data$log_W2
log_y = data$log_y

W1.const = data$W1.const
W2.const = data$W2.const
y.const = data$y.const
out.const = data$out.const

Knots.driect = data$Knots.driect
mul_c = 1

p0=0.1
p0_list = c(0.1,0.25,0.5,0.75,0.9,0.95,0.99)
for(sim_idx in 1:10){
  for(p0 in p0_list){
    # f_name = sprintf('../debugging/income_naive_%s_%s.RData',p0,sim_idx)
    # if(!file.exists(file=f_name)){
    #   NQR_res = NQR_w_MME(y,W1,W2,p0,inp.min = min(W2),inp.max = max(W2),multiply_c = 10,inp.version = 1,N.Knots = 30, Knots.direct = Knots.driect)
    #   NQR_res$X_trace = colMeans(NQR_res$X_trace)
    #   save(NQR_res, file=f_name)
    # }
    
    # f_name = sprintf('../debugging/income_log_y_%s_%s.RData',p0,sim_idx)
    f_name = sprintf('../debugging/income_log_y_out%s_W2%s_%s_%s.RData',out.const,W2.const,p0,sim_idx)
    if(!file.exists(file=f_name)){
      NQR_Ylog_res = NQR_w_MME(log_y,W1,W2,p0,inp.min = min(W2),inp.max = max(W2),multiply_c = mul_c,inp.version = 1,N.Knots = 30, Knots.direct = Knots.driect)
      NQR_Ylog_res$X_trace = colMeans(NQR_Ylog_res$X_trace)
      save(NQR_Ylog_res, file=f_name)
    }
    
    
    # f_name = sprintf('../debugging/income_log_log_%s_%s.RData',p0,sim_idx)
    f_name = sprintf('../debugging/income_log_log_out%s_%s_%s_%s.RData',out.const,W2.const,p0,sim_idx)
    if(!file.exists(file=f_name)){
      NQR_YLog_WLog_res = NQR_w_MME(log_y,log_W1,log_W2,p0,inp.min = min(log_W2),inp.max = max(log_W2),multiply_c = mul_c,inp.version = 1,N.Knots = 30)
      NQR_YLog_WLog_res$X_trace = colMeans(NQR_YLog_WLog_res$X_trace)
      save(NQR_YLog_WLog_res, file=f_name)
    }
    
    # f_name = sprintf('../debugging/income_log_%s_%s.RData',p0,sim_idx)
    # if(!file.exists(file=f_name)){
    #   NQR_WLog_res = NQR_w_MME(y,log_W1,log_W2,p0,inp.min = min(log_W2),inp.max = max(log_W2),multiply_c = 10,inp.version = 1,N.Knots = 30)
    #   NQR_WLog_res$X_trace = colMeans(NQR_WLog_res$X_trace)
    #   save(NQR_WLog_res, file=f_name)
    # }
  }
}
# 
# p0=0.9
# Knots.driect = exp(seq(min(log_W2),max(log_W2),length.out = 30))-W2.const
# NQR_naive_Knots.direct = NQR_w_MME(y,W1,W2,p0,inp.min = min(W2),inp.max = max(W2),multiply_c = 1,inp.version = 1,N.Knots = 30,Knots.direct=Knots.driect)
# NQR_Log_Log_res = NQR_w_MME(log_y,log_W1,log_W2,p0,inp.min = min(log_W2),inp.max = max(log_W2),multiply_c = 1,inp.version = 1,N.Knots = 30)
# NQR_Log_y_res = NQR_w_MME(log_y,W1,W2,p0,inp.min = min(W2),inp.max = max(W2),multiply_c = 1,inp.version = 1,N.Knots = 30,Knots.direct=Knots.driect)
# 
# 
# 
# tau.i = Knots.driect
# W2.noise = W2+rnorm(length(W2),0,0.1)
# order_idx = order(W2.noise)
# W2.noise.ord = W2.noise[order_idx]
# y.ord = y[order_idx]
# 
# mspline=spline(x = W2.noise,y = y,xout = tau.i)
# y.est=mspline$y
# points(tau.i,y.est,type='o',col=2,lwd=3)
# 
# fit.ns <- lm(y.ord~ ns(x = W2.noise.ord, knots = tau.i[-c(1,length(tau.i))]) )
# y.est=predict(fit.ns, data.frame(W2.noise.ord=tau.i))
# points(tau.i,y.est,type='o',col=3,lwd=3)
# 
# BQR_res = mBayesQR(y,cbind(1,W2),p0)
# beta.est=colMeans(BQR_res$beta_trace)
# 
# g0=rep(0,30)
# for(param.idx in 1:length(beta.est)){
#   g0 = g0 + beta.est[param.idx]*tau.i^(param.idx-1) # for quadratic g0
# }
# 
# 
# 
# plot(W2,y)
# 
# X.Knots = exp(seq(min(log_W2),max(log_W2),length.out = 30))-W2.const
# xout = seq(min(W2),max(W2),length.out = 1000)
# 
# ####
# model = NQR_naive_Knots.direct
# g.est = colMeans(model$g_trace)
# fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
# y.est=predict(fit.ns, data.frame(X.Knots=xout))
# 
# points(xout,y.est,type='l',col=2,lwd=3)
# 
# ####
# model = NQR_Log_Log_res
# g.est = colMeans(model$g_trace)
# fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
# y.est=predict(fit.ns, data.frame(X.Knots=xout))
# 
# points(xout,exp(y.est)-y.const,type='l',col=3,lwd=3)
# 
# ####
# model = NQR_Log_y_res
# g.est = colMeans(model$g_trace)
# points(exp(tau.i)-W2.const,g0,type='o',col=2)
# fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
# y.est=predict(fit.ns, data.frame(X.Knots=xout))
# 
# points(xout,exp(y.est)-y.const,type='l',col=4,lwd=3)

####