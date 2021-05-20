rm(list = ls())

setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')
source('./income_analysis/income_preprocess.R')


data = load_data(out.const = 0.0001,W2.const = 1000)
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

sim_idx = 1
p0 = 0.1
n_xout = 1000
plot(W2,y)


plot(W2,y,main='y log')
p0_list = c(0.1,0.25,0.5,0.75,0.9,0.95,0.99)
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_log_y_out%s_W2%s_%s_%s.RData',out.const,W2.const,p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_Ylog_res$g_trace)
        
        X.Knots = NQR_Ylog_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(xout,exp(y.est)-y.const,type='l',col=2,lwd=3)
}


plot(W2,y,main='W2 log, y log')
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_log_log_out$s_%s_%s_%s.RData',out.const,W2.const,p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_YLog_WLog_res$g_trace)
        
        X.Knots = NQR_YLog_WLog_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(exp(xout)-W2.const,exp(y.est)-y.const,type='l',col=2,lwd=3)
}
par(mfrow=c(1,1))




######################## Previous Analysis ###################################################################

data = load_data(out.const = 0.001,W2.const = 5000)
W1 = data$W1
W2 = data$W2
y = data$y
log_W1 = data$log_W1
log_W2 = data$log_W2
log_y = data$log_y

W1.const = data$W1.const
W2.const = data$W2.const
y.const= data$y.const

Knots.driect = data$Knots.driect

# plot(W2,y,main='W2 log, y log',xlim=c(0,40000))
# p0_list = c(0.1,0.25,0.5,0.75,0.9)
# for(p0 in p0_list){
#   f_name = sprintf('../debugging/income_log_log_%s.RData',p0)
#   load(f_name)
# 
#   g.est = colMedians(NQR_Log_Log_res$g_trace)
# 
#   # points(exp(NQR_Log_Log_res$Knots),exp(g.est),type='o',col=3,lwd=3)
# 
#   X.Knots = NQR_Log_Log_res$Knots
#   fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
#   xout = seq(1,12,length.out = 1000)
#   y.est=predict(fit.ns, data.frame(X.Knots=xout))
# 
#   points(exp(xout),exp(y.est),type='l',col=2,lwd=3)
# }

# log,log 버젼만 우선 해서 y가 제대로 truncated 되었는지를 확인하는게 더 좋겠다.

sim_idx = 1
p0 = 0.1
n_xout = 1000

#### Income Naive needs more iteration!!!------------------------------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(W2,y,main='naive')
p0_list = c(0.1,0.25,0.5,0.75,0.9)
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_naive_%s_%s.RData',p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_res$g_trace)
        cat(g.est[1:10],'\n')
        
        X.Knots = NQR_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(xout,y.est,type='l',col=2,lwd=3)
}

plot(W2,y,main='W log')
p0_list = c(0.1,0.25,0.5,0.75,0.9)
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_log_%s_%s.RData',p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_WLog_res$g_trace)
        cat(g.est[1:10],'\n')
        
        X.Knots = NQR_WLog_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(exp(xout)-W2.const,y.est,type='o',col=2,lwd=3)
}


plot(W2,y,main='y log')
p0_list = c(0.1,0.25,0.5,0.75,0.9)
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_log_y_%s_%s.RData',p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_Ylog_res$g_trace)
        
        X.Knots = NQR_Ylog_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(xout,exp(y.est)-y.const,type='l',col=2,lwd=3)
}


plot(W2,y,main='W2 log, y log')
for(p0 in p0_list){
  y.est_save = matrix(NA,ncol=n_xout,nrow=10)
  for(sim_idx in 1:10){
    tryCatch(
      {
        f_name = sprintf('../debugging/income_log_log_%s_%s.RData',p0,sim_idx)
        load(f_name)
        
        g.est = colMedians(NQR_YLog_WLog_res$g_trace)
        
        X.Knots = NQR_YLog_WLog_res$Knots
        fit.ns <- lm(g.est ~ ns(x = X.Knots, knots = X.Knots[-c(1,length(X.Knots))]))
        xout = seq(min(X.Knots),max(X.Knots),length.out = n_xout)
        y.est.tmp = predict(fit.ns, data.frame(X.Knots=xout))
        y.est_save[sim_idx,] = y.est.tmp
      },
      error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
  }
  y.est = colMeans(y.est_save,na.rm = T)
  points(exp(xout)-W2.const,exp(y.est)-y.const,type='l',col=2,lwd=3)
}
par(mfrow=c(1,1))
