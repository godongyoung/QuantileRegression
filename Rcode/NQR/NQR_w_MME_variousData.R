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
beta=-c(3,5)
theta=c(5,-5)
X_demean = T

# sigma2_11=1
# sigma2_22=1
# sigma2_11=0.8**2
# sigma2_22=1.2**2
sigma2_11=1
sigma2_22=2
sigma2_33=1.5
sigma2_44=1

n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)

make_data = function(X,W1.v2=F,W2.v2=F){
  set.seed(sim_idx)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  delta3=rnorm(n,0,sd=sqrt(sigma2_33))
  delta4=rnorm(n,0,sd=sqrt(sigma2_44))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  W3=X%*%beta+delta3
  W4=X%*%theta+delta3
  if(W1.v2){
    W1 = 1/5*X[,1] + 9/5*X[,2] + 1/5*(X[,2])**2 + delta1
  }
  if(W2.v2){
    W2 = -21 + 9.3*log(X[,2]+10) + delta2
  }
  
  return(list('W1'=W1, 'W2'=W2,'W3'=W3,'W4'=W4))
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
p0 = 0.25
is.t=''
data.type=3
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
          # {if(if.W1.scaled){
          #   f_name = sprintf('../debugging/NQR_data%s_%swME_W1scaled_%s_%s.RData',data.type,is.t,p0,sim_idx)
          #   W_list$W1 = (scale(W_list$W1)*sd(W_list$W2))
          # }
          # else{
          #   f_name = sprintf('../debugging/NQR_data%s_%swME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # }}
          # if(sigma2_22!=sigma2_11){
          #   f_name = sprintf('../debugging/NQR_data%s_%swME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
          # }
          # if(!file.exists(file=f_name)){
          #   NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
          #   NQR_res$X_trace = colMeans(NQR_res$X_trace)
          #   NQR_res$mux_trace = mean(NQR_res$mux_trace)
          #   NQR_res$alpha_trace = colMeans(NQR_res$alpha_trace)
          #   save(NQR_res, file=f_name)
          # }
          # 
          # # woME model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%swoME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(!file.exists(file=f_name)){
          #   NQR_wo_ME_res=NQR(y,X,p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
          #   save(NQR_wo_ME_res, file=f_name)
          # }
          # 
          # # naive model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%sW2_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(sigma2_22!=sigma2_11){
          #   f_name = sprintf('../debugging/NQR_data%s_%sW2_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
          # }
          # if(!file.exists(file=f_name)){
          #   NQR_W2_ME_res=NQR(y,cbind(1,W_list$W2),p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots,inp.sigma.g = 0.0025)
          #   save(NQR_W2_ME_res, file=f_name)
          #   NQR_W2_ME_res$g_accept_ratio
          # }

          # # 3ME model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%s3ME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(sigma2_22!=sigma2_11){
          #   f_name = sprintf('../debugging/NQR_data%s_%s3ME_sig1%s_sig2%s_sig3%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,sigma2_33,p0,sim_idx)
          # }
          # if(!file.exists(file=f_name)){
          #   NQR_3ME_res = NQR_w_3ME(y = y,W1 = W_list$W1,W2 = W_list$W2, W3 = W_list$W3, p0 = p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots,Knots.direct = NA,inp.sigma.g = 0.005)
          #   NQR_3ME_res$X_trace = colMeans(NQR_3ME_res$X_trace)
          #   save(NQR_3ME_res, file=f_name)
          # }
          
          # 4ME model --------------------------
          f_name = sprintf('../debugging/NQR_data%s_%s4ME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          if(sigma2_22!=sigma2_11){
            f_name = sprintf('../debugging/NQR_data%s_%s4ME_sig1%s_sig2%s_sig3%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,sigma2_33,p0,sim_idx)
          }
          if(!file.exists(file=f_name)){
            NQR_4ME_res = NQR_w_4ME(y = y,W1 = W_list$W1,W2 = W_list$W2, W3 = W_list$W3, W4 = W_list$W4, p0 = p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots,Knots.direct = NA,inp.sigma.g = 0.005)
            NQR_4ME_res$X_trace = colMeans(NQR_4ME_res$X_trace)
            save(NQR_4ME_res, file=f_name)
          }
          
          # # SME model --------------------------
          # f_name = sprintf('../debugging/NQR_data%s_%swSME_cvar_%s_%s.RData',data.type,is.t,p0,sim_idx) # with inp.sigma.g = 0.01
          # f_name = sprintf('../debugging/NQR_data%s_%swSME_%s_%s.RData',data.type,is.t,p0,sim_idx) # with inp.sigma.g = 0.0025, inp.mul = 10
          # f_name = sprintf('../debugging/NQR_data%s_%swSME_mul%s_%s_%s.RData',data.type,is.t,inp.mul,p0,sim_idx)
          # cheat.fname = sprintf('../debugging/NQR_data%s_%swME_%s_%s.RData',data.type,is.t,p0,sim_idx)
          # if(sigma2_22!=sigma2_11){
          #   f_name = sprintf('../debugging/NQR_data%s_%swSME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx) # with inp.sigma.g = 0.0025, inp.mul = 10
          #   f_name = sprintf('../debugging/NQR_data%s_%swSME_sig1%s_sig2%s_mul%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,inp.mul,p0,sim_idx)
          #   cheat.fname = sprintf('../debugging/NQR_data%s_%swME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
          # }          
          # if(!file.exists(file=f_name)){
          #   NQR_SME_res = list()
          #   NQR_SME_res$g_accept_ratio = 0
          #   # Repeat Fitting until it have higher convergence ratio
          #   while(NQR_SME_res$g_accept_ratio<0.05){
          #     
          #     NQR_SME_res = NQR_w_SME(y = y,W1 = W_list$W1,W2 = W_list$W2,p0 = p0,inp.min = -5,inp.max = 5,inp.version = 1,multiply_c = inp.mul,
          #                             N.Knots = inp.N.Knots,Knots.direct = NA,
          #                             inp.sigma.g = 0.0025, g0.cheat.fname = cheat.fname)
          #   }
          #   NQR_SME_res$X_trace = colMeans(NQR_SME_res$X_trace)
          #   NQR_SME_res$mux_trace = mean(NQR_SME_res$mux_trace)
          #   save(NQR_SME_res, file=f_name)
          # }
          
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


# Function --------------------------------------------------------------------------------
smooth.y=function(knots,g.tau,xout,version=1){
  if(version==1){
    mspline=spline(x = knots,y = g.tau,xout = xout)
    y.est=mspline$y
  }
  if(version==2){
    msmooth.spline=smooth.spline(x = knots,y = g.tau,cv = NA,lambda = lambda.t)
    mspline=predict(msmooth.spline,xout)
    y.est=mspline$y
  }
  if(version==3){
    g.tau=as.numeric(g.tau)
    fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
    y.est=predict(fit.ns, data.frame(knots=xout))
  }
  return(y.est)
}

tmp.p0 = 0.9
gen_y.p0 = function(data.type,tmp.p0){
  if(data.type==1){
    x1range = seq(0, 1,length.out = N)
    y.p0 = sin(12*(x1range+0.2))/(x1range+0.2) 
    res = y.p0 + qt(tmp.p0,df=3)
  }
  if(data.type==2){
    x1range = seq(0, 1,length.out = N)
    y.p0 = sin(10*x1range)
    res = y.p0 + (x1range+0.25)/(0.1)*qnorm(tmp.p0,0,inp.sd)
  }
  if(data.type==3){
    x1range = seq(0, 1,length.out = N)
    y.p0 = x1range*sin(2.5*pi*x1range) 
    res = y.p0 + qnorm(p = tmp.p0,mean = 0,sd = inp.sd)
  }
  return(res)
}
save_data = g_save;type='wME'
m.boxplot=function(save_data,p0,type,data.type,inp.sub=NA){
  valid_cnt = sum(!(is.na(save_data[,1])))
  
  colnames(save_data)=Knots
  ds=cbind(rep(Knots,each=dim(save_data)[1]),as.numeric(save_data))
  colnames(ds)=c('Knots','value')
  
  if(data.type==1){inp.ylim = c(-5,5)}
  if(data.type==2){inp.ylim = c(-3,3)}
  if(data.type==3){inp.ylim = c(-1.5,2)}
  
  
  boxplot(value~Knots, data=ds,ylim=inp.ylim,main=sprintf('%s\'s Box plot of %s with Knots %s for %s simulation',type,p0,inp.N.Knots,valid_cnt),names = round(Knots,1),
          sub=inp.sub) 
  for(tmp.p0 in p0_list){
    y.p0 = gen_y.p0(data.type,tmp.p0)
    points(seq(1:length(Knots)),y.p0,type='l',lwd=2,col=3)
  }
  y.p0 = gen_y.p0(data.type,p0)
  points(seq(1:length(Knots)),y.p0,type='l',lwd=3,col=2)
}

library(ggplot2)
library(MASS)

# df = data.frame(cbind(matrix(save_data,ncol=1),as.factor(rep(Knots[1:30],each = 500))))
# colnames(df) = c('g','Knots')
# df$Knots = as.factor(df$Knots)
# str(df)
# ggplot(df, aes(x = Knots, y = g)) + 
#   geom_boxplot(width=0.8) +
#   ggtitle("Box Plot by Car Type")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

save_data = g_save
save_data_list = g_save_list
m.boxplo.v2=function(save_data_list,p0,data.type){
  par(mfrow = c(1,1))
  fname = sprintf('../Figure/data%s_%s.png',data.type,p0)
  tiff(fname, units="in", width=6*1.5, height=4*1.5, res=600)
  
  if(data.type==1){inp.ylim = c(-5,5)}
  if(data.type==2){inp.ylim = c(-3,3)}
  if(data.type==3){inp.ylim = c(-1.5,2)}
  
  
  plot(NA,NA,xlim=c(-5,5),ylim=inp.ylim,xlab = 'X', ylab = 'y')
  for(tmp.p0 in p0_list){
    y.p0 = gen_y.p0(data.type,tmp.p0)
    points(Knots,y.p0,type='l',lwd=1,lty='dotted',col=alpha(1, alpha = 0.5))
  }
  y.p0 = gen_y.p0(data.type,p0)
  points(Knots,y.p0,type='l',lwd=2,col=1)
  
  for(each in c('woME','wME','W2')){
    g_mean = colMeans(save_data_list[[each]],na.rm = T)
    if(each=='woME'){tmp.lty='twodash';tmp.col=gg_color_hue(3)[1]}
    if(each=='wME'){tmp.lty='longdash';tmp.col=gg_color_hue(3)[2]}
    if(each=='W2'){tmp.lty='dotdash';tmp.col=gg_color_hue(3)[3]}
    points(Knots,g_mean,type='l',col=tmp.col,lty = tmp.lty,lwd=2)
    
    # g_quantile = apply(save_data_list[[each]],MARGIN = 2,FUN = function(x) quantile(x ,probs = c(0.025,0.975),na.rm = T))
    # points(Knots,g_quantile[1,],type='l',col=2,lty='dotted')
    # points(Knots,g_quantile[2,],type='l',col=2,lty='dotted')
  }
  # legend(4, 2.8, legend=c('woME','wME','W2'),
  #        col=2, lty=c('twodash','longdash','dotdash'), cex=1,lwd=2)
  
  dev.off()
}


HPD_ratio = function(g_trace,lb=0.025,ub=0.975){
  g_quantile = apply(g_trace,MARGIN = 2,FUN = function(x) quantile(x ,probs = c(lb,ub)))
  y.p0 = gen_y.p0(data.type,tmp.p0 = p0 )
  diff_mat = (g_quantile)-matrix(rep(y.p0,each = 2),nrow=2)
  HPD_include = (sign(diff_mat[1,]) * sign(diff_mat[2,]) ) == rep(-1,N)
  return(mean(HPD_include))
}

make_X_shit = function(Mu_x, Xmul){
  if(X_demean){
    X_shit = Mu_x*Xmul
  }
  else{X_shit = 0}
  return(X_shit)
}

true.g.func = function(inp.x,tmp.p0 = p0){
  if(data.type==1){
    inp.x = (inp.x+X_shit)/Xmul
    y.p0 = sin(12*(inp.x+0.2))/(inp.x+0.2) 
    res = y.p0 + qt(tmp.p0,df=3)
  }
  if(data.type==2){
    inp.x = (inp.x+X_shit)/Xmul
    y.p0 = sin(10*inp.x) 
    res = y.p0 + (inp.x+0.25)/(0.1)*qnorm(tmp.p0,0,inp.sd)
  }
  if(data.type==3){
    inp.x = (inp.x+X_shit)/Xmul
    y.p0 = inp.x*sin(2.5*pi*inp.x) 
    res = y.p0 + qnorm(p = tmp.p0,mean = 0,sd = inp.sd)
  }
  return(res)
}

ise_func = function(x){
  # return (((true.g.func(x)-est.g.func(x))/sd(true.g.func(integral_range)))**2)
  return ((true.g.func(x)-est.g.func(x))**2)
}

# g.tau=as.numeric(g.est);knots = Knots
# fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
# y.est=predict(fit.ns, data.frame(knots=integral_range))
est.g.func = function(inp.x){
  predict(fit.ns, data.frame(knots=inp.x))
}
# save_data = g_save;type='wME';data.type=1

# Define True parameter--------------------------------------------------------------------------------
