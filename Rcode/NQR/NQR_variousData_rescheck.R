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
library(MASS)

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
n=1000
alpha=c(4,3)
beta=-c(3,5)
theta=c(5,-5)

# sigma2_11=1
# sigma2_22=1
# sigma2_11=0.8**2
# sigma2_22=1.2**2
sigma2_11=1
sigma2_22=2
sigma2_33=1.5
sigma2_44=1

# Simulation start --------------------------------------------------------------------------------
nmax=100
if(sigma2_22!=sigma2_11){
  nmax=100
}

is.plot=F
to_see_list = c('wME','woME','W2','SME')
# to_see_list = c('wME')
# to_see_list = c('wME','woME')
# to_see_list = c('wME','SME','cvar','woME','W2')
to_see = to_see_list[1]
inp.N.Knots = 30
inp.mul = 50
N = inp.N.Knots


X_demean = T
Xmul = 10
X_shit = make_X_shit(0.5,Xmul)
integral_range = seq(-5,5,length.out = 1000)
cal_MSE = T

# Simulation check --------------------------------------------------------------------------------




n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)
p0=0.1
sim_idx=1

data.type = 3
is.t=''
inp.sd = 1
if(data.type==1){
  is.t='t_'
}
if(data.type==2){
  inp.sd = 0.07
}
if(data.type==3){
  inp.sd = 0.5  
}

# to_see_list = c('SME')
to_see_list = c('wME')
# to_see_list = c('W2')
to_see = to_see_list[1]

summary_list = list()
par(mfrow=c(length(p0_list),1))
if(length(to_see_list)>1){par(mfrow=c(length(to_see_list),1))}

for(p0 in p0_list){
# for(p0 in c(0.1)){
  g_save_list = list()
  for(to_see in to_see_list){
    # Define save matrix for current setting------------------------------------------------------------------------------------------
    
    accept_g_save = rep(NA,nmax)
    g_save=matrix(NA,ncol=inp.N.Knots,nrow=nmax)
    HPD_save = rep(NA,nmax)
    ISE_save = rep(NA,nmax)
    MSE_save = rep(NA,nmax)
    g_var_save = rep(NA,nmax)
    
    alpha_save = matrix(NA,ncol=2,nrow=nmax)
    beta_save = matrix(NA,ncol=2,nrow=nmax)
    theta_save = matrix(NA,ncol=2,nrow=nmax)
    X_save=matrix(NA,ncol=n,nrow=nmax)
    mux_save=rep(NA,nmax)
    sigma2_11_save=rep(NA,nmax)
    sigma2_22_save=rep(NA,nmax)
    sigma2_33_save=rep(NA,nmax)
    sigma2_44_save=rep(NA,nmax)
    sigma2_xx_save=rep(NA,nmax)
    
    
    summary_list[[as.character(p0)]][[to_see]][['ISE']]=list()
    summary_list[[as.character(p0)]][[to_see]][['HPD']]=list()
    summary_list[[as.character(p0)]][[to_see]][['MSE']]=list()
    
    # Load each iteration result------------------------------------------------------------------------------------------
    tic()
    for(sim_idx in 1:nmax){  
      tryCatch(
        {
          
          if(to_see == 'wME'){
            f_name = sprintf('../debugging/NQR_data%s_%swME_%s_%s.RData',data.type,is.t,p0,sim_idx)
            if(sigma2_22!=sigma2_11){
              f_name = sprintf('../debugging/NQR_data%s_%swME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
            }
            load(file= f_name)
            g.est=colMeans(NQR_res$g_trace)
            Knots = NQR_res$Knots
            if(is.null(dim(NQR_res$alpha_trace))){
              alpha_save[sim_idx,]=(NQR_res$alpha_trace)
            }
            else{
              alpha_save[sim_idx,]=colMeans(NQR_res$alpha_trace)
            }
            # X.est=colMeans(NQR_res$X_trace)
            mux.est=mean(NQR_res$mux_trace)
            sigma2_11.est=mean(NQR_res$sigma2_11_trace)
            sigma2_22.est=mean(NQR_res$sigma2_22_trace)
            sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
            sigma2_11_save[sim_idx] = sigma2_11.est
            sigma2_22_save[sim_idx] = sigma2_22.est
            sigma2_xx_save[sim_idx] = sigma2_xx.est
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_res$g_trace)
          }
          
          if(to_see == 'SME'){
            f_name = sprintf('../debugging/NQR_data%s_%swSME_cvar_%s_%s.RData',data.type,is.t,p0,sim_idx) # with inp.sigma.g = 0.01
            f_name = sprintf('../debugging/NQR_data%s_%swSME_%s_%s.RData',data.type,is.t,p0,sim_idx)
            if(sigma2_22!=sigma2_11){
              f_name = sprintf('../debugging/NQR_data%s_%swSME_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
              if(inp.mul!=10){
                f_name = sprintf('../debugging/NQR_data%s_%swSME_sig1%s_sig2%s_mul%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,inp.mul,p0,sim_idx)
              }
            }
            load(file=f_name)

            g.est = colMeans(NQR_SME_res$g_trace) 
            Knots = NQR_SME_res$Knots
            mux.est=mean(NQR_SME_res$mux_trace)
            sigma2_22.est=mean(NQR_SME_res$sigma2_22_trace)
            sigma2_xx.est=mean(NQR_SME_res$sigma2_xx_trace)
            sigma2_22_save[sim_idx] = sigma2_22.est
            sigma2_xx_save[sim_idx] = sigma2_xx.est
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_SME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_SME_res$g_trace)
            g_var_save[sim_idx] = median(colVars(NQR_SME_res$g_trace))
          }
          
          if(to_see == '3ME'){
            f_name = sprintf('../debugging/NQR_data%s_%s3ME_%s_%s.RData',data.type,is.t,p0,sim_idx)
            if(sigma2_22!=sigma2_11){
              f_name = sprintf('../debugging/NQR_data%s_%s3ME_sig1%s_sig2%s_sig3%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,sigma2_33,p0,sim_idx)
            }
            load(file=f_name)
            
            g.est = colMeans(NQR_3ME_res$g_trace) 
            Knots = NQR_3ME_res$Knots
            mux.est=mean(NQR_3ME_res$mux_trace)
            alpha_save[sim_idx,] = colMeans(NQR_3ME_res$alpha_trace)
            beta_save[sim_idx,] = colMeans(NQR_3ME_res$beta_trace)

            sigma2_11_save[sim_idx] = mean(NQR_3ME_res$sigma2_11_trace)
            sigma2_22_save[sim_idx] = mean(NQR_3ME_res$sigma2_22_trace)
            sigma2_33_save[sim_idx] = mean(NQR_3ME_res$sigma2_33_trace)
            sigma2_xx_save[sim_idx] = mean(NQR_3ME_res$sigma2_xx_trace)
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_3ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_3ME_res$g_trace)
            g_var_save[sim_idx] = median(colVars(NQR_3ME_res$g_trace))
          }
          
          if(to_see == '4ME'){
            f_name = sprintf('../debugging/NQR_data%s_%s4ME_%s_%s.RData',data.type,is.t,p0,sim_idx)
            if(sigma2_22!=sigma2_11){
              f_name = sprintf('../debugging/NQR_data%s_%s4ME_sig1%s_sig2%s_sig3%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,sigma2_33,p0,sim_idx)
            }
            load(file=f_name)
            
            g.est = colMeans(NQR_4ME_res$g_trace) 
            Knots = NQR_4ME_res$Knots
            mux.est=mean(NQR_4ME_res$mux_trace)
            alpha_save[sim_idx,] = colMeans(NQR_4ME_res$alpha_trace)
            beta_save[sim_idx,] = colMeans(NQR_4ME_res$beta_trace)
            theta_save[sim_idx,] = colMeans(NQR_4ME_res$theta_trace)
            
            sigma2_11_save[sim_idx] = mean(NQR_4ME_res$sigma2_11_trace)
            sigma2_22_save[sim_idx] = mean(NQR_4ME_res$sigma2_22_trace)
            sigma2_33_save[sim_idx] = mean(NQR_4ME_res$sigma2_33_trace)
            sigma2_44_save[sim_idx] = mean(NQR_4ME_res$sigma2_44_trace)
            sigma2_xx_save[sim_idx] = mean(NQR_4ME_res$sigma2_xx_trace)
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_4ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_4ME_res$g_trace)
            g_var_save[sim_idx] = median(colVars(NQR_4ME_res$g_trace))
          }
          # {
          #   par(mfrow=c(1,1))
          #   is.t=''
          #   Xmul = 10
          #   X_shit = make_X_shit(0.5,Xmul)
          #   # # An Introduction to Kernel and Nearest-Neighbor Nonparametric Regression #############################################################################################
          #   set.seed(sim_idx)
          #   x1 = runif(n,0,1)
          #   inp.sd = 0.5
          #   y = x1*sin(2.5*pi*x1) + rnorm(n,0,inp.sd)
          #   X=cbind(1,x1*Xmul-X_shit)
          #   Xrange = seq(-5, 5,length.out = 100)
          #   x1range = seq(0, 1,length.out = 100)
          #   plot(X[,2],y)
          #   for(p0.tmp in p0_list){
          #     y.p0 = x1range*sin(2.5*pi*x1range) + qnorm(p0.tmp,0,inp.sd)
          #     points(Xrange,y.p0,col=2,lwd=2,type='l')
          #   }
          #   points(NQR_SME_res$Knots,colMeans(NQR_SME_res$g_trace),type = 'l',col=2,lwd=3)
          # }
          # plot(X[,2],NQR_SME_res$X_trace)
          # plot(NQR_SME_res$X_trace,y)
          # ts.plot(NQR_SME_res$g_trace[,1])
          # ts.plot(NQR_SME_res$sigma2_22_trace)
          # ts.plot(NQR_SME_res$sigma2_xx_trace)
          
          
          if(to_see == 'cvar'){
            f_name = sprintf('../debugging/NQR_data%s_%swME_cvar_%s_%s.RData',data.type,is.t,p0,sim_idx)
            load(file=f_name)
            
            g.est = colMeans(NQR_cvar_res$g_trace) 
            Knots = NQR_cvar_res$Knots
            
            alpha_save[sim_idx,]=NQR_cvar_res$alpha_trace
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_cvar_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_cvar_res$g_trace)
          }
          
          if(to_see == 'woME'){
            load(file=sprintf('../debugging/NQR_data%s_%swoME_%s_%s.RData',data.type,is.t,p0,sim_idx))
            
            
            g.est = colMeans(NQR_wo_ME_res$g_trace) 
            Knots = NQR_wo_ME_res$Knots
            # stopifnot( sum(Knots!=NQR_wo_ME_res$Knots)==0 )
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_wo_ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_wo_ME_res$g_trace)
          }
          
          
          if(to_see == 'W2'){
            f_name = sprintf('../debugging/NQR_data%s_%sW2_%s_%s.RData',data.type,is.t,p0,sim_idx)
            if(sigma2_22!=sigma2_11){
              f_name = sprintf('../debugging/NQR_data%s_%sW2_sig1%s_sig2%s_%s_%s.RData',data.type,is.t,sigma2_11,sigma2_22,p0,sim_idx)
            }
            load(file = f_name)
            
            g.est = colMeans(NQR_W2_ME_res$g_trace) 
            Knots = NQR_W2_ME_res$Knots
            # stopifnot( sum(Knots!=NQR_W2_ME_res$Knots)==0 )
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_W2_ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_W2_ME_res$g_trace)
            g_var_save[sim_idx] = median(colVars(NQR_W2_ME_res$g_trace))
          }
          
          g.tau=as.numeric(g.est);knots = Knots
          fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
          if(is.plot){
            plot(integral_range,(est.g.func(integral_range)-mean(true.g.func(integral_range)))/sd(true.g.func(integral_range)))
            points(integral_range,(true.g.func(integral_range)-mean(true.g.func(integral_range)))/sd(true.g.func(integral_range)),type='l')
            points(integral_range,tmp.int(integral_range),type='l')
          }
          
          ISE = integrate(ise_func,-5,5)$value
          ISE_save[sim_idx] = ISE          
          
          if(cal_MSE){
            set.seed(sim_idx)
            x1 = runif(n,0,1)
            if(data.type==2){
              x1 = (seq(1,n)-1)/n
            }
            x1 = x1*10-5
            MSE = mean(ise_func(x1))
            (ise_func(x1[1]))
            (true.g.func(x1[1])-est.g.func(x1[1]))**2
            MSE_save[sim_idx] = MSE
          }
          
          
          if(is.plot){
            par(mfrow=c(2,2))
            plot(X.est,X[,2]);abline(0,1)
            plot(X.est,W1);abline(alpha.est)
            plot(X.est,y);points(tau.i,g.est,type='l')
            par(mfrow=c(1,1))
          }
        },
        error = function(e) print(e))
    }
    toc()
    
    colMeans(alpha_save)
    mean(sigma2_11_save) # 0.2, true value is 1
    mean(sigma2_22_save) # 2.0, which is true value
    mean(sigma2_xx_save) 
    summary(accept_g_save)
    
    # Calculate the summary statistics from iterated result------------------------------------------------------------------------------------------
    nconverge_idx=which(accept_g_save<0.1)
    # if((to_see=='SME')|(to_see=='cvar')){
    # if((to_see=='cvar')){
    #   nconverge_idx=which(accept_g_save<0.01)
    # }
    if(length(nconverge_idx)==0){
      nconverge_idx = nmax+1
    }
    if(to_see=='W2'){
      nconverge_idx = nmax+1
    }
    
    
    inp.sub = NA
    if(to_see=='wME'){
      inp.sub = sprintf('alpha: %s,%s, sig1^2: %s, sig2^2: %s',round(colMeans(alpha_save),3)[1],round(colMeans(alpha_save),3)[2],round(mean(sigma2_11_save),3),round(mean(sigma2_22_save),3))
    }
    if((to_see=='W2')|(to_see=='SME')){
      inp.sub = sprintf('g accept: %s, g_var: %s,sig2: %s',round(median(accept_g_save,na.rm = T),3),round(mean(g_var_save,na.rm = T),6),round(mean(sigma2_22_save,na.rm = T),3))
    }
    if(to_see=='3ME'){
      inp.sub = sprintf('alpha: %s,%s,beta: %s,%s, sig1^2: %s, sig2^2: %s, sig3^2: %s',
                        round(colMeans(alpha_save,na.rm = T),3)[1],
                        round(colMeans(alpha_save,na.rm = T),3)[2],
                        round(colMeans(beta_save,na.rm = T),3)[1],
                        round(colMeans(beta_save,na.rm = T),3)[2],
                        round(mean(sigma2_11_save,na.rm = T),3),round(mean(sigma2_22_save,na.rm = T),3),round(mean(sigma2_33_save,na.rm = T),3))
    }
    
    if(to_see=='4ME'){
      inp.sub = sprintf('alpha: %s,%s,beta: %s,%s, theta: %s,%s, sig1^2: %s, sig2^2: %s, sig3^2: %s, sig4^2: %s',
                        round(colMeans(alpha_save,na.rm = T),3)[1],
                        round(colMeans(alpha_save,na.rm = T),3)[2],
                        round(colMeans(beta_save,na.rm = T),3)[1],
                        round(colMeans(beta_save,na.rm = T),3)[2],
                        round(colMeans(theta_save,na.rm = T),3)[1],
                        round(colMeans(theta_save,na.rm = T),3)[2],
                        round(mean(sigma2_11_save,na.rm = T),3),round(mean(sigma2_22_save,na.rm = T),3),
                        round(mean(sigma2_33_save,na.rm = T),3),round(mean(sigma2_44_save,na.rm = T),3))
    }
    
    m.boxplot(g_save[-nconverge_idx,],p0,type=paste(to_see,'HPD:',round(mean(HPD_save[-nconverge_idx],na.rm = T),3),',MISE:',round(mean(ISE_save[-nconverge_idx],na.rm = T),3)),data.type = data.type,inp.sub = inp.sub)  
    g_save_list[[to_see]]=g_save[-nconverge_idx,]
    
    summary_list[[as.character(p0)]][[to_see]][['ISE']][['mean']] = mean(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['ISE']][['sd']] = sd(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['mean']] = mean(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['sd']] = sd(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['mean']] = mean(MSE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['sd']] = sd(MSE_save[-nconverge_idx],na.rm = T)
  }
  
  # m.boxplo.v2(g_save_list,p0,data.type)
}
par(mfrow=c(1,1))
cat(sum(is.na(g_save[,1]))/nmax*100,'% is not yout done\n')
