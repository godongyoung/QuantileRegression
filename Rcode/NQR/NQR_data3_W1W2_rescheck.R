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

m.boxplot=function(save_data,p0,type,data.type){
  valid_cnt = sum(!(is.na(save_data[,1])))
  
  colnames(save_data)=Knots
  ds=cbind(rep(Knots,each=dim(save_data)[1]),as.numeric(save_data))
  colnames(ds)=c('Knots','value')
  
  if(data.type==1){inp.ylim = c(-5,5)}
  if(data.type==2){inp.ylim = c(-3,3)}
  if(data.type==3){inp.ylim = c(-1.5,2)}
  
  
  boxplot(value~Knots, data=ds,ylim=inp.ylim,main=sprintf('%s\'s Box plot of %s with Knots %s for %s simulation',type,p0,inp.N.Knots,valid_cnt),names = round(Knots,1)) 
  for(tmp.p0 in p0_list){
    y.p0 = gen_y.p0(data.type,tmp.p0)
    points(seq(1:length(Knots)),y.p0,type='l',lwd=2,col=3)
  }
  y.p0 = gen_y.p0(data.type,p0)
  points(seq(1:length(Knots)),y.p0,type='l',lwd=3,col=2)
}


m.boxplo.v2t=function(save_data,p0,type,data.type){
  valid_cnt = sum(!(is.na(save_data[,1])))
  
  colnames(save_data)=Knots
  ds=cbind(rep(Knots,each=dim(save_data)[1]),as.numeric(save_data))
  colnames(ds)=c('Knots','value')
  
  if(data.type==1){inp.ylim = c(-5,5)}
  if(data.type==2){inp.ylim = c(-3,3)}
  if(data.type==3){inp.ylim = c(-1.5,2)}
  
  boxplot(value~Knots, data=ds,ylim=inp.ylim,names = round(Knots,1)) 
  plot(NA,NA,xlim=c(-5,5),ylim=inp.ylim)
  for(tmp.p0 in p0_list){
    y.p0 = gen_y.p0(data.type,tmp.p0)
    points(Knots,y.p0,type='l',lwd=2,lty='dotted')
  }
  y.p0 = gen_y.p0(data.type,p0)
  points(Knots,y.p0,type='l',lwd=3,col=1)
  
  for(each in names(g_save_list)){
    g_mean = colMeans(g_save_list[[each]])
    points(Knots,g_mean,type='l',col=2)
    g_quantile = apply(g_save_list[[each]],MARGIN = 2,FUN = function(x) quantile(x ,probs = c(0.025,0.975)))
    
    points(Knots,g_quantile[1,],type='l',col=2,lty='dotted')
    points(Knots,g_quantile[2,],type='l',col=2,lty='dotted')
  }
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
  return ((true.g.func(x)-est.g.func(x))**2)
}

# g.tau=as.numeric(g.est);knots = Knots
# fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
# y.est=predict(fit.ns, data.frame(knots=integral_range))
est.g.func = function(inp.x){
  predict(fit.ns, data.frame(knots=inp.x))
}
# save_data = g_save;type='wME';data.type=1

library(ggplot2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

m.boxplo.v2=function(save_data_list,p0,data.type,save=F){
  if(W1.V2&W2.V2){
    f_name = sprintf('../Figure/data%s_W1W2_%s.png',data.type,p0)
  }
  else{
    f_name = sprintf('../Figure/data%s_W1%sW2%s_ver%s_%s.png',data.type,W1.V2,W2.V2,inp.version,p0)
  }
  
  # tiff(fname, units="in", width=6*1.5, height=4*1.5, res=600)
  if(save){
    png( f_name, width = 6*1.5, height=4*1.5, units = "in", res = 600, pointsize = 13)  
  }
  par(mfrow = c(1,1), mar = c(4,4,1,4))
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
  if(save){
    dev.off()
  }
  
}


# Define True parameter--------------------------------------------------------------------------------
n=1000
alpha=c(4,3)

sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
nmax=500

is.plot=F
to_see_list = c('wME','woME','W2')
to_see = to_see_list[1]
inp.N.Knots = 30
inp.mul = 10
N = inp.N.Knots


X_demean = T
Xmul = 10
X_shit = make_X_shit(0.5,Xmul)
integral_range = seq(-5,5,length.out = 1000)
cal_MSE = T

# Simulation check --------------------------------------------------------------------------------
alpha_save = matrix(NA,ncol=2,nrow=nmax)
X_save=matrix(NA,ncol=n,nrow=nmax)
X_save2=matrix(NA,ncol=n,nrow=nmax)
mux_save=rep(NA,nmax)
sigma2_11_save=rep(NA,nmax)
sigma2_22_save=rep(NA,nmax)
sigma2_xx_save=rep(NA,nmax)


n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)
p0=0.9
sim_idx=1

data.type = 3
is.t=''
inp.sd = 1
if(data.type==2){
  inp.sd = 0.07
}
if(data.type==3){
  inp.sd = 0.5  
}

# W1.W2.case = 3
# if(W1.W2.case==1){W1.V2=T;W2.V2=F}
# if(W1.W2.case==2){W1.V2=F;W2.V2=T}
# if(W1.W2.case==3){W1.V2=T;W2.V2=T}
W1.V2=T;W2.V2=F
inp.version = 1

summary_list = list()
par(mfrow=c(length(p0_list),1))
par(mfrow=c(length(to_see_list),1))
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
    
    summary_list[[as.character(p0)]][[to_see]][['ISE']]=list()
    summary_list[[as.character(p0)]][[to_see]][['HPD']]=list()
    summary_list[[as.character(p0)]][[to_see]][['MSE']]=list()
    
    # Load each iteration result------------------------------------------------------------------------------------------
    tic()
    for(sim_idx in 1:nmax){  
      tryCatch(
        {
          if(to_see == 'wME'){
            #####
            # set.seed(sim_idx)
            # x1 = runif(n,0,1)
            # X=cbind(1,x1*Xmul-X_shit)
            #####
            if(W1.V2&W2.V2){
              f_name = sprintf('../debugging/NQR_data%s_W1W2_ver%s_%swME_%s_%s.RData',data.type,inp.version,is.t,p0,sim_idx)
            }
            else{
              f_name = sprintf('../debugging/NQR_data%s_W1%sW2%s_ver%s_%swME_%s_%s.RData',data.type,W1.V2,W2.V2,inp.version,is.t,p0,sim_idx)
            }
            
            load(file=f_name)
            g.est=colMeans(NQR_res$g_trace)
            Knots = NQR_res$Knots
            # alpha.est=colMeans(NQR_res$alpha_trace)
            # X.est=colMeans(NQR_res$X_trace)
            # mux.est=mean(NQR_res$mux_trace)
            # alpha_save[sim_idx,] = alpha.est
            sigma2_11.est=mean(NQR_res$sigma2_11_trace)
            sigma2_22.est=mean(NQR_res$sigma2_22_trace)
            sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
            #####
            sigma2_11_save[sim_idx]=sigma2_11.est
            sigma2_22_save[sim_idx]=sigma2_22.est
            
            set.seed(sim_idx)
            x1 = runif(n,0,1)
            X_save2[sim_idx,] = x1*10-5
            
            {if(is.null(dim(NQR_res$X_trace)[2])){
              X_save[sim_idx,] = NQR_res$X_trace
            }
            else{
              X_save[sim_idx,] = colMeans(NQR_res$X_trace)
            }}
            
            #####
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_res$g_trace)
            
            
          }
          
          if(to_see == 'woME'){
              f_name = sprintf('../debugging/NQR_data%s_%swoME_%s_%s.RData',data.type,is.t,p0,sim_idx)  
            load(file = f_name)
            
            g.est = colMeans(NQR_wo_ME_res$g_trace) 
            Knots = NQR_wo_ME_res$Knots
            # stopifnot( sum(Knots!=NQR_wo_ME_res$Knots)==0 )
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_wo_ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_wo_ME_res$g_trace)
          }
          
          
          if(to_see == 'W2'){
            if(W1.V2&W2.V2){
              f_name = sprintf('../debugging/NQR_data%s_W1W2_ver%s_%sW2_%s_%s.RData',data.type,inp.version,is.t,p0,sim_idx)   
            }
            if(W2.V2==FALSE){
              f_name = sprintf('../debugging/NQR_data%s_%sW2_%s_%s.RData',data.type,is.t,p0,sim_idx)   
            }
            else{
              f_name = sprintf('../debugging/NQR_data%s_W1%sW2%s_ver%s_%sW2_%s_%s.RData',data.type,W1.V2,W2.V2,inp.version,is.t,p0,sim_idx)   
            }
            load(file = f_name)
            
            g.est = colMeans(NQR_W2_ME_res$g_trace) 
            Knots = NQR_W2_ME_res$Knots
            # stopifnot( sum(Knots!=NQR_W2_ME_res$Knots)==0 )
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_W2_ME_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_W2_ME_res$g_trace)
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
        error = function(e) cat(sim_idx,'of',p0,'is not done yet. \n'))
    }
    toc()

    # Calculate the summary statistics from iterated result------------------------------------------------------------------------------------------
    nconverge_idx=which(accept_g_save<0.01)
    # nconverge_idx = nmax+1
    if(length(nconverge_idx)==0){
      nconverge_idx = nmax+1
    }
    if((to_see=='W2')){
      nconverge_idx = nmax+1
    }
    m.boxplot(g_save[-nconverge_idx,],p0,type=paste(to_see,'HPD:',round(mean(HPD_save[-nconverge_idx],na.rm = T),3),',MISE:',round(mean(ISE_save[-nconverge_idx],na.rm = T),3)),data.type = data.type)
    g_save_list[[to_see]]=g_save[-nconverge_idx,]
    
    summary_list[[as.character(p0)]][[to_see]][['ISE']][['mean']] = mean(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['ISE']][['sd']] = sd(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['mean']] = mean(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['sd']] = sd(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['mean']] = mean(MSE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['sd']] = sd(MSE_save[-nconverge_idx],na.rm = T)
  }
  m.boxplo.v2(g_save_list,p0,data.type,save = T)
}
par(mfrow=c(1,1))
cat(sum(is.na(g_save[,1]))/nmax*100,'% is not yout done\n')


# Check X.true & X.est for 'wME'------------------------------------------------------------------------------------------
sim_idx = 1
X.true = X_save2[sim_idx,]
X.est = X_save[sim_idx,]
plot(X.true,X.est,main=sprintf('X.true Vs X.est Ver%s',inp.version));abline(0,1)


# Check plot of X & W1 or X & W2------------------------------------------------------------------------------------------

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




set.seed(1)
x1 = runif(n,0,1)
X=cbind(1,x1*Xmul-X_shit)
inp.version = 3
W_list = make_data(X,W1.V2,W2.V2,version = inp.version)
plot(X[,2],W_list$W1,main = sprintf('X = a0 + a1*W2 + e in Ver%s',inp.version));abline(lm(W_list$W1~X[,2])$coeff)
plot(X[,2],W_list$W2,main = sprintf('X = W1 + e in Ver%s',inp.version));abline(0,1)

