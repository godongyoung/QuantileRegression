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


save_data = g_save[-nconverge_idx,]
type='wo'
m.boxplot=function(save_data,p0,type){
  y.p0=2+sin(Knots)+qnorm(p = p0,mean = 0,sd = 0.1)
  
  colnames(save_data)=Knots
  ds=cbind(rep(Knots,each=dim(save_data)[1]),as.numeric(save_data))
  colnames(ds)=c('Knots','value')
  boxplot(value~Knots, data=ds,ylim=c(0,4),main=sprintf('%s\'s Box plot of %s',type,p0),names = round(Knots,1)) 
  points(seq(1:length(Knots)),y.p0,type='l',lwd=3,col=2)
}


N=30

# Define True parameter--------------------------------------------------------------------------------
n=1000
alpha=c(4,3)

Mu_x=5
sigma2=1
sigma2_xx=3
sigma2_11=1
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
sim_idx=4
p0=0.25
nmax=500
is.plot=F
is.cal_quantile=F
if.short = T
if.NQR_wo_ME = T


# Simulation check --------------------------------------------------------------------------------


alpha_save = matrix(NA,ncol=2,nrow=nmax)
g_save=matrix(NA,ncol=30,nrow=nmax)
X_save=matrix(NA,ncol=n,nrow=nmax)
mux_save=rep(NA,nmax)
sigma2_11_save=rep(NA,nmax)
sigma2_22_save=rep(NA,nmax)
sigma2_xx_save=rep(NA,nmax)
accept_g_save = rep(NA,nmax)


accept_g_woME_save = rep(NA,nmax)
g_woME_save=matrix(NA,ncol=30,nrow=nmax)

load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
Knots=NQR_res$Knots
Knots=seq(from = 0.1,to = 2*Mu_x,length.out = 30)


#Data type test start #########################################################################################
set.seed(20210401)
n=1e4
x1i=runif(n=n,min=0,max=2*Mu_x)
x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sqrt(sigma2_xx))
X=cbind(1,x1i)
X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)

inp.sd=1
y=2+sin(x1i)+rnorm(n,0,inp.sd)
plot(x1i,y,xlim=c(0,10))
p0=0.25
for(p0 in c(0.1,0.25,0.5,0.75,0.9)){
  y.p0=2+sin(Knots)+qnorm(p = p0,mean = 0,sd = inp.sd)
  points(Knots,y.p0,type='l',lwd=3,col=2)
}



inp.sd=1
y=2+3*log(x1i)+rnorm(n,0,x1i*inp.sd)
plot(x1i,y,xlim=c(0,10))
p0=0.25
for(p0 in c(0.1,0.25,0.5,0.75,0.9)){
  y.p0=2+3*log(Knots)+Knots*qnorm(p = p0,mean = 0,sd = inp.sd)
  points(Knots,y.p0,type='l',lwd=3,col=2)
  
  points(8,quantile(y[(7<x1i)&(x1i<9)],p0),col=3,lwd=5)
}

#Data type test end#########################################################################################
p0_list=c(0.1,0.25,0.5,0.75,0.9)
length(p0_list)

par(mfrow=c(length(p0_list),1))
for(p0 in p0_list){
  tic()
  for(sim_idx in 1:nmax){
    tryCatch(
      {
        # Load MCMC result--------------------------------------------------------------------------------
        {if(!if.short){
          load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx)) # this is done by inp_version=1
          alpha.est=colMeans(NQR_res$alpha_trace)
          g.est=colMeans(NQR_res$g_trace)
          X.est=colMeans(NQR_res$X_trace)
          mux.est=mean(NQR_res$mux_trace)
          sigma2_11.est=mean(NQR_res$sigma2_11_trace)
          sigma2_22.est=mean(NQR_res$sigma2_22_trace)
          sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
          Knots = NQR_res$Knots
        }
          if(if.short){
            load(file=sprintf('../debugging/NQR_short_%s_%s.RData',p0,sim_idx)) 
            alpha.est = NQR_res_short$alpha.est
            g.est = NQR_res_short$g.est
            X.est = NQR_res_short$X.est
            mux.est = NQR_res_short$mux.est
            sigma2_11.est = NQR_res_short$sigma2_11.est
            sigma2_22.est = NQR_res_short$sigma2_22.est
            sigma2_xx.est = NQR_res_short$sigma2_xx.est
            Knots = NQR_res_short$Knots
          }}
        
        if(if.NQR_wo_ME){
          load(file=sprintf('../debugging/NQR_woME_short_%s_%s.RData',p0,sim_idx)) 
          g.est_woME = NQR_res_woME_short$g.est 
          
          stopifnot( sum(Knots!=NQR_res_woME_short$Knots)==0 )
        }
        
        alpha_save[sim_idx,]=alpha.est
        g_save[sim_idx,]=g.est
        X_save[sim_idx,]=X.est
        mux_save[sim_idx]=mux.est
        sigma2_11_save[sim_idx]=sigma2_11.est
        sigma2_22_save[sim_idx]=sigma2_22.est
        sigma2_xx_save[sim_idx]=sigma2_xx.est
        accept_g_save[sim_idx]=NQR_res$g_accept_ratio
        
        g_woME_save[sim_idx,]=g.est_woME
        accept_g_woME_save[sim_idx]=NQR_res_woME_short$g_accept_ratio
        
        if(is.plot){
          par(mfrow=c(2,2))
          plot(X.est,X[,2]);abline(0,1)
          plot(X.est,W1);abline(alpha.est)
          plot(X.est,y);points(tau.i,g.est,type='l')
          par(mfrow=c(1,1))
          
          
          # mean(NQR_res$sigma2_11_trace)
          # mean(NQR_res$sigma2_22_trace)
          # mean(NQR_res$sigma2_xx_trace)
          # 
          # ts.plot(NQR_res$sigma2_11_trace)
          # ts.plot(NQR_res$sigma2_22_trace)
          # ts.plot(NQR_res$sigma2_xx_trace)
          # 
        }
      },
      error = function(e) cat(sim_idx,'is not done yet. \n'))
  }
  toc()
  
  
  # nconverge_idx=which(accept_g_woME_save<0.1)
  # {if(length(nconverge_idx)==0){m.boxplot(g_woME_save,p0,type='woME')}
  # else {m.boxplot(g_woME_save[-nconverge_idx,],p0,type='woME')}}
  
  
  nconverge_idx=which(accept_g_save<0.1)
  {if(length(nconverge_idx)==0){m.boxplot(g_save,p0,type='wME')}
    else {m.boxplot(g_save[-nconverge_idx,],p0,type='wME')}}
  
}
par(mfrow=c(1,1))
cat(sum(is.na(g_woME_save[,1]))/nmax,'% is not yout done\n')


g_save

# plot(X[,2],y)
# points(X[,2],y.est,type='l')
# boxplot(y,y.est);abline(h=mean(y.quantile_save[-nconverge_idx]))
# boxplot(y.est.quantile_save[-nconverge_idx]);abline(h=mean(y.quantile_save[-nconverge_idx]))
# boxplot(y.est.quantile_save2[-nconverge_idx]);abline(h=mean(y.quantile_save[-nconverge_idx]))
# hist(y.est.quantile_save[-nconverge_idx],nclass=100);abline(v=mean(y.quantile_save[-nconverge_idx]))
# boxplot(bias.quantile_save[-nconverge_idx]);abline(h=0)
# Make larger data for Ground Truth--------------------------------------------------------------------------------
set.seed(sim_idx)
n=1e4
# x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sigma2_xx)
x1i=runif(n=n,min=0,max=2*Mu_x)
x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
X=cbind(1,x1i)
X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
y=2+sin(x1i)+rnorm(n,0,0.1)

#generate W1,W2
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

W1=X%*%alpha+delta1
W2=X[,2]+delta2


# Check result --------------------------------------------------------------------------------
nconverge_idx=which(accept_g_save<0.1)
hist(accept_g_save[-nconverge_idx])


par(mfrow=c(4,2))
plot(X[,2],y,main='X vs Y, with g.est');points(tau.i,colMedians(g_save[-nconverge_idx,],na.rm = T),type='l',col=2,lwd=3)
plot(X[,2],W1,main="X vs W1 with alpha.est");abline(colMedians(alpha_save[-nconverge_idx,],na.rm=T),col=2,lwd=3)

hist(alpha_save[-nconverge_idx,1],nclass=100);abline(v=alpha[1],col=2,lwd=3)
hist(alpha_save[-nconverge_idx,2],nclass=100);abline(v=alpha[2],col=2,lwd=3)

hist(mux_save[-nconverge_idx],nclass=100);abline(v=Mu_x,col=2,lwd=3)
hist(sigma2_11_save[-nconverge_idx],nclass=100);abline(v=sigma2_11,col=2,lwd=3)
hist(sigma2_22_save[-nconverge_idx],nclass=100);abline(v=sigma2_22,col=2,lwd=3)
hist(sigma2_xx_save[-nconverge_idx],nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
par(mfrow=c(1,1))


# mean(mux_save,na.rm = T)
# median(sigma2_11_save,na.rm = T)
# hist(sigma2_11_save,nclass=100)
# mean(sigma2_22_save,na.rm = T)
# mean(sigma2_xx_save,na.rm = T)

# Debugging weired result --------------------------------------------------------------------------------
sim_idx=254
sim_idx=3
condition = which(abs(alpha_save[,2])>10)
condition = which(abs(alpha_save[,2])<10)

for(idx in 1:length(condition)){
  sim_idx=condition[idx]
  # sim_idx=idx
  if(sim_idx %in% nconverge_idx){next}
  
  load(file=sprintf('../debugging/NQR_%s_%s.RData',p0,sim_idx))
  alpha.est=colMeans(NQR_res$alpha_trace)
  g.est=colMeans(NQR_res$g_trace)
  X.est=colMeans(NQR_res$X_trace)
  mux.est=mean(NQR_res$mux_trace)
  sigma2_11.est=mean(NQR_res$sigma2_11_trace)
  sigma2_22.est=mean(NQR_res$sigma2_22_trace)
  sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
  
  
  set.seed(sim_idx)
  n=1000
  # x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sigma2_xx)
  x1i=runif(n=n,min=0,max=2*Mu_x)
  x1i=rnorm(n,Mu_x,sqrt(sigma2_xx))
  X=cbind(1,x1i)
  X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
  y=2+sin(x1i)+rnorm(n,0,0.1)
  
  #generate W1,W2
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  # par(mfrow=c(4,2))
  # plot(X.est,X[,2],main='X.est Vs X with y=x line');abline(0,1)
  # plot(X.est,W1,main='X.est Vs W1 with alpha.est');abline(alpha.est)
  # plot(X.est,W2,main='X.est Vs W2 with y=x line');abline(0,1)
  # plot(X.est,y,main='X.est Vs Y with beta.est');points(tau.i,g.est,type='l',col=2,lwd=3)
  # 
  # hist(NQR_res$sigma2_11_trace,nclass=100);abline(v=sigma2_11,col=2,lwd=3)
  # hist(NQR_res$sigma2_22_trace,nclass=100);abline(v=sigma2_22,col=2,lwd=3)
  # hist(NQR_res$sigma2_xx_trace,nclass=100);abline(v=sigma2_xx,col=2,lwd=3)
  # 
  # ts.plot(NQR_res$sigma2_22_trace)
  # par(mfrow=c(1,1))
  
  NQR_res$g_accept_ratio
  NQR_res$l_accept_ratio
  NQR_res$x_accept_ratio
  sigma2_11.est
  sigma2_22.est
  sigma2_xx.est
  
  par(mfrow=c(3,2))
  ts.plot(NQR_res$g_trace[,1],main=sim_idx)
  ts.plot(NQR_res$alpha_trace[,1])
  ts.plot(NQR_res$X_trace[,1])
  ts.plot(NQR_res$sigma2_11_trace)
  ts.plot(NQR_res$sigma2_22_trace)
  ts.plot(NQR_res$sigma2_xx_trace)
  par(mfrow=c(1,1))
  
  accept_g_save[condition]
  
  median(accept_g_save,na.rm = T)

}
# sim_idx=condition[1]

