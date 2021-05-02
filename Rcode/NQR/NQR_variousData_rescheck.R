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

HPD_ratio = function(g_trace,lb=0.05,ub=0.95){
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
  return (true.g.func(x)-est.g.func(x))**2
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
mux_save=rep(NA,nmax)
sigma2_11_save=rep(NA,nmax)
sigma2_22_save=rep(NA,nmax)
sigma2_xx_save=rep(NA,nmax)


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

summary_list = list()
par(mfrow=c(length(p0_list),1))
par(mfrow=c(length(to_see_list),1))

for(p0 in p0_list){
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
            load(file=sprintf('../debugging/NQR_data%s_%swME_%s_%s.RData',data.type,is.t,p0,sim_idx))
            g.est=colMeans(NQR_res$g_trace)
            Knots = NQR_res$Knots
            alpha.est=colMeans(NQR_res$alpha_trace)
            X.est=colMeans(NQR_res$X_trace)
            mux.est=mean(NQR_res$mux_trace)
            sigma2_11.est=mean(NQR_res$sigma2_11_trace)
            sigma2_22.est=mean(NQR_res$sigma2_22_trace)
            sigma2_xx.est=mean(NQR_res$sigma2_xx_trace)
            
            g_save[sim_idx,]=g.est
            accept_g_save[sim_idx]=NQR_res$g_accept_ratio
            HPD_save[sim_idx] = HPD_ratio(NQR_res$g_trace)
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
            load(file=sprintf('../debugging/NQR_data%s_%sW2_%s_%s.RData',data.type,is.t,p0,sim_idx))
            
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
    nconverge_idx=which(accept_g_save<0.1)
    if(length(nconverge_idx)==0){
      nconverge_idx = nmax+1
    }
    
    m.boxplot(g_save[-nconverge_idx,],p0,type=paste(to_see,'HPD:',round(mean(HPD_save[-nconverge_idx],na.rm = T),3),',MISE:',round(mean(ISE_save[-nconverge_idx],na.rm = T),3)),data.type = data.type)

    summary_list[[as.character(p0)]][[to_see]][['ISE']][['mean']] = mean(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['ISE']][['sd']] = sd(ISE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['mean']] = mean(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['HPD']][['sd']] = sd(HPD_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['mean']] = mean(MSE_save[-nconverge_idx],na.rm = T)
    summary_list[[as.character(p0)]][[to_see]][['MSE']][['sd']] = sd(MSE_save[-nconverge_idx],na.rm = T)
  }
}
par(mfrow=c(1,1))
cat(sum(is.na(g_save[,1]))/nmax*100,'% is not yout done\n')


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
  inp.sd=1
  y=2+sin(x1i)+rnorm(n,0,inp.sd)
  
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


par(mfrow=c(1,1))
NQR_res$g_accept_ratio
NQR_res$x_accept_ratio
ts.plot(NQR_res$g_trace[,1])
acf(NQR_res$g_trace[,1])
ts.plot(NQR_res$X_trace[,1])
acf(NQR_res$X_trace[,1])
NQR_res$l_accept_ratio
mean(NQR_res$lambda_trace)
mean(NQR_res$mux_trace)
mean(NQR_res$sigma2_11_trace)
mean(NQR_res$sigma2_22_trace)
mean(NQR_res$sigma2_xx_trace)
X.est = colMeans(NQR_res$X_trace)
colMeans(NQR_res$alpha_trace)
alpha
set.seed(sim_idx)
x1 = runif(n,-4,2.5)
plot(X.est,x1)




HPD_ratio(NQR_res$g_trace)

boxplot(g_save)
points(seq(1,30),(g_quantile)[1,],type='l')
points(seq(1,30),(g_quantile)[2,],type='l')
points(seq(1,30),y.p0,type='l')
