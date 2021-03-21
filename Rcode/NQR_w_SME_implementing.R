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

########
set.seed(20210317)
# Make data--------------------------------------------------------------------------------
p0=0.25
n=1000
inp.version=1
multiply_c=5 
N=30

Mu_x=5
sigma2_xx=1
sigma2_22=1

# x1i=rtruncnorm(n = n,a = 0,b = 2*Mu_x,mean=Mu_x,sd=sigma2_xx)
x1i=runif(n=n,min=0,max=2*Mu_x)
X=cbind(1,x1i)
y=2+sin(x1i)+rnorm(n,0,0.1)

#generate w1,w2
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

W2=X[,2]+delta2



#Start#################

# Function --------------------------------------------------------------------------------
qloss=function(u,p0){
  return (u*(p0-(u<0)))
}

smooth.y=function(knots,g.tau,xout,version=1){
  if(version==1){
    mspline=spline(x = knots,y = g.tau,xout = xout)
    
  }
  if(version==2){
    msmooth.spline=smooth.spline(x = knots,y = g.tau,cv = NA,lambda = lambda.t)
    mspline=predict(msmooth.spline,xout)
  }
  
  return(mspline$y)
}


log.likeli.g=function(g.t, lambda.t,X.t){
  y.est=smooth.y(tau.i,g.t,X.t[,2],version=inp.version)
  ui=y-y.est
  term1 = -sum(qloss(ui,p0))
  
  term2 = -0.5 * lambda.t * (g.t) %*% K %*% t(g.t)
  return(term1+term2)
}

log.likeli.l=function(lambda.f,lambda.b,g.t){
  term1 = (N-2)/2*log(lambda.f) 
  term2 = -0.5 * lambda.f * (g.t) %*% K %*% t(g.t)
  term3 = dgamma(x = lambda.f,shape = labmda.a, scale = labmda.b, log = T)
  term4 = dlnorm(x = lambda.b, meanlog = log(lambda.f), sdlog = sqrt(jump_lambda),log = T)  
  return (term1 + term2 + term3)
}

log.likeli.x=function(X.1t,lambda.t){ ### This can be changed with individual updates!
  y.est=smooth.y(tau.i,g.t,X.1t,version=inp.version)
  ui=y-y.est
  term1 = -sum(qloss(ui,p0))
  term3 = sum(dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T))
  term4 = sum(dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T))
  return (sum(term1 + term3 + term4))
}

# Set Defaults --------------------------------------------------------------------------------

niter    <- 30000*multiply_c
nburn    <- 5000*multiply_c
nthin    <- 5*multiply_c
nprint   <- 10000*multiply_c
nmcmc=(niter-nburn)/nthin

n = length(y)
inp.min = min(X[,2]); inp.max = max(X[,2])
# tau.i=seq(from = min(X[,2]),to = max(X[,2]),length.out = N)
tau.i=seq(from = inp.min,to = inp.max,length.out = N)

hi=diff(tau.i) # because we define tau.i as equally spaced, hi has same value for all i.
# Make Q matrix
Q=matrix(0,nrow=N,ncol=(N-2))
for(j in seq(2,N-2)){
  Q[(j-1),j] = 1/hi[(j-1)]
  Q[j,j] = -1/hi[(j-1)] - 1/hi[j]
  Q[(j+1),j] = 1/hi[j]
}
## This is suspectable!!!
Q[1,1] = -1/hi[(j-1)] - 1/hi[j]
Q[2,1] = 1/hi[j]

# Make R matrix
R=matrix(0,nrow=(N-2),ncol=(N-2))
for(i in seq(2,N-2)){
  R[i,i] = (hi[i-1] + hi[i])/3
  if(i<(N-2)){
    R[(i+1),i] = hi[i]/6
    R[i,(i+1)] = hi[i]/6
  }
}
## This is suspectable!!!
R[1,1] = (hi[i-1] + hi[i])/3
R[1,2] = hi[i]/6
R[2,1] = hi[i]/6

# Make K matrix
K=Q%*%solve(R)%*%t(Q)
stopifnot(N-2==qr(K)$rank)
ev=eigen(K)
mus=ev$values[1:(N-2)]

# initial values --------------------------------------------------------------------------------
alpha.t=rnorm(2,0,1)
sigma2_xx.t=10#runif(1)
sigma2_22.t=10#runif(1)
mux.t=rnorm(1,0,sqrt(100))
# X.1t=W2
X.1t=rtruncnorm(n = n,a = 0,b = 10,mean=5,sd=2)
X.t=cbind(1,X.1t)

BQR_res = mBayesQR(y,X.t,p0)
beta.est=colMeans(BQR_res$beta_trace)
g0=rep(0,N)
for(param.idx in 1:length(beta.est)){
  g0 = g0 + beta.est[param.idx]*tau.i^(param.idx-1) # for quadratic g0
}


msmooth.spline=smooth.spline(x = tau.i,y = g0,control.spar = list('maxit'=1,'trace'=F))
lambda0=msmooth.spline$lambda

g.t=matrix(g0,nrow = 1)
lambda.t=lambda0


# Check with NQR --------------------------------------------------------------------------------
plot(X[,2],y, main="True Data + BayesQR init")
X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
# for(p0 in c(.1,.25,.5,.75,.9)){
#   NQR_res=NQR(y,X,p0)
#   g.est=colMeans(NQR_res$g_trace)
#   lambda.est=mean(NQR_res$lambda_trace)
#   mspline=spline(x = NQR_res$Knots,y = g.est,xout = X_range)
#   
#   points(X_range,mspline$y,type = 'l')    
# }
points(tau.i,g0,type = 'l')

# Set Priors --------------------------------------------------------------------------------

#Gamma for lambda
labmda.b=0.1/lambda0
labmda.a=lambda0/labmda.b

aaa=0.1
#IG for sigma2_xx
prior_ax=aaa;prior_bx=aaa
#IG for sigma2_22
prior_a2=aaa;prior_b2=aaa

#N for mu_x
var_param=100
pr_mean_mux=0;pr_var_mux=var_param


# Make trace --------------------------------------------------------------------------------
g_trace=matrix(NA,ncol=N,nrow=nmcmc)
lambda_trace=rep(NA,nmcmc)
sigma2_xx_trace=rep(NA,nmcmc)
sigma2_22_trace=rep(NA,nmcmc)
mux_trace=rep(NA,nmcmc)
X_trace=matrix(NA,nrow=nmcmc,ncol=n)   

# Set jumping rules --------------------------------------------------------------------------------

pr.sigma.g=0.01
jump_g=pr.sigma.g^2*diag(N)
jump_lambda = 0.1
pr.sigma.x=0.08
jump_x=pr.sigma.x^2*diag(N)

accept_g = 0
accept_l = 0
accept_x = 0



##### Debugging part
# plot(X.t[,2],y)
# mspline=spline(x = tau.i,y = g.t,xout = tau.i)
# points(tau.i,mspline$y,type = 'l')
# 
# mspline=spline(x = tau.i,y = g.t,xout = X.t[,2])
# mspline=spline(x = tau.i,y = g.star,xout = X.t[,2])
# ui=y-mspline$y
# term1 = -sum(qloss(ui,p0));term1
# 
# X.1.star = rnorm(n = n,mean = X.1t,sd = pr.sigma.x)
# plot(X.1t,W2);abline(0,1)
# log.likeli.x=function(X.1t){ ### This can be changed with individual updates!
# mspline=spline(x = tau.i,y = g.t,xout = X.1t) # this is where X.t contributes
# ui=y-mspline$y
# term1 = -sum(qloss(ui,p0))
# term2 = dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t,sd = sqrt(sigma2_11.t),log = T)
# term3 = dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T)
# term4 = dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T)
#   return (sum(term1 + term2 + term3 + term4))
# }
sigma2_22.t=1

g.true=smooth.y(X[,2],y,tau.i,version=2)
ui=y-smooth.y(tau.i,g.t,X.t[,2],version=inp.version)
term1 = -sum(qloss(ui,p0));term1

ui=y-smooth.y(tau.i,g.true,X.t[,2],version=inp.version)
term1 = -sum(qloss(ui,p0));term1

plot(X[,2],y)
points(tau.i,smooth.y(tau.i,g.true,tau.i,version=inp.version),type = 'l')

###########Debugging part end

# iter start --------------------------------------------------------------------------------
tic()
for(iter in 1:niter){
  if(iter%%nprint==1){cat(iter,'th iter is going\n')}
  # Sample g-----------------------------------------------------------------
  g.star = rmvnorm(n = 1,mu = g.t,sigma = jump_g)
  log_accept.r = log.likeli.g(g.star,lambda.t,X.t = X.t) - log.likeli.g(g.t,lambda.t,X.t = X.t)
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    g.t = g.star
    accept_g = accept_g + 1
  }
  
  # Sample lambda-----------------------------------------------------------------
  # lamda.star = exp(rnorm(1,mean = log(lambda.t),sd = sqrt(jump_lambda)))
  lambda.star = rlnorm(1,meanlog = log(lambda.t),sdlog = sqrt(jump_lambda))
  log_accept.r = log.likeli.l(lambda.star,lambda.t,g.t) - log.likeli.l(lambda.t,lambda.star,g.t)  
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    lambda.t = lambda.star
    accept_l = accept_l + 1
  }
  
  #sample X.1t-----------------------------------------------------------------
  ### Update at all X1i.t at once
  # X.1.star = rmvnorm(n = 1,mu = X.1t,sigma = jump_x)
  X.1.star = rnorm(n = n,mean = X.1t,sd = pr.sigma.x)
  log_accept.r = log.likeli.x(X.1.star) - log.likeli.x(X.1t)
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    X.1t = X.1.star
    X.t=cbind(1,X.1t)
    accept_x = accept_x + 1
  }
  ### Individual update
  # for(idx in 1:n){
  #   X.1.star.each = rnorm(n = 1,mean = X.1t[idx],sd = pr.sigma.x)
  #   log_accept.r = log.likeli.x.each(X.1.star.each,idx) - log.likeli.x.each(X.1t[idx],idx)
  #   log.u = log(runif(1))
  #   if(log.u < log_accept.r){
  #     X.1t[idx] = X.1.star.each
  #     accept_x = accept_x + 1
  #   }
  # }
  # X.t=cbind(1,X.1t)
  
  #sample sigma2_xx-----------------------------------------------------------------
  post_ax=0.5*n+prior_ax
  post_bx=prior_bx+1/2*sum((X.1t-mux.t)^2)
  sigma2_xx.t=rinvgamma(1,post_ax,post_bx)
  stopifnot(sigma2_xx.t>0)
  #cat('accept2\n')
  
  #sample sigma2_22-----------------------------------------------------------------
  # post_a2=0.5*n+prior_a2
  # post_b2=prior_b2+1/2*sum((W2-X.1t)^2)
  # sigma2_22.t=rinvgamma(1,post_a2,post_b2)
  # stopifnot(sigma2_22.t>0)
  #cat('accept4\n')
  
  #sample mux-----------------------------------------------------------------
  post_V_mux=1/(n/sigma2_xx.t+1/pr_var_mux)
  post_M_mux=(sum(X.1t)/sigma2_xx.t+pr_mean_mux/pr_var_mux)/(n/sigma2_xx.t+1/pr_var_mux)
  mux.t=rnorm(n = 1,mean = post_M_mux,sd = sqrt(post_V_mux))
  #cat('accept5\n')
  
  #cat('------------------------------------------------------------------\n')
  
  
  # Save samples in each thin-----------------------------------------------------------------
  if((iter > nburn) & (iter %% nthin== 0) ){
    thinned_idx=(iter-nburn)/nthin
    g_trace[thinned_idx,]=g.t
    lambda_trace[thinned_idx]=lambda.t
    sigma2_xx_trace[thinned_idx]=sigma2_xx.t
    sigma2_22_trace[thinned_idx]=sigma2_22.t
    mux_trace[thinned_idx]=mux.t
    X_trace[thinned_idx,] = X.1t
    
    
    
    #### Plotting ###############################
    if((iter > nburn) & (iter %% (nthin*1000)== 0) ){
      par(mfrow=c(3,2))
      plot(colMeans(X_trace,na.rm = T),X[,2]);abline(0,1)
      
      idx=100
      ts.plot(X_trace[,idx])
      abline(h=X[idx,2])
      
      ts.plot(g_trace[,1])
      
      # plot(colMeans(alpha_trace)[1]+colMeans(alpha_trace)[2]*colMeans(X_trace),W1);abline(0,1)
      plot(colMeans(X_trace,na.rm = T),W2);abline(0,1)
      
      # ts.plot(mux_trace)
      
      g.est=colMeans(g_trace,na.rm = T)
      mspline=spline(x = tau.i,y = g.est,xout = tau.i)
      plot(X[,2],y)
      points(tau.i,mspline$y,type = 'l')    
      
      plot(colMeans(X_trace,na.rm = T),y)
      points(tau.i,mspline$y,type = 'l')    
      abline(v=tau.i)
      par(mfrow=c(1,1))
      #### Plotting ###############################
    }
  }
}
toc()



accept_g/niter
accept_x/niter

mean(mux_trace)
mean(colMeans(X_trace)[(W2>6)&(W2<8)])
mean(sigma2_22_trace)
mean(sigma2_xx_trace)