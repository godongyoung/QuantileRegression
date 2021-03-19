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
type='BayesQR'
p0=0.9

n=200
P=2

# Make data--------------------------------------------------------------------------------

Mu_x=10
sigma2=1
sigma2_xx=1
sigma2_11=1
sigma2_22=1

x1i=rnorm(n = n,Mu_x,sqrt(sigma2_xx))
X=cbind(1,x1i)

beta=c(1,1)
alpha=c(1,1)

if(type=='BayesQR'){
  x1i=runif(n=n,min=0,max=10)
  X=cbind(1,x1i)
  beta=c(1,2)
  alpha=c(1,1)
  ei=rnorm(n=n, mean=0, sd=.6*x1i)
  beta.true=c(beta[1]+0*qnorm(p0,0,1),beta[2]+0.6*qnorm(p0,0,1))
}

#generate w1,w2
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

w1=X%*%alpha+delta1
w2=X[,2]+delta2

# Generate error based on input types & p0 ---------------------------------------
{if (type=='N'){
  ei=rnorm(n,0,1)
}
  if (type=='t'){
    ei=rt(n,df = 3)
  }
  if (type=='locshit'){
    ei=rnorm(n,0,1)*(1+x1i)
  }  
  if (type=='MixN'){
    if(p0==0.1){to_subtract=1.43}
    if(p0==0.25){to_subtract=0.72}
    if(p0==0.5){to_subtract=-0.02}
    if(p0==0.75){to_subtract=-0.77}
    if(p0==0.9){to_subtract=-1.555}
    
    g_indx=rbinom(n = n,size = 1,prob = c(0.1))
    n1=sum(g_indx==0)
    n2=sum(g_indx==1)
    samp1=rnorm(n=n1,mean=to_subtract,sd=1)
    samp2=rnorm(n=n2,mean=(to_subtract+1),sd=5)
    samp=c(samp1,samp2)
    ei=samp
  }}

y=cbind(1,x1i,x1i^2)%*%c(1,-5,1)+5*ei
plot(X[,2],y,main='True data plot')
# Calculate True beta based on input types & p0
{if (type=='N'){
  beta.true=beta+c(qnorm(p0,0,1),0)
}
  if (type=='t'){
    beta.true=beta+c(qt(p0,df = 3),0)
  }
  if (type=='locshit'){
    beta.true=c(beta[1]+qnorm(p0,0,1),beta[2]+qnorm(p0,0,1))
  }  
  if (type=='MixN'){
    beta.true=beta
  }}

inp.min=min(X[,2]);inp.max=max(X[,2])
W1=w1;W2=w2




#Start#################

# Function --------------------------------------------------------------------------------
qloss=function(u,p0){
  return (u*(p0-(u<0)))
}

log.likeli.g=function(g.t, lambda.t,X.t){
  mspline=spline(x = tau.i,y = g.t,xout = X.t[,2])
  ui=y-mspline$y
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


log.likeli.x=function(X.1t){ ### This can be changed with individual updates!
  mspline=spline(x = tau.i,y = g.t,xout = X.1t) # this is where X.t contributes
  ui=y-mspline$y
  term1 = -sum(qloss(ui,p0))
  term2 = dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t,sd = sqrt(sigma2_11.t),log = T)
  term3 = dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T)
  term4 = dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T)
  return (sum(term1 + term2 + term3 + term4))
}

log.likeli.x.each=function(X.1t.each,idx){ ### This can be changed with individual updates!
  mspline=spline(x = tau.i,y = g.t,xout = X.1t.each) # this is where X.t contributes
  ui=y[idx]-mspline$y
  term1 = -(qloss(ui,p0))
  term2 = dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t.each,sd = sqrt(sigma2_11.t),log = T)
  term2 = -1/(2*sigma2_11.t)*(W1[idx]-(alpha.t[1]+alpha.t[2]*X.1t.each))^2
  term3 = -1/(2*sigma2_22.t)*(W2[idx]-X.1t.each)^2
  term4 = -1/(2*sigma2_xx.t)*(X.1t.each-mux.t)^2
  return (sum(term1 + term2 + term3 + term4))
}
# Set Defaults --------------------------------------------------------------------------------
niter    <- 30000
nburn    <- 5000
nthin    <- 5
nprint   <- 10000
nmcmc=(niter-nburn)/nthin

n = length(y)
N=30
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
sigma2_11.t=10#runif(1)
sigma2_22.t=10#runif(1)
mux.t=rnorm(1,0,sqrt(100))
X.1t=cbind(1,W1)%*%alpha.t #Let's start with 21, obs data
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

# Set Priors --------------------------------------------------------------------------------

#Gamma for lambda
labmda.b=0.1/lambda0
labmda.a=lambda0/labmda.b

aaa=0.01
#IG for sigma2_xx
prior_ax=aaa;prior_bx=aaa
#IG for sigma2_11
prior_a1=aaa;prior_b1=aaa
#IG for sigma2_22
prior_a2=aaa;prior_b2=aaa

#N for mu_x
var_param=100
pr_mean_mux=0;pr_var_mux=var_param

#N for alpha
pr_mean_alpha = rep(0,2)
pr_sd_alpha = var_param*diag(2)

# Make trace --------------------------------------------------------------------------------
g_trace=matrix(NA,ncol=N,nrow=nmcmc)
lambda_trace=rep(NA,nmcmc)
sigma2_xx_trace=rep(NA,nmcmc)
sigma2_11_trace=rep(NA,nmcmc)
sigma2_22_trace=rep(NA,nmcmc)
mux_trace=rep(NA,nmcmc)
alpha_trace=matrix(NA,nrow=nmcmc,ncol=2)
X_trace=matrix(NA,nrow=nmcmc,ncol=n)   

# Set jumping rules --------------------------------------------------------------------------------

pr.sigma.g=0.01
jump_g=pr.sigma.g^2*diag(N)
jump_lambda = 0.1
pr.sigma.x=0.01
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
# log.likeli.x=function(X.1t){ ### This can be changed with individual updates!
#   mspline=spline(x = tau.i,y = g.t,xout = X.1t) # this is where X.t contributes
#   ui=y-mspline$y
#   term1 = -sum(qloss(ui,p0))
#   term2 = dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t,sd = sqrt(sigma2_11.t),log = T)
#   term3 = dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T)
#   term4 = dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T)
#   return (sum(term1 + term2 + term3 + term4))
# }
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
  
  #sample sigma2_11-----------------------------------------------------------------
  post_a1=0.5*n+prior_a1
  post_b1=prior_b1+1/2*sum((W1-(alpha.t[1]+alpha.t[2]*X.1t))^2)
  sigma2_11.t=rinvgamma(1,post_a1,post_b1)
  stopifnot(sigma2_11.t>0)
  #cat('accept3\n')
  
  #sample sigma2_22-----------------------------------------------------------------
  post_a2=0.5*n+prior_a2
  post_b2=prior_b2+1/2*sum((W2-X.1t)^2)
  sigma2_22.t=rinvgamma(1,post_a2,post_b2)
  stopifnot(sigma2_22.t>0)
  #cat('accept4\n')
  
  #sample mux-----------------------------------------------------------------
  post_V_mux=1/(n/sigma2_xx.t+1/pr_var_mux)
  post_M_mux=(sum(X.1t)/sigma2_xx.t+pr_mean_mux/pr_var_mux)/(n/sigma2_xx.t+1/pr_var_mux)
  mux.t=rnorm(n = 1,mean = post_M_mux,sd = sqrt(post_V_mux))
  #cat('accept5\n')
  
  #sample alpha-----------------------------------------------------------------
  ### Can I sample 2 alphsa in Gibbs at the same time? i.e. Collapsed Gibbs is okay? Maybe it is okay! it satisfy the condition
  post_M_alpha=solve(1/sigma2_11.t*t(X.t)%*%X.t+solve(pr_sd_alpha))%*%(1/sigma2_11.t*t(X.t)%*%W1+solve(pr_sd_alpha)%*%pr_mean_alpha)
  post_V_alpha=solve(1/sigma2_11.t*t(X.t)%*%X.t+solve(pr_sd_alpha))
  alpha.t=mvrnorm(n = 1,mu = post_M_alpha,Sigma = post_V_alpha)
  #cat('accept6\n')
  #cat('alpha.t',alpha.t,'\n')
  #cat('------------------------------------------------------------------\n')
  
  
  # Save samples in each thin-----------------------------------------------------------------
  if((iter > nburn) & (iter %% nthin== 0) ){
    thinned_idx=(iter-nburn)/nthin
    g_trace[thinned_idx,]=g.t
    lambda_trace[thinned_idx]=lambda.t
    sigma2_xx_trace[thinned_idx]=sigma2_xx.t
    sigma2_11_trace[thinned_idx]=sigma2_11.t
    sigma2_22_trace[thinned_idx]=sigma2_22.t
    mux_trace[thinned_idx]=mux.t
    alpha_trace[thinned_idx,]=alpha.t
    X_trace[thinned_idx,] = X.1t
  }
}
toc()
plot(colMeans(X_trace),X[,2])
idx=100
ts.plot(X_trace[,idx])
abline(h=X[idx,2])

ts.plot(g_trace[,1])

plot(colMeans(alpha_trace)[1]+colMeans(alpha_trace)[2]*colMeans(X_trace),W1);abline(0,1)
plot(colMeans(X_trace),W2);abline(0,1)

mean(sigma2_11_trace)
mean(sigma2_22_trace)
mean(mux_trace)
accept_g/niter
accept_x/niter
#End#################
# NQR_res=NQR_w_MME(y,W1,W2,p0,inp.min,inp.max )
# plot(colMeans(NQR_res$X_trace),X[,2])
# NQR_res$lambda_trace
# NQR_res$g_accept_ratio
# NQR_res$g_trace
ts.plot(mux_trace)

g.est=colMeans(g_trace)
mspline=spline(x = tau.i,y = g.est,xout = tau.i)
plot(X[,2],y)
points(tau.i,mspline$y,type = 'l')    
plot(colMeans(X_trace),y)
points(tau.i,mspline$y,type = 'l')    
