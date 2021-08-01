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

X_demean = T
make_X_shit = function(Mu_x, Xmul){
  if(X_demean){
    X_shit = Mu_x*Xmul
  }
  else{X_shit = 0}
  
  return(X_shit)
}

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
    fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,length(knots))]) )
    y.est=predict(fit.ns, data.frame(knots=xout))
  }
  return(y.est)
}


log.likeli.x.i=function(X.1t){ ### This can be changed with individual updates!
  y.est = smooth.y(tau.i,g.y.t,X.1t,version=inp.version)
  term1 = (dnorm(y,mean = y.est,sd = sqrt(sigma2.t),log = T))
  
  w.est = smooth.y(tau.i,g.w.t,X.1t,version=inp.version)
  term2 = (dnorm(W1,mean = w.est,sd = sqrt(sigma2_11.t),log = T))
  
  term3 = (dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T))
  term4 = (dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T))
  return (term1 + term2 + term3 + term4)
}

log.likeli.g.y=function(g.y.t, lambda.t,X.1t){
  y.est.t = smooth.y(tau.i,g.y.t,X.1t,version=inp.version)
  term1 = sum(dnorm(y,mean = y.est.t,sd = sqrt(sigma2.t),log = T))
  term2 = -0.5 * lambda.t * (g.y.t) %*% K %*% t(g.y.t)
  return(term1+term2)
}
log.likeli.g.w=function(g.w.t, lambda.t,X.1t){
  w.est.t = smooth.y(tau.i,g.w.t,X.1t,version=inp.version)
  term1 = sum(dnorm(W1,mean = w.est.t,sd = sqrt(sigma2_11.t),log = T))
  term2 = -0.5 * lambda.t * (g.w.t) %*% K %*% t(g.w.t)
  return(term1+term2)
}


log.likeli.l=function(lambda.f,lambda.b,g.t){
  term1 = (N-2)/2*log(lambda.f) 
  term2 = -0.5 * lambda.f * (g.t) %*% K %*% t(g.t)
  term3 = dgamma(x = lambda.f,shape = labmda.a, scale = labmda.b, log = T)
  term4 = dlnorm(x = lambda.b, meanlog = log(lambda.f), sdlog = sqrt(jump_lambda),log = T)  
  return (term1 + term2 + term3)
}
# Define True parameter--------------------------------------------------------------------------------
inp.version = 1
Knots.direct = NA
N.Knots = 30
inp.min = -5;inp.max = 5
n=1000
multiply_c=2
beta=c(3,5)

Mu_x=5
sigma2=1
sigma2_xx=1
sigma2_11=0.5
sigma2_22=1

# Simulation start --------------------------------------------------------------------------------
sim_idx=1
nmax=500
is.plot=F

# Make data--------------------------------------------------------------------------------
Xmul = 10
X_shit = make_X_shit(0.5,Xmul)
set.seed(sim_idx)
x1 = runif(n,0,1)
inp.sd = 1
y = sin(12*(x1+0.2))/(x1+0.2) + rnorm(n,0,inp.sd)
X=cbind(1,x1*Xmul-X_shit)

#generate W1,W2
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))


W1=x1*sin(2.5*pi*x1)+delta1
W2=X[,2]+delta2

plot(W2,y)
plot(X[,2],y)
plot(X[,2],W1)


# set default--------------------------------------------------------------------------------
niter    <- 30000*multiply_c
nburn    <- 5000*multiply_c
nthin    <- 5*multiply_c
nprint   <- 10000*multiply_c
nmcmc=(niter-nburn)/nthin


n=length(y)
N=N.Knots
# tau.i=seq(from = min(X[,2]),to = max(X[,2]),length.out = N)
{if(all(is.na(Knots.direct))){
  tau.i=seq(from = inp.min,to = inp.max,length.out = N)  
}
  else{
    tau.i = Knots.direct
  }}

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
sigma2.t=10#runif(1)
sigma2_xx.t=10#runif(1)
sigma2_11.t=10#runif(1)
sigma2_22.t=10#runif(1)
mux.t=rnorm(1,0,sqrt(100))

accept_n=0


X.1t=W2
X.t=cbind(1,X.1t)

# For starting with NCS
fit.ns <- lm(W1~ ns(x = X.1t, knots = tau.i[-c(1,length(tau.i))]) )
g.w0 = predict(fit.ns, data.frame(X.1t=tau.i))


msmooth.spline=smooth.spline(x = tau.i,y = g.w0,control.spar = list('maxit'=1,'trace'=F))
lambda.w0=msmooth.spline$lambda

g.w.t=matrix(g.w0,nrow = 1)
lambda.w.t=lambda.w0

fit.ns <- lm(y~ ns(x = X.1t, knots = tau.i[-c(1,length(tau.i))]) )
g.y0 = predict(fit.ns, data.frame(X.1t=tau.i))


msmooth.spline=smooth.spline(x = tau.i,y = g.y0,control.spar = list('maxit'=1,'trace'=F))
lambda.y0=msmooth.spline$lambda

g.y.t=matrix(g.y0,nrow = 1)
lambda.y.t=lambda.y0
# Prior Setting ----------------------------------------------------------------------  
#N for beta
pr_mean_beta=rep(0,2)
pr_var_beta=100*diag(2)

#IG for sigma
aaa=0.01
prior_a=aaa;prior_b=aaa
#IG for sigma2_xx
prior_ax=aaa;prior_bx=aaa
#IG for sigma2_11
prior_a1=aaa;prior_b1=aaa
#IG for sigma2_22
prior_a2=aaa;prior_b2=aaa

#N for mu_x
var_param=100
pr_mean_mux=0;pr_var_mux=var_param

#Gamma for lambda
labmda.b=0.1/lambda.y0
labmda.a=lambda.y0/labmda.b
# labmda.a=0.1/lambda0
# labmda.b=lambda0/labmda.a

# Make Trace ----------------------------------------------------------------------
X_trace=matrix(NA,nrow=nmcmc,ncol=n) 
sigma2_trace=rep(NA,nmcmc)
sigma2_xx_trace=rep(NA,nmcmc)
sigma2_11_trace=rep(NA,nmcmc)
sigma2_22_trace=rep(NA,nmcmc)
mux_trace=rep(NA,nmcmc)
g.y_trace=matrix(NA,ncol=N,nrow=nmcmc)
lambda.y_trace=rep(NA,nmcmc)
g.w_trace=matrix(NA,ncol=N,nrow=nmcmc)
lambda.w_trace=rep(NA,nmcmc)

# Set jumping rules --------------------------------------------------------------------------------

sigma.g=0.01
sigma.x=0.1
jump_g=sigma.g^2*diag(N)
jump_lambda = 1
jump_x=sigma.x^2*diag(n)

accept_g.y = 0
accept_l.y = 0
accept_g.w = 0
accept_l.w = 0
accept_x = 0

# Iteration--------------------------------------------------------------------------------
tic()
for(iter in 1:niter){
  if(iter%%nprint==1){cat(iter,'th iter is going\n')}
  
  # Sample g.y-----------------------------------------------------------------
  g.y.star = rmvnorm(n = 1,mu = g.y.t,sigma = jump_g)
  log_accept.r = log.likeli.g.y(g.y.star,lambda.y.t,X.1t = X.1t) - log.likeli.g.y(g.y.t,lambda.y.t,X.1t = X.1t)
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    g.y.t = g.y.star
    accept_g.y = accept_g.y + 1
  }
  y.est.t = smooth.y(tau.i,g.y.t,X.1t,version=inp.version)
  
  # Sample lambda.y-----------------------------------------------------------------
  # lamda.star = exp(rnorm(1,mean = log(lambda.t),sd = sqrt(jump_lambda)))
  lambda.y.star = rlnorm(1,meanlog = log(lambda.y.t),sdlog = sqrt(jump_lambda))
  log_accept.r = log.likeli.l(lambda.y.star,lambda.y.t,g.y.t) - log.likeli.l(lambda.y.t,lambda.y.star,g.y.t)  
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    lambda.y.t = lambda.y.star
    accept_l.y = accept_l.y + 1
  }
  
  
  #sample sigma-----------------------------------------------------------------
  post_a=0.5*n + prior_a
  post_b=prior_b + 1/2*sum((y-y.est.t)^2)
  sigma2.t=rinvgamma(1,post_a,post_b)
  stopifnot(sigma2.t>0)
  
  
  #sample X.1t-----------------------------------------------------------------
  # Individual version
  X.1.star = rnorm(n = n,mean = X.1t,sd = sigma.x)
  log_accept.r.vec = log.likeli.x.i(X.1.star) - log.likeli.x.i(X.1t)
  log.u.vec = log(runif(n))
  X.1t[log.u.vec < log_accept.r.vec] = X.1.star[log.u.vec < log_accept.r.vec]
  accept_x = accept_x + sum(log.u.vec < log_accept.r.vec)
  X.t=cbind(1,X.1t)
  
  w.est.t = smooth.y(tau.i,g.w.t,X.1t,version=inp.version)
  
  #sample sigma2_11-----------------------------------------------------------------
  post_a1 = 0.5*n+prior_a1
  post_b1 = prior_b1+1/2*sum((W1-w.est.t)^2)
  sigma2_11.t = rinvgamma(1,post_a1,post_b1)
  stopifnot(sigma2_11.t>0)
  #cat('accept1\n')
  
  #sample sigma2_22-----------------------------------------------------------------
  post_a2=0.5*n+prior_a2
  post_b2=prior_b2+1/2*sum((W2-X.1t)^2)
  sigma2_22.t = rinvgamma(1,post_a2,post_b2)
  stopifnot(sigma2_22.t>0)
  #cat('accept2\n')
  
  #sample sigma2_xx-----------------------------------------------------------------
  post_ax=0.5*n+prior_ax
  post_bx=prior_bx+1/2*sum((X.1t-mux.t)^2)
  sigma2_xx.t=rinvgamma(1,post_ax,post_bx)
  stopifnot(sigma2_xx.t>0)
  #cat('accept3\n')
  
  #sample mux-----------------------------------------------------------------
  post_V_mux = 1/(n/sigma2_xx.t+1/pr_var_mux)
  post_M_mux = post_V_mux * (sum(X.1t)/sigma2_xx.t+pr_mean_mux/pr_var_mux)
  mux.t=rnorm(n = 1,mean = post_M_mux,sd = sqrt(post_V_mux))
  #cat('accept5\n')
  
  
  # Sample g.w-----------------------------------------------------------------
  g.w.star = rmvnorm(n = 1,mu = g.w.t,sigma = jump_g)
  log_accept.r = log.likeli.g.w(g.w.star,lambda.w.t,X.1t = X.1t) - log.likeli.g.w(g.w.t,lambda.w.t,X.1t = X.1t)
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    g.w.t = g.w.star
    accept_g.w = accept_g.w + 1
  }
  
  # Sample lambda.w-----------------------------------------------------------------
  # lamda.star = exp(rnorm(1,mean = log(lambda.w.t),sd = sqrt(jump_lambda)))
  lambda.w.star = rlnorm(1,meanlog = log(lambda.w.t),sdlog = sqrt(jump_lambda))
  log_accept.r = log.likeli.l(lambda.w.star,lambda.w.t,g.w.t) - log.likeli.l(lambda.w.t,lambda.w.star,g.w.t)  
  
  log.u = log(runif(1))
  if(log.u < log_accept.r){
    lambda.w.t = lambda.w.star
    accept_l.w = accept_l.w + 1
  }
  
  
  # Save samples in each thin-----------------------------------------------------------------
  if((iter > nburn) & (iter %% nthin== 0) ){
    thinned_idx=(iter-nburn)/nthin
    X_trace[thinned_idx,]=X.1t
    sigma2_trace[thinned_idx]=sigma2.t
    sigma2_11_trace[thinned_idx]=sigma2_11.t
    sigma2_22_trace[thinned_idx]=sigma2_22.t
    sigma2_xx_trace[thinned_idx]=sigma2_xx.t
    mux_trace[thinned_idx]=mux.t
    g.y_trace[thinned_idx,]=g.y.t
    g.w_trace[thinned_idx,]=g.w.t
    lambda.y_trace[thinned_idx]=lambda.y.t
    lambda.w_trace[thinned_idx]=lambda.w.t
  }
}
toc()
res_list=list()
res_list[['X_trace']]=X_trace
res_list[['sigma2_trace']]=sigma2_trace
res_list[['sigma2_11_trace']]=sigma2_11_trace
res_list[['sigma2_22_trace']]=sigma2_22_trace
res_list[['sigma2_xx_trace']]=sigma2_xx_trace
res_list[['mux_trace']]=mux_trace
res_list[['g.y_trace']]=g.y_trace
res_list[['g.w_trace']]=g.w_trace
res_list[['lambda.w_trace']]=lambda.w_trace
res_list[['lambda.y_trace']]=lambda.y_trace
res_list[['g.y_accept_ratio']]=accept_g.y/niter
res_list[['l.y_accept_ratio']]=accept_l.y/niter
res_list[['g.w_accept_ratio']]=accept_g.w/niter
res_list[['l.w_accept_ratio']]=accept_l.w/niter
res_list[['Knots']]=tau.i
res_list[['inp.version']]=inp.version


res_list$l.w_accept_ratio
res_list$l.y_accept_ratio
res_list$g.w_accept_ratio
res_list$g.y_accept_ratio
X.est = colMeans(res_list$X_trace)
plot(X[,2],X.est);abline(0,1)
g.y.est = colMeans(res_list$g.y_trace)
g.w.est = colMeans(res_list$g.w_trace)
plot(X[,2],W1)
points(res_list$Knots,g.w.est,type='l')
plot(X[,2],y)
points(res_list$Knots,g.y.est,type='l')


mean(res_list$sigma2_11_trace)
mean(res_list$sigma2_22_trace)
mean(res_list$sigma2_xx_trace)
mean(res_list$sigma2_trace)
ts.plot(res_list$sigma2_trace)
