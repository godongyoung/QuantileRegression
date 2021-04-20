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

sigma2_11=1
sigma2_22=1

n=1000
p0_list=c(0.1,0.25,0.5,0.75,0.9)

make_data = function(X){
  set.seed(sim_idx)
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  W1=X%*%alpha+delta1
  W2=X[,2]+delta2
  
  return(list('W1'=W1, 'W2'=W2))
}



# Loop start #############################################################################################

sim_idx=1
p0=0.5
# ESL data1 #############################################################################################

set.seed(sim_idx)

inp.sd = 1
X = runif(n,0,1)
y = sin(12*(X+0.2))/(X+0.2) + rnorm(n,0,inp.sd)
Xrange = seq(0,1,length.out = 100)
plot(X,y)
for(p0 in p0_list){
  y.p0 = sin(12*(Xrange+0.2))/(Xrange+0.2) + qnorm(p0,0,inp.sd)
  points(Xrange,y.p0,col=2,lwd=2,type='l')
}
X = 10*X
X = cbind(1,X)
W_list = make_data(X)
plot(X[,2],W_list$W2);abline(0,1)
plot(X[,]%*%alpha,W_list$W1);abline(0,1)
# NQR_res = NQR_w_MME(y,W_list$W1,W_list$W2,p0,inp.min = 0,inp.max = 1,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)

W1 = W_list$W1
W2 = W_list$W2

inp.min = 0
inp.max = 1*10
inp.version = 1
multiply_c = 10
N.Knots = inp.N.Knots




# Function --------------------------------------------------------------------------------
log.likeli.g=function(g.t, lambda.t,X.1t){
  y.est=smooth.y(tau.i,g.t,X.1t,version=inp.version)
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

log.likeli.x=function(X.1t){ ### This can be changed with individual updates!
  y.est=smooth.y(tau.i,g.t,X.1t,version=inp.version)
  ui = y-y.est
  term1 = -sum(qloss(ui,p0))
  term2 = sum(dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t,sd = sqrt(sigma2_11.t),log = T))
  term3 = sum(dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T))
  term4 = sum(dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T))
  return (term1 + term2 + term3 + term4)
}

log.likeli.x.i=function(X.1t){ ### This can be changed with individual updates!
  y.est=smooth.y(tau.i,g.t,X.1t,version=inp.version)
  ui=y-y.est
  term1 = -(qloss(ui,p0))
  term2 = (dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.1t,sd = sqrt(sigma2_11.t),log = T))
  term3 = (dnorm(W2,mean = X.1t,sd = sqrt(sigma2_22.t),log = T))
  term4 = (dnorm(X.1t,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T))
  return (term1 + term2 + term3 + term4)
}
# Set Defaults --------------------------------------------------------------------------------
niter    <- 30000*multiply_c
nburn    <- 5000*multiply_c
nthin    <- 5*multiply_c
nprint   <- 10000*multiply_c
nmcmc=(niter-nburn)/nthin


n = length(y)
N=N.Knots
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
alpha.t=rmvnorm(1,c(0,0),diag(2))
sigma2_xx.t=10#runif(1)
sigma2_11.t=10#runif(1)
sigma2_22.t=10#runif(1)
mux.t=rnorm(1,0,sqrt(100))
# X.1t=cbind(1,W1)%*%t(alpha.t) #Let's start with 21, obs data
X.1t=W2
# X.1t=(cbind(1,W1)%*%t(alpha.t) +W2)/2
# X.1t=rnorm(n,0,100)
X.t=cbind(1,X.1t)

BQR_res = mBayesQR(y,X.t,p0)
beta.est=colMeans(BQR_res$beta_trace)
if (sum(is.nan(beta.est))>0){ # for unknown error
  QR_res=mQR(mcycle$accel,cbind(1,mcycle$times),p0)
  beta.est=QR_res$beta_est
}

g0=rep(0,N)
for(param.idx in 1:length(beta.est)){
  g0 = g0 + beta.est[param.idx]*tau.i^(param.idx-1) # for quadratic g0
}
g0 = smooth.y(tau.i,sin(12*(tau.i+0.2))/(tau.i+0.2),tau.i,version=3) # Starts with cheating value

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
pr_mean_alpha=rep(0,2)
pr_var_alpha=var_param*diag(2)

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

sigma.g=0.01
sigma.x=0.1
jump_g=sigma.g^2*diag(N)
jump_lambda = 1
jump_x=sigma.x^2*diag(n)

accept_g = 0
accept_l = 0
accept_x = 0

# iter start --------------------------------------------------------------------------------
tic()
for(iter in 1:niter){
  if(iter%%nprint==1){cat(iter,'th iter is going\n')}
  # Sample g-----------------------------------------------------------------
  g.star = rmvnorm(n = 1,mu = g.t,sigma = jump_g)
  log_accept.r = log.likeli.g(g.star,lambda.t,X.1t = X.1t) - log.likeli.g(g.t,lambda.t,X.1t = X.1t)
  
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
  ### This can be changed with individual updates!
  # X.1.star = t(rmvnorm(n = 1,mu = X.1t,sigma = jump_x))
  # X.1.star = rnorm(n = n,mean = X.1t,sd = sigma.x)
  # log_accept.r = log.likeli.x(X.1.star) - log.likeli.x(X.1t)
  # log.u = log(runif(1))
  # if(log.u < log_accept.r){
  #   X.1t = X.1.star
  #   X.t=cbind(1,X.1t)
  #   accept_x = accept_x + 1
  # }
  # Individual version
  X.1.star = rnorm(n = n,mean = X.1t,sd = sigma.x)
  log_accept.r.vec = log.likeli.x.i(X.1.star) - log.likeli.x.i(X.1t)
  log.u.vec = log(runif(n))
  X.1t[log.u.vec < log_accept.r.vec] = X.1.star[log.u.vec < log_accept.r.vec]
  accept_x = accept_x + sum(log.u.vec < log_accept.r.vec)
  X.t=cbind(1,X.1t)
  
  #sample sigma2_xx-----------------------------------------------------------------
  post_ax=0.5*n+prior_ax
  post_bx=prior_bx+1/2*sum((X.1t-mux.t)^2)
  sigma2_xx.t=rinvgamma(1,post_ax,post_bx)
  stopifnot(sigma2_xx.t>0)
  #cat('accept2\n')
  
  #sample sigma2_11-----------------------------------------------------------------
  post_a1=0.5*n+prior_a1
  post_b1=prior_b1+1/2*sum((W1-(X.t%*%t(alpha.t)))^2)
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
  post_V_mux = 1/(n/sigma2_xx.t+1/pr_var_mux)
  post_M_mux = post_V_mux * (sum(X.1t)/sigma2_xx.t+pr_mean_mux/pr_var_mux)
  mux.t = rnorm(n = 1,mean = post_M_mux,sd = sqrt(post_V_mux))
  #cat('accept5\n')
  
  #sample alpha-----------------------------------------------------------------
  post_V_alpha = solve(t(X.t)%*%X.t/sigma2_11.t + solve(pr_var_alpha))
  post_M_alpha = post_V_alpha %*% (t(X.t)%*%W1/sigma2_11.t + solve(pr_var_alpha)%*%pr_mean_alpha)
  alpha.t = rmvnorm(n = 1, mu = post_M_alpha,sigma = post_V_alpha)
  #cat('accept6\n')
  
  
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


res_list=list()
res_list[['g_trace']]=g_trace
res_list[['lambda_trace']]=lambda_trace
res_list[['g_accept_ratio']]=accept_g/niter
res_list[['l_accept_ratio']]=accept_l/niter
res_list[['X_trace']]=X_trace
res_list[['x_accept_ratio']]=accept_x/(niter*n)
res_list[['alpha_trace']]=alpha_trace
res_list[['mux_trace']]=mux_trace
res_list[['sigma2_22_trace']]=sigma2_22_trace
res_list[['sigma2_11_trace']]=sigma2_11_trace
res_list[['sigma2_xx_trace']]=sigma2_xx_trace
res_list[['Knots']]=tau.i
res_list[['inp.version']]=inp.version
9071.406/60

p0
plot(X[,2],y)
# NQR_res = NQR(y,X,0.1,inp.min = 0,inp.max = 1,inp.version = 1,multiply_c = inp.mul,N.Knots = inp.N.Knots)
# load('../debugging/NQR_data1_0.1_1.RData')
# points(res_list$Knots,colMeans(NQR_res$g_trace),type='o',col=2,lwd=3)
# plot(X[,2]/10,colMeans(NQR_res$X_trace));abline(0,1)
points(res_list$Knots,colMeans(res_list$g_trace),type='o',col=2,lwd=3)
points(res_list$Knots,g0,type='l',col=2,lwd=3)

res_list$g_accept_ratio
res_list$x_accept_ratio

ts.plot(g_trace[,1])
alpha.est = colMeans(alpha_trace)
ts.plot(alpha_trace[,1])
ts.plot(alpha_trace[,2])
alpha

plot(X[,2],colMeans(res_list$X_trace));abline(0,1)
plot(X[,2],W2);abline(0,1)

X.est = colMeans(X_trace)
plot(X.est,W2);abline(0,1)
mean(sigma2_22_trace)
plot(cbind(1,X.est)%*%alpha.est,W1);abline(0,1)
mean(sigma2_11_trace)

mean(sigma2_xx_trace)
plot(X.est,X[,2]);abline(0,1)
mean(mux_trace)
which(X.est<(-1))
ts.plot(X_trace[,8])

X.1t[8]
X.1.star[8]
(log.u.vec < log_accept.r.vec)[8]
X.1.star = rnorm(n = n,mean = X.1t,sd = sigma.x)



X.tmp = X.1.star
X.1t
y.est = smooth.y(tau.i,g.t,X.tmp,version=1)
plot(X.1t,y.est)
y.est[8]
plot(y,y.est)
plot(X.tmp,y.est)
plot(X.1t,y.est)

ui=y-y.est
term1 = -(qloss(ui,p0))
term2 = (dnorm(W1,mean = alpha.t[1]+alpha.t[2]*X.tmp,sd = sqrt(sigma2_11.t),log = T))
term3 = (dnorm(W2,mean = X.tmp,sd = sqrt(sigma2_22.t),log = T))
term4 = (dnorm(X.tmp,mean = mux.t,sd = sqrt(sigma2_xx.t),log = T))
c(term1[8],term2[8],term3[8],term4[8])

c(term1[8],term2[8],term3[8],term4[8])
log_accept.r.vec = log.likeli.x.i(X.1.star) - log.likeli.x.i(X.1t)
log.u.vec = log(runif(n))
X.1t[log.u.vec < log_accept.r.vec] = X.1.star[log.u.vec < log_accept.r.vec]
accept_x = accept_x + sum(log.u.vec < log_accept.r.vec)
X.t=cbind(1,X.1t)

plot(X.est,y)
points(res_list$Knots,colMeans(res_list$g_trace),type='o',col=2,lwd=3)

g.tau = g.t
knots = tau.i
xout = X.tmp

?ns
tic()
for(idx in 1:1000){
  g.tau=as.numeric(g.tau)
  fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,length(knots))]) )
  y.est=predict(fit.ns, data.frame(knots=xout))
}
toc()

tic()
for(idx in 1:1000){
  g.tau=as.numeric(g.tau)
  fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,length(knots))],Boundary.knots = knots[c(1,length(knots))]) )
  y.est=predict(fit.ns, data.frame(knots=xout))
}
toc()