rm(list = ls())
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)
library(Brq)
library(bayesQR)
library(SuppDists)
library(statmod)
source('fn_wo_ME.R')
# Function --------------------------------------------------------------------------------
# These are functions for defining the domain range of parameter gamma.
tmp_func=function(gamma_t){
  tmp.y=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)
  tmp.y
}


tmp_func1=function(gamma_t){
  tmp.y=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)-(p0)
  # tmp.y=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)-(1-p0)
  tmp.y
}

tmp_func2=function(gamma_t){
  # tmp.y=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)-(p0)
  tmp.y=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)-(1-p0)
  tmp.y
}

GAL_w_SME=function(y,W1,p0){
  total_X.P=2
  
  times=1
  niter    <- 30000*times
  nburn    <- 5000*times
  nthin    <- 5*times
  nprint   <- 10000*times
  nmcmc=(niter-nburn)/nthin
  
  #unif for gamma
  U=nleqslv(0,tmp_func1)$x
  L=-nleqslv(0,tmp_func2)$x
  
  #N for beta
  pr_mean_beta=rep(0,total_X.P)
  pr_var_beta=100*diag(total_X.P)
  
  #IG for sigma
  prior_a=0.01
  prior_b=0.01
  
  #IG for sigma_delta
  prior_c=0.01
  prior_d=0.01
  
  #N for alpha
  # pr_mean_alpha=rep(0,total_T.P)
  # pr_var_alpha=100*diag(total_T.P)
  
  
  beta_trace=matrix(NA,nrow=nmcmc,ncol=total_X.P)
  sigma_trace=rep(NA,nmcmc)
  gamma_trace=rep(NA,nmcmc)
  X_trace=matrix(NA,nrow=nmcmc,ncol=n) #X,S,K is auxilary variable
  # alpha_trace=matrix(NA,nrow=nmcmc,ncol=total_T.P)
  sigma_delta_trace=rep(NA,nmcmc)
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  gamma_t=runif(1,L,U)
  # gamma_t=0.1937465
  beta_t=rnorm(total_X.P,0,1)
  alpha_t=alpha
  sigma_t=10#rinvgamma(1,prior_a,prior_b)
  sigma_delta_t=10#rinvgamma(1,prior_c,prior_d)
  s_t=rexp(n)
  k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = 0,sd = 1)
  X_1t=W1%*%alpha_t
  X_t=cbind(1,X_1t)
  
  accept_n=0
  tic()
  for(iter in 1:niter){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Define A,B,C with current gamma_t
    g_gamma=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)
    p.gamma_t=ifelse(gamma_t>0,p0/g_gamma,1-(1-p0)/g_gamma)
    stopifnot(abs(p.gamma_t)<1)
    A=(1-2*p.gamma_t)/(p.gamma_t*(1-p.gamma_t));B=2/(p.gamma_t*(1-p.gamma_t));C=1/(ifelse(gamma_t>0,1,0)-p.gamma_t)
    v_t=s_t*sigma_t
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_var_beta) + t(X_t/sqrt(v_t))%*%(X_t/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_var_beta)%*%pr_mean_beta + colSums(X_t*((y-(sigma_t*C*abs(gamma_t)*k_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
    beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    stopifnot(sum(is.nan(beta_t))==0)
    
    #sample ki-----------------------------------------------------------------
    post_var_ki=1/((C*gamma_t)**2*sigma_t/(B*v_t)+1)
    post_mean_ki=post_var_ki*C*abs(gamma_t)*(y - (X_t%*%t(beta_t)+A*v_t))/(B*v_t)
    k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = post_mean_ki,sd = sqrt(post_var_ki)) #rtruncnorm can sample multiple case, when input param is vector
    
    #sample_vi (or si)-----------------------------------------------------------------
    a_t=(y - (X_t%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t))**2/(B*sigma_t)
    b_t=2/sigma_t+A**2/(B*sigma_t)
    v_t = 1/rinvgauss(n, mean = sqrt(b_t/a_t), dispersion = 1/ b_t) # which is same for rgig
    stopifnot(sum(is.na(v_t))==0)
    # for(k in 1:n){
    #     v_t[k]=rgig(n = 1,lambda = 0.5,chi =a_t[k],psi = b_t)
    # }
    s_t=v_t/sigma_t
    
    #sample X_1t-----------------------------------------------------------------
    V_X = (sigma_t**2*B*s_t*sigma_delta_t**2)/((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
    M_X = (sigma_delta_t**2*(y-beta_t[1]-sigma_t*C*abs(gamma_t)*k_t-sigma_t*A*s_t)*beta_t[2]+sigma_t**2*B*s_t*W1%*%alpha_t) /((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
    X_1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X_t=cbind(1,X_1t)
    
    #sample sigma-----------------------------------------------------------------
    nu=-(prior_a+1.5*n)
    c=2*prior_b + 2*sum(v_t) + sum((y - (X_t%*%t(beta_t)+A*v_t))**2/(B*v_t))
    d=sum((C*gamma_t*k_t)**2/(B*v_t))
    sigma_t=rgig(n = 1,lambda = nu,chi = c,psi = d)
    
    #sample_gamma using MH wtn Gibbs-----------------------------------------------------------------
    u=rnorm(1,0,1) #sd2:0.01467065 #sd1:0.03030813 #sd0.5:0.05372212
    # because it is not a randomwalk sampler, sd should be at least 1 to fully explore the uniform interval.
    gamma_star=1/(1+exp(-u))*(U-L)-abs(L) #normal proposal on logit scale #Is it right....?
    # hist(gamma_star,nclass=100)
    # plot(1/(1+exp(-seq(-10,10)))*(U-L)-abs(L))
    # max(1/(1+exp(-seq(-10,10)))*(U-L)-abs(L))
    
    #Define A.star,B.star,C.star with gamma_star
    g_gamma.star=2*pnorm(-abs(gamma_star))*exp(gamma_star**2/2)
    p.star=ifelse(gamma_star>0,p0/g_gamma.star,1-(1-p0)/g_gamma.star)
    stopifnot(abs(p.star)<1)
    A.star=(1-2*p.star)/(p.star*(1-p.star));    B.star=2/(p.star*(1-p.star));   C.star=1/(ifelse(gamma_star>0,1,0)-p.star)
    
    numer=sum(dnorm(x = y,mean = X_t%*%t(beta_t)+sigma_t*C.star*abs(gamma_star)*k_t+sigma_t*A.star*s_t,sd = sigma_t*sqrt(B.star*s_t),log = T))
    denom=sum(dnorm(x = y,mean = X_t%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t+sigma_t*A*s_t,sd = sigma_t*sqrt(B*s_t),log = T))
    
    u2=log(runif(1))
    ratio=numer-denom
    u2<ratio
    if(u2<ratio){
      gamma_t=gamma_star
      accept_n=accept_n+1
    } 
    else{
      gamma_t=gamma_t
    }
    # gamma_t=ifelse(u2<ratio,gamma_star,gamma_t) #This code won't able to count accept counts
    stopifnot((L<gamma_t)&(gamma_t<U))
    
    #sample sigma_delta-----------------------------------------------------------------
    post_c=prior_c+1.5*n
    post_d=0.5*sum((X_1t-W1%*%alpha_t)**2)+prior_d
    sigma_delta_t=sqrt(rinvgamma(1,post_c,post_d))
    #cat('sigma_delta_t',sigma_delta_t,'\n')
    
    #sample alpha-----------------------------------------------------------------
    # V_delta=solve(1/sigma_delta_t**2*(t(W1)%*%W1+sigma_delta_t**2*solve(pr_var_alpha)))
    # M_delta=solve(t(W1)%*%W1+sigma_delta_t**2*solve(pr_var_alpha))%*%t(X_1t%*%W1+sigma_delta_t**2*t(pr_mean_alpha)%*%solve(pr_var_alpha))
    # alpha_t=mvrnorm(n = 1,mu = M_delta,Sigma = V_delta)
    #cat('alpha_t',alpha_t,'\n')
    #cat('------------------------------------------------------------------\n')
    stopifnot(sum(alpha_t==alpha)==2)
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta_t
      sigma_trace[thinned_idx]=sigma_t
      gamma_trace[thinned_idx]=gamma_t
      X_trace[thinned_idx,]=X_1t
      # alpha_trace[thinned_idx,]=alpha_t
      sigma_delta_trace[thinned_idx]=sigma_delta_t
      # cat(beta_t,gamma_t,sigma_t,'for',iter,'th iter\n')
    }
  }
  toc()
  res_list=list()
  res_list[['beta_trace']] = beta_trace
  res_list[['sigma_trace']] = sigma_trace
  res_list[['gamma_trace']] = gamma_trace
  res_list[['X_trace']] = X_trace
  res_list[['sigma_delta_trace']] = sigma_delta_trace
  
  return(res_list)
}

GAL_w_SME_nalpha=function(y,W1,p0){
  n=length(y)
  W1=W2
  W1=cbind(1,W1)
  total_X.P=2
  
  times=1
  niter    <- 30000*times
  nburn    <- 5000*times
  nthin    <- 5*times
  nprint   <- 10000*times
  nmcmc=(niter-nburn)/nthin
  
  #unif for gamma
  U=nleqslv(0,tmp_func1)$x
  L=-nleqslv(0,tmp_func2)$x
  
  #N for beta
  pr_mean_beta=rep(0,total_X.P)
  pr_var_beta=100*diag(total_X.P)
  
  #IG for sigma
  prior_a=0.01
  prior_b=0.01
  
  #IG for sigma_delta
  prior_c=0.01
  prior_d=0.01
  
  #N for alpha
  # pr_mean_alpha=rep(0,total_T.P)
  # pr_var_alpha=100*diag(total_T.P)
  
  
  beta_trace=matrix(NA,nrow=nmcmc,ncol=total_X.P)
  sigma_trace=rep(NA,nmcmc)
  gamma_trace=rep(NA,nmcmc)
  X_trace=matrix(NA,nrow=nmcmc,ncol=n) #X,S,K is auxilary variable
  # alpha_trace=matrix(NA,nrow=nmcmc,ncol=total_T.P)
  sigma_delta_trace=rep(NA,nmcmc)
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  gamma_t=runif(1,L,U)
  # gamma_t=0.1937465
  beta_t=rnorm(total_X.P,0,1)
  alpha_t=c(0,1)
  sigma_t=10#rinvgamma(1,prior_a,prior_b)
  sigma_delta_t=10#rinvgamma(1,prior_c,prior_d)
  s_t=rexp(n)
  k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = 0,sd = 1)
  X.1t=W1%*%alpha_t
  X.t=cbind(1,X.1t)
  
  accept_n=0
  tic()
  for(iter in 1:niter){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Define A,B,C with current gamma_t
    g_gamma=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)
    p.gamma_t=ifelse(gamma_t>0,p0/g_gamma,1-(1-p0)/g_gamma)
    stopifnot(abs(p.gamma_t)<1)
    A=(1-2*p.gamma_t)/(p.gamma_t*(1-p.gamma_t));B=2/(p.gamma_t*(1-p.gamma_t));C=1/(ifelse(gamma_t>0,1,0)-p.gamma_t)
    v_t=s_t*sigma_t
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_var_beta) + t(X.t/sqrt(v_t))%*%(X.t/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_var_beta)%*%pr_mean_beta + colSums(X.t*((y-(sigma_t*C*abs(gamma_t)*k_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
    beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    stopifnot(sum(is.nan(beta_t))==0)
    
    #sample ki-----------------------------------------------------------------
    post_var_ki=1/((C*gamma_t)**2*sigma_t/(B*v_t)+1)
    post_mean_ki=post_var_ki*C*abs(gamma_t)*(y - (X.t%*%t(beta_t)+A*v_t))/(B*v_t)
    k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = post_mean_ki,sd = sqrt(post_var_ki)) #rtruncnorm can sample multiple case, when input param is vector
    
    #sample_vi (or si)-----------------------------------------------------------------
    a_t=(y - (X.t%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t))**2/(B*sigma_t)
    b_t=2/sigma_t+A**2/(B*sigma_t)
    v_t = 1/rinvgauss(n, mean = sqrt(b_t/a_t), dispersion = 1/ b_t) # which is same for rgig
    stopifnot(sum(is.na(v_t))==0)
    s_t=v_t/sigma_t
    
    #sample X.1t-----------------------------------------------------------------
    V_X = (sigma_t**2*B*s_t*sigma_delta_t**2)/((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
    M_X = (sigma_delta_t**2*(y-beta_t[1]-sigma_t*C*abs(gamma_t)*k_t-sigma_t*A*s_t)*beta_t[2]+sigma_t**2*B*s_t*W1%*%alpha_t) /((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
    X.1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X.t=cbind(1,X.1t)
    
    #sample sigma-----------------------------------------------------------------
    nu=-(prior_a+1.5*n)
    c=2*prior_b + 2*sum(v_t) + sum((y - (X.t%*%t(beta_t)+A*v_t))**2/(B*v_t))
    d=sum((C*gamma_t*k_t)**2/(B*v_t))
    sigma_t=rgig(n = 1,lambda = nu,chi = c,psi = d)
    
    #sample_gamma using MH wtn Gibbs-----------------------------------------------------------------
    u=rnorm(1,0,1) #sd2:0.01467065 #sd1:0.03030813 #sd0.5:0.05372212
    gamma_star=1/(1+exp(-u))*(U-L)-abs(L) #normal proposal on logit scale #Is it right....?
    
    #Define A.star,B.star,C.star with gamma_star
    g_gamma.star=2*pnorm(-abs(gamma_star))*exp(gamma_star**2/2)
    p.star=ifelse(gamma_star>0,p0/g_gamma.star,1-(1-p0)/g_gamma.star)
    stopifnot(abs(p.star)<1)
    A.star=(1-2*p.star)/(p.star*(1-p.star));    B.star=2/(p.star*(1-p.star));   C.star=1/(ifelse(gamma_star>0,1,0)-p.star)
    
    numer=sum(dnorm(x = y,mean = X.t%*%t(beta_t)+sigma_t*C.star*abs(gamma_star)*k_t+sigma_t*A.star*s_t,sd = sigma_t*sqrt(B.star*s_t),log = T))
    denom=sum(dnorm(x = y,mean = X.t%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t+sigma_t*A*s_t,sd = sigma_t*sqrt(B*s_t),log = T))
    
    u2=log(runif(1))
    ratio=numer-denom
    u2<ratio
    if(u2<ratio){
      gamma_t=gamma_star
      accept_n=accept_n+1
    } 
    else{
      gamma_t=gamma_t
    }
    # gamma_t=ifelse(u2<ratio,gamma_star,gamma_t) #This code won't able to count accept counts
    stopifnot((L<gamma_t)&(gamma_t<U))
    
    #sample sigma_delta-----------------------------------------------------------------
    post_c=prior_c+1.5*n
    post_d=0.5*sum((X.1t-W1%*%alpha_t)**2)+prior_d
    sigma_delta_t=sqrt(rinvgamma(1,post_c,post_d))
    #cat('sigma_delta_t',sigma_delta_t,'\n')
    
    #cat('------------------------------------------------------------------\n')
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta_t
      sigma_trace[thinned_idx]=sigma_t
      gamma_trace[thinned_idx]=gamma_t
      X_trace[thinned_idx,]=X.1t
      sigma_delta_trace[thinned_idx]=sigma_delta_t
      # cat(beta_t,gamma_t,sigma_t,'for',iter,'th iter\n')
    }
  }
  toc()
  res_list=list()
  res_list[['beta_trace']] = beta_trace
  res_list[['sigma_trace']] = sigma_trace
  res_list[['gamma_trace']] = gamma_trace
  res_list[['X_trace']] = X_trace
  res_list[['sigma_delta_trace']] = sigma_delta_trace
  
  return(res_list)
}
NQR_w_MME=function(y,W1,W2,p0,inp.min,inp.max,multiply_c=2,inp.version=3){
  # Function --------------------------------------------------------------------------------
  qloss=function(u,p0){
    return (u*(p0-(u<0)))
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
      fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
      y.est=predict(fit.ns, data.frame(knots=xout))
    }
    return(y.est)
  }
  
  
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
  alpha.t=rmvnorm(1,c(0,0),diag(2))
  sigma2_xx.t=10#runif(1)
  sigma2_11.t=10#runif(1)
  sigma2_22.t=10#runif(1)
  mux.t=rnorm(1,0,sqrt(100))
  X.1t=cbind(1,W1)%*%t(alpha.t) #Let's start with 21, obs data
  X.1t=W2
  # X.1t=(cbind(1,W1)%*%t(alpha.t) +W2)/2
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
  return(res_list)
}


NQR_w_SME=function(y,W2,p0,is.plot=F){
  print('Current version is wo alpha, & fix sigma2_22 for simplification')
  n=length(y)
  inp.version=3
  multiply_c=1
  N=30
  
  # Function --------------------------------------------------------------------------------
  qloss=function(u,p0){
    return (u*(p0-(u<0)))
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
      fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
      y.est=predict(fit.ns, data.frame(knots=xout))
    }
    return(y.est)
  }
  
  
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
  if(is.plot){
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
  }
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
  sigma2_22.t=1
  
  g.true=smooth.y(X[,2],y,tau.i,version=2)
  # ui=y-smooth.y(tau.i,g.t,X.t[,2],version=inp.version)
  # term1 = -sum(qloss(ui,p0));term1
  # 
  # ui=y-smooth.y(tau.i,g.true,X.t[,2],version=inp.version)
  # term1 = -sum(qloss(ui,p0));term1
  # 
  # plot(X[,2],y)
  # points(tau.i,smooth.y(tau.i,g.true,tau.i,version=inp.version),type = 'l')
  ###########Debugging part end
  
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
      
      
      if(is.plot){
        #### Plotting ###############################
        if((iter > nburn) & (iter %% (nthin*1000)== 0) ){
          par(mfrow=c(3,2))
          plot(colMeans(X_trace,na.rm = T),X[,2],main='X.est Vs X');abline(0,1)
          
          idx=100
          ts.plot(X_trace[,idx],main='trace for x_100')
          abline(h=X[idx,2])
          
          ts.plot(g_trace[,1],main='trace for g_1')
          
          # plot(colMeans(alpha_trace)[1]+colMeans(alpha_trace)[2]*colMeans(X_trace),W1);abline(0,1)
          plot(colMeans(X_trace,na.rm = T),W2,main='X.est vs W2');abline(0,1)
          
          # ts.plot(mux_trace)
          
          g.est=colMeans(g_trace,na.rm = T)
          mspline=spline(x = tau.i,y = g.est,xout = tau.i)
          plot(X[,2],y,main='g.est in X')
          points(tau.i,mspline$y,type = 'l')    
          
          plot(colMeans(X_trace,na.rm = T),y,main='g.est in X.est')
          points(tau.i,mspline$y,type = 'l')    
          abline(v=tau.i)
          par(mfrow=c(1,1))
          #### Plotting ###############################
        }        
      }

    }
  }
  toc()
  res_list=list()
  res_list[['g_trace']]=g_trace
  res_list[['lambda_trace']]=lambda_trace
  res_list[['sigma2_xx_trace']]=sigma2_xx_trace
  res_list[['mux_trace']]=mux_trace
  res_list[['X_trace']]=X_trace
  res_list[['l_accept_ratio']]=accept_l/niter
  res_list[['g_accept_ratio']]=accept_g/niter
  res_list[['x_accept_ratio']]=accept_x/niter
  res_list[['Knots']]=tau.i
  return(res_list)
}

ALD_w_MME=function(y,W1,W2,p0){
  n=length(y)
  
  # set default--------------------------------------------------------------------------------
  niter    <- 30000
  nburn    <- 5000
  nthin    <- 5
  nprint   <- 10000
  nmcmc=(niter-nburn)/nthin
  
  
  # Prior Setting ----------------------------------------------------------------------  
  #Define A,B
  stopifnot(abs(p0)<1)
  A=(1-2*p0)/(p0*(1-p0));B=2/(p0*(1-p0));
  
  #N for beta
  pr_mean_beta=rep(0,2)
  pr_var_beta=100*diag(2)
  
  #IG for sigma
  aaa=0.01
  prior_a=aaa
  prior_b=aaa
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
  
  # Make Trace ----------------------------------------------------------------------
  beta_trace=matrix(NA,nrow=nmcmc,ncol=2)
  sigma_trace=rep(NA,nmcmc)
  s_trace=matrix(NA,nrow=nmcmc,ncol=n)
  sigma2_xx_trace=rep(NA,nmcmc)
  sigma2_11_trace=rep(NA,nmcmc)
  sigma2_22_trace=rep(NA,nmcmc)
  mux_trace=rep(NA,nmcmc)
  alpha_trace=matrix(NA,nrow=nmcmc,ncol=2)
  X_trace=matrix(NA,nrow=nmcmc,ncol=n) #X,S,K is auxilary variable, but let's keep X only.
  ## z and s are auxilary variable
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  beta.t=rmvnorm(1,c(0,0),diag(2))
  alpha.t=rmvnorm(1,c(0,0),diag(2))
  sigma.t=10#runif(1)
  sigma2_xx.t=10#runif(1)
  sigma2_11.t=10#runif(1)
  sigma2_22.t=10#runif(1)
  mux.t=rnorm(1,pr_mean_mux,sqrt(pr_var_mux))
  s.t=rexp(n)
  
  accept_n=0
  
  X.1t=cbind(1,W1)%*%t(alpha.t) #Let's start with 21, obs data
  X.1t=W2
  X.t=cbind(1,X.1t)
  
  v.t=s.t*sigma.t
  
  
  tic()
  for(iter in 1:niter){
    # for(iter in 1:10){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_var_beta) + t(X.t/sqrt(v.t))%*%(X.t/sqrt(v.t))/(B*sigma.t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_var_beta)%*%pr_mean_beta + colSums(X.t*((y-(A*v.t))/(B*sigma.t*v.t))[1:n]))
    beta.t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    
    #sample_vi (or zi)-----------------------------------------------------------------
    a_t=(y - (X.t%*%t(beta.t)))**2/(B*sigma.t)
    b_t=2/sigma.t+A**2/(B*sigma.t)
    # v.t=sapply(X = a_t, FUN = function(x) rgig(n = 1,lambda = 0.5, chi = x, psi = b_t))
    for(ii in 1:n){
      v.t[ii]=rgig(n = 1,lambda = 0.5,chi =a_t[ii],psi = b_t)
    }
    # v.t=rgig(n = 1,lambda = 0.5,chi =a_t,psi = b_t) ### error version
    s.t=v.t/sigma.t
    
    #sample sigma-----------------------------------------------------------------
    nu=(prior_a+1.5*n)
    c=0.5 * (2*prior_b + 2*sum(v.t) + sum((y - (X.t%*%t(beta.t)+A*v.t))**2/(B*v.t)))
    sigma.t=rinvgamma(n = 1,shape = nu ,scale = c)
    s.t=v.t/sigma.t #Is this right??? or should I just leave s.t as previous(i.e.unupdated) one?
    
    
    #sample X.1t-----------------------------------------------------------------
    V_X= 1 / (beta.t[2]^2/(sigma.t*B*v.t)+alpha.t[2]^2/sigma2_11.t+1/sigma2_22.t+1/sigma2_xx.t)
    M_X=V_X * ((beta.t[2]*(y-(beta.t[1]+A*v.t)))/(sigma.t*B*v.t)+
                 alpha.t[2]*(W1-alpha.t[1])/sigma2_11.t+
                 W2/sigma2_22.t+
                 mux.t/sigma2_xx.t)
    X.1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X.t=cbind(1,X.1t)
    
    #sample sigma2_xx-----------------------------------------------------------------
    post_ax=0.5*n+prior_ax
    post_bx=prior_bx+1/2*sum((X.1t-mux.t)^2)
    sigma2_xx.t=rinvgamma(1,post_ax,post_bx)
    stopifnot(sigma2_xx.t>0)
    #cat('accept2\n')
    
    #sample sigma2_11-----------------------------------------------------------------
    post_a1=0.5*n+prior_a1
    post_b1=prior_b1+1/2*sum((W1-X.t%*%t(alpha.t))^2)
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
    post_V_alpha = solve(t(X.t)%*%X.t/sigma2_11.t + solve(pr_var_alpha))
    post_M_alpha = post_V_alpha %*% (t(X.t)%*%W1/sigma2_11.t + solve(pr_var_alpha)%*%pr_mean_alpha)
    alpha.t = rmvnorm(n = 1, mu = post_M_alpha,sigma = post_V_alpha)
    #cat('accept6\n')
    #cat('alpha.t',alpha.t,'\n')
    #cat('------------------------------------------------------------------\n')
    
    
    
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta.t
      sigma_trace[thinned_idx]=sigma.t
      sigma2_xx_trace[thinned_idx]=sigma2_xx.t
      sigma2_11_trace[thinned_idx]=sigma2_11.t
      sigma2_22_trace[thinned_idx]=sigma2_22.t
      mux_trace[thinned_idx]=mux.t
      alpha_trace[thinned_idx,]=alpha.t
      X_trace[thinned_idx,]=X.1t
      s_trace[thinned_idx,]=s.t
      
    }
  }
  toc()
  
  res_list=list()
  res_list[['beta_trace']]=beta_trace
  res_list[['sigma_trace']]=sigma_trace
  res_list[['sigma2_xx_trace']]=sigma2_xx_trace
  res_list[['sigma2_11_trace']]=sigma2_11_trace
  res_list[['sigma2_22_trace']]=sigma2_22_trace
  res_list[['mux_trace']]=mux_trace
  res_list[['alpha_trace']]=alpha_trace
  res_list[['X_trace']]=X_trace
  
  return(res_list)  
}

BLR_w_MME=function(y,W1,W2){
  n=length(y)
  
  # set default--------------------------------------------------------------------------------
  niter    <- 30000*multiply_c
  nburn    <- 5000*multiply_c
  nthin    <- 5*multiply_c
  nprint   <- 10000*multiply_c
  nmcmc=(niter-nburn)/nthin
  
  
  # Prior Setting ----------------------------------------------------------------------  
  #N for beta
  pr_mean_beta=rep(0,2)
  pr_var_beta=100*diag(2)
  
  #IG for sigma
  aaa=0.01
  prior_a=aaa
  prior_b=aaa
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
  
  # Make Trace ----------------------------------------------------------------------
  beta_trace=matrix(NA,nrow=nmcmc,ncol=2)
  X_trace=matrix(NA,nrow=nmcmc,ncol=n) 
  sigma2_trace=rep(NA,nmcmc)
  sigma2_xx_trace=rep(NA,nmcmc)
  sigma2_11_trace=rep(NA,nmcmc)
  sigma2_22_trace=rep(NA,nmcmc)
  mux_trace=rep(NA,nmcmc)
  alpha_trace=matrix(NA,nrow=nmcmc,ncol=2)
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  beta.t=rmvnorm(1,c(0,0),diag(2))
  alpha.t=rmvnorm(1,c(0,0),diag(2))
  sigma2.t=10#runif(1)
  sigma2_xx.t=10#runif(1)
  sigma2_11.t=10#runif(1)
  sigma2_22.t=10#runif(1)
  mux.t=rnorm(1,pr_mean_mux,sqrt(pr_var_mux))
  
  accept_n=0
  
  
  X.1t=cbind(1,W1)%*%t(alpha.t) #Let's start with 21, obs data
  X.1t=W2
  X.t=cbind(1,X.1t)
  
  tic()
  for(iter in 1:niter){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve(t(X.t)%*%(X.t)/sigma2.t + solve(pr_var_beta))
    post_mean_beta = post_sd_beta %*% (solve(pr_var_beta)%*%pr_mean_beta + t(X.t)%*%y/sigma2.t)
    beta.t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)# tmp=c(0,0);idx=1
    
    #sample sigma-----------------------------------------------------------------
    post_a=0.5*n + prior_a
    post_b=prior_b + 1/2*sum((y-X.t%*%t(beta.t))^2)
    sigma2.t=rinvgamma(1,post_a,post_b)
    stopifnot(sigma2.t>0)
    
    #sample X.1t-----------------------------------------------------------------
    V_X = 1/(beta.t[2]^2/sigma2.t + alpha.t[2]^2/sigma2_11.t + 1/sigma2_22.t + 1/sigma2_xx.t)
    M_X = V_X * (beta.t[2]*(y-beta.t[1])/sigma2.t +
                   alpha.t[2]*(W1-alpha.t[1])/sigma2_11.t +
                   W2/sigma2_22.t + 
                   mux.t/sigma2_xx.t
    )
    X.1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X.t=cbind(1,X.1t)
    
    #sample sigma2_11-----------------------------------------------------------------
    post_a1 = 0.5*n+prior_a1
    post_b1 = prior_b1+1/2*sum((W1-X.t%*%t(alpha.t))^2)
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
    
    #sample alpha-----------------------------------------------------------------
    post_V_alpha = solve(t(X.t)%*%X.t/sigma2_11.t + solve(pr_var_alpha))
    post_M_alpha = post_V_alpha %*% (t(X.t)%*%W1/sigma2_11.t + solve(pr_var_alpha)%*%pr_mean_alpha)
    alpha.t = rmvnorm(n = 1, mu = post_M_alpha,sigma = post_V_alpha)
    #cat('accept6\n')
    #cat('alpha.t',alpha.t,'\n')
    #cat('------------------------------------------------------------------\n')
    
    
    
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta.t
      X_trace[thinned_idx,]=X.1t
      sigma2_trace[thinned_idx]=sigma2.t
      sigma2_11_trace[thinned_idx]=sigma2_11.t
      sigma2_22_trace[thinned_idx]=sigma2_22.t
      sigma2_xx_trace[thinned_idx]=sigma2_xx.t
      mux_trace[thinned_idx]=mux.t
      alpha_trace[thinned_idx,]=alpha.t
    }
  }
  toc()
  res_list=list()
  res_list[['beta_trace']]=beta_trace
  res_list[['X_trace']]=X_trace
  res_list[['sigma2_trace']]=sigma2_trace
  res_list[['sigma2_11_trace']]=sigma2_11_trace
  res_list[['sigma2_22_trace']]=sigma2_22_trace
  res_list[['sigma2_xx_trace']]=sigma2_xx_trace
  res_list[['mux_trace']]=mux_trace
  res_list[['alpha_trace']]=alpha_trace
  
  return(res_list)
}

