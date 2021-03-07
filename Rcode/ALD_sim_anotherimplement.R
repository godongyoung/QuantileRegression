rm(list = ls())
setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
# library(mvtnorm)
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)

# Function --------------------------------------------------------------------------------
nsim=200
beta_save=matrix(NA,nrow=nsim,ncol=2)
alpha_save=matrix(NA,nrow=nsim,ncol=2)
for(sim_idx in 1:nsim){
  tic()
  set.seed(sim_idx)
  # Gen data --------------------------------------------------------------------------------
  n=200
  mux.true=1
  sigma2_xx.true=1
  sigma2_11.true=1
  sigma2_22.true=1
  
  alpha.true=c(1,1)
  beta.true=c(1,1)
  P.beta=length(beta.true)
  P.alpha=length(alpha.true)
  sigma2.true=1
  
  xi.true=rnorm(n,mean = mux.true,sd = sqrt(sigma2_xx.true))
  e2i=rnorm(n,0,sd = sqrt(sigma2_22.true))
  e1i=rnorm(n,0,sd = sqrt(sigma2_11.true))
  ei=rnorm(n,0,sd = sqrt(sigma2.true))
  
  w1i=alpha.true[1]+alpha.true[2]*xi.true+e1i
  w2i=xi.true+e2i
  X=cbind(1,xi.true)
  yi=X%*%beta.true+ei
  
  # Settings for MCMC --------------------------------------------------------------------------------
  niter    <- 30000
  nburn    <- 5000
  nthin    <- 5
  nprint   <- 10000
  nmcmc=(niter-nburn)/nthin
  
  xi.trace=matrix(NA,nrow = nmcmc,ncol=n)
  si.trace=matrix(NA,nrow = nmcmc,ncol=n)
  beta.trace=matrix(NA,nrow = nmcmc,ncol=P.beta)
  sigma2.trace=rep(NA,nmcmc)
  sigma2_11.trace=rep(NA,nmcmc)
  sigma2_22.trace=rep(NA,nmcmc)
  sigma2_xx.trace=rep(NA,nmcmc)
  mux.trace=rep(NA,nmcmc)
  alpha.trace=matrix(NA,nrow = nmcmc,ncol=P.alpha)
  
  # Priors --------------------------------------------------------------------------------
  pr_mean_beta=rep(0,P.beta)
  pr_var_beta=100*diag(P.beta)
  
  sigma.a=0.1
  sigma.b=0.1
  
  sigma2_xx.a=0.1
  sigma2_xx.b=0.1
  
  sigma2_11.a=0.1
  sigma2_11.b=0.1
  
  sigma2_22.a=0.1
  sigma2_22.b=0.1
  
  pr_mean_alpha=rep(0,P.alpha)
  pr_var_alpha=100*diag(P.alpha)
  
  pr_mean_mux=0
  pr_var_mux=100
  
  # init values --------------------------------------------------------------------------------
  xi1.t=rnorm(n,0,1)
  X.t=cbind(1,xi1.t)
  si.t=rexp(n,1)
  beta.t=rnorm(P.beta,0,1)
  sigma.t=runif(1,0,3)
  vi.t=sigma.t*si.t
  sigma2_11.t=runif(1,0,3)
  sigma2_22.t=runif(1,0,3)
  sigma2_xx.t=runif(1,0,3)
  mux.t=rnorm(1,0,1)
  alpha.t=rmvnorm(1,rep(0,P.alpha),diag(P.alpha))
  
  # Start Gibbs --------------------------------------------------------------------------------
  p0=0.5
  A=(1-2*p0)
  B=2/(p0*(1-p0))
  for(iter in 1:niter){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    # Sample beta ----------------------------------------
    
    V.b=solve(t(X.t/sqrt((B*sigma.t*vi.t)))%*%(X.t/sqrt((B*sigma.t*vi.t))) + solve(pr_var_beta))
    M.b=V.b %*% ((t(X.t/(B*sigma.t*vi.t))%*%(yi-A*vi.t)) + solve(pr_var_beta)%*%pr_mean_beta)
    beta.t=rmvnorm(n = 1,mu = M.b,sigma = V.b)
    
    # Sample vi ---------------------------------------- 
    term1 = (yi- X.t%*%t(beta.t))^2/(B*sigma.t)
    term2 = 2/sigma.t + A^2/(B*sigma.t)
    for(ii in 1:n){
      vi.t[ii]=rgig(n = 1,lambda = 0.5,chi = term1[ii],psi = term2)
    }
    si.t=vi.t/sigma.t
    
    # Sample sigma ---------------------------------------- 
    ig.sigma.a=1.5*n+sigma.a
    ig.sigma.b=0.5*sum((yi-(X.t%*%t(beta.t)+A*vi.t))^2/(B*vi.t))+sum(vi.t)+sigma.b
    # tmp=0
    # for(ii in 1:n){
    #   tmp=tmp+(yi[ii]- (X.t[ii,])%*%t(beta.t) + A*vi.t[ii])^2/(B*vi.t[ii])
    # }
    sigma.t=rinvgamma(1,ig.sigma.a,ig.sigma.b)
    
    # Sample xi1 ---------------------------------------- 
    xi1.t
    V.x=1/(beta.t[2]^2/(sigma.t*B*vi.t) + alpha.t[2]^2/sigma2_11.t + 1/sigma2_22.t + 1/sigma2_xx.t)
    M.x=(((beta.t[2]*(yi-(beta.t[1]+A*vi.t)))/(sigma.t*B*vi.t) + (alpha.t[2]*(w1i-alpha.t[1]))/(sigma2_11.t) + w2i/sigma2_22.t + mux.t/sigma2_xx.t)
         /(beta.t[2]^2/(sigma.t*B*vi.t) + alpha.t[2]^2/sigma2_11.t + 1/sigma2_22.t + 1/sigma2_xx.t))
    xi1.t=rnorm(n = n,mean = M.x,sd = sqrt(V.x))
    X.t=cbind(1,xi1.t)
    
    # Sample sigma2_11 ---------------------------------------- 
    ig.sigma2_11.a=0.5*n+sigma2_11.a
    ig.sigma2_11.b=sigma2_11.b  + 0.5*sum((w1i-(X.t%*%t(alpha.t)))^2) 
    sigma2_11.t=rinvgamma(1,ig.sigma2_11.a,ig.sigma2_11.b) 
    
    # Sample sigma2_22 ---------------------------------------- 
    ig.sigma2_22.a=0.5*n+sigma2_22.a
    ig.sigma2_22.b=sigma2_22.b  + 0.5*sum((w2i - xi1.t)^2) 
    sigma2_22.t=rinvgamma(1,ig.sigma2_22.a,ig.sigma2_22.b) 
    
    # Sample sigma2_xx ---------------------------------------- 
    ig.sigma2_xx.a=0.5*n+sigma2_xx.a
    ig.sigma2_xx.b=sigma2_xx.b  + 0.5*sum((xi1.t-mux.t)^2) 
    sigma2_xx.t=rinvgamma(1,ig.sigma2_xx.a,ig.sigma2_xx.b) 
    
    # Sample mux ---------------------------------------- 
    V.mx=1/(n/sigma2_xx.t + 1/pr_var_mux)
    M.mx=(sum(xi1.t)/sigma2_xx.t + pr_mean_mux/pr_var_mux) / (n/sigma2_xx.t + 1/pr_var_mux)
    mux.t=rnorm(n = 1,mean = M.mx,sd = sqrt(V.mx))
    
    # Sample alpha ---------------------------------------- 
    V.alpha = solve(t(X.t)%*%X.t/sigma2_11.t + solve(pr_var_alpha))
    M.alpha = V.alpha %*% ((t(X.t)%*%w1i)/sigma2_11.t + solve(pr_var_alpha)%*%pr_mean_alpha )
    alpha.t = rmvnorm(n = 1,mu = M.alpha,sigma = V.alpha)
    
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      s_idx=(iter-nburn)/nthin
      xi.trace[s_idx,]=xi1.t
      si.trace[s_idx,]=si.t
      beta.trace[s_idx,]=beta.t
      sigma2.trace[s_idx]=sigma.t^2
      sigma2_11.trace[s_idx]=sigma2_11.t
      sigma2_22.trace[s_idx]=sigma2_22.t
      sigma2_xx.trace[s_idx]=sigma2_xx.t
      mux.trace[s_idx]=mux.t
      alpha.trace[s_idx,]=alpha.t
    }
  }
  toc()
  save.image(file=sprintf('../debugging/ALD_sim_anotherimplement_%s.RData',sim_idx))
  alpha_save[sim_idx,]=colMeans(alpha.trace)
  beta_save[sim_idx,]=colMeans(beta.trace)
}


hist(beta_save[,1],nclass=10)
colMedians(alpha_save)
colMedians(beta_save)
colMeans(alpha_save)
colMeans(beta_save)

beta.true
alpha.true
p0
