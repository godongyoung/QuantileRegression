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
  pr_sd_beta=100*diag(total_X.P)
  
  #IG for sigma
  prior_a=0.01
  prior_b=0.01
  
  #IG for sigma_delta
  prior_c=0.01
  prior_d=0.01
  
  #N for alpha
  # pr_mean_alpha=rep(0,total_T.P)
  # pr_sd_alpha=100*diag(total_T.P)
  
  
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
    post_sd_beta = solve( solve(pr_sd_beta) + t(X_t/sqrt(v_t))%*%(X_t/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_sd_beta)%*%pr_mean_beta + colSums(X_t*((y-(sigma_t*C*abs(gamma_t)*k_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
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
    # V_delta=solve(1/sigma_delta_t**2*(t(W1)%*%W1+sigma_delta_t**2*solve(pr_sd_alpha)))
    # M_delta=solve(t(W1)%*%W1+sigma_delta_t**2*solve(pr_sd_alpha))%*%t(X_1t%*%W1+sigma_delta_t**2*t(pr_mean_alpha)%*%solve(pr_sd_alpha))
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
