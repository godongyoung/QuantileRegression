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

my_hist=function(inp_data,true_value,inp_text){
  hist(inp_data,nclass=100,main=inp_text,xlab=inp_text)
  mtext(sprintf('Median : %s',round(median(inp_data),2)),side=3)
  abline(v=true_value,col=2,lwd=3)
}
nmax=200
for(sim_idx in 1:nmax){
  # Make data--------------------------------------------------------------------------------
  set.seed(sim_idx)
  p0=0.25
  n=200
  P=2
  
  P=P-1
  covmat=matrix(NA,P,P)
  for(i in 1:P){
    for(j in 1:P){
      covmat[i,j]=0.5**(abs(i-j))
    }
  }
  
  Mu_x=10
  sigma2=1
  sigma2_xx=2
  sigma2_11=2
  sigma2_22=2
  
  
  X=rmvnorm(n = n,mu = rep(Mu_x,P),sigma = covmat);dim(X)
  colMeans(X);cor(X)
  X=cbind(1,X)
  P=P+1
  beta=c(3,1.5,2)
  beta=c(3,1.5,2,1,2,3,4,5,6,7,8,9,10)[1:P]
  # beta=sample(1:100,P)
  # beta=c(10,2)
  beta=c(50,50)
  alpha=sample(1:100,2)
  
  
  #generate w1,w2
  delta1=rnorm(n,0,sd=sqrt(sigma2_11))
  delta2=rnorm(n,0,sd=sqrt(sigma2_22))
  
  w2=X[,2]+delta2
  w1=X%*%alpha+delta1
  
  # to_subtract=qnorm(p = p0,mean = -0,sd = 3)
  # ei=rnorm(n = n,mean = to_subtract,sd = 3)
  # 
  # ### For laplace distn-------------
  # library(rmutil)
  # to_subtract=qlaplace(p = p0, m=0, s=3)
  # ei=rlaplace(n = n, m = to_subtract, s = 3)
  # ###
  
  ##### For Mixture Normal distn-------------
  to_subt_func=function(tmp.mu){
    # tmp.mu=0
    g_indx=rbinom(n = 1e5,size = 1,prob = c(0.1))
    n1=sum(g_indx==0)
    n2=sum(g_indx==1)
    
    samp1=rnorm(n=n1,mean=tmp.mu,sd=1)
    samp2=rnorm(n=n2,mean=(tmp.mu+1),sd=5)
    
    samp=c(samp1,samp2)
    should.zero=as.numeric(quantile(samp,p0))
    should.zero
  }
  to_subtract=nleqslv(0.76,to_subt_func,method = 'Broyden')$x;to_subtract
  
  g_indx=rbinom(n = n,size = 1,prob = c(0.1))
  n1=sum(g_indx==0)
  n2=sum(g_indx==1)
  samp1=rnorm(n=n1,mean=to_subtract,sd=1)
  samp2=rnorm(n=n2,mean=(to_subtract+1),sd=5)
  samp=c(samp1,samp2)
  ei=samp
  #####
  
  y=(X%*%beta)+ei
  
  # set default--------------------------------------------------------------------------------
  niter    <- 30000
  nburn    <- 5000
  nthin    <- 5
  nprint   <- 10000
  nmcmc=(niter-nburn)/nthin
  
  #unif for gamma
  # L=-sqrt(1-p0)
  # U=sqrt(p0)
  
  U=nleqslv(0,tmp_func1)$x
  L=-nleqslv(0,tmp_func2)$x
  L;U
  
  #N for beta
  pr_mean_beta=rep(0,P)
  pr_sd_beta=100*diag(P)
  
  #IG for sigma
  aaa=2
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
  pr_sd_alpha=var_param*diag(2)
  
  
  beta_trace=matrix(NA,nrow=nmcmc,ncol=P)
  sigma_trace=rep(NA,nmcmc)
  gamma_trace=rep(NA,nmcmc)
  s_trace=matrix(NA,nrow=nmcmc,ncol=n)
  k_trace=matrix(NA,nrow=nmcmc,ncol=n)
  sigma2_xx_trace=rep(NA,nmcmc)
  sigma2_11_trace=rep(NA,nmcmc)
  sigma2_22_trace=rep(NA,nmcmc)
  mux_trace=rep(NA,nmcmc)
  alpha_trace=matrix(NA,nrow=nmcmc,ncol=2)
  X_trace=matrix(NA,nrow=nmcmc,ncol=n) #X,S,K is auxilary variable, but let's keep X only.
  
  ## z and s are auxilary variable
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  gamma_t=runif(1,L,U)
  # gamma_t=0.1937465
  beta_t=rmvnorm(1,c(0,0),diag(2))
  alpha_t=rnorm(2,0,1)
  sigma_t=runif(1)
  sigma2_xx_t=runif(1)
  sigma2_11_t=runif(1)
  sigma2_22_t=runif(1)
  mux_t=rnorm(1,pr_mean_mux,sqrt(pr_var_mux))
  s_t=rexp(n)
  k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = 0,sd = 1)
  
  accept_n=0
  
  
  X_1t=cbind(1,w1)%*%alpha_t #Let's start with 21, obs data
  X_t=cbind(1,X_1t)
  
  jump_gamma_sd=0.5
  
  tic()
  for(iter in 1:niter){
    # for(iter in 1:10){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Define A,B,C with current gamma_t
    g_gamma=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)
    p.gamma_t=ifelse(gamma_t>0,p0/g_gamma,1-(1-p0)/g_gamma)
    stopifnot(abs(p.gamma_t)<1)
    A=(1-2*p.gamma_t)/(p.gamma_t*(1-p.gamma_t));B=2/(p.gamma_t*(1-p.gamma_t));C=1/(ifelse(gamma_t>0,1,0)-p.gamma_t)
    
    v_t=s_t*sigma_t
    # alpha_t=abs(gamma_t)/(ifelse(gamma_t>0,1,0)-p)
    
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_sd_beta) + t(X/sqrt(v_t))%*%(X/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_sd_beta)%*%pr_mean_beta + colSums(X*((y-(sigma_t*C*abs(gamma_t)*k_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
    beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    # cat(beta_t,'\n')
    
    #sample_vi (or zi)-----------------------------------------------------------------
    a_t=(y - (X%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t))**2/(B*sigma_t)
    b_t=2/sigma_t+A**2/(B*sigma_t)
    for(ii in 1:n){
      v_t[ii]=rgig(n = 1,lambda = 0.5,chi =a_t[ii],psi = b_t)
    }
    # v_t=rgig(n = 1,lambda = 0.5,chi =a_t,psi = b_t)
    s_t=v_t/sigma_t
    
    #sample ki-----------------------------------------------------------------
    post_var_ki=1/((C*gamma_t)**2*sigma_t/(B*v_t)+1)
    post_mean_ki=post_var_ki*C*abs(gamma_t)*(y - (X%*%t(beta_t)+A*v_t))/(B*v_t)
    k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = post_mean_ki,sd = sqrt(post_var_ki)) #rtruncnorm can sample multiple case, when input param is vector
    
    #sample sigma-----------------------------------------------------------------
    nu=-(prior_a+1.5*n)
    c=2*prior_b + 2*sum(v_t) + sum((y - (X%*%t(beta_t)+A*v_t))**2/(B*v_t))
    d=sum((C*gamma_t*k_t)**2/(B*v_t))
    sigma_t=rgig(n = 1,lambda = nu,chi = c,psi = d)
    s_t=v_t/sigma_t #Is this right??? or should I just leave s_t as previous(i.e.unupdated) one?
    
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
    
    numer=sum(dnorm(x = y,mean = X%*%t(beta_t)+sigma_t*C.star*abs(gamma_star)*k_t+sigma_t*A.star*s_t,sd = sigma_t*sqrt(B.star*s_t),log = T))
    denom=sum(dnorm(x = y,mean = X%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t+sigma_t*A*s_t,sd = sigma_t*sqrt(B*s_t),log = T))
    
    ###########################
    # p.tmp=ifelse(gamma_t<0,1,0) + (p0-ifelse(gamma_t<0,1,0))/g_gamma
    # p.plus=p.tmp-ifelse(gamma_t>0,1,0)
    # p.minus=p-ifelse(gamma_t<0,1,0)
    # y.star=(y-X%*%t(beta_t))/sigma_t
    # 
    # (2*p.tmp*(1-p.tmp)/sigma_t*(pnorm(-p.plus*y.star[1]/abs(gamma_t)+p.minus/p.plus *abs(gamma_t))-pnorm(p.minus/p.plus *abs(gamma_t)))*exp(-p.minus*y.star[1]+gamma_t**2/2*(p.minus/p.plus)**2)*ifelse(y.star[1]/gamma_t>0,1,0)
    #     + pnorm(-abs(gamma_t)+p.plus*y.star[1]/abs(gamma_t)*ifelse(y.star[1]/gamma_t>0,1,0))*exp(-p.plus*y.star[1]+gamma_t**2/2))
    ###########################
    u2=log(runif(1))
    ratio=numer-denom
    u2<ratio
    if(u2<ratio){
      gamma_t=gamma_star
      A=A.star;B=B.star;C=C.star
      accept_n=accept_n+1
    } 
    else{
      gamma_t=gamma_t
    }
    # gamma_t=ifelse(u2<ratio,gamma_star,gamma_t) #This code won't able to count accept counts
    stopifnot((L<gamma_t)&(gamma_t<U))
    
    #sample X_1t-----------------------------------------------------------------
    V_X= 1 / (beta_t[2]^2/(sigma_t*B*v_t)+alpha_t[2]^2/sigma2_11_t+1/sigma2_22_t+1/sigma2_xx_t)
    #mean((1/V_X-(sigma2_11_t*sigma2_22_t*sigma2_xx_t*beta_t[2]^2+sigma_t*B*v_t*sigma2_22_t*sigma2_xx_t*alpha_t[2]^2+sigma_t*B*v_t*sigma2_11_t*sigma2_xx_t+sigma_t*B*v_t*sigma2_11_t*sigma2_22_t)/(sigma_t*B*v_t*sigma2_11_t*sigma2_22_t*sigma2_xx_t))<1e-10) # just for debugging
    M_X=((beta_t[2]*(y-(beta_t[1]+sigma_t*C*abs(gamma_t)*k_t+A*v_t)))/(sigma_t*B*v_t)+
           alpha_t[2]*(w1-alpha_t[1])/sigma2_11_t+
           w2/sigma2_22_t+
           mux_t/sigma2_xx_t)/
      (beta_t[2]^2/(sigma_t*B*v_t)+alpha_t[2]^2/sigma2_11_t+1/sigma2_22_t+1/sigma2_xx_t)
    # mean((M_X-(sigma2_11_t*sigma2_22_t*sigma2_xx_t*beta_t[2]*(y-(beta_t[1]+sigma_t*C*abs(gamma_t)*k_t+A*v_t))+
    #     sigma_t*B*v_t*sigma2_22_t*sigma2_xx_t*alpha_t[2]*(w1-alpha_t[1])+
    #     sigma_t*B*v_t*sigma2_11_t*sigma2_xx_t*w2+
    #     sigma_t*B*v_t*sigma2_11_t*sigma2_22_t*mux_t)/(sigma2_11_t*sigma2_22_t*sigma2_xx_t*beta_t[2]^2+sigma_t*B*v_t*sigma2_22_t*sigma2_xx_t*alpha_t[2]^2+sigma_t*B*v_t*sigma2_11_t*sigma2_xx_t+sigma_t*B*v_t*sigma2_11_t*sigma2_22_t))<1e-10) # just for debugging
    X_1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X_t=cbind(1,X_1t)
    
    #sample sigma2_xx-----------------------------------------------------------------
    post_ax=0.5*n+prior_ax
    # post_ax=1.5*n+prior_ax
    post_bx=prior_bx+1/2*sum((X_1t-mux_t)^2)
    sigma2_xx_t=rinvgamma(1,post_ax,post_bx)
    stopifnot(sigma2_xx_t>0)
    #cat('accept2\n')
    
    #sample sigma2_11-----------------------------------------------------------------
    post_a1=0.5*n+prior_a1
    # post_a1=1.5*n+prior_a1
    post_b1=prior_b1+1/2*sum((w1-(alpha_t[1]+alpha_t[2]*X_1t))^2)
    sigma2_11_t=rinvgamma(1,post_a1,post_b1)
    stopifnot(sigma2_11_t>0)
    #cat('accept3\n')
    
    #sample sigma2_22-----------------------------------------------------------------
    post_a2=0.5*n+prior_a2
    # post_a2=1.5*n+prior_a2
    post_b2=prior_b2+1/2*sum((w2-X_1t)^2)
    sigma2_22_t=rinvgamma(1,post_a2,post_b2)
    stopifnot(sigma2_22_t>0)
    #cat('accept4\n')
    
    #sample mux-----------------------------------------------------------------
    post_V_mux=1/(n/sigma2_xx_t+1/pr_var_mux)
    post_M_mux=(sum(X_1t)/sigma2_xx_t+pr_mean_mux/pr_var_mux)/(n/sigma2_xx_t+1/pr_var_mux)
    mux_t=rnorm(n = 1,mean = post_M_mux,sd = sqrt(post_V_mux))
    #cat('accept5\n')
    
    #sample alpha-----------------------------------------------------------------
    ### Can I sample 2 alphsa in Gibbs at the same time? i.e. Collapsed Gibbs is okay? Maybe it is okay! it satisfy the condition
    post_M_alpha=solve(1/sigma2_11_t*t(X_t)%*%X_t+solve(pr_sd_alpha))%*%(1/sigma2_11_t*t(X_t)%*%w1+solve(pr_sd_alpha)%*%pr_mean_alpha)
    post_V_alpha=solve(1/sigma2_11_t*t(X_t)%*%X_t+solve(pr_sd_alpha))
    alpha_t=mvrnorm(n = 1,mu = post_M_alpha,Sigma = post_V_alpha)
    #cat('accept6\n')
    #cat('alpha_t',alpha_t,'\n')
    #cat('------------------------------------------------------------------\n')
    
    
    
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta_t
      sigma_trace[thinned_idx]=sigma_t
      gamma_trace[thinned_idx]=gamma_t
      sigma2_xx_trace[thinned_idx]=sigma2_xx_t
      sigma2_11_trace[thinned_idx]=sigma2_11_t
      sigma2_22_trace[thinned_idx]=sigma2_22_t
      mux_trace[thinned_idx]=mux_t
      alpha_trace[thinned_idx,]=alpha_t
      X_trace[thinned_idx,]=X_1t
      s_trace[thinned_idx,]=s_t
      k_trace[thinned_idx,]=k_t
      
    }
  }
  toc()
  # save.image(file=sprintf('./debugging/QR_alpharandom_ii_%s.RData',sim_idx))
}

tic()
mean_bias=matrix(NA,ncol=6,nrow=nmax)
median_bias=matrix(NA,ncol=6,nrow=nmax)
for(sim_idx in 1:nmax){
  load(file=sprintf('../debugging/QR_alpharandom_ii_%s.RData',sim_idx))
  # load(file=sprintf('../debugging/QR_betarandom_%s.RData',sim_idx))
  aa=colMedians(beta_trace)
  bb=colMeans(beta_trace)
  median_bias[sim_idx,1:2]=(aa-beta)
  mean_bias[sim_idx,1:2]=(bb-beta)
  
  aa=colMedians(alpha_trace)
  bb=colMeans(alpha_trace)
  median_bias[sim_idx,3:4]=(aa-alpha)
  mean_bias[sim_idx,3:4]=(bb-alpha)
  
  median_bias[sim_idx,5]=median(sigma2_11_trace)-sigma2_11
  mean_bias[sim_idx,5]=mean(sigma2_11_trace)-sigma2_11
  
  median_bias[sim_idx,6]=median(sigma2_22_trace)-sigma2_22
  mean_bias[sim_idx,6]=mean(sigma2_22_trace)-sigma2_22
}
toc()

beta_hist=function(inp_data){
  label_list=c('beta0','beta1','alpha0','alpha1','sigma2_11','sigma2_22')
  par(mfrow=c(3,2))
  for(ii in 1:dim(inp_data)[2]){
    tmp_data=inp_data[,ii]
    tmp_data=tmp_data[!is.na(tmp_data)]
    quant=quantile(tmp_data,c(0.025,0.975))
    quant_data=tmp_data[(tmp_data>quant[1])&(tmp_data<quant[2])]
    p.v=round(t.test(quant_data, mu=0)$p.value,10)
    hist(quant_data,nclass=100,main=sprintf('%s bias\n t.test : %s\n mean : %s',label_list[ii],p.v,round(mean(quant_data),2)),xlab=label_list[ii])
    abline(v=0,col=2,lwd=3)  
  }
  par(mfrow=c(1,1))
}

# beta_hist(median_bias)
beta_hist(mean_bias)
