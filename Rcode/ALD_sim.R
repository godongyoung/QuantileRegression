#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
type = args[1]
p0= as.numeric(args[2])
start.idx = as.numeric(args[3])
end.idx = as.numeric(args[4])
print(sprintf("type:%s / p0:%s / start:%s / end:%s",type,p0,start.idx,end.idx))

# rm(list = ls())
# tryCatch( # set wording directory to current source file folder
#   {setwd(getSrcDirectory()[1])}
#   ,error=function(e){setwd(dirname(rstudioapi::getActiveDocumentContext()$path))}
# )
# setwd('C:/Users/admin/Documents/dongyoung/Rcode')
setwd('D:/dy_result/QuantileRegression/Rcode')
# setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
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

# Simulation start --------------------------------------------------------------------------------
for(sim_idx in start.idx:end.idx){
  set.seed(sim_idx)
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
  
  y=(X%*%beta)+ei
  
  # Calculate True beta based on input types & p0
  {if (type=='N'){
    beta.true=beta+c(qnorm(p0,0,1),0)
  }
    if (type=='t'){
      beta.true=beta+c(qt(p0,df = 3),0)
    }
    if (type=='locshit'){
      beta.true=c(beta[1]+qnorm(p0,0,1),beta[2]*qnorm(p0,0,1))
    }  
    if (type=='MixN'){
      beta.true=beta
    }}
  
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
  pr_mean_beta=rep(0,P)
  pr_sd_beta=100*diag(P)
  
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
  pr_sd_alpha=var_param*diag(2)
  
  # Make Trace ----------------------------------------------------------------------
  beta_trace=matrix(NA,nrow=nmcmc,ncol=P)
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
  beta_t=rmvnorm(1,c(0,0),diag(2))
  alpha_t=rnorm(2,0,1)
  sigma_t=10#runif(1)
  sigma2_xx_t=10#runif(1)
  sigma2_11_t=10#runif(1)
  sigma2_22_t=10#runif(1)
  mux_t=rnorm(1,pr_mean_mux,sqrt(pr_var_mux))
  s_t=rexp(n)
  
  accept_n=0
  
  
  X_1t=cbind(1,w1)%*%alpha_t #Let's start with 21, obs data
  X_t=cbind(1,X_1t)
  
  
  

  
  
  v_t=s_t*sigma_t
  tic()
  for(iter in 1:niter){
    # for(iter in 1:10){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    # alpha_t=abs(gamma_t)/(ifelse(gamma_t>0,1,0)-p)
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_sd_beta) + t(X_t/sqrt(v_t))%*%(X_t/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_sd_beta)%*%pr_mean_beta + colSums(X_t*((y-(A*v_t))/(B*sigma_t*v_t))[1:n]))
    beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    
    #sample_vi (or zi)-----------------------------------------------------------------
    a_t=(y - (X_t%*%t(beta_t)))**2/(B*sigma_t)
    b_t=2/sigma_t+A**2/(B*sigma_t)
    # v_t=sapply(X = a_t, FUN = function(x) rgig(n = 1,lambda = 0.5, chi = x, psi = b_t))
    for(ii in 1:n){
      v_t[ii]=rgig(n = 1,lambda = 0.5,chi =a_t[ii],psi = b_t)
    }
    # v_t=rgig(n = 1,lambda = 0.5,chi =a_t,psi = b_t) ### error version
    s_t=v_t/sigma_t
    
    #sample sigma-----------------------------------------------------------------
    nu=(prior_a+1.5*n)
    c=0.5 * (2*prior_b + 2*sum(v_t) + sum((y - (X_t%*%t(beta_t)+A*v_t))**2/(B*v_t)))
    sigma_t=rinvgamma(n = 1,shape = nu ,scale = c)
    s_t=v_t/sigma_t #Is this right??? or should I just leave s_t as previous(i.e.unupdated) one?
    

    #sample X_1t-----------------------------------------------------------------
    V_X= 1 / (beta_t[2]^2/(sigma_t*B*v_t)+alpha_t[2]^2/sigma2_11_t+1/sigma2_22_t+1/sigma2_xx_t)
    M_X=((beta_t[2]*(y-(beta_t[1]+A*v_t)))/(sigma_t*B*v_t)+
           alpha_t[2]*(w1-alpha_t[1])/sigma2_11_t+
           w2/sigma2_22_t+
           mux_t/sigma2_xx_t)/
      (beta_t[2]^2/(sigma_t*B*v_t)+alpha_t[2]^2/sigma2_11_t+1/sigma2_22_t+1/sigma2_xx_t)
    X_1t=rnorm(n = n,mean = M_X,sd = sqrt(V_X))
    X_t=cbind(1,X_1t)
    
    #sample sigma2_xx-----------------------------------------------------------------
    post_ax=0.5*n+prior_ax
    post_bx=prior_bx+1/2*sum((X_1t-mux_t)^2)
    sigma2_xx_t=rinvgamma(1,post_ax,post_bx)
    stopifnot(sigma2_xx_t>0)
    #cat('accept2\n')
    
    #sample sigma2_11-----------------------------------------------------------------
    post_a1=0.5*n+prior_a1
    post_b1=prior_b1+1/2*sum((w1-(alpha_t[1]+alpha_t[2]*X_1t))^2)
    sigma2_11_t=rinvgamma(1,post_a1,post_b1)
    stopifnot(sigma2_11_t>0)
    #cat('accept3\n')
    
    #sample sigma2_22-----------------------------------------------------------------
    post_a2=0.5*n+prior_a2
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
      sigma2_xx_trace[thinned_idx]=sigma2_xx_t
      sigma2_11_trace[thinned_idx]=sigma2_11_t
      sigma2_22_trace[thinned_idx]=sigma2_22_t
      mux_trace[thinned_idx]=mux_t
      alpha_trace[thinned_idx,]=alpha_t
      X_trace[thinned_idx,]=X_1t
      s_trace[thinned_idx,]=s_t
      
    }
  }
  toc()
  save.image(file=sprintf('../debugging/ALD_%s_%s_%s.RData',type,p0,sim_idx))
}

nmax=10
sim_idx=7
tic()
mean_bias=matrix(NA,ncol=6,nrow=nmax)
median_bias=matrix(NA,ncol=6,nrow=nmax)
summary_mat=matrix(NA,ncol=6,nrow=nmax)
for(sim_idx in 1:nmax){
  load(file=sprintf('../debugging/ALD_%s_%s_%s.RData',type,p0,sim_idx))
  aa=colMedians(beta_trace)
  bb=colMeans(beta_trace)
  median_bias[sim_idx,1:2]=(aa-beta.true)
  mean_bias[sim_idx,1:2]=(bb- beta.true)
  summary_mat[sim_idx,1:2]=bb

  aa=colMedians(alpha_trace)
  bb=colMeans(alpha_trace)
  median_bias[sim_idx,3:4]=(aa-alpha)
  mean_bias[sim_idx,3:4]=(bb-alpha)
  summary_mat[sim_idx,3:4]=bb

  median_bias[sim_idx,5]=median(sigma2_11_trace)-sigma2_11
  mean_bias[sim_idx,5]=mean(sigma2_11_trace)-sigma2_11

  median_bias[sim_idx,6]=median(sigma2_22_trace)-sigma2_22
  mean_bias[sim_idx,6]=mean(sigma2_22_trace)-sigma2_22
}
toc()
colMeans(summary_mat[1:10,1:4])
colMedians(summary_mat[1:10,1:4])
colVars(summary_mat[1:10,1:4])
beta.true
# 
# beta_hist=function(inp_data){
#   label_list=c('beta0','beta1','alpha0','alpha1','sigma2_11','sigma2_22')
#   par(mfrow=c(3,2))
#   for(ii in 1:dim(inp_data)[2]){
#     tmp_data=inp_data[,ii]
#     tmp_data=tmp_data[!is.na(tmp_data)]
#     quant=quantile(tmp_data,c(0.025,0.975))
#     quant_data=tmp_data[(tmp_data>quant[1])&(tmp_data<quant[2])]
#     p.v=round(t.test(quant_data, mu=0)$p.value,10)
#     hist(quant_data,nclass=100,main=sprintf('%s bias\n t.test : %s\n mean : %s',label_list[ii],p.v,round(mean(quant_data),2)),xlab=label_list[ii])
#     abline(v=0,col=2,lwd=3)
#   }
#   par(mfrow=c(1,1))
# }
# 
# beta_hist(median_bias)
# beta_hist(mean_bias)
