rm(list = ls())
setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
# library(mvtnorm)
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)
library(Brq)
library(quantreg)
library(bayesQR)

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

GAL_QR=function(y,X,p0){
  # set.seed(20210305)
  
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
  P=dim(X)[2]
  n=dim(X)[1]
  pr_mean_beta=rep(0,P)
  pr_sd_beta=100*diag(P)
  
  #IG for sigma
  prior_a=2
  prior_b=2
  # prior_a=3/2
  # prior_b=0.1/2
  
  beta_trace=matrix(NA,nrow=nmcmc,ncol=P)
  sigma_trace=rep(NA,nmcmc)
  gamma_trace=rep(NA,nmcmc)
  ## z and s are auxilary variable
  
  # Iteration--------------------------------------------------------------------------------
  
  #init value
  gamma_t=runif(1,L,U)
  # gamma_t=0.1937465
  beta_t=rnorm(3,0,1)
  # sigma_t=rinvgamma(1,prior_a,prior_b)
  sigma_t=1
  z_t=rexp(n)
  
  
  s_t=rtruncnorm(n = n,a = 0,b = Inf,mean = 0,sd = 1)
  
  accept_n=0
  tic()
  for(iter in 1:niter){
    # for(iter in 1:10){
    if(iter%%nprint==1){cat(iter,'th iter is going\n')}
    
    #Define A,B,C with current gamma_t
    g_gamma=2*pnorm(-abs(gamma_t))*exp(gamma_t**2/2)
    p.gamma_t=ifelse(gamma_t>0,p0/g_gamma,1-(1-p0)/g_gamma)
    stopifnot(abs(p.gamma_t)<1)
    A=(1-2*p.gamma_t)/(p.gamma_t*(1-p.gamma_t));B=2/(p.gamma_t*(1-p.gamma_t));C=1/(ifelse(gamma_t>0,1,0)-p.gamma_t)
    
    v_t=z_t*sigma_t
    # alpha_t=abs(gamma_t)/(ifelse(gamma_t>0,1,0)-p)
    
    
    #Sample beta-----------------------------------------------------------------
    post_sd_beta = solve( solve(pr_sd_beta) + t(X/sqrt(v_t))%*%(X/sqrt(v_t))/(B*sigma_t) )
    post_mean_beta = post_sd_beta %*% (solve(pr_sd_beta)%*%pr_mean_beta + colSums(X*((y-(sigma_t*C*abs(gamma_t)*s_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
    beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
    # cat(beta_t,'\n')
    
    #sample_vi (or zi)-----------------------------------------------------------------
    a_t=(y - (X%*%t(beta_t)+sigma_t*C*abs(gamma_t)*s_t))**2/(B*sigma_t)
    b_t=2/sigma_t+A**2/(B*sigma_t)
    for(ii in 1:n){
      v_t[ii]=rgig(n = 1,lambda = 0.5,chi =a_t[ii],psi = b_t)
    }
    # v_t=rgig(n = 1,lambda = 0.5,chi =a_t,psi = b_t)
    z_t=v_t/sigma_t
    
    #sample si-----------------------------------------------------------------
    post_var_si=1/((C*gamma_t)**2*sigma_t/(B*v_t)+1)
    post_mean_si=post_var_si*C*abs(gamma_t)*(y - (X%*%t(beta_t)+A*v_t))/(B*v_t)
    s_t=rtruncnorm(n = n,a = 0,b = Inf,mean = post_mean_si,sd = sqrt(post_var_si)) #rtruncnorm can sample multiple case, when input param is vector
    
    #sample sigma-----------------------------------------------------------------
    nu=-(prior_a+1.5*n)
    c=2*prior_b + 2*sum(v_t) + sum((y - (X%*%t(beta_t)+A*v_t))**2/(B*v_t))
    d=sum((C*gamma_t*s_t)**2/(B*v_t))
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
    
    numer=sum(dnorm(x = y,mean = X%*%t(beta_t)+sigma_t*C.star*abs(gamma_star)*s_t+sigma_t*A.star*z_t,sd = sigma_t*sqrt(B.star*z_t),log = T))
    denom=sum(dnorm(x = y,mean = X%*%t(beta_t)+sigma_t*C*abs(gamma_t)*s_t+sigma_t*A*z_t,sd = sigma_t*sqrt(B*z_t),log = T))
    
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
    
    # Save samples in each thin-----------------------------------------------------------------
    if((iter > nburn) & (iter %% nthin== 0) ){
      thinned_idx=(iter-nburn)/nthin
      beta_trace[thinned_idx,]=beta_t
      sigma_trace[thinned_idx]=sigma_t
      gamma_trace[thinned_idx]=gamma_t
      # cat(beta_t,gamma_t,sigma_t,'for',iter,'th iter\n')
    }
  }
  toc()
  return(colMeans(beta_trace))
}

# ####Check about location shift model
# set.seed(1209)
# beta0=0
# beta1=2
# gamma0=0
# gamma1=0.5
# n=1000
# xi=rnorm(n,4,1)
# ei=rnorm(n)
# y=beta0+beta1*xi+(gamma0+gamma1*xi)*ei
# X=cbind(1,xi)
# plot(xi,y)

# Make Data----------------------------------
# ####
# data("ImmunogG")
# head(ImmunogG)
# 
# y=ImmunogG$IgG
# X=cbind(1,ImmunogG$Age,ImmunogG$Age^2)

### BMQR for linear models data-----------------------------------
tic()
nmax=1000
beta_save1=matrix(NA,ncol=2,nrow=nmax)
beta_save2=matrix(NA,ncol=2,nrow=nmax)
beta_save3=matrix(NA,ncol=2,nrow=nmax)
for(idx in 1:nmax){
  # make data ------------------------------
  set.seed(idx)
  n=200
  
  xi=seq(0,20,length.out = n)
  ei=rnorm(n,0,sd=(1+0.5*xi))
  yi=5+2*xi+ei
  # plot(xi,yi)
  
  X=cbind(1,xi)
  y=yi
  # set quantile ------------------------------
  
  p0=0.5
  # hist(ei,nclass = 100);
  # Q.tau=quantile(ei,p0)
  # true_beta=beta+c(Q.tau,0,0);true_beta
  
  # estimation starts ------------------------------
  ###
  beta.est=GAL_QR(y,X,p0)
  beta_save1[idx,]=beta.est
  
  ###
  # res=rq(y~X[,-1],p0)
  # beta_save2[idx,]=res$coefficients
  # 
  # ###
  # out <- bayesQR(y~X[,-1], quantile=c(p0), ndraw=500)
  # sum <- summary(out, burnin=50)
  # beta_save3[idx,]=sum[[1]]$betadraw[,1]
}
toc()

true_beta=c(5,2)

# beta_save1.adj=beta_save1[-which(abs(beta_save1[,1])>35),]

df0.5=data.frame(matrix(NA,nrow=3,ncol=((dim(X)[2]*3+2))),row.names = c('GAL','QR','BayesQR'))
# df0.5[1,]=c((colMeans(beta_save1.adj)-true_beta)/1, NA, (colVars(beta_save1.adj)), NA, (colMeans(beta_save1.adj)-true_beta)^2+ colVars(beta_save1.adj))
df0.5[1,]=c((colMeans(beta_save1)-true_beta)/1, NA, (colVars(beta_save1)), NA, (colMeans(beta_save1)-true_beta)^2+ colVars(beta_save1))
df0.5[2,]=c((colMeans(beta_save2)-true_beta)/1, NA, (colVars(beta_save2)), NA, (colMeans(beta_save2)-true_beta)^2+ colVars(beta_save2))
df0.5[3,]=c((colMeans(beta_save3)-true_beta)/1, NA, (colVars(beta_save3)), NA, (colMeans(beta_save3)-true_beta)^2+ colVars(beta_save3))

colnames(df0.5)=c(paste0('bias',1:dim(X)[2]),NA,paste0('var',1:dim(X)[2]),NA,paste0('mse',1:dim(X)[2]))
# save.image(file='../debugging/Debugging_w_BMQR_data.RData')
# load(file='../debugging/Debugging_w_BMQR_data.RData')
df0.5[2,]

par(mfrow=c(2,2))
hist(beta_save1[,1],nclass=30,main='GAL')
# hist(beta_save1.adj[,1],nclass=30)
hist(beta_save2[,1],nclass=30,main='quantreg')
hist(beta_save3[,1],nclass=30,main='BayesQR')
par(mfrow=c(1,1))