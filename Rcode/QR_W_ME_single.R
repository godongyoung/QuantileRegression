rm(list = ls())
setwd('C:/MyFiles/Research/quantile regression/')
# library(mvtnorm)
library(MCMCpack)
library(truncnorm)    
library(nleqslv)
library(tictoc)
library(GIGrvg)
library(Rfast)

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


# Make case1 data--------------------------------------------------------------------------------
# set.seed(20200827)
# p0=0.25
# n=1000
# P=3
# 
# covmat=matrix(NA,P,P)
# for(i in 1:P){
#   for(j in 1:P){
#     covmat[i,j]=0.5**(abs(i-j))
#   }
# }
# 
# 
# X=rmvnorm(n = n,mu = rep(0,(P)),sigma = covmat);dim(X)
# colMeans(X);cor(X)
# X=cbind(rep(1,n),X)
# 
# total_P=P+1
# beta=c(3,1.5,2,1,2,3,4,5,6,7,8,9,10)[1:total_P]
# 
# to_subtract=qnorm(p = p0,mean = -0,sd = 3)
# ei=rnorm(n = n,mean = to_subtract,sd = 3)
# y=(X%*%beta)+ei
# 
# plot(X%*%beta,y)

# Make Case2 data--------------------------------------------------------------------------------
set.seed(20200827)
n=1000
T.P=1
total_T.P=T.P+1
# t <- as.matrix(mvrnorm(n=n,mu = c(3,5),Sigma = 1*diag(T.P)));sum(t)
t <- matrix(runif(T.P*n,min=0, max=10),ncol=T.P)
t <- cbind(1,t)
alpha=as.matrix(sample(1:10,total_T.P,replace = T)/2);alpha
sigma_delta=2
X<-t%*%alpha+rnorm(n,0,sd=sigma_delta);summary(X)
if(min(X)<0){
    alpha[1]=alpha[1]+abs(min(X)) # shift using true intercept
    X<-t%*%alpha+rnorm(n,0,sd=sigma_delta);summary(X)
}

# X <- runif(n=n,min=0,max=10)
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)
plot(X,y)
X=cbind(rep(1,n),X)
# rm(X)

total_X.P=2
p0=0.25
# set default--------------------------------------------------------------------------------
# p0=0.95
for(p0 in c(0.05,0.25,0.5,0.75,0.95)){
    # for(p0 in c(0.5)){
    times=4
    niter    <- 30000*times
    nburn    <- 5000*times
    nthin    <- 5*times
    nprint   <- 1000*times
    nmcmc=(niter-nburn)/nthin
    
    #unif for gamma
    U=nleqslv(0,tmp_func1)$x
    L=-nleqslv(0,tmp_func2)$x
    L;U
    
    #N for beta
    pr_mean_beta=rep(0,total_X.P)
    pr_sd_beta=3*diag(total_X.P)
    
    #IG for sigma
    prior_a=3
    prior_b=2
    
    #IG for sigma_delta
    prior_c=3
    prior_d=2
    
    #N for alpha
    pr_mean_alpha=rep(0,total_T.P)
    pr_sd_alpha=3*diag(total_T.P)
    
    
    beta_trace=matrix(NA,nrow=nmcmc,ncol=total_X.P)
    sigma_trace=rep(NA,nmcmc)
    gamma_trace=rep(NA,nmcmc)
    X_trace=matrix(NA,nrow=nmcmc,ncol=n) #X,S,K is auxilary variable
    alpha_trace=matrix(NA,nrow=nmcmc,ncol=total_T.P)
    sigma_delta_trace=rep(NA,nmcmc)
    
    # Iteration--------------------------------------------------------------------------------
    
    #init value
    gamma_t=runif(1,L,U)
    # gamma_t=0.1937465
    beta_t=rnorm(total_X.P,0,1)
    alpha_t=rnorm(total_T.P,0,1)
    sigma_t=rinvgamma(1,prior_a,prior_b)
    sigma_delta_t=rinvgamma(1,prior_c,prior_d)
    s_t=rexp(n)
    k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = 0,sd = 1)
    X_1t=t%*%alpha_t
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
        # alpha_t=abs(gamma_t)/(ifelse(gamma_t>0,1,0)-p)
        
        
        #Sample beta-----------------------------------------------------------------
        post_sd_beta = solve( solve(pr_sd_beta) + t(X_t/sqrt(v_t))%*%(X_t/sqrt(v_t))/(B*sigma_t) )
        post_mean_beta = post_sd_beta %*% (solve(pr_sd_beta)%*%pr_mean_beta + colSums(X_t*((y-(sigma_t*C*abs(gamma_t)*k_t+A*v_t))/(B*sigma_t*v_t))[1:n]))
        beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_sd_beta)
        
        #sample ki-----------------------------------------------------------------
        post_var_ki=1/((C*gamma_t)**2*sigma_t/(B*v_t)+1)
        post_mean_ki=post_var_ki*C*abs(gamma_t)*(y - (X_t%*%t(beta_t)+A*v_t))/(B*v_t)
        k_t=rtruncnorm(n = n,a = 0,b = Inf,mean = post_mean_ki,sd = sqrt(post_var_ki)) #rtruncnorm can sample multiple case, when input param is vector
        
        #sample_vi (or si)-----------------------------------------------------------------
        a_t=(y - (X_t%*%t(beta_t)+sigma_t*C*abs(gamma_t)*k_t))**2/(B*sigma_t)
        b_t=2/sigma_t+A**2/(B*sigma_t)
        for(k in 1:n){
            v_t[k]=rgig(n = 1,lambda = 0.5,chi =a_t[k],psi = b_t)
        }
        s_t=v_t/sigma_t
        #cat('v_t',v_t[1:3],'\n')
        
        #sample X_1t-----------------------------------------------------------------
        V_X = (sigma_t**2*B*s_t*sigma_delta_t**2)/((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
        # 1/(beta_t[2]**2/(sigma_t**2*B*s_t)+1/sigma_delta_t**2)
        M_X = (sigma_delta_t**2*(y-beta_t[1]-sigma_t*C*abs(gamma_t)*k_t-sigma_t*A*s_t)*beta_t[2]+sigma_t**2*B*s_t*t%*%alpha_t) /((sigma_delta_t*beta_t[2])**2+sigma_t**2*B*s_t)
        #cat('M_X',M_X[1:3],'\n')
        X_1t=rnorm(n = n,mean = M_X,sd = V_X)
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
        # post_c=prior_c+0.5*n
        post_c=prior_c+1.5*n
        post_d=0.5*sum((X_1t-t%*%alpha_t)**2)+prior_d
        sigma_delta_t=sqrt(rinvgamma(1,post_c,post_d))
        #cat('sigma_delta_t',sigma_delta_t,'\n')
        
        #sample alpha-----------------------------------------------------------------
        V_delta=solve(1/sigma_delta_t**2*(t(t)%*%t+sigma_delta_t**2*solve(pr_sd_alpha)))
        M_delta=solve(t(t)%*%t+sigma_delta_t**2*solve(pr_sd_alpha))%*%t(X_1t%*%t+sigma_delta_t**2*t(pr_mean_alpha)%*%solve(pr_sd_alpha))
        alpha_t=mvrnorm(n = 1,mu = M_delta,Sigma = V_delta)
        #cat('alpha_t',alpha_t,'\n')
        #cat('------------------------------------------------------------------\n')
        
        # Save samples in each thin-----------------------------------------------------------------
        if((iter > nburn) & (iter %% nthin== 0) ){
            thinned_idx=(iter-nburn)/nthin
            beta_trace[thinned_idx,]=beta_t
            sigma_trace[thinned_idx]=sigma_t
            gamma_trace[thinned_idx]=gamma_t
            X_trace[thinned_idx,]=X_1t
            alpha_trace[thinned_idx,]=alpha_t
            sigma_delta_trace[thinned_idx]=sigma_delta_t
            # cat(beta_t,gamma_t,sigma_t,'for',iter,'th iter\n')
        }
    }
    toc()
    
    # ts.plot(beta_trace[,1])
    # plot(X[,2],y)
    # plot(X_t%*%colMeans(beta_trace[1000:which(is.na(beta_trace[,1]))[1],],na.rm = T),y)
    # abline(a=0,b=1)
    # 
    # plot(t[,2],y)
    # plot(t[,2],X[,2])
    # abline(a=0,b=1)
    
    # save(beta_trace,X_trace,alpha_trace,file=sprintf('Measurement_Trace_Single_%s.RData',p0))
    # save(beta_trace,X_trace,alpha_trace,file=sprintf('Measurement_Trace_Single_%s_iterlong.RData',p0))
}



plot(t[,2], y, main="", cex=.6, xlab="t")
for(p0 in c(0.05,0.25,0.5,0.75,0.95)){
    load(file=sprintf('Measurement_Trace_Single_%s.RData',p0))
    # beta.est=colMeans(beta_trace[1:(which(is.na(beta_trace[,1]))[1]-1),])
    # alpha.est=colMeans(alpha_trace[1:(which(is.na(beta_trace[,1]))[1]-1),])
    
    alpha.est=colMeans(alpha_trace[,])
    beta.est=colMeans(beta_trace[,])
    abline(a=beta.est[1]+beta.est[2]*alpha.est[1],b=beta.est[2]*alpha.est[2])
    # abline(a=beta.est[1],beta.est[2])
    Sys.sleep(1)
}
for(p0 in c(0.05,0.25,0.5,0.75,0.95)){
    load(file=sprintf('Measurement_Trace_Single_%s.RData',p0))
    par(mfrow=c(2,2))
    ts.plot(alpha_trace[,1])
    ts.plot(alpha_trace[,2])
    ts.plot(beta_trace[,1])
    ts.plot(beta_trace[,2])
    par(mfrow=c(1,1))
}

par(mfrow=c(3,2))
for(p0 in c(0.05,0.25,0.5,0.75,0.95)){
    load(file=sprintf('Measurement_Trace_Single_%s.RData',p0))
    # X.est=colMeans(X_trace[1:(which(is.na(beta_trace[,1]))[1]-1),])
    X.est=colMeans(X_trace[,])
    plot(X.est,X[,2],main=p0)
    abline(a=0,b=1)
    
}
par(mfrow=c(1,1))
