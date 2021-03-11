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

GAL_wo_ME=function(y,X,p0,seed=T){
    if(seed){
        set.seed(20210310)    
    }
    
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
    beta_t=rnorm(P,0,1)
    sigma_t=10#rinvgamma(1,prior_a,prior_b)
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
        v_t = 1/rinvgauss(n, mean = sqrt(b_t/a_t), dispersion = 1/ b_t) # which is same for rgig
        # v_t = 1/rinvGauss(n, nu = sqrt(b_t/a_t), lambda = b_t) # which is same for rgig!
        stopifnot(sum(is.na(v_t))==0)
        # for(ii in 1:n){
        #     v_t[ii]=rgig(n = 1,lambda = 0.5,chi =a_t[ii],psi = b_t)
        # }
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
    
    res_list=list()
    res_list[['beta_trace']]=beta_trace
    res_list[['sigma_trace']]=sigma_trace
    res_list[['gamma_trace']]=gamma_trace
    return(res_list)
}


ALD_wo_ME = function(y,X,p0){
    # set default--------------------------------------------------------------------------------
    niter    <- 30000
    nburn    <- 5000
    nthin    <- 5
    nprint   <- 10000
    nmcmc=(niter-nburn)/nthin
    
    P=dim(X)[2]
    n=dim(X)[1]
    
    pr_mean_beta=rep(0,P)
    pr_var_beta=100*diag(P)
    
    sigma.a=0.01
    sigma.b=0.01
    
    # Trace --------------------------------------------------------------------------------
    beta_trace=matrix(NA,nrow=nmcmc,ncol=P)
    sigma_trace=rep(NA,nmcmc)
    s_trace=matrix(NA,nrow=nmcmc,ncol=n)
    # Init val --------------------------------------------------------------------------------
    beta_t=rnorm(2,0,100)
    sigma_t=10
    s_t=rexp(n) # si is auxliary variable
    v_t=s_t*sigma_t
    
    #Define A,B
    stopifnot(abs(p0)<1)
    A=(1-2*p0)/(p0*(1-p0));B=2/(p0*(1-p0));
    
    # Start iter --------------------------------------------------------------------------------
    tic()
    for(iter in 1:niter){
        if(iter%%nprint==1){cat(iter,'th iter is going\n')}
        
        #Sample beta-----------------------------------------------------------------
        post_var_beta = solve(t(X/sqrt(v_t))%*%(X/sqrt(v_t))/(B*sigma_t) + solve(pr_var_beta))
        post_mean_beta = post_var_beta %*% (colSums(X*((y-(A*v_t))/(B*sigma_t*v_t))[1:n]) + solve(pr_var_beta)%*%pr_mean_beta)
        beta_t=rmvnorm(n = 1,mu = post_mean_beta,sigma = post_var_beta)
        
        #Sample vi----------------------------------------------------------------- 
        v_t.a = (y-X%*%t(beta_t))^2/(B*sigma_t)
        v_t.b = 2/sigma_t + A^2/(B*sigma_t)
        v_t = 1/rinvgauss(n, mean = sqrt(v_t.b/v_t.a), dispersion = 1/ v_t.b) # which is same for rgig
        # v_t = 1/rinvGauss(n, nu = sqrt(v_t.b/v_t.a), lambda = v_t.b) # which is same for rgig
        # for(ii in 1:n){
        #   v_t[ii]=rgig(n = 1,lambda = 0.5,chi = v_t.a[ii],psi = v_t.b)
        # }
        
        #Sample sigma-----------------------------------------------------------------  
        sigma.post.a = 1.5*n + sigma.a
        sigma.post.b = 0.5*sum((y - (X%*%t(beta_t)+A*v_t))^2/(B*v_t)) + sum(v_t) + sigma.b
        sigma_t = rinvgamma(n = 1,shape = sigma.post.a ,scale = sigma.post.b)
        
        # Save samples in each thin-----------------------------------------------------------------
        if((iter > nburn) & (iter %% nthin== 0) ){
            thinned_idx=(iter-nburn)/nthin
            beta_trace[thinned_idx,]=beta_t
            sigma_trace[thinned_idx]=sigma_t
            s_trace[thinned_idx,]=s_t
        }
    }
    toc()
    
    res_list=list()
    res_list[['beta_trace']]=beta_trace
    res_list[['sigma_trace']]=sigma_trace
    res_list[['s_trace']]=s_trace
    return(res_list)  
}


mBayesQR=function(y,X,p0,norm.approx=T){
    out <- bayesQR(y~X[,-1], quantile=c(p0), ndraw=500,normal.approx = norm.approx)
    
    res_list=list()
    res_list[['beta_trace']]=out[[1]]$betadraw
    res_list[['sigma_trace']]=out[[1]]$sigmadraw
    return(res_list)
}

# res=GAL_wo_ME(y,X,p0)
# bayesQR_res=mBayesQR(y,X,p0)
# 
# 
# colMeans(res[['beta_trace']])
# colMeans(bayesQR_res[['beta_trace']])
# beta.true
# 
# 
# mean(res[['sigma_trace']])
# mean(bayesQR_res[['sigma_trace']])

# # Check for inverse Gaussian vs Generalized inverse gamma---------------------------
# nmax=10000
# tic()
# chk_mat=matrix(NA,ncol=n,nrow=nmax)
# for(ii in 1:nmax){
#     chk_mat[ii,]=1/rinvGauss(n, nu = sqrt(b_t/a_t), lambda = b_t)
# }
# toc()
# 
# tic()
# chk_mat2=matrix(NA,ncol=n,nrow=nmax)
# for(ii in 1:nmax){
#     tmp=rep(NA,n)
#     for(jj in 1:n){
#         tmp[jj]=rgig(n = 1,lambda = 0.5,chi = a_t[jj],psi = b_t)
#     }
#     chk_mat2[ii,]=tmp
# }
# toc()
# 
# library(statmod)
# tic()
# chk_mat3=matrix(NA,ncol=n,nrow=nmax)
# for(ii in 1:nmax){
#     chk_mat3[ii,]=1/rinvgauss(n, mean = sqrt(b_t/a_t), dispersion = 1/ b_t)
# }
# toc()
# 
# colMeans(chk_mat)
# colMeans(chk_mat2)
# colMeans(chk_mat3)
# 
# colVars(chk_mat)
# colVars(chk_mat2)
# colVars(chk_mat3)