rm(list = ls())
setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
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
get_validmat=function(mat){
  mat=mat[!is.na(mat[,1]),]  
  return (mat)
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
gen_error=function(type,n,x2){
  if(type=='N'){
    ei=rnorm(n,0,1)
  }
  if(type=='t'){
    ei=rt(n = n,df = 3)
  }
  if(type=='locshift'){
    ei=(1+x2)*rnorm(n,0,1)
  }
  return (ei)
}

gen_truebeta=function(type,p0){
  if(type=='N'){
    to_add=c(qnorm(p0),0,0)
  }
  if(type=='t'){
    to_add=c(qt(p = p0,df = 3),0,0)
  }
  if(type=='locshift'){
    to_add=c(qnorm(p0),qnorm(p0),0)
  }
  return (to_add)  
}


total_df=list()
P=3
# type='N'
tic()
for (type in c('N','t','locshift')){
  df=list()
  for(p0 in c(0.1,0.5,0.9)){
    # nmax=1000
    nmax=100
    beta_save1=matrix(NA,ncol=P,nrow=nmax)
    beta_save2=matrix(NA,ncol=P,nrow=nmax)
    beta_save3=matrix(NA,ncol=P,nrow=nmax)
    for(idx in 1:nmax){
      tryCatch({
        # make data ------------------------------
        set.seed(idx)
        n=100
        
        x2=rnorm(n,0,1)
        x3=rnorm(n,0,1)
        X=cbind(1,x2,x3)
        stopifnot(P==dim(X)[2])
        ei=gen_error(type=type,n=n,x2=x2)
        beta=c(1,1,1)
        y=X%*%beta+ei
        
        # set quantile ------------------------------
        # p0=0.5
        # hist(ei,nclass = 100);
        # Q.tau=quantile(ei,p0)
        # true_beta=beta+c(Q.tau,0,0);true_beta
        
        # estimation starts ------------------------------
        ###
        GAL_res=GAL_wo_ME(y,X,p0)
        beta_save1[idx,]=colMeans(GAL_res[['beta_trace']])
        
        ###
        res=rq(y~X[,-1],p0)
        beta_save2[idx,]=res$coefficients

        ###
        BayesQR_res=mBayesQR(y,X,p0)
        beta_save3[idx,]=colMeans(BayesQR_res[['beta_trace']])
      }, error= function(e) {cat("Error", "\n")},
      warning= function(w) {cat("Warning", "\n")})      
    }
    
    true_beta=beta+gen_truebeta(type = type,p0 = p0)

    # hist(beta_save1[,1],nclass=30)
    # hist(beta_save2[,1],nclass=30)
    # hist(beta_save3[,1],nclass=30)
    name=paste0('df',p0)
    df[[name]]=list()
    df[[name]][['GAL']]=beta_save1
    df[[name]][['QR']]=beta_save2
    df[[name]][['BayesQR']]=beta_save3
    
    # df[[name]]=data.frame(matrix(NA,nrow=3,ncol=((dim(X)[2]*3+2))),row.names = c('GAL','QR','BayesQR'))
    # df[[name]][1,]=c((colMeans(beta_save1)-true_beta)/1, NA, sqrt(colVars(beta_save1)), NA, sqrt((colMeans(beta_save1)-true_beta)^2+ colVars(beta_save1)))
    # df[[name]][2,]=c((colMeans(beta_save2)-true_beta)/1, NA, sqrt(colVars(beta_save2)), NA, sqrt((colMeans(beta_save2)-true_beta)^2+ colVars(beta_save2)))
    # df[[name]][3,]=c((colMeans(beta_save3)-true_beta)/1, NA, sqrt(colVars(beta_save3)), NA, sqrt((colMeans(beta_save3)-true_beta)^2+ colVars(beta_save3)))
    # colnames(df[[name]])=c(paste0('bias',1:dim(X)[2]),NA,paste0('sd',1:dim(X)[2]),NA,paste0('rmse',1:dim(X)[2]))
  }  
  total_df[[type]]=df
}
toc()

# save.image(file='../debugging/Debugging_w_GibbsBQR_data.RData')
# load(file='../debugging/Debugging_w_GibbsBQR_data.RData')
beta_svae.model=type_df[['df0.1']][['BayesQR']]

var(beta_svae.model[,1],na.rm = T)
colVars(beta_svae.model[-216,],na.rm = T)
colVars(beta_svae.model[,],na.rm = T)

for (type in names(total_df)){
  type_df=total_df[[type]]
  for (qversion in names(type_df)){
    p0=as.numeric(substr(qversion,start = (nchar(qversion)-2),stop = nchar(qversion)))
    true_beta=beta+gen_truebeta(type = type,p0 = p0)
    
    q_type_df=type_df[[qversion]]
    df=data.frame(matrix(NA,nrow=3,ncol=((dim(X)[2]*3+2))),row.names = c('GAL','QR','BayesQR'))
    cnt=1
    for (mversion in names(q_type_df)){
      beta_svae.model = q_type_df[[mversion]]
      beta_svae.model[which(is.nan(beta_svae.model[,1])),]=rep(NA,3)
      df[cnt,]=c((colMeans(beta_svae.model,na.rm = T)-true_beta)/1, NA, sqrt(colVars(beta_svae.model,na.rm = T)), NA, sqrt((colMeans(beta_svae.model,na.rm = T)-true_beta)^2+ colVars(beta_svae.model,na.rm = T)))
      cnt=cnt+1
    }
    colnames(df)=c(paste0('bias',1:dim(X)[2]),NA,paste0('sd',1:dim(X)[2]),NA,paste0('rmse',1:dim(X)[2]))
    print(c(type,p0))
    print(df)
  }
}
beta_svae.model


type='locshift'
n
p0=0.1
# set.seed(10)
x2=rnorm(n,0,1)
# x3=rnorm(n,0,1)
# X=cbind(1,x2,x3)
X=cbind(1,x2,x2^2)
ei=gen_error(type=type,n=n,x2=x2)
beta=c(1,1,1)
y=X%*%beta+ei
p0=0.55
true_beta=beta+gen_truebeta(type = type,p0 = p0)
plot(X[,2],y)
X_range=seq(min(X[,2]),max(X[,2]),length.out = 100)
y.hat=true_beta[1]+true_beta[2]*X_range+true_beta[3]*X_range^2
lines(X_range,y.hat)

abline(true_beta)
