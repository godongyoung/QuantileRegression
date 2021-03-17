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
beta_save3.2=matrix(NA,ncol=2,nrow=nmax)

for(idx in 1:nmax){
  # make data ------------------------------
  set.seed(idx)
  n=100
  
  xi=seq(0,20,length.out = n)
  ei=rnorm(n,0,sd=(1+0.5*xi))
  yi=5+2*xi+ei
  
  X=cbind(1,xi)
  y=yi
  
  # set quantile ------------------------------
  p0=0.5
  
  # estimation starts ------------------------------
  ###
  GAL_res=GAL_wo_ME(y,X,p0)
  beta_save1[idx,]=colMeans(GAL_res[['beta_trace']])
  
  ###
  res=rq(y~X[,-1],p0)
  beta_save2[idx,]=res$coefficients

  ###
  BQR_res = mBayesQR(y,X,p0)
  beta_save3[idx,]=colMeans(BQR_res[['beta_trace']])
  
  ###
  ALD_res = ALD_wo_ME(y,X,p0)
  beta_save3.2[idx,]=colMeans(ALD_res[['beta_trace']])
}
toc()

true_beta=c(5,2) + c(1*qnorm(p0),0.5*qnorm(p0))

# beta_save1.adj=beta_save1[-which(abs(beta_save1[,1])>35),]

df0.5=data.frame(matrix(NA,nrow=4,ncol=((dim(X)[2]*3+2))),row.names = c('GAL','QR','BayesQR','ALD'))
# df0.5[1,]=c((colMeans(beta_save1.adj)-true_beta)/1, NA, (colVars(beta_save1.adj)), NA, (colMeans(beta_save1.adj)-true_beta)^2+ colVars(beta_save1.adj))
df0.5[1,]=c((colMeans(beta_save1,na.rm = T)-true_beta)/1, NA, (colVars(beta_save1,na.rm = T)), NA, (colMeans(beta_save1,na.rm = T)-true_beta)^2+ colVars(beta_save1,na.rm = T))
df0.5[2,]=c((colMeans(beta_save2,na.rm = T)-true_beta)/1, NA, (colVars(beta_save2,na.rm = T)), NA, (colMeans(beta_save2,na.rm = T)-true_beta)^2+ colVars(beta_save2,na.rm = T))
df0.5[3,]=c((colMeans(beta_save3,na.rm = T)-true_beta)/1, NA, (colVars(beta_save3,na.rm = T)), NA, (colMeans(beta_save3,na.rm = T)-true_beta)^2+ colVars(beta_save3,na.rm = T))
df0.5[4,]=c((colMeans(beta_save3.2,na.rm = T)-true_beta)/1, NA, (colVars(beta_save3.2,na.rm = T)), NA, (colMeans(beta_save3.2,na.rm = T)-true_beta)^2+ colVars(beta_save3.2,na.rm = T))

colnames(df0.5)=c(paste0('bias',1:dim(X)[2]),NA,paste0('var',1:dim(X)[2]),NA,paste0('mse',1:dim(X)[2]))
df0.5

par(mfrow=c(2,2))
hist(beta_save1[,1],nclass=30,main='GAL')
hist(beta_save2[,1],nclass=30,main='quantreg')
hist(beta_save3[,1],nclass=30,main='BayesQR')
hist(beta_save3.2[,1],nclass=30,main='ALD')
par(mfrow=c(1,1))

# save.image(file='../debugging/Debugging_w_BMQR_data2.RData')
# save.image(file='../debugging/Debugging_w_BMQR_data.RData')
load((file='../debugging/Debugging_w_BMQR_data.RData'))
