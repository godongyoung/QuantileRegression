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
library(bayesQR)
library(quantreg)




# Data --------------------------------------------------------------------------------
data("ImmunogG")
head(ImmunogG)

y=ImmunogG$IgG
X=cbind(1,ImmunogG$Age,ImmunogG$Age^2)

X_range=seq(from = min(X[,2]),to = max(X[,2]),length.out = 1000)
plot(X[,2],y)
for (p0 in c(.05,.25,.5,.75,.95)){
  # # QR
  # mQR_res = mQR(y,X,p0)
  # beta.est=mQR_res[['beta_est']]
  # y_est=beta.est[1]+beta.est[2]*X_range+beta.est[3]*X_range^2
  # points(X_range,y_est,type = 'l')
  
  # BQR
  BQR_res = mBayesQR(y,X,p0)
  beta.est=colMeans(BQR_res$beta_trace)
  y_est=beta.est[1]+beta.est[2]*X_range+beta.est[3]*X_range^2
  points(X_range,y_est,type = 'l')  
}

plot(X[,2],y)
for(p0 in c(.05,.25,.5,.75,.95)){
  NQR_res=NQR(y,X,p0)
  g.est=colMeans(NQR_res$g_trace)
  lambda.est=mean(NQR_res$lambda_trace)
  mspline=spline(x = NQR_res$Knots,y = g.est,xout = X_range)
  
  points(X_range,mspline$y,type = 'l')    
}
