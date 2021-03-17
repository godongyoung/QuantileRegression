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

# Check with BQR --------------------------------------------------------------------------------
plot(X[,2],y, main="BQR fitted lines")
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

# Check with NQR --------------------------------------------------------------------------------
plot(X[,2],y, main="NQR fitted lines")
for(p0 in c(.05,.25,.5,.75,.95)){
  NQR_res=NQR(y,X,p0)
  g.est=colMeans(NQR_res$g_trace)
  lambda.est=mean(NQR_res$lambda_trace)
  mspline=spline(x = NQR_res$Knots,y = g.est,xout = X_range)
  
  points(X_range,mspline$y,type = 'l')    
}

# Check with empirical QR --------------------------------------------------------------------------------
N=10
plot(X[,2],y, main="Empirical QR fitted lines")
for(p0 in c(.05,.25,.5,.75,.95)){
  empirical_y.est=rep(NA,(N-1))
  for(idx in 1:(N-1)){
    Knots_idx=30/N*idx
    target=(NQR_res$Knots[Knots_idx]<X[,2])&(X[,2]<NQR_res$Knots[(Knots_idx+30/N)])
    empirical_y.est[idx]=quantile(y[target],probs = p0)
  }
  points(NQR_res$Knots,c(rep(empirical_y.est[1],30/N),rep(empirical_y.est,each=30/N)),type = 'l')    
}
