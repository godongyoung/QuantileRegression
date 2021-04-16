setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')

# ESL------------------------------------------------------------------------------------------
library(ElemStatLearn)

data("bone") #bone mineral
head(bone);dim(bone)
plot(bone$age,bone$spnbmd)

# write.table(bone,'./DATA/bone.txt')
n=1000
X = runif(n,0,1)
y = sin(12*(X+0.2))/(X+0.2) + rnorm(n,0,1)
plot(X,y)

X = runif(n,0,1)
y = sin(4*X) + rnorm(n,0,1)
plot(X,y)

# Design Adaptive Nonparametric Regression------------------------------------------------------------------------------------------

# X = runif(n,0,1)
X = runif(n,-3,3)
X = cbind(rnorm(n/2,-1,1), rnorm(n/2,1.75,0.25));hist(X,nclass=100)
y = sin(2.5*X) + 0.4*rnorm(n,0,1)
plot(X,y)

# An Introduction to Kernel and Nearest-Neighbor Nonparametric Regression ------------------------------------------------------------------------------------------
X = runif(n,0,1)
y = X*sin(2.5*pi*X) + rnorm(n,0,0.1)
plot(X,y)

chicago = read.table('./DATA/T67.1.txt')[,-c(1,2,3)]
colnames(chicago) = c('Zip Code', 'Racial comp', 'Fire', 'Theft' ,'Age', 'Voluntary market activity', 'Involuntary market activity', 'Income')
head(chicago)
X = chicago$Theft
y = chicago$'Voluntary market activity'

adj.const = 10
X = rep(X,adj.const)
y = rep(y,adj.const)
n = length(y)
X = X + rnorm(n,0,4)
y = y + rnorm(n)
X = cbind(1,X)
plot(X[,2],y)
# p0 = 0.1
# NQR_res = NQR(y,X,p0,inp.min = 0,inp.max = 150,inp.version = 1,multiply_c = 3)
# g.est = colMeans(NQR_res$g_trace)
# Knots = NQR_res$Knots
# points(Knots,g.est,type='l')

# Bayesian nonparametric quantile regression ------------------------------------------------------------------------------------------
# I will follow the procedure they described. 
# They firstly fit the ncs on mcycle data (didn't specified the knots)
# And then based on the fitted value of ncs at 30 equally spaced on X, 
# they generated 1000 replicated with gaussian error on each spaces
data("mcycle")
y=mcycle$accel
x1=mcycle$times
n=length(y)

# x1=x1+rnorm(n,0,0.1)
# order_idx=order(x1)
# x1=x1[order_idx]
# y=y[order_idx]

N=10
Basis.Knots = seq(min(x1),max(x1),length.out = N)
plot.Knots = seq(min(x1),max(x1),length.out = 1000)
plot(x1,y)

data = data.frame(y,x1)
colnames(data) = c('y','x1')
model = lm(y~
             ns(x1,knots=Basis.Knots[-c(1,length(Basis.Knots))],
                # Boundary.knots = Basis.Knots[c(1,length(Basis.Knots))]
                ),
           data = data )
y.est = predict(model,newdata = list(x1=plot.Knots))
points(plot.Knots,y.est,type='l')

# start generating data
N=30
Knots = seq(min(x1),max(x1),length.out = N)
y.est.smooth = predict(model,newdata = list(x1=Knots))
y.generate = sapply(X = y.est.smooth,FUN = function(x) x+rnorm(1000,0,sd = 20))
y.generate = as.numeric(y.generate)
x.generate = rep(Knots,each = 1000)
plot(x.generate,y.generate)

y.est.smooth.true = predict(model,newdata = list(x1=plot.Knots)) + qnorm(0.95,0,20)
points(plot.Knots,y.est.smooth.true,type='l',col=4,lwd=3)

X.generate = cbind(1,x.generate)
# NQR_res = NQR(y.generate,X.generate,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 1,multiply_c = 6)
# save(NQR_res, file=sprintf('./NQR_v1_fitted_%s.RData',0.95))

# NQR_res = NQR(y.generate,X.generate,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 3,multiply_c = 6)
# save(NQR_res, file=sprintf('./NQR_v3_fitted_%s.RData',0.95))

load(file=sprintf('./NQR_v1_fitted_%s.RData',0.95))
# ts.plot(NQR_res$g_trace[,10])
g.est = colmeans(NQR_res$g_trace[3000:5000,])
points(NQR_res$Knots,g.est,type='l',col=2,lwd=3)

n=length(y.generate)
sigma2_11 = 1
sigma2_22 = 1
delta1=rnorm(n,0,sd=sqrt(sigma2_11))
delta2=rnorm(n,0,sd=sqrt(sigma2_22))

alpha=c(4,3)
W1=X.generate%*%alpha+delta1
W2=X.generate[,1]+delta2

# NQR_w_MME_res = NQR_w_MME(y.generate,W1,W2,0.95,inp.min = min(x1),inp.max = max(x1),inp.version = 1,multiply_c = 6)
# save(NQR_w_MME_res, file=sprintf('./NQRwMME_v1_fitted_%s.RData',0.95))
load(file=sprintf('./NQRwMME_v1_fitted_%s.RData',0.95))

# ts.plot(NQR_w_MME_res$g_trace[,10])
g.est = colmeans(NQR_w_MME_res$g_trace[3000:5000,])
points(NQR_w_MME_res$Knots,g.est,type='o',col=3,lwd=3)
plot(x.generate,colMeans(NQR_w_MME_res$X_trace));abline(0,1)

# DPpackage: Bayesian Semi- and Nonparametric Modeling in R ------------------------------------------------------------------------------------------
set.seed(0)
n <- 500
x1i <- runif(n)
X = cbind(1,x1i)
y1 <- x1i + rnorm(n, 0, sqrt(0.01))
y2 <- x1i^4 + rnorm(n, 0, sqrt(0.04))
u <- runif(n)
prob <- exp(-2 * x1i)
y <- ifelse(u < prob, y1, y2)
plot(X[,2],y)