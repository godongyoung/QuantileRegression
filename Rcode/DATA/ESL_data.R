setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')

# ESL------------------------------------------------------------------------------------------
library(ElemStatLearn)

data("bone") #bone mineral
head(bone);dim(bone)
plot(bone$age,bone$spnbmd)

data("ozone")
dim(ozone)


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
data("mcycle")
y=mcycle$accel
x1=mcycle$times
X = cbind(1,x1)
N = 30
Knots = seq(min(X[,2]),max(X[,2]),length.out = N)
plot(X[,2],y)
y.est = smooth.y(X[,2],y,Knots,version = 3)
points(Knots,y.est,type='l')


g.tau=as.numeric(y)
fit.ns <- lm(g.tau~ ns(x = knots, knots = knots[-c(1,N)]) )
y.est=predict(fit.ns, data.frame(knots=xout))