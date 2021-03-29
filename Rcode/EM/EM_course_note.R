# Code for showing what happens in measurement error model.----------------------------------------------------------------------------------------

### Set up
B = 1000
n = 200
rst1 = matrix(NA,B,4)
rst2 = matrix(NA,B,4)

### MC Simulation
set.seed(1105)
for(i in 1:B){
  e = rnorm(n)
  u = rnorm(n)
  x = rgamma(n,1,1)
  z = x + u
  y = 1 + 0.5*x + e
  dat = data.frame(y,x,z)
  
  fit1 = lm(y ~ x, data=dat)
  fit2 = lm(y ~ z, data=dat)
  
  rst1[i,] = c(fit1$coeff, summary(fit1)$sigma^2 , sqrt(mean((y - fit1$fitted)^2)))
  rst2[i,] = c(fit2$coeff, summary(fit2)$sigma^2 , sqrt(mean((y - fit2$fitted)^2)))
  
  if(i%%100==0){
    cat("_____Iteratin___",'\n')
    print(i)
    cat("_____True Model___",'\n')
    print(round(apply(rst1[1:i,],2,mean),3))
    print(round(apply(rst1[1:i,],2,sd),3))
    cat("_____Measurement Error Model___",'\n')
    print(round(apply(rst2[1:i,],2,mean),3))
    print(round(apply(rst2[1:i,],2,sd),3))    
  } ## print
} ## MC simulation




# Gauss Newton ----------------------------------------------------------------------------------------
rm(list=ls()) 

n=400 
B=1000 
rst=matrix(NA,B,9) 
colnames(rst)=c('beta0','beta1','alpha0','alpha1','Mux','sigma2_xx','sigma2','sigma2_11','sigma2_22')

set.seed(1105)
for(i in 1:B){ 
  ### Make data ------------------------
  x=rgamma(n,1,1) 
  e=rnorm(n) 
  u1=rnorm(n) 
  u2=rnorm(n) 
  y=0.3+1.2*x+e 
  W1=0.5+0.9*x+u2
  W2=x+u1 
  
  
  dat=data.frame(y,W2,W1) 
  vdat=var(dat) 
  eta.hat=c(colMeans(dat),diag(vdat),vdat[1,2:3],vdat[2,3])

  theta=rep(1,9) 
  eta.theta=vector(length=9) 
  A=matrix(0,9,9)
  
  iter=0 
  repeat{
    iter=iter+1 
    theta0=theta
    
    eta.theta[1]=theta[1]+theta[2]*theta[5]
    eta.theta[2]=theta[5]
    eta.theta[3]=theta[3]+theta[4]*theta[5]
    eta.theta[4]=theta[2]^2*theta[6]+theta[7]
    eta.theta[5]=theta[6]+theta[8]
    eta.theta[6]=theta[4]^2*theta[6]+theta[9]
    eta.theta[7]=theta[2]*theta[6]
    eta.theta[8]=theta[2]*theta[4]*theta[6]
    eta.theta[9]=theta[4]*theta[6]
    
    Z=eta.hat-eta.theta
    A[1,c(1,2,5)]=c(1,theta[5],theta[2])
    A[2,5]=1
    A[3,3:5]=c(1,theta[5],theta[4])
    A[4,c(2,6,7)]=c(2*theta[2]*theta[6],theta[2]^2,1)
    A[5,c(6,8)]=c(1,1)
    A[6,c(4,6,9)]=c(2*theta[4]*theta[6],theta[4]^2,1)
    A[7,c(2,6)]=c(theta[6],theta[2])
    A[8,c(2,4,6)]=c(theta[4]*theta[6],theta[2]*theta[6],theta[2]*theta[4])
    A[9,c(4,6)]=c(theta[6],theta[4])
    
    theta=theta0+solve(t(A)%*%A)%*%t(A)%*%Z
    dif=sum((theta-theta0)^2)
    
    if(dif<1e-6){break}
    }
  rst[i,]=theta
  if(i%%500==0){
    cat("__________________MonteCarlo Mean and SE__________________","\n")
    # cat('Mean:',round(apply(rst[1:i,],2,mean),2),'\n')
    # cat('Sd:',round(apply(rst[1:i,],2,sd),2),'\n')
    print(round(apply(rst[1:i,],2,mean),2))
    print(round(apply(rst[1:i,],2,sd),2))
  }
}


