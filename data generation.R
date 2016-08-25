library(data.table)
#library(tictoc)
library(numDeriv)
library(survival)

#---DATA GENERATION---

n=1000
m=2
eta.true=1
alpha0.true=1.8
theta.true=0.4
beta1.true=1.3
beta2.true=2
beta3.true=0.5
par.true=c(alpha0.true,eta.true,theta.true,beta1.true, beta2.true, beta3.true)
id=rep(1:n,each=m)
u=rep(rgamma(n,shape=1/theta.true,scale=theta.true),each=m)
x=rnorm(n*m)
z=rnorm(n*m)
alpha=alpha0.true/(u*exp(beta1.true*x + beta2.true*z + beta3.true*z*x))^(1/eta.true)
t=rweibull(n*m,shape=eta.true,scale=alpha)

#right-censoring
#c <- rep(Inf,n*m)
c <- runif(n*m,0,10)
delta <- as.numeric(t<c)
t <- pmin(t,c)

#strong left-truncation
#tstar <- rep(0,n*m)
tstar <- runif(n*m,0,2)
incl <- t>tstar
incl <- ave(x=incl,id,FUN=sum)==m
id <- id[incl]
x <- x[incl]
z <- z[incl]
t <- t[incl]
delta <- delta[incl]
tstar <- tstar[incl]
data <- data.frame(t, tstar, delta, x, z, id)

d <- as.vector(tapply(delta,id,FUN=sum))
D <- max(d)
j <- 0:(D-1)

#---LIKELIHOOD FUNCTION---
logp <- c(log(alpha0.true),log(eta.true),log(theta.true),beta1.true, beta2.true, beta3.true)
formula=Surv(tstar, t, delta) ~ x + z + z*x

clusterid<-"id"


#---DATA GENERATION 2---

n=1000
m=3
N=100
eta.true=1
alpha0.true=1.8
theta.true=0.4
beta1.true=1.3
beta2.true=2
beta3.true=0.5


logp <- c(log(alpha0.true),log(eta.true),log(theta.true),beta1.true, beta2.true, beta3.true)
formula=Surv(tstar, t, delta) ~ x + z + z*x

est.all <- matrix(nrow=N,ncol=length(logp))
se.all <-  matrix(nrow=N,ncol=length(logp))

for(i in 1:N){
  print(i)
  id=rep(1:n,each=m)
  u=rep(rgamma(n,shape=1/theta.true,scale=theta.true),each=m)
  x=rnorm(n*m)
  z=rnorm(n*m)
  alpha=alpha0.true/(u*exp(beta1.true*x + beta2.true*z + beta3.true*z*x))^(1/eta.true)
  t=rweibull(n*m,shape=eta.true,scale=alpha)
  
  #right-censoring
  #c <- rep(Inf,n*m)
  c <- runif(n*m,0,10)
  delta <- as.numeric(t<c)
  t <- pmin(t,c)
  
  #strong left-truncation
  #tstar <- rep(0,n*m)
  tstar <- runif(n*m,0,2)
  incl <- t>tstar
  incl <- ave(x=incl,id,FUN=sum)==m
  
  id <- id[incl]
  x <- x[incl]
  z <- z[incl]
  t <- t[incl]
  delta <- delta[incl]
  tstar <- tstar[incl]
  data <- data.frame(t, tstar, delta, x, z, id)
  
  fit <-frailty(formula, data, logp, clusterid="id") 
  est.all[i, ] <- fit$par
  nclust <- length(unique(id))
  se.all[i, ] <- sqrt(-diag(solve(fit$hessian)))
  
}

print("logp")
print(logp)
print("colMeans(est.all)")
print(colMeans(est.all))

print("apply(X=est.all,MARGIN=2,FUN=sd)")
print(apply(X=est.all,MARGIN=2,FUN=sd))
print("colMeans(se.all)")
print(colMeans(se.all))
