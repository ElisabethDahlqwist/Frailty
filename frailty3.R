rm(list=ls())
library(data.table)
library(tictoc)
library(numDeriv)

#---DATA GENERATION---

n=20
m=2
eta.true=1
alpha0.true=1.8
theta.true=0.4
beta.true=1.3
par.true=c(alpha0.true,eta.true,theta.true,beta.true)
id=rep(1:n,each=m)
u=rep(rgamma(n,shape=1/theta.true,scale=theta.true),each=m)
x=rnorm(n*m)
alpha=alpha0.true/(u*exp(beta.true*x))^(1/eta.true)
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
t <- t[incl]
delta <- delta[incl]
tstar <- tstar[incl]

d <- tapply(delta,id,FUN=sum)
D <- max(d)
j <- 0:(D-1)

#---LIKELIHOOD FUNCTION---

like <- function(logp){
 
  alpha <- exp(logp[1]) 
  eta <- exp(logp[2])
  theta <- exp(logp[3]) 
  beta <- logp[4]
  
  h <- delta*log(eta*t^(eta-1)/alpha^eta*exp(beta*x))  
  H <- (t/alpha)^eta*exp(beta*x) 
  Hstar <- (tstar/alpha)^eta*exp(beta*x)
  temp <- data.table(h,H,Hstar)
  temp <- temp[, j = lapply(.SD, sum), by = id]
  h <- temp$h
  H <- temp$H
  Hstar <- temp$Hstar
  G <- d*log(theta)+cumsum(c(0,log(1/theta+j)))[d+1]
  ll <- -mean(G+h+1/theta*log(1+theta*Hstar)-(1/theta+d)*log(1+theta*H))
  
  return(ll)
  
}

gradientfunc <- function(logp){
  
  alpha <- exp(logp[1]) 
  eta <- exp(logp[2])
  theta <- exp(logp[3]) 
  beta <- logp[4]
  
  h.eta <- delta*(1+eta*(log(t)-log(alpha)))
  h.beta <- delta*x
  H <- (t/alpha)^eta*exp(beta*x) 
  Hstar <- (tstar/alpha)^eta*exp(beta*x)
  H.eta <- eta*log(t/alpha)*H 
  Hstar.eta <- eta*log(tstar/alpha)*Hstar
  Hstar.eta[tstar==0] <- 0
  H.beta <- x*H
  Hstar.beta <- x*Hstar
  
  temp <- data.table(h.eta,h.beta,
    H,Hstar,
    H.eta,Hstar.eta,
    H.beta,Hstar.beta)
  temp <- temp[, j = lapply(.SD, sum), by = id]
 
  h.alpha <- -d*eta 
  h.eta <- temp$h.eta
  h.beta <- temp$h.beta
  H <- temp$H
  Hstar <- temp$Hstar
  H.alpha <- -eta*temp$H
  H.eta <- temp$H.eta
  Hstar.eta <- temp$Hstar.eta
  H.beta <- temp$H.beta
  Hstar.beta <- temp$Hstar.beta
  
  Hstar.alpha <- -eta*Hstar  
  K <- H/(1+theta*H)
  Kstar <- Hstar/(1+theta*Hstar)
  G.theta <- d-cumsum(c(0,1/(1+theta*j)))[d+1]

  dl.dalpha <- h.alpha+Hstar.alpha/(1+theta*Hstar)-(1+theta*d)*H.alpha/(1+theta*H)
  dl.deta <- h.eta+Hstar.eta/(1+theta*Hstar)-(1+theta*d)*H.eta/(1+theta*H)
  dl.dtheta <- G.theta+1/theta*(log(1+theta*H)-log(1+theta*Hstar))+Kstar-
    (1+d*theta)*K
  dl.dbeta <- h.beta+Hstar.beta/(1+theta*Hstar)-(1+theta*d)*H.beta/(1+theta*H)
  
  gradient <- -c(mean(dl.dalpha),mean(dl.deta),mean(dl.dtheta),mean(dl.dbeta)) 
  names(gradient) <- c("logalpha","logeta","logtheta","beta") 
  
  return(gradient)

}

hessianfunc <- function(logp){

  alpha <- exp(logp[1]) 
  eta <- exp(logp[2])
  theta <- exp(logp[3]) 
  beta <- logp[4]
  
  h.eta <- delta*(1+eta*(log(t)-log(alpha)))
  H <- (t/alpha)^eta*exp(beta*x) 
  Hstar <- (tstar/alpha)^eta*exp(beta*x)
  H.eta <- eta*log(t/alpha)*H 
  Hstar.eta <- eta*log(tstar/alpha)*Hstar
  Hstar.eta[tstar==0] <- 0
  H.eta.eta <- H.eta+eta^2*(log(t/alpha))^2*H 
  Hstar.eta.eta <- Hstar.eta+eta^2*(log(tstar/alpha))^2*Hstar
  Hstar.eta.eta[tstar==0] <- 0
  H.eta.beta <- eta*log(t/alpha)*x*H
  Hstar.eta.beta <- eta*log(tstar/alpha)*x*Hstar
  Hstar.eta.beta[tstar==0] <- 0
  H.beta <- x*H
  Hstar.beta <- x*Hstar
  H.beta.beta <- x^2*H
  Hstar.beta.beta <- x^2*Hstar 
 
  temp <- data.table(h.eta,
    H,Hstar,
    H.eta,Hstar.eta,H.eta.eta,Hstar.eta.eta,H.eta.beta,Hstar.eta.beta,
    H.beta,Hstar.beta,H.beta.beta,Hstar.beta.beta)
  temp <- temp[, j = lapply(.SD, sum), by = id]
  
  h.eta <- temp$h.eta
  H <- temp$H
  Hstar <- temp$Hstar
  H.eta <- temp$H.eta
  Hstar.eta <- temp$Hstar.eta
  H.eta.eta <- temp$H.eta.eta
  Hstar.eta.eta <- temp$Hstar.eta.eta
  H.eta.beta <- temp$H.eta.beta
  Hstar.eta.beta <- temp$Hstar.eta.beta
  H.beta <- temp$H.beta
  Hstar.beta <- temp$Hstar.beta
  H.beta.beta <- temp$H.beta.beta
  Hstar.beta.beta <- temp$Hstar.beta.beta
  
  h.alpha.alpha <- 0
  h.alpha.eta <- -d*eta
  h.eta.eta <- h.eta-d
  H.alpha <- -eta*H
  Hstar.alpha <- -eta*Hstar
  H.alpha.alpha <- eta^2*H
  Hstar.alpha.alpha <- eta^2*Hstar  
  H.alpha.eta <- -eta*(H+H.eta)
  Hstar.alpha.eta <- -eta*(Hstar+Hstar.eta)
  H.alpha.beta <- -eta*H.beta
  Hstar.alpha.beta <- -eta*Hstar.beta
  
  K <- H/(1+theta*H)
  Kstar <- Hstar/(1+theta*Hstar)
  G.theta.theta <- cumsum(c(0,theta*j/(1+theta*j)^2))[d+1]
  
  hessian <- matrix(nrow=4,ncol=4)
  dl.dalpha.dalpha <- mean(h.alpha.alpha+Hstar.alpha.alpha/(1+theta*Hstar)-
    theta*(Hstar.alpha/(1+theta*Hstar))^2-
    (1+theta*d)*(H.alpha.alpha/(1+theta*H)-
    theta*(H.alpha/(1+theta*H))^2))  
  dl.dalpha.deta <- mean(h.alpha.eta+Hstar.alpha.eta/(1+theta*Hstar)-
    theta*Hstar.alpha*Hstar.eta/(1+theta*Hstar)^2-
    (1+theta*d)*(H.alpha.eta/(1+theta*H)-
    theta*H.alpha*H.eta/(1+theta*H)^2))  
  dl.dalpha.dtheta <- mean(theta*(-Hstar.alpha*Hstar/(1+theta*Hstar)^2+
    H.alpha*H/(1+theta*H)^2-d*(H.alpha/(1+theta*H)-theta*H.alpha*H/(1+theta*H)^2)))
  dl.dalpha.dbeta <- mean(Hstar.alpha.beta/(1+theta*Hstar)-
    theta*Hstar.alpha*Hstar.beta/(1+theta*Hstar)^2-
    (1+theta*d)*(H.alpha.beta/(1+theta*H)-
    theta*H.alpha*H.beta/(1+theta*H)^2))
  dl.deta.deta <- mean(h.eta.eta+Hstar.eta.eta/(1+theta*Hstar)-
    theta*(Hstar.eta/(1+theta*Hstar))^2-(1+theta*d)*
    (H.eta.eta/(1+theta*H)-theta*(H.eta/(1+theta*H))^2))
  dl.deta.dtheta <- mean(theta*(-Hstar.eta*Hstar/(1+theta*Hstar)^2+
    H.eta*H/(1+theta*H)^2-d*(H.eta/(1+theta*H)-
    theta*H.eta*H/(1+theta*H)^2)))
  dl.deta.dbeta <- mean(Hstar.eta.beta/(1+theta*Hstar)-
    theta*Hstar.eta*Hstar.beta/(1+theta*Hstar)^2-
    (1+theta*d)*(H.eta.beta/(1+theta*H)-
    theta*H.eta*H.beta/(1+theta*H)^2))
  dl.dtheta.dtheta <- mean(G.theta.theta+1/theta*(log(1+theta*Hstar)-log(1+theta*H))+
    K-Kstar+theta*(K^2-Kstar^2)+d*theta*K*(theta*K-1))
  dl.dtheta.dbeta <- mean(theta*(-Hstar.beta*Hstar/(1+theta*Hstar)^2+
    H.beta*H/(1+theta*H)^2-d*(H.beta/(1+theta*H)-
    theta*H.beta*H/(1+theta*H)^2)))
  dl.dbeta.dbeta <- mean(Hstar.beta.beta/(1+theta*Hstar)-
    theta*(Hstar.beta/(1+theta*Hstar))^2-(1+theta*d)*
    (H.beta.beta/(1+theta*H)-theta*(H.beta/(1+theta*H))^2))
  
  hessian[1,1] <- -dl.dalpha.dalpha  
  hessian[1,2] <- -dl.dalpha.deta 
  hessian[1,3] <- -dl.dalpha.dtheta
  hessian[1,4] <- -dl.dalpha.dbeta 
  hessian[2,2] <- -dl.deta.deta 
  hessian[2,3] <- -dl.deta.dtheta 
  hessian[2,4] <- -dl.deta.dbeta
  hessian[3,3] <- -dl.dtheta.dtheta 
  hessian[3,4] <- -dl.dtheta.dbeta
  hessian[4,4] <- -dl.dbeta.dbeta
  
  hessian[lower.tri(hessian)] <- t(hessian)[lower.tri(hessian)]
  
  colnames(hessian) <- c("logalpha","logeta","logtheta","beta")
  rownames(hessian) <- colnames(hessian)
  
  return(hessian)
}

#---Old code

logp <- c(log(alpha0.true),log(eta.true),log(theta.true),beta.true)

print("true parameter values")
print(logp)

fit=optim(par=logp,fn=like,gr=gradientfunc,method="BFGS",hessian=FALSE)  
print("optim with gradient")
print(fit$par)

gradient <- gradientfunc(logp) 
print("gradient")       
print(gradient)

numDeriv.hess <- hessian(func=like,x=logp)
print("numerical hessian")
print(numDeriv.hess)


hessian <- hessianfunc(logp) 
print("analytic hessian")       
print(hessian)

#-- New code
logp <- c(log(alpha0.true),log(eta.true),log(theta.true),beta.true)

formula <- Surv(tstar, t, delta) ~ x
clusterid<-"id"

print("true")
print(logp)
print("with grad")
estimates<-frailty_model(formula, data, logp, clusterid="id")

estimates





