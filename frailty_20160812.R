rm(list=ls())
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

frailty_model <- function(formula, data, logp, clusterid){
  ## Renaming input arguments
  X <- model.matrix(formula, data)[, -1]
  clusterid <- data[, clusterid]
  nbeta <- length(attr(terms(formula), "term.labels"))
  npar <- 3 + length(attr(terms(formula), "term.labels"))
  if(missing(logp)) logp <- c(rep(0, npar))
  tstar <- data[, as.character(terms(formula)[[2]][2])]
  t <- data[, as.character(terms(formula)[[2]][3])]
  delta <- data[, as.character(terms(formula)[[2]][4])]
  d <- as.vector(tapply(delta,clusterid,FUN=sum))
  
  #### Likelihood #### 
  like <- function(logp){
    
    # Defining parameter names
    alpha <- exp(logp[1]) 
    eta <- exp(logp[2])
    theta <- exp(logp[3]) 
    beta <- logp[4:npar]
    
    # Constructing elements in the log-likelihood
    B <- as.vector(X%*%beta)
    h <- delta*log(eta*t^(eta-1)/alpha^eta*exp(B))
    H <- as.vector((t/alpha)^eta*exp(B))
    Hstar <- (tstar/alpha)^eta*exp(B)
    temp <- data.table(h,H,Hstar)
    temp <- temp[, j = lapply(.SD, sum), by = clusterid]
    h <- temp$h
    H <- temp$H
    Hstar <- temp$Hstar
    G <- d*log(theta)+cumsum(c(0,log(1/theta+j)))[d+1]
    
    # Log-likelihood
    ll <- -mean(G+h+1/theta*log(1+theta*Hstar)-(1/theta+d)*log(1+theta*H))
    return(ll)
  }
  
  # the function temp aggregate the data table x by cluster id
  temp <- function(x){
    temp <- data.table(x)
    temp <- as.matrix(temp[, j = lapply(.SD, sum), by = clusterid])[, -1]
  }
  
  #### Gradient ####
  gradientfunc <- function(logp, score){
    if(missing(score)) score = FALSE
    alpha <- exp(logp[1]) 
    eta <- exp(logp[2])
    theta <- exp(logp[3]) 
    beta <- logp[4:npar]
    
    # Constructing elements for gradient
    B <- as.vector(X%*%beta)
    h.eta <- delta*(1+eta*(log(t)-log(alpha)))
    H <- (t/alpha)^eta*exp(B) 
    Hstar <- (tstar/alpha)^eta*exp(B)
    H.eta <- eta*log(t/alpha)*H 
    Hstar.eta <- eta*log(tstar/alpha)*Hstar
    Hstar.eta[tstar==0] <- 0
    H.beta <- X * H
    Hstar.beta <- X* Hstar 
    h.beta <- X*delta 
    
    # Aggregate all elements that are sums over cluster id
    h.alpha <- -d*eta 
    h.eta <- temp(h.eta)
    h.beta <- temp(h.beta)
    H <- temp(H)
    Hstar <- temp(Hstar)
    H.alpha <- -eta*H
    H.eta <- temp(H.eta)
    Hstar.eta <- temp(Hstar.eta)
    H.beta <- temp(H.beta)
    Hstar.beta <- temp(Hstar.beta)
    
    Hstar.alpha <- -eta*Hstar  
    K <- H/(1+theta*H)
    Kstar <- Hstar/(1+theta*Hstar)
    G.theta <- d-cumsum(c(0,1/(1+theta*j)))[d+1]
    
    # Second derivatives of the log-likelihood
    dl.dalpha <- h.alpha+Hstar.alpha/(1+theta*Hstar)-(1+theta*d)*H.alpha/(1+theta*H)
    dl.deta <- h.eta+Hstar.eta/(1+theta*Hstar)-(1+theta*d)*H.eta/(1+theta*H)
    dl.dtheta <- G.theta+1/theta*(log(1+theta*H)-log(1+theta*Hstar))+Kstar-(1+d*theta)*K
    dl.dbeta <- h.beta+Hstar.beta/(1+theta*Hstar)-(1+theta*d)*H.beta/(1+theta*H)
    
    # Gradient for all parameters
    gradient <- -c(mean(dl.dalpha), mean(dl.deta), mean(dl.dtheta), colMeans(dl.dbeta)) 
    
    # Why not always use score?
    if(score == TRUE){
      score <- -cbind(dl.dalpha, dl.deta, dl.dtheta, dl.dbeta) 
      return(score)
    }
    else {return(gradient)}
    #names(gradient) <- c("logalpha","logeta","logtheta","beta") 
  }
  
  #### Hessian ####
  hessianfunc <- function(logp){
    
    alpha <- exp(logp[1]) 
    eta <- exp(logp[2])
    theta <- exp(logp[3]) 
    beta <- logp[4:npar]
    
    B <- as.vector(X%*%beta)
    XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
    h.eta <- delta*(1+eta*(log(t)-log(alpha)))
    H <- (t/alpha)^eta*exp(B) 
    Hstar <- (tstar/alpha)^eta*exp(B)
    H.eta <- eta*log(t/alpha)*H 
    Hstar.eta <- eta*log(tstar/alpha)*Hstar
    Hstar.eta[tstar==0] <- 0
    H.eta.eta <- H.eta+eta^2*(log(t/alpha))^2*H 
    Hstar.eta.eta <- Hstar.eta+eta^2*(log(tstar/alpha))^2*Hstar
    Hstar.eta.eta[tstar==0] <- 0
    H.eta.beta <- eta*log(t/alpha)*(H*X)
    Hstar.eta.beta <- eta*log(tstar/alpha)*(Hstar*X)
    Hstar.eta.beta[tstar==0] <- 0
    H.beta <- cbind(H[rep(1:length(H), nbeta)]*XX, clusterid)
    #dim(H.beta) <- c(nrow(X), nbeta+1, nbeta)
    Hstar.beta <- cbind(Hstar[rep(1:length(H), nbeta)]*XX, clusterid)
    H.beta.beta <- H*X^2
    Hstar.beta.beta <- Hstar*X #was "Hstar*X" but shouldnt it be "Hstar*X^2"?
    
    # Aggregate
    h.eta <- temp(h.eta)
    H <- temp(H)
    Hstar <- temp(Hstar)
    H.eta <- temp(H.eta)
    Hstar.eta <- temp(Hstar.eta)
    H.eta.eta <- temp(H.eta.eta)
    Hstar.eta.eta <- temp(Hstar.eta.eta)
    H.eta.beta <- temp(H.eta.beta)
    Hstar.eta.beta <- temp(Hstar.eta.beta)
    H.beta <- temp(H.beta)
    Hstar.beta <- temp(Hstar.beta)
    H.beta.beta <- temp(H.beta.beta)
    Hstar.beta.beta <- temp(Hstar.beta.beta)
    
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
    
    # Hessian, derivative of gradient of alpha with respect to all parameters
    dl.dalpha.dalpha <- -mean(h.alpha.alpha+Hstar.alpha.alpha/(1+theta*Hstar)-
                                theta*(Hstar.alpha/(1+theta*Hstar))^2-
                                (1+theta*d)*(H.alpha.alpha/(1+theta*H)-
                                               theta*(H.alpha/(1+theta*H))^2))  
    dl.dalpha.deta <- -mean(h.alpha.eta+Hstar.alpha.eta/(1+theta*Hstar)-
                              theta*Hstar.alpha*Hstar.eta/(1+theta*Hstar)^2-
                              (1+theta*d)*(H.alpha.eta/(1+theta*H)-
                                             theta*H.alpha*H.eta/(1+theta*H)^2))  
    dl.dalpha.dtheta <- -mean(theta*(-Hstar.alpha*Hstar/(1+theta*Hstar)^2+
                                       H.alpha*H/(1+theta*H)^2-d*(H.alpha/(1+theta*H)-theta*H.alpha*H/(1+theta*H)^2)))
    dl.dalpha.dbeta <- -t(colMeans(Hstar.alpha.beta/(1+theta*Hstar)-
                                     theta*Hstar.alpha*Hstar.beta/(1+theta*Hstar)^2-
                                     (1+theta*d)*(H.alpha.beta/(1+theta*H)-
                                                    theta*H.alpha*H.beta/(1+theta*H)^2)))
    
    dl.dalpha <- cbind(dl.dalpha.dalpha, dl.dalpha.deta, dl.dalpha.dtheta, dl.dalpha.dbeta)
    
    # Hessian, derivative of gradient of eta with respect to all parameters
    dl.deta.deta <- -mean(h.eta.eta+Hstar.eta.eta/(1+theta*Hstar)-
                            theta*(Hstar.eta/(1+theta*Hstar))^2-(1+theta*d)*
                            (H.eta.eta/(1+theta*H)-theta*(H.eta/(1+theta*H))^2))
    dl.deta.dtheta <- -mean(theta*(-Hstar.eta*Hstar/(1+theta*Hstar)^2+
                                     H.eta*H/(1+theta*H)^2-d*(H.eta/(1+theta*H)-
                                                                theta*H.eta*H/(1+theta*H)^2)))
    dl.deta.dbeta <- -t(colMeans(Hstar.eta.beta/(1+theta*Hstar)-
                                   theta*Hstar.eta*Hstar.beta/(1+theta*Hstar)^2-
                                   (1+theta*d)*(H.eta.beta/(1+theta*H)-
                                                  theta*H.eta*H.beta/(1+theta*H)^2)))
    
    dl.deta <- cbind(dl.dalpha.deta, dl.deta.deta, dl.deta.dtheta, dl.deta.dbeta)
    
    # Hessian, derivative of gradient of theta with respect to all parameters
    dl.dtheta.dtheta <- -mean(G.theta.theta+1/theta*(log(1+theta*Hstar)-log(1+theta*H))+
                                K-Kstar+theta*(K^2-Kstar^2)+d*theta*K*(theta*K-1))
    dl.dtheta.dbeta <- -t(colMeans(theta*(-Hstar.beta*Hstar/(1+theta*Hstar)^2+
                                            H.beta*H/(1+theta*H)^2-d*(H.beta/(1+theta*H)-
                                                                        theta*H.beta*H/(1+theta*H)^2))))
    
    dl.dtheta <- cbind(dl.dalpha.dtheta, dl.deta.dtheta, dl.dtheta.dtheta, dl.dtheta.dbeta)
    
    # Hessian, derivative of gradient of beta with respect to all parameters
    
    ## Wrong dimension?
    dl.dbeta.dbeta <- -t(colMeans(Hstar.beta.beta/(1+theta*Hstar)-
                                    theta*(Hstar.beta/(1+theta*Hstar))^2-(1+theta*d)*
                                    (H.beta.beta/(1+theta*H)-theta*(H.beta/(1+theta*H))^2)))
    
    dl.dbeta <- cbind(t(dl.dalpha.dbeta), t(dl.deta.dbeta), t(dl.dtheta.dbeta))
    
    hessian <- rbind(dl.dalpha, dl.deta, dl.dtheta, dl.dbeta)
    
    #colnames(hessian) <- c("logalpha","logeta","logtheta","beta")
    #rownames(hessian) <- colnames(hessian)
    
    return(hessian)
  }
  
  ########################################
  fit=optim(par=logp,fn=like,gr=gradientfunc,method="BFGS",hessian=FALSE)
  par <- fit$par
  score <- gradientfunc(par, score=TRUE)
  
  #### Output ####
  out <- c(list(par = par, score = score))
  class(out) <- "frailty_model"
  return(out)
  
}


print("true")
print(logp)
print("with grad")
estimates<-frailty_model(formula, data, logp, clusterid="id")
estimates







