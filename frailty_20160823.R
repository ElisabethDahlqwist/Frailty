rm(list=ls())
library(data.table)
#library(tictoc)
library(numDeriv)
library(survival)

#---DATA GENERATION---

n=20
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

frailty <- function(formula, data, logp, clusterid){
  call <- match.call()
  ## Defining input arguments
  X <- as.matrix(model.matrix(formula, data)[, -1])
  clusterid <- data[, clusterid]
  n <- nrow(X)
  ncluster <- length(unique(clusterid))
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
    beta <- as.matrix(logp[4:npar])
    
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
  
  # The function "temp" aggregate the data table x by cluster id
  temp <- function(x){
    temp <- data.table(x)
    temp <- as.matrix(temp[, j = lapply(.SD, sum), by = clusterid])[, -1]
  }
  
  #### Gradient ####
  gradientfunc <- function(logp, score=F){
    if(missing(score)) score = FALSE
    alpha <- exp(logp[1]) 
    eta <- exp(logp[2])
    theta <- exp(logp[3]) 
    beta <- as.matrix(logp[4:npar])
    
    # Constructing elements for gradient
    B <- as.vector(X%*%beta)
    h.eta <- delta*(1+eta*(log(t)-log(alpha)))
    h.beta <- X*delta 
    H <- (t/alpha)^eta*exp(B) 
    Hstar <- (tstar/alpha)^eta*exp(B)
    H.eta <- eta*log(t/alpha)*H 
    Hstar.eta <- eta*log(tstar/alpha)*Hstar
    Hstar.eta[tstar==0] <- 0
    H.beta <- X * H
    Hstar.beta <- X* Hstar 
    
    
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
    dl.dbeta <- as.matrix(h.beta+Hstar.beta/(1+theta*Hstar)-(1+theta*d)*H.beta/(1+theta*H))
    
    # Gradient for all parameters
    gradient <- -c(mean(dl.dalpha), mean(dl.deta), mean(dl.dtheta), colMeans(dl.dbeta)) 
    
    # Score function for all parameters and individuals
    if(score == TRUE){
      score <- -cbind(dl.dalpha, dl.deta, dl.dtheta, dl.dbeta) 
      return(score=score)
    }
    else {return(gradient)}
    #names(gradient) <- c("logalpha","logeta","logtheta","beta") 
  }
  
  #### Hessian ####
  hessianfunc <- function(logp){
    
    alpha <- exp(logp[1]) 
    eta <- exp(logp[2])
    theta <- exp(logp[3]) 
    beta <- as.matrix(logp[4:npar])
    
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
    Hstar.beta <- cbind(Hstar[rep(1:length(H), nbeta)]*XX, clusterid)
    
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
    
    h.alpha.alpha <- 0
    h.alpha.eta <- -d*eta
    h.eta.eta <- h.eta-d
    H.alpha <- -eta*H
    Hstar.alpha <- -eta*Hstar
    H.alpha.alpha <- eta^2*H
    Hstar.alpha.alpha <- eta^2*Hstar  
    H.alpha.eta <- -eta*(H+H.eta)
    Hstar.alpha.eta <- -eta*(Hstar+Hstar.eta)

    K <- H/(1+theta*H)
    Kstar <- Hstar/(1+theta*Hstar)
    G.theta.theta <- cumsum(c(0,theta*j/(1+theta*j)^2))[d+1]
    
    # Hessian, derivative of gradient of alpha with respect to all parameters except beta
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
    
    dl.dalpha <- cbind(dl.dalpha.dalpha, dl.dalpha.deta, dl.dalpha.dtheta)
    
    # Hessian, derivative of gradient of eta with respect to all parameters except beta
    dl.deta.deta <- -mean(h.eta.eta+Hstar.eta.eta/(1+theta*Hstar)-
                            theta*(Hstar.eta/(1+theta*Hstar))^2-(1+theta*d)*
                            (H.eta.eta/(1+theta*H)-theta*(H.eta/(1+theta*H))^2))
    dl.deta.dtheta <- -mean(theta*(-Hstar.eta*Hstar/(1+theta*Hstar)^2+
                                     H.eta*H/(1+theta*H)^2-d*(H.eta/(1+theta*H)-
                                                  theta*H.eta*H/(1+theta*H)^2)))
    
    dl.deta <- cbind(dl.dalpha.deta, dl.deta.deta, dl.deta.dtheta)
    
    # Hessian, derivative of gradient of theta with respect to all parameters except beta
    dl.dtheta.dtheta <- -mean(G.theta.theta+1/theta*(log(1+theta*Hstar)-log(1+theta*H))+
                                K-Kstar+theta*(K^2-Kstar^2)+d*theta*K*(theta*K-1))
    dl.dtheta <- cbind(dl.dalpha.dtheta, dl.deta.dtheta, dl.dtheta.dtheta)
    
    # Hessian, derivative of gradient of beta with respect to all parameters
    H <- (t/alpha)^eta*exp(B) 
    Hstar <- (tstar/alpha)^eta*exp(B)
    XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
    nbeta_rep <- rep(1:nbeta, each = nrow(X))
    
    ### Creating squared cluster sums of H.beta and Hstar.beta
    H.beta <- as.matrix(temp(H * X))
    H.beta2 <- H.beta[rep(1:nrow(H.beta), nbeta), ] * c(H.beta)
    Hstar.beta <- as.matrix(temp(Hstar * X))
    Hstar.beta2 <- Hstar.beta[rep(1:nrow(Hstar.beta), nbeta), ] * c(Hstar.beta)
    
    ### Creating Cross products of covariates multiplied with H and Hstar
    Hstar.beta.beta <- data.table(nbeta_rep, clusterid, Hstar * XX)
    H.beta.beta <- data.table(nbeta_rep, clusterid, H * XX)
    
    ### Aggregate H and Hstar over clusters
    H <- temp(H)
    Hstar <- temp(Hstar)
    
    ### Calculating Hstar2 <- theta*(sum(H*X))^2/(1+theta*sum(H))^2 and H2 <- (1+d*theta)*(sum(H*X))^2/(1+theta*H)^2
    Hstar2 <- theta * Hstar.beta2 / (1 + theta * Hstar)^2
    H2 <- theta * (1+d*theta) * H.beta2 / (1 + theta * H)^2
    
    ### Aggregate Hstar.beta.beta and H.beta.beta over cluster
    Hstar.beta.beta <- data.table(clusterid, nbeta_rep, Hstar.beta.beta)
    Hstar.beta.beta <- as.matrix(Hstar.beta.beta[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])[, -1:-2, drop=FALSE] # because columns are droped it is no longer a matrix
    
    H.beta.beta <- data.table(clusterid, nbeta_rep, H.beta.beta)
    H.beta.beta <- as.matrix(H.beta.beta[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])[, -1:-2, drop=FALSE] 
    
    ### Calculate Hstar1 <- Hstar.beta.beta/(1+theta*Hstar) and H1 <- theta * (1 + d * theta)*H.beta.beta/(1+theta*H)
    Hstar1 <- Hstar.beta.beta / (1 + theta * Hstar)
    H1 <- (1 + d * theta) * H.beta.beta / (1 + theta * H)
    
    dl.dbeta.dbeta <- -(Hstar.beta.beta/(1+theta*Hstar)-
                        theta*Hstar.beta2/(1+theta*Hstar)^2-(1+theta*d)*
                       (H.beta.beta/(1+theta*H)-theta*H.beta2/(1+theta*H)^2))
    
    ## Aggregate over clusters
    nbeta_rep2 <- rep(1:nbeta, each = length(H))
    dl.dbeta.dbeta <- data.table(nbeta_rep2, dl.dbeta.dbeta)
    dl.dbeta.dbeta <- as.matrix(dl.dbeta.dbeta[, j = lapply(.SD, mean), by = .(nbeta_rep2)])[, -1]
   
    # Derivative of gradient of alpha with respect to beta
    H.alpha.beta <- -eta*H.beta
    Hstar.alpha.beta <- -eta*Hstar.beta
    
    dl.dalpha.dbeta <- -(colMeans(as.matrix(Hstar.alpha.beta/(1+theta*Hstar)-
                                               theta*Hstar.alpha*Hstar.beta/(1+theta*Hstar)^2-
                                               (1+theta*d)*(H.alpha.beta/(1+theta*H)-
                                                theta*H.alpha*H.beta/(1+theta*H)^2))))
    dl.dalpha <- cbind(dl.dalpha, t(dl.dalpha.dbeta))
    
    # Derivative of gradient of eta with respect to beta
    dl.deta.dbeta <- -t(colMeans(as.matrix(Hstar.eta.beta/(1+theta*Hstar)-
                                             theta*Hstar.eta*Hstar.beta/(1+theta*Hstar)^2-
                                             (1+theta*d)*(H.eta.beta/(1+theta*H)-
                                                            theta*H.eta*H.beta/(1+theta*H)^2))))
    dl.deta <- cbind(dl.deta, dl.deta.dbeta)
    
    # Derivative of gradient of theta with respect to beta
    dl.dtheta.dbeta <- -t(colMeans(as.matrix(theta*(-Hstar.beta*Hstar/(1+theta*Hstar)^2+ H.beta*H/(1+theta*H)^2
                                                    -d*(H.beta/(1+theta*H)-theta*H.beta*H/(1+theta*H)^2)))))
    dl.dtheta <- cbind(dl.dtheta, dl.dtheta.dbeta)
    
    ### Derivative of gradient of beta with respect to all parameters
    dl.dbeta <- rbind(dl.dalpha.dbeta, dl.deta.dbeta, dl.dtheta.dbeta, dl.dbeta.dbeta)
    
    hessian <- cbind(t(dl.dalpha), t(dl.deta), t(dl.dtheta), dl.dbeta)
    #colnames(hessian) <- c("logalpha","logeta","logtheta","beta")
    #rownames(hessian) <- colnames(hessian)
    
    return(hessian)
  }
  
  ########################################
  ### Optimize the likelihood function
  fit <- optim(par=logp, fn=like, gr=gradientfunc, method="BFGS", hessian=FALSE)
  par <- fit$par
  
  ### Get the score function
  score <- gradientfunc(par, score=TRUE)
  
  ### Get the hessian
  hessian <- hessianfunc(par)
  
  #### Output ####
  out <- c(list(X=X, fit = fit, par = par, score = score, hessian = hessian, call = call,
                formula = formula, data = data, logp = logp, clusterid =clusterid,
                ncluster = ncluster, n = n))
  class(out) <- "frailty"
  return(out)
  
}

print("true")
print(logp)
estimates<-frailty(formula, data, logp, clusterid="id")
estimates


print.frailty<-function(x, ...){
  cat("Call:", "\n")
  print.default(x$call)
  cat("\nEstimated parameters in the Gamma-Weibull frailty model:", "\n")
  cat("\n")
  table.est <- t(as.matrix(x$par))
  colnames(table.est) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(x$X)))
  r <- rep("", , 1)
  rownames(table.est) <- c(r)
  print.default(table.est)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster)
  }

summary.frailty <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  
    ans <- list(AF = AF, CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n, n.cases = object$n.cases,
                n.cluster = object$n.cluster, modelcall = modelcall, call = object$call, method = method,
                formula = object$formula, exposure = object$exposure, outcome = object$outcome,
                fit = fit, sandwich = object$sandwich)
  }
  class(ans) <- "summary.AF"
  return(ans)
}



