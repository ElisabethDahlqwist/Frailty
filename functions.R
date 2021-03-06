b <- matrix(c(2:9), ncol=2,nrow=4)
a <- c(1,2)

bb <- b[rep(1:nrow(b), ncol(b)), ]
bbb <- c(b)*bb






idd <- c(1,1,2,2)
iddd <- rep(idd, ncol(b))
par_rep <- sort(rep(1:ncol(b), nrow(b)))
temp <- data.table(bbb,iddd, par_rep)
temp <- as.matrix(temp[, j = lapply(.SD, sum), by = .(par_rep, iddd)])[, -1]
data1 <- matrix(c(temp[, 1], a*temp[,2:(ncol(b)+1)]), ncol=3, nrow=nrow(b))

###############
H <- (t/alpha)^eta*exp(B) 
Hstar <- (tstar/alpha)^eta*exp(B)
XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
nbeta_rep <- rep(1:nbeta, each = nrow(X))
#clusterid_rep <- rep(clusterid, nbeta)

H <- (t/alpha)^eta*exp(B) 
Hstar <- (tstar/alpha)^eta*exp(B)
H_X <- temp(H * X)
H_Xsqr <- H_X[rep(1:nrow(H_X), nbeta), ] * c(H_X)
Hstar_X <- temp(Hstar * X)
Hstar_Xsqr <- Hstar_X[rep(1:nrow(H_X), nbeta), ] * c(Hstar_X)

Hstar_first <- data.table(clusterid, nbeta_rep, Hstar * XX)
H_first <- data.table(clusterid, nbeta_rep, H * XX)

H <- temp(H)
Hstar <- temp(Hstar)
Hstar_second <- theta * Hstar_Xsqr / (1 + theta * Hstar)^2
H_second <- theta * (1+d*theta) * H_Xsqr / (1 + theta * H)^2

data_Hstar <- data.table(clusterid, nbeta_rep, Hstar_first)
data_Hstar_aggr <- as.matrix(data_Hstar[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])

data_H <- data.table(clusterid, nbeta_rep, H_first)
data_H_aggr <- as.matrix(data_H[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])

Hstar_first <- data_Hstar_aggr[, -1:-2] / (1 + theta * Hstar)
H_first <- theta * (1 + d * theta) * data_H_aggr[, -1:-2] / (1 + theta * H)


dl.dbeta.dbeta <- data_Hstar_aggr[, -1:-2] / (1 + theta * Hstar) - theta * Hstar_Xsqr / (1 + theta * Hstar)^2
- theta * (1 + d * theta) * data_H_aggr[, -1:-2] / (1 + theta * H) + theta * (1+d*theta) * H_Xsqr / (1 + theta * H)^2

## aggregate over clusters

nbeta_rep2 <- rep(1:nbeta, each = length(H))
dl.dbeta.dbeta <- data.table(nbeta_rep2, dl.dbeta.dbeta)
dl.dbeta.dbeta <- as.matrix(dl.dbeta.dbeta[, j = lapply(.SD, sum), by = .(nbeta_rep2)])[, -1]
dl.dbeta.dbeta

######################################

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

dl.dbeta.dbeta <- -t(colMeans(Hstar.beta.beta/(1+theta*Hstar)-
                                theta*(Hstar.beta/(1+theta*Hstar))^2-(1+theta*d)*
                                (H.beta.beta/(1+theta*H)-theta*(H.beta/(1+theta*H))^2)))

dl.dbeta <- cbind(t(dl.dalpha.dbeta), t(dl.deta.dbeta), t(dl.dtheta.dbeta))

hessian <- rbind(dl.dalpha, dl.deta, dl.dtheta, dl.dbeta)

################ old stuff

apply(b, 1, prod) 

nbeta <- 3

sqpart <- b*b

rbind(b*b[,1], b*b[,2], b*b[,3])

mprod <- function(b, i) b*b[,i]
apply(b, 1, mprod, i="1")

x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
cave <- function(x, c1, c2) c(mean(x[c1]), mean(x[c2]))
apply(x, 1, cave,  c1 = "x1", c2 = c("x1","x2"))

kk<- matrix(data=NA, 3,3)
hfunc <- function(X, nbeta){
  for(i in 1:nbeta){
    kk <- X*X[,i]
    return(kk)
  }
}

kronecker(b[,1], b)