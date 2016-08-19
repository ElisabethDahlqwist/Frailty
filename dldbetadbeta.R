# Hessian, derivative of gradient of beta with respect to all parameters
H <- (t/alpha)^eta*exp(B) 
Hstar <- (tstar/alpha)^eta*exp(B)
XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
nbeta_rep <- rep(1:nbeta, each = nrow(X))
#clusterid_rep <- rep(clusterid, nbeta)

### Creating squared cluster sums of H.beta and Hstar.beta
H.beta <- as.matrix(temp(H * X))
H.beta2 <- H_X[rep(1:nrow(H_X), nbeta), ] * c(H_X)
Hstar.beta <- as.matrix(temp(Hstar * X))
Hstar.beta2 <- Hstar_X[rep(1:nrow(Hstar_X), nbeta), ] * c(Hstar_X)

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
Hstar.beta.beta <- as.matrix(data_Hstar[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])[, -1:-2, drop=FALSE] # because columns are droped it is no longer a matrix

H.beta.beta <- data.table(clusterid, nbeta_rep, H.beta.beta)
H.beta.beta <- as.matrix(H.beta.beta[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid)])[, -1:-2, drop=FALSE] 

### Calculate Hstar1 <- Hstar.beta.beta/(1+theta*Hstar) and H1 <- theta * (1 + d * theta)*H.beta.beta/(1+theta*H)
Hstar1 <- Hstar.beta.beta / (1 + theta * Hstar)
H1 <- (1 + d * theta) * H.beta.beta / (1 + theta * H)

dl.dbeta.dbeta <- colMeans(Hstar.beta.beta/(1+theta*Hstar)-
                         theta*Hstar.beta2/(1+theta*Hstar)^2-(1+theta*d)*
                         (H.beta.beta/(1+theta*H)-theta*H.beta2/(1+theta*H)^2))


test7 <- matrix(c(
        Hstar.beta.beta,
        Hstar.beta2,
        Hstar,
        H.beta.beta,
        H.beta2,
        H,
        theta*d), ncol=7)

## aggregate over clusters

nbeta_rep2 <- rep(1:nbeta, each = length(H))
dl.dbeta.dbeta <- data.table(nbeta_rep2, dl.dbeta.dbeta)
dl.dbeta.dbeta <- as.matrix(dl.dbeta.dbeta[, j = lapply(.SD, sum), by = .(nbeta_rep2)])[, -1]

### Second derivative of log-likelihood with respect to beta
dl.dbeta <- cbind(t(dl.dalpha.dbeta), t(dl.deta.dbeta), t(dl.dtheta.dbeta), dl.dbeta.dbeta)

hessian <- rbind(dl.dalpha, dl.deta, dl.dtheta, dl.dbeta)
