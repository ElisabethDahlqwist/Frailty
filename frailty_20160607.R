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



#aggregate

XX <- X[rep(1:nrow(X), nbeta), ]
XXX <- c(X)*XX
clusterid_rep <- rep(clusterid, nbeta)
nbeta_rep <- sort(rep(1:nbeta, nrow(X)))
data_rep <- data.table(XXX,clusterid_rep, nbeta_rep)
data_rep <- as.matrix(data_rep[, j = lapply(.SD, sum), by = .(nbeta_rep, clusterid_rep)])
data_new <- matrix(c(data_rep[, 1], data_rep[, 2], H * data_rep[, 3:(nbeta+2)]), ncol = (nbeta+2), byrow=F)


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