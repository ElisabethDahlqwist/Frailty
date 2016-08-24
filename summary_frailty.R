CI <- function(Est, Std.Error, confidence.level, CI.transform){
  if(CI.transform == "untransformed"){
    lower <- Est - abs(qnorm((1 - confidence.level) / 2)) * Std.Error
    upper <- Est + abs(qnorm((1 - confidence.level) / 2)) * Std.Error
  }
  if(CI.transform == "log"){
    lower <- Est * exp( - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / Est)
    upper <- Est * exp(abs(qnorm((1 - confidence.level) / 2)) * Std.Error / Est)
  }
  if(CI.transform == "logit"){
    logit <- function(x) log(x / (1 - x))
    lower <- exp(logit(Est) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (Est * (1 - Est))) / (1 + exp(logit(Est) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (Est * (1 - Est))))
    upper <- exp(logit(Est) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (Est * (1 - Est))) / (1 + exp(logit(Est) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (Est * (1 - Est))))
  }
  CI <- cbind(lower, upper)
  return(CI)
}


print.frailty<-function(x, ...){
  cat("\nEstimated parameters in the Gamma-Weibull frailty model", "\n")
  cat("\n")
  cat("Call:", "\n")
  print.default(x$call)
  cat("\nCoefficients:")
  cat("\n")
  table.est <- as.matrix(cbind(x$par, x$se))
  rownames(table.est) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(x$X)))
  colnames(table.est) <- c("Estimate", "Std. Error")
  print.default(table.est)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster)
}

summary.frailty <- function(object, CI.transform, confidence.level, digits = max(3L, getOption("digits") - 3L), ...){
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  
  ## Inference
  ### Standard errors
  se <- sqrt(-diag(solve(object$hessian)))
  zvalue <- object$par / se
  pvalue <- 2 * pnorm( - abs(zvalue))
  confidence.interval <- CI(Est = object$par, Std.Error = se,
                            confidence.level = confidence.level,
                            CI.transform = CI.transform)
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")
  
  ## Score functions for each individual
  colnames(object$score) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(object$X)))
  cat("Score", "\n")
  print.default(object$score)
  cat("\n")
  
  ## Hessian
  colnames(object$hessian) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(object$X)))
  rownames(object$hessian) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(object$X)))
  cat("Hessian", "\n")
  print.default(object$hessian)
  
  ans <- list(par = object$par, se = se, zvalue = zvalue, pvalue = pvalue, score = object$score, X = object$X,
              hessian = object$hessian, call = object$call,
              formula = object$formula, modelmatrix = object$X, data = object$data,
              logp = object$logp, clusterid =object$clusterid, ncluster = object$ncluster,
              n = object$n, CI.transform = CI.transform, confidence.level = confidence.level,
              confidence.interval = confidence.interval)
  class(ans) <- "summary.frailty"
  return(ans)
}


print.summary.frailty <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...){
  ## Function call
  cat("Call:", "\n")
  print.default(x$call)
  cat("\nEstimated parameters in the Gamma-Weibull frailty model:", x$CI.transform, CI.text,  "Wald CI:", "\n")
  cat("\n")
  level <- x$confidence.level * 100
  CI.text <- paste0(as.character(level),"%")
  table.est <- cbind(x$par, x$se, x$zvalue, x$pvalue, x$confidence.interval)
  rownames(table.est) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003B8)", names(as.data.frame(x$X)))
  colnames(table.est) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Lower limit", "Upper limit")
  print.default(table.est)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster, "\n")
  cat("\n")
}

summary(estimates)
