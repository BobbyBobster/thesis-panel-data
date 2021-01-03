library(plm)
library(data.table)
library(tidyverse)

consistency.within.est.paper <- function(
  panel.set, 
  dep.var, 
  indep.vars, 
  group.index,
  time.index,
  alpha = 0.05,
  no.gen = FALSE) {
  
  # within estimator
  form <- as.formula(paste(dep.var, paste(indep.vars, collapse=" + "), sep=" ~ "))
  within.est <- plm(
    form,
    data = panel.set,
    model = "within")
  
  # function for setup
  diffMatrix <- function (timespan, a) {
    mat <- (-1) * diag(timespan)[ , 1:(timespan - a)]
    dif <- diag(timespan)[ , (1 + a):timespan]
    mat + dif
  }
  
  # setup
  if (is.null(panel.set$id)) panel.set["id"] <- panel.set[group.index]
  panel.set <- as.data.table(panel.set)
  
  n <- as.numeric(count(unique(panel.set[, ..group.index])))
  Tt <- as.numeric(count(unique(panel.set[, ..time.index])))
  k <- length(indep.vars)
  B <- diffMatrix(Tt - 1, 1)
  R <- kronecker(B, diag(k))
  
  
  beta.hat <- c()
  beta.hat.unstacked <- c()
  for (jdx in 1:(Tt - 1)) {
    Diff <- diffMatrix(Tt, jdx)
    Delta <- Diff %*% t(Diff)
    
    first.term.paper <- matrix(0, nrow = k, ncol = k)
    second.term.paper <- matrix(0, nrow = k, ncol = 1)
    
    for (idx in 1:n) {
      x.i <- data.matrix(panel.set[id == idx][, ..indep.vars])
      y.i <- data.matrix(panel.set[id == idx][, ..dep.var])
      first.term.paper <- first.term.paper + t(x.i) %*% Delta %*% x.i
      second.term.paper <- second.term.paper + t(x.i) %*% Delta %*% y.i
    }
    beta.j <- solve(first.term.paper) %*% second.term.paper
    
    beta.hat.unstacked <- cbind(beta.hat.unstacked, beta.j)
    beta.hat <- rbind(beta.hat, beta.j)
  }
  
  
  middle.term.wald <- matrix(0, nrow = (Tt - 2) * k, ncol = (Tt - 2) * k)
  middle.uu <- matrix(0, nrow = (Tt - 1) * k, ncol = (Tt - 1) * k)
  for (idx in 1:n) {
    x.i <- data.matrix(panel.set[id == idx][, ..indep.vars])
    y.i <- data.matrix(panel.set[id == idx][, ..dep.var])
    
    u.i <- c()
    for (jdx in 1:(Tt - 1)) {
      Diff <- diffMatrix(Tt, jdx)
      Delta <- Diff %*% t(Diff)
      
      first.term.u <- matrix(0, nrow = k, ncol = k)
      for (ldx in 1:n) {
        x.l <- data.matrix(panel.set[id == ldx][, ..indep.vars])
        first.term.u <- first.term.u + t(x.l) %*% Delta %*% x.l
      }
      
      e.ij <- y.i - x.i %*% beta.hat.unstacked[, jdx]
      
      row.u <- solve(first.term.u) %*% t(x.i) %*% Delta %*% e.ij
      u.i <- rbind(u.i, row.u)
    }
    
    middle.uu <- middle.uu + u.i %*% t(u.i)
    middle.term.wald <- middle.term.wald + t(R) %*% u.i %*% t(u.i) %*% R
  }
  
  Sigma.hat <-  middle.uu
  beta.std.err <- sqrt(diag(Sigma.hat) / n)
  
  
  # Wald statistic 
  wald <- t(beta.hat) %*% R %*% MASS::ginv(middle.term.wald) %*% t(R) %*% beta.hat
  chisq.bound <- qchisq(1 - alpha, df = (Tt - 2) * k)
  p.value <- 1 - pchisq(wald, df = (Tt - 2) * k)
  
  
  # Gen Wald statistic
  if (!no.gen) {
    gen.wald <- t(beta.hat) %*% R %*% t(R) %*% beta.hat
    
    Sigma.R <- t(R) %*% Sigma.hat %*% R
    evs <- eigen(Sigma.R)$values
    
    gen.chisq.bound <- NA
    gen.p.value <- CompQuadForm::imhof(q = gen.wald, lambda = evs)$Qq
  } else {
    gen.wald <- NA
    gen.chisq.bound <- NA
    gen.p.value <- NA
  }
  
  
  
  
  # necessary for plotting
  beta.hat.unstacked <- matrix(nrow = k, ncol = Tt - 1)
  beta.std.err.unstacked <- matrix(nrow = k, ncol = Tt - 1)
  rownames(beta.hat.unstacked) <- indep.vars
  for (idx in 1:(Tt - 1)) {
    rs <- 1 + (idx - 1) * k
    re <- rs + (k - 1)
    beta.hat.unstacked[, idx] <- beta.hat[rs:re, ]
    beta.std.err.unstacked[, idx] <- beta.std.err[rs:re]
  }
  
  # return object
  consist.test <- list(
    k = k, n = n, Tt = Tt, R = R, B = B, alpha = alpha,
    within.est = within.est,
    beta.hat = beta.hat,
    beta.std.err = beta.std.err,
    beta.hat.unstacked = beta.hat.unstacked,
    beta.std.err.unstacked = beta.std.err.unstacked,
    Sigma.hat = Sigma.hat,
    #rss = rss, tss = tss, R.sq = R.sq,
    wald = wald,
    chisq.bound = chisq.bound,
    p.value = p.value,
    gen.wald = gen.wald,
    gen.chisq.bound = gen.chisq.bound,
    gen.p.value = gen.p.value)
  class(consist.test) <- "ConsistencyWithinEstimator"
  consist.test
}
