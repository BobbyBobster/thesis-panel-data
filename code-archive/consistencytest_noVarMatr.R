library(plm)
library(data.table)
library(tidyverse)

consistency.novar <- function(
  panel.set, 
  dep.var, 
  indep.vars, 
  group.index,
  time.index,
  alpha = 0.05) {
  
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
  B <- diffMatrix(Tt, 1)
  R <- kronecker(B, diag(k))
  
  # calcs
  beta.hat <- rep(NA, times = k * (Tt - 1))
  first.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = (Tt - 1) * k)
  second.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = 1)
  for (idx in 1:n) {
    x <- data.matrix(panel.set[id == idx][, ..indep.vars])
    stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
    for (timespan in 1:(Tt - 1)) {
      rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
      rowend <- rowstart + (Tt - timespan - 1)
      colstart <- 1 + (timespan - 1) * k
      colend <- colstart + (k - 1)
      stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
    }
    
    y <- data.matrix(panel.set[id == idx][, ..dep.var])
    stacked.y <- c()
    for (timespan in 1:(Tt - 1)) {
      stacked.y <- rbind(stacked.y, t(diffMatrix(Tt, timespan)) %*% y)
    }
    
    first.term.appC <- first.term.appC + t(stacked.x) %*% stacked.x
    second.term.appC <- second.term.appC + t(stacked.x) %*% stacked.y
  }
  beta.hat <- solve((1 / n) * first.term.appC) %*% ((1 / n) * second.term.appC)
  
  # vcov
  inner.vcov <- matrix(0, nrow = (Tt - 1) * k, ncol = (Tt - 1) * k)
  for (idx in 1:n) {
    x <- data.matrix(panel.set[id == idx][, ..indep.vars])
    stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
    for (timespan in 1:(Tt - 1)) {
      rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
      rowend <- rowstart + (Tt - timespan - 1)
      colstart <- 1 + (timespan - 1) * k
      colend <- colstart + (k - 1)
      stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
    }
    y <- data.matrix(panel.set[id == idx][, ..dep.var])
    stacked.y <- c()
    for (timespan in 1:(Tt - 1)) {
      stacked.y <- rbind(stacked.y, t(diffMatrix(Tt, timespan)) %*% y)
    }
    
    u.i <- stacked.y - stacked.x %*% beta.hat
    
    inner.vcov <- inner.vcov + t(stacked.x) %*% u.i %*% t(u.i) %*% stacked.x
  }
  
  Sigma.hat <- solve((1 / n) * first.term.appC) %*% ((1 / n) * inner.vcov) %*% solve((1 / n) * first.term.appC)
  
  beta.std.err <- sqrt(diag(Sigma.hat) / n)
  
  
  
  # plotting
  # Keep current graphics settings
  oldpar <- par(no.readonly = TRUE)
  
  beta.hat.unstacked <- matrix(nrow = k, ncol = Tt - 1)
  beta.std.err.unstacked <- matrix(nrow = k, ncol = Tt - 1)
  rownames(beta.hat.unstacked) <- indep.vars
  for (idx in 1:(Tt - 1)) {
    rs <- 1 + (idx - 1) * k
    re <- rs + (k - 1)
    beta.hat.unstacked[, idx] <- beta.hat[rs:re, ]
    beta.std.err.unstacked[, idx] <- beta.std.err[rs:re]
  }
  
  
  wald <- t(beta.hat) %*% t(R) %*% R %*% beta.hat
  chisq.bound <- qchisq(1 - alpha, df = (Tt - 2) * k)
  
  
  consist.test <- list(
    within.est = within.est,
    beta.hat = beta.hat,
    beta.std.err = beta.std.err,
    beta.hat.unstacked = beta.hat.unstacked,
    beta.std.err.unstacked = beta.std.err.unstacked,
    Sigma.hat = Sigma.hat,
    wald = wald,
    chisq.bound = chisq.bound,
    k = k, n = n, Tt = Tt, R = R, B = B)
  class(consist.test) <- "ConsistencyNoVar"
  consist.test
}

# Plotting ====
plot.ConsistencyNoVar <- function (object)
{
  oldpar <- par(no.readonly = TRUE)
  
  # Create nice rectangle of plots
  plt.row <- round(sqrt(object$k))
  plt.col <- ceiling(sqrt(object$k))
  par(
    mfrow = c(plt.row, plt.col),
    mar = c(5.1, 4.1, 0.1, 2.1),
    oma = c(0, 0, 4, 0)
  )
  
  for (row in 1:object$k) {
    plot(object$beta.hat.unstacked[row, ], 
      type = "b", 
      xlab = "timespan",
      ylab = rownames(object$beta.hat.unstacked)[row],
      ylim = c(
        min(object$beta.hat.unstacked[row, ]) - abs(max(object$beta.std.err.unstacked[row, ]) * 2),
        max(object$beta.hat.unstacked[row, ]) + abs(max(object$beta.std.err.unstacked[row, ]) * 2)))
    abline(h = object$within.est$coefficients[row], lty = 3)
    abline(h = 0, lty = 2)
    arrows(
      x0 = 1:(object$Tt - 1), 
      y0 = object$beta.hat.unstacked[row, ] - 1.96 * object$beta.std.err.unstacked[row, ],
      x1 = 1:(object$Tt - 1),
      y1 = object$beta.hat.unstacked[row, ] + 1.96 * object$beta.std.err.unstacked[row, ],
      code = 3,
      angle = 90,
      length = 0.05)
  }
  # reset old graphics params
  par(oldpar)
}
