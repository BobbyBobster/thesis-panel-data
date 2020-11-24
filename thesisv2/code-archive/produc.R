library(tidyverse)
library(plm)
data("Produc", package="plm")

Produc <- Produc %>% 
  mutate(
    lgsp = log(gsp),
    lpcap = log(pcap),
    lpc = log(pc),
    lemp = log(emp),
    lunemp = log(unemp))

produc <- pdata.frame(Produc, index = c("state", "year"))


## Setup ====
diffMatrix <- function (timespan, a) {
  mat <- (-1) * diag(timespan)[ , 1:(timespan - a)]
  dif <- diag(timespan)[ , (1+a):timespan]
  mat + dif
}

## Setup ====
panel.set <- produc

Tt <- 17
n <- 48

indep.var <- quo(lgsp)
dep.vars <- quos(lpcap, lpc, lemp, lunemp)
group.index <- quo(state)
time.index <- quo(year)

k <- length(dep.vars)
#B <- rbind(diag(Tt - 1), 0) - rbind(0, diag(Tt - 1))
B <- diffMatrix(Tt, 1)
R <- kronecker(B, diag(k))




beta.hat <- rep(NA, times = k * (Tt - 1))
first.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = (Tt - 1) * k)
second.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = 1)
for (idx in levels(produc$state)) {
  x <- panel.set %>% filter(!!group.index == idx) %>% select(!!!dep.vars) %>%  data.matrix()
  stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
  for (timespan in 1:(Tt - 1)) {
    rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
    rowend <- rowstart + (Tt - timespan - 1)
    colstart <- 1 + (timespan - 1) * k
    colend <- colstart + (k - 1)
    stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
  }
  
  y <- panel.set %>% filter(!!group.index == idx) %>% select(!!indep.var) %>% data.matrix()
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
for (idx in levels(produc$state)) {
  x <- panel.set %>% filter(!!group.index == idx) %>% select(!!!dep.vars) %>%  data.matrix()
  stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
  for (timespan in 1:(Tt - 1)) {
    rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
    rowend <- rowstart + (Tt - timespan - 1)
    colstart <- 1 + (timespan - 1) * k
    colend <- colstart + (k - 1)
    stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
  }
  y <- panel.set %>% filter(!!group.index == idx) %>% select(!!indep.var) %>% data.matrix()
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
rownames(beta.hat.unstacked) <- c("lpcap", "lpc", "lemp", "lunemp")
for (tdx in 1:(Tt - 1)) {
  rs <- 1 + (tdx - 1) * k
  re <- rs + (k - 1)
  beta.hat.unstacked[, tdx] <- beta.hat[rs:re, ]
  beta.std.err.unstacked[, tdx] <- beta.std.err[rs:re]
}

# Plotting ====
{
  par(
    mfrow = c(2, 2),
    mar = c(5.1, 4.1, 0.1, 2.1),
    oma = c(0, 0, 4, 0)
  )
  for (row in 1:k) {
    plot(beta.hat.unstacked[row, ], 
      type = "b", 
      xlab = "timespan",
      ylab = rownames(beta.hat.unstacked)[row],
      ylim = c(
        min(beta.hat.unstacked[row, ]) - abs(max(beta.std.err.unstacked[row, ]) * 2),
        max(beta.hat.unstacked[row, ]) + abs(max(beta.std.err.unstacked[row, ]) * 2)))
    abline(h = mean(beta.hat.unstacked[row, ]), lty = 3)
    abline(h = 0, lty = 2)
    arrows(
      x0 = 1:(Tt - 1), 
      y0 = beta.hat.unstacked[row, ] - 1.96 * beta.std.err.unstacked[row, ],
      x1 = 1:(Tt - 1),
      y1 = beta.hat.unstacked[row, ] + 1.96 * beta.std.err.unstacked[row, ],
      code = 3,
      angle = 90,
      length = 0.05)
  }
  par(oldpar)
}

B.1 <- R
wald <- n * (t(beta.hat) %*% t(B.1) %*% MASS::ginv(B.1 %*% Sigma.hat %*% t(B.1)) %*% B.1 %*% beta.hat)
wald
pchisq(wald, df = Tt - 2)



