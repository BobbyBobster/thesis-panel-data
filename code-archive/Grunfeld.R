library(plm)
library(data.table)
library(tidyverse)

data("Grunfeld", package = "plm")

grunfeld <- pdata.frame(Grunfeld, 
  index = c("firm", "year"), 
  drop.index = FALSE)
grunfeld <- as.data.table(grunfeld)

form <- inv ~ value + capital

#plot(form, data = grunfeld, col = firm)

grun.fe <- plm(form,
  data = grunfeld,
  model = "within")
fe.effects <- fixef(grun.fe, type = "level")
grun.re <- update(grun.fe, model = "random")
phtest(grun.fe, grun.re)

## Setup ====
diffMatrix <- function (timespan, a) {
  mat <- (-1) * diag(timespan)[ , 1:(timespan - a)]
  dif <- diag(timespan)[ , (1+a):timespan]
  mat + dif
}


panel.set <- grunfeld

iv <- c("inv")
dv <- c("value", "capital")
gid <- c("firm")
tid <- c("year")

n <- as.numeric(count(unique(panel.set[, ..gid])))
Tt <- as.numeric(count(unique(panel.set[, ..tid])))
k <- length(dv)
B <- diffMatrix(Tt, 1)
R <- kronecker(B, diag(k))


beta.hat <- rep(NA, times = k * (Tt - 1))
first.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = (Tt - 1) * k)
second.term.appC <- matrix(0, nrow = (Tt - 1) * k, ncol = 1)
for (idx in 1:n) {
  #x <- panel.set %>% filter(!!group.index == idx) %>% select(!!!dep.vars) %>%  data.matrix()
  x <- data.matrix(panel.set[firm == idx][, ..dv])
  stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
  for (timespan in 1:(Tt - 1)) {
    rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
    rowend <- rowstart + (Tt - timespan - 1)
    colstart <- 1 + (timespan - 1) * k
    colend <- colstart + (k - 1)
    stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
  }
  
  #y <- panel.set %>% filter(!!group.index == idx) %>% select(!!indep.var) %>% data.matrix()
  y <- data.matrix(panel.set[firm == idx][, ..iv])
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
  #x <- panel.set %>% filter(!!group.index == idx) %>% select(!!!dep.vars) %>%  data.matrix()
  x <- data.matrix(panel.set[firm == idx][, ..dv])
  stacked.x <- matrix(0, nrow = (Tt - 1) * Tt / 2, ncol = (Tt - 1) * k)
  for (timespan in 1:(Tt - 1)) {
    rowstart <- 1 + ifelse(timespan != 1, sum((Tt - 1):(Tt - timespan + 1)), 0)
    rowend <- rowstart + (Tt - timespan - 1)
    colstart <- 1 + (timespan - 1) * k
    colend <- colstart + (k - 1)
    stacked.x[rowstart:rowend, colstart:colend] <- t(diffMatrix(Tt, timespan)) %*% x
  }
  #y <- panel.set %>% filter(!!group.index == idx) %>% select(!!indep.var) %>% data.matrix()
  y <- data.matrix(panel.set[firm == idx][, ..iv])
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
rownames(beta.hat.unstacked) <- c("value", "capital")
for (idx in 1:(Tt - 1)) {
  rs <- 1 + (idx - 1) * k
  re <- rs + (k - 1)
  beta.hat.unstacked[, idx] <- beta.hat[rs:re, ]
  beta.std.err.unstacked[, idx] <- beta.std.err[rs:re]
}

# Plotting ====
{

}

wald <- n * (t(beta.hat) %*% t(R) %*% MASS::ginv(R %*% Sigma.hat %*% t(R)) %*% R %*% beta.hat)
wald
pchisq(wald, df = Tt - 2)


