library(CompQuadForm)
library(MASS)
library(tidyverse)

# Simulate some normals
y.mean <- c(5, 6, 7)
y.vcov <- rbind(
  c(1, 0.5, 0.3),
  c(0.5, 1, 0.5),
  c(0.3, 0.5, 1))  

A <- rbind(
  c(1, 2, 3),
  c(2, 1, 2),
  c(3, 2, 1))
quads <- rep(NA, 10000)
for (ix in seq(1, 10000)) {
  ys <- mvrnorm(n = 1, mu = y.mean, Sigma = y.vcov)
  
  quads[ix] <- t(ys) %*% A %*% ys
}
plot(density(quads))

K <- eigen(y.vcov)$vectors
Delta <- diag(eigen(y.vcov)$values)
K %*% Delta %*% t(K)

Gamma <- sqrt(Delta) %*% t(K) %*% A %*% K %*% sqrt(Delta)
U <- eigen(Gamma)$vectors
Lambda <- diag(eigen(Gamma)$values)
U %*% Lambda %*% t(U)

P <- K %*% sqrt(Delta) %*% U + (diag(3) - y.vcov %*% solve(y.vcov)) %*% A %*% K %*% sqrt(Delta) %*% U %*% solve(Lambda)
delta <- t(y.mean) %*% (diag(3) - solve(y.vcov) %*% y.vcov) %*%
        (A - A %*% y.vcov %*% solve(y.vcov %*% A %*% y.vcov) %*% y.vcov %*% A) %*%
        (diag(3) - y.vcov %*% solve(y.vcov)) %*% y.mean


# https://math.stackexchange.com/questions/3700360/whats-the-distribution-of-xyxzyz-where-x-y-z-are-independent-standard-no
B <- rbind(
  c(0, 1, 1),
  c(1, 0, 1),
  c(1, 1, 0))
ws <- t(mvrnorm(n = 1000, mu = c(0, 0, 0), Sigma = diag(3)))
w.quad <- (t(ws) %*% B %*% ws)
plot(density(w.quad))

# 
{
  xs <- seq(-20, 20, 0.5)
  n <- length(xs)
  points <- rep(NA, n)
  for (idx in 1:n) {
    points[idx] <- 
      (1 - 
          imhof(xs[idx], 
            lambda = diag(Lambda),
            )$Qq)
  }
  plot(y = points, x = xs, type = "l", main = "Gen. Chi^2 cdf")
  pdf <- diff(points) / diff(xs)
  plot(pdf, x = xs[2:n], type = "l", main = "Gen. Chi^2 pdf")
}

