library(mvtnorm)
library(CompQuadForm)
library(doParallel)

# run coru_data first

cwe <- cnv <- cwe.wages

n <- cnv$n
R <- cnv$R
Tt <- cnv$Tt
Sigma.hat <- cnv$Sigma.hat
  
Sigma.R <- t(R) %*% Sigma.hat %*% R
evs <- eigen(Sigma.R)$values

CompQuadForm::imhof(q = cwe$gen.wald, lambda = evs)$Qq




evs <- c(4,1,0.5)

find.alpha.bound <- function (evs, alpha = 0.05) {
  imhof.value <- function (point) { 
    (CompQuadForm::imhof(q = point, lambda = evs)$Qq - alpha)**2
  }
  
  optim(0.1, 
    imhof.value,
    method = "BFGS")
}



# Simulate gen. chi^2 vars with estimated vcov matrix from wages dataset
N <- 2000
k <- rep(NA, N)
x <- rep(NA, N)
x.o <- rep(NA, N)
for (idx in 1:N) {
  obs <- mvtnorm::rmvnorm(n = 1, sigma = Sigma.R)
  x[idx] <- obs %*% t(obs)
  
  obs.o <- rnorm(n = 3)
  x.o[idx] <- sum(evs * obs.o**2)
}
plot(density(x), col = "red")
lines(density(x.o), col = "blue")
#





cl <- makeCluster(3, type="FORK")
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = 2020)
{
  N <- 10
  xs <- rep(NA, N)
  xs <- foreach(1:N) %dopar% {
    obs <- mvtnorm::rmvnorm(n = 1, sigma = Sigma.R)
    obs %*% t(obs)
  }
  #xs <- as.numeric(xs)
}
#{
#  N <- 5000
#  sumchisim <- rep(NA, N)
#  sumchisim <- foreach(1:N) %dopar% {
#    chi <- rchisq(n = 63, df = 1)
#    t(evs) %*% chi
#  }
#  sumchisim <- as.numeric(sumchisim)
#}
stopCluster(cl)
rm(cl)


# package CompQuadForm for the "true" distribution
{
  step.size <- 0.1
  step.end <- 20
  pts.imhof <- rep(NA, step.end / step.size)
  xs.imhof <- seq(0, step.end - step.size, by = step.size)
  for (point in xs.imhof) {
    pts.imhof[point * (1 / step.size) + 1] <- 
      (1 - CompQuadForm::imhof(q = point, lambda = evs)$Qq)
  }
  
  grad.imhof <- pracma::gradient(pts.imhof, h1 = step.size)
}
#pts.me <- rep(NA, end * step.size)
#for (point in sequence) {
#  pts.me[point * (1 / step.size)] <- 
#    sum(evs * (1 - pchisq(point, 1)))
#}


# 
Tt <- cnv$Tt
k <- cnv$k

V <- Sigma.hat
evd.V <- eigen(Sigma.hat)
Delta <- diag(evd.V$values)
K <- evd.V$vectors
Gamma <- sqrt(Delta) %*% t(K) %*% t(R) %*% R %*% K %*% sqrt(Delta)
evd.Gamma <- eigen(Gamma)
Lambda <- diag(evd.Gamma$values)
U <- evd.Gamma$vectors
P <- K %*% sqrt(Delta) %*% U +
  (diag((Tt - 1) * k) - V %*% solve(V)) %*% 
  t(R) %*% R %*% K %*% sqrt(Delta) %*% U %*% solve(Lambda)
zs <- rmvnorm(
  n = 10000, 
  mean = t(P) %*% beta.0)



#

plot(
  x = xs.imhof,
  y = pracma::gradient(pts.imhof, h1 = step.size), 
  col = "green",
  type = "l",
  lty = 2)
lines(
  density(xs, bw = "SJ"), 
  col = "blue",
  main = "Gen. Chi^2: Simulation vs package",
  xlab = "x" )
{
  if (step.size >= 1) {
    crit.val <- length(pts.imhof) - length(pts.imhof[pts.imhof > 0.95])
  } else if (step.size < 1 & step.size > 0) {
    crit.val <- (length(pts.imhof) - length(pts.imhof[pts.imhof > 0.95])) * step.size
  }
  abline( # crit value of test stat
    v = crit.val,
    col = "blue")
}
abline(
  v = cnv$wald,
  col = "green")
legend(
  x = 150,
  y = 0.025,
  legend = c("Simulated", "imhof"),
  col = c("blue", "red"), 
  lty = c(1, 2))

