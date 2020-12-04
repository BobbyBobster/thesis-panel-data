data.gen.ovb.power <- function (nobs, rho.lag, nperson = 5000) {
#nperson <- 5000
#nobs <- 10

fe.sd <- x.sd <- 1
u.sd <- 0.5
z.sd <- 1
theta.sd <- 0.6
eta.sd <- 0.6
theta.eta.cov <- rbind(
  c(0.6, -0.216),
  c(-0.216, 0.6))

beta.lag <- 0.05
beta1 <- 0 # intercept
beta2 <- 2
beta3 <- 1
beta.fe <- 1.5

gamma <- 1
#rho.lag <- 0.6
delta.lag <- 0.3

delta1 <- 0
delta2 <- 0.5

# Add a time index, first define a group apply function
# that applies by group index.
gapply <- function(x, group, fun) {
  returner <- numeric(length(group))
  for (idx in unique(group)) 
    returner[idx == group] <- get("fun")(x[idx == group])
  returner
}

# Data setup
constantdata <- data.frame(id = 1:nperson, fe = rnorm(nperson, sd = fe.sd))
fulldata <- constantdata[rep(1:nperson, each = nobs), ]
fulldata$t <- gapply(rep(1, length(fulldata$id)), 
  group = fulldata$id, 
  fun = cumsum)

for (idx in 0:(nperson - 1)) {
  pos <- 1 + idx * nobs
  ar.errs <- mvtnorm::rmvnorm(nobs, sigma = theta.eta.cov)
  fulldata$x2[pos:(pos + nobs - 1)] <- arima.sim(n = nobs, model = list(ar = rho.lag), innov = ar.errs[, 1])
  fulldata$z[pos:(pos + nobs - 1)] <- arima.sim(n  = nobs, model = list(ar = delta.lag), innov = ar.errs[, 2])
}

u <- rnorm(nobs * nperson, mean = 0, sd = u.sd)

fulldata$y <- 
  beta1 + 
  beta2 * fulldata$x2 + 
  gamma * fulldata$z +
  u

fulldata <- as.data.table(fulldata)
}
