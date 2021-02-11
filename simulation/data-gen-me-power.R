data.gen.me.power <- function (nobs, rho.lag, nperson = 5000, delta.lag = 0.3) {
  fe.sd <- 1
  eps.sd <- 1
  theta.sd <- 1.2
  eta.sd <- 0.8
  theta.eta.cov <- rbind(
    c(theta.sd**2, 0),
    c(0, eta.sd**2))
  
  beta <- 1
  #delta.lag <- 0.3
  
  
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
    fulldata$ksi[pos:(pos + nobs - 1)] <- 
      arima.sim(n = nobs, model = list(ar = rho.lag), innov = ar.errs[, 1])
    fulldata$nu[pos:(pos + nobs - 1)] <- 
      arima.sim(n  = nobs, model = list(ar = delta.lag), innov = ar.errs[, 2])
  }
  
  fulldata$x2 <- fulldata$ksi + fulldata$nu
  
  fulldata$eps <- rnorm(nobs * nperson, mean = 0, sd = eps.sd)
  
  fulldata$y <- 
    fulldata$fe + 
    beta * fulldata$ksi + 
    fulldata$eps
  
  fulldata <- as.data.table(fulldata)
}