data.gen.ovb.size <- function (nobs, rho.lag, nperson = 5000) {
  fe.sd <- 1
  eps.sd <- 0.5
  theta.sd <- 0.6
  eta.sd <- 0.6
  # this should be changed to use theta.sd**2 and eta.sd**2 (ONLY IF THERE IS TIME)
  theta.eta.cov <- rbind(
    c(0.6, -0.216),
    c(-0.216, 0.6))
  
  beta1 <- 0 # intercept
  beta2 <- 2
  
  gamma <- 0 # no z term, so no ommitted variable -> to calc size
  delta.lag <- 0.3
  
  
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
    fulldata$x2[pos:(pos + nobs - 1)] <- 
      arima.sim(n = nobs, model = list(ar = rho.lag), innov = ar.errs[, 1])
    #fulldata$z[pos:(pos + nobs - 1)] <- 
    #  arima.sim(n  = nobs, model = list(ar = delta.lag), innov = ar.errs[, 2])
  }
  
  eps <- rnorm(nobs * nperson, mean = 0, sd = eps.sd)
  
  fulldata$y <- 
    beta1 + 
    beta2 * fulldata$x2 + 
    #gamma * fulldata$z +
    eps
  
  fulldata <- as.data.table(fulldata)
}
