data.gen.si.size <- function (nobs, rho.lag, nperson = 5000) {
  fe.x.sd <- 1
  fe.y.sd <- 1
  eps.sd <- 2
  theta.sd <- 1
  
  beta <- 1
  alpha <- 0 # alpha zero so no simulaneity
  
  
  gapply <- function(x, group, fun) {
    returner <- numeric(length(group))
    for (idx in unique(group)) 
      returner[idx == group] <- get("fun")(x[idx == group])
    returner
  }
  
  # Data setup
  constantdata <- data.frame(id = 1:nperson, 
    fe.x = rnorm(nperson, sd = fe.x.sd), 
    fe.y = rnorm(nperson, sd = fe.y.sd))
  fulldata <- constantdata[rep(1:nperson, each = nobs), ]
  fulldata$t <- gapply(rep(1, length(fulldata$id)), 
    group = fulldata$id, 
    fun = cumsum)
  
  for (idx in 0:(nperson - 1)) {
    pos <- 1 + idx * nobs
    fulldata$u[pos:(pos + nobs - 1)] <- 
      arima.sim(n = nobs, model = list(ar = rho.lag), sd = theta.sd)
  }
  
  eps <- rnorm(nobs * nperson, sd = eps.sd) 
  
  fulldata$y <- 
    (fulldata$fe.y + beta * fulldata$fe.x) / (1 - alpha * beta) +
    (beta * fulldata$u + eps) / (1 - alpha * beta)
    
  fulldata$x2 <-
    (fulldata$fe.x + alpha * fulldata$fe.y) / (1 - alpha * beta) +
    (fulldata$u + alpha * eps) / (1 - alpha * beta)
  
  fulldata <- as.data.table(fulldata)
}
