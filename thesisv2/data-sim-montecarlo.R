

## econometricsbysimulation.com
# and http://online.sfsu.edu/mbar/ECON312_files/OmittedSimul.html
nperson <- 100
nobs <- 5

# Main equation 
#fe.sd <- 0.6019
#x.sd <- 0.3600
#u.sd <- 0.3982
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
rho.lag <- 0.6
delta.lag <- 0.3


# Auxiliary equation
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



{
B <- 100
accepts <- 0
p.vals <- rep(0, B)
gen.accepts <- 0
gen.p.vals <- rep(0, B)
within.ests <- matrix(0, nrow = B, ncol = 2)
pooled.ests <- matrix(0, nrow = B, ncol = 2)
gen.walds <- rep(0, B)

  # First generate our data using the time constant effects
  constantdata <- data.frame(id = 1:nperson, fe = rnorm(nperson, sd = fe.sd))
  
  # We expand our data by nobs
  fulldata <- constantdata[rep(1:nperson, each = nobs), ]
  
  # Using the generalized apply function coded above
  fulldata$t <- gapply(rep(1, length(fulldata$id)), 
    group = fulldata$id, 
    fun = cumsum)
  

  for (bdx in 1:B) {
    # Now we are ready to calculate the time variant xs
    # And our unobservable error
    # included regressor
    #fulldata$x2 <- fulldata$fe +rnorm(nobs * nperson, sd = x.sd)
    u <- rnorm(nobs * nperson, mean = 0, sd = u.sd)
    
    # omitted regressor
    #eps <- rnorm(nobs * nperson, mean = 0, sd = 1)
    #x3 <- delta1 + delta2 * fulldata$x2 + eps
    
    
for (idx in 0:(nperson - 1)) {
  pos <- 1 + idx * nobs
  ar.errs <- mvtnorm::rmvnorm(nobs, sigma = theta.eta.cov)
  fulldata$x2[pos:(pos + nobs - 1)] <- arima.sim(n = nobs, model = list(ar = rho.lag), innov = ar.errs[, 1])
  fulldata$z[pos:(pos + nobs - 1)] <- arima.sim(n  = nobs, model = list(ar = delta.lag), innov = ar.errs[, 2])
}
    # Finally we are ready to simulate our y variables
    fulldata$y <- 
      beta1 + 
      beta2 * fulldata$x2 + 
      gamma * fulldata$z +
      u
    
    
    fulldata.plm <- pdata.frame(fulldata, index = c("id", "t"))
    cwe.sim <- consistency.within.est(panel.set = fulldata.plm, 
      dep.var = c("y"), indep.vars = c("x2"), 
      group.index = c("id"), time.index = "t",
      no.gen = TRUE)
    
    if (cwe.sim$wald < cwe.sim$chisq.bound) { 
      accepts <- accepts + 1 
    }
    p.vals[bdx] <- cwe.sim$p.value
    within.ests[bdx, ] <- cwe.sim$within.est$coefficients
    pooled.ests[bdx, ] <- plm(y ~ x2, data = fulldata.plm, model = "pooling")$coefficients
    
    #gen.walds[bdx] <- cwe.sim$gen.wald
    #if (gen.walds[bdx] < cwe.sim$gen.chisq.bound) {
    #  gen.accepts <- gen.accepts + 1
    #}
    #gen.p.vals[bdx] <- cwe.sim$gen.p.value
  }
}
{
withins <- hist(within.ests[, 1], breaks = 10, plot = FALSE)
#wwww <- hist(within.ests[, 2], breaks = 10, plot = FALSE)
pooleds <- hist(pooled.ests[, 2], breaks = 10, plot = FALSE)
plot(withins, col=rgb(0,0,1,1/4), xlim = c(0,5), main = "")
#plot(wwww, col=rgb(0,1,0,1/4), xlim = c(0,5), main = "")
#plot(pooleds, col=rgb(1,0,0,1/4), xlim = c(0,5), add = TRUE)
abline(v = beta2, col = "green", lty = 1)
legend(
  x = 2.5, y = 20,
  legend = c("withins", "pooleds", "true beta2 == 2.0"),
  col = c("blue", "red", "green"),
  lty = c(1, 1))
title(paste(
  "model y ~ x2 \n",
  #"DGP: y = b1 + b2 * x2 + u, \n",
  "DGP: y = b1 + b2 * x2 + b3 * x3 + fe + u, \n",
  "where: x3 = d1 + d2 * x2 + e \n",
  "test accepts: ", accepts / B, " gen test accepts: ", gen.accepts / B))
}
