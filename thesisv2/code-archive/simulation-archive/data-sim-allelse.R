


# error correlated with indep vars
#mu <- c(2, 1, 0)
#Sigma <- matrix(1, nrow = 3, ncol = 3)
#Sigma <- rbind(
#  c(2, 1.2, 2),
#  c(1.2, 2, 1.4),
#  c(2, 1.4, 3))
#Sigma <- diag(3)

mu <- c(1, 0)
Sigma <- diag(c(2, 1))

B <- 100
hits <- 0
p.vals <- rep(NA, B)
for (bdx in 1:B) {
  rawvars <- MASS::mvrnorm(n = nobs * nperson, mu=mu, Sigma=Sigma)
  
  fulldata$raw1 <- rawvars[, 1]
  fulldata$rawerror <- rawvars[, 2]
  #fulldata$rawerror <- rawvars[, 3]
  
  fulldata$ycor <- beta * fulldata$raw1 + 0.5 * fulldata$fe + fulldata$rawerror
  
  fulldata.plm <- pdata.frame(fulldata, index = c("id", "t"))
  dv.sim.cor <- c("ycor")
  iv.sim.cor <- c("raw1")
  gid.sim.cor <- c("id")
  tid.sim.cor <- c("t")
  cwe.cor <- consistency.within.est(panel.set = fulldata.plm, 
    dep.var = dv.sim.cor, indep.vars = iv.sim.cor, 
    group.index = gid.sim.cor, time.index = tid.sim.cor)
  
  if (cwe.cor$wald < cwe.cor$chisq.bound) { 
    hits <- hits + 1 
  }
  p.vals[bdx] <- cwe.cor$p.value
}
#plot(cwe.cor, main = "error corr with indep")




## Practical Test for Strict Exogeneity ... by Liangjun Su, Yonghui Zhang et al.
N <- 100
Tt <- 9
alpha <- rnorm(mean = 1, sd = 0.25)
eps <- rnorm(mean = 0, sd = 1)
err <- rnorm(mean = 0, sd = 1)

x.init <- 0.5
beta <- 1
deltas <- c(0, 0.1, 0.2)
#y[i][t] <- x[i][t] + x[






#








## Wansbeek Spierdijk paper
# Sim vars
N <- 100
Tt <- 10


# Omitted variables 
# Params
beta <- 1
gamma <- 1
var.eps <- 0.25
var.theta <- 0.36
var.eta <- 0.36
corr.theta.eta <- -0.6

rho <- 0.6

y[idx, tdx] <- alpha[idx] + beta * x[idx, tdx] + gamma * z[idx, tdx] + eps[idx, tdx]





#### DOESNT WORK ####


nperson <- 200
nobs <- 5

# Main equation 
#fe.sd <- 0.6019
#x.sd <- 0.3600
#u.sd <- 0.3982
fe.sd <- x.sd <- u.sd <- 1

beta1 <- -0.5 # intercept
beta2 <- 2
beta3 <- 1

# Auxiliary equation
delta1 <- 0.5
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
  within.ests <- rep(0, B)
  pooled.ests <- matrix(0, nrow = B, ncol = 2)
  gen.walds <- rep(0, B)
  
  # First generate our data using the time constant effects
  constantdata <- data.frame(id = 1:nperson, fe = rnorm(nperson, sd = fe.sd), start.y = rnorm(nperson, sd = x.sd))
  
  # We expand our data by nobs
  fulldata <- constantdata[rep(1:nperson, each = nobs), ]
  
  # Using the generalized apply function coded above
  fulldata$t <- gapply(rep(1, length(fulldata$id)), 
    group = fulldata$id, 
    fun = cumsum)
  
  # lagged variables as AR(1) process
  ys <- rep(0, nperson * nobs)
  for (idx in 0:(nperson-1)) {
    ars <- arima.sim(model = list(ar = 0.8), n = nobs, sd=0.5)
    ys[(1 + idx * nobs):(5 + idx * nobs)] <- ars[1:5]
  }
  
  fulldata$ys <- ys
  fulldata$y.lag <- 0
  Hmisc::Lag(fulldata$ys, -1)
  
  for (bdx in 1:B) {
    
    # Finally we are ready to simulate our y variables
    fulldata$y <- 
      beta1 + 
      beta2 * fulldata$x2 + 
      #beta3 * x3 +
      #.5 * fulldata$fe + 
      u
    
    fulldata.plm <- pdata.frame(fulldata, index = c("id", "t"))
    cwe.sim <- consistency.within.est(panel.set = fulldata.plm, 
      dep.var = c("y"), indep.vars = c("x2"), 
      group.index = c("id"), time.index = "t")
    
    if (cwe.sim$wald < cwe.sim$chisq.bound) { 
      accepts <- accepts + 1 
    }
    p.vals[bdx] <- cwe.sim$p.value
    within.ests[bdx] <- cwe.sim$within.est$coefficients
    pooled.ests[bdx, ] <- plm(y ~ x2, data = fulldata.plm, model = "pooling")$coefficients
    
    gen.walds[bdx] <- cwe.sim$gen.wald
    if (gen.walds[bdx] < cwe.sim$gen.chisq.bound) {
      gen.accepts <- gen.accepts + 1
    }
    gen.p.vals[bdx] <- cwe.sim$gen.p.value
  }
}
{
  withins <- hist(within.ests, breaks = 10, plot = FALSE)
  pooleds <- hist(pooled.ests[, 2], breaks = 10, plot = FALSE)
  plot(withins, col=rgb(0,0,1,1/4), xlim = c(1, 4), main = "")
  plot(pooleds, col=rgb(1,0,0,1/4), xlim = c(1, 4), add = TRUE)
  abline(v = beta2, col = "green", lty = 1)
  legend(
    x = 3, y = 20,
    legend = c("withins", "pooleds", "true beta2 == 2.0"),
    col = c("blue", "red", "green"),
    lty = c(1, 1))
  title(paste(
    "model y ~ x2 \n",
    "DGP: y = b1 + b2 * x2 + u, \n",
    #"DGP: y = b1 + b2 * x2 + b3 * x3 + fe + u, \n",
    #"where: x3 = d1 + d2 * x2 + e \n",
    "test accepts: ", accepts / B, " gen test accepts: ", gen.accepts / B))
}