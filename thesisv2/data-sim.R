

## econometricsbysimulation.com
nperson <- 20
nobs <- 10

fe.sd <- 0.01 
x.sd <- 1.5
u.sd <- 0.01
beta <- 2

# Add a time index, first define a group apply function
# that applies by group index.
gapply <- function(x, group, fun) {
  returner <- numeric(length(group))
  for (idx in unique(group)) 
    returner[idx == group] <- get("fun")(x[idx == group])
  returner
}


B <- 200
C <- 20
p.vals.total <- matrix(NA, nrow = B, ncol = C)
for (cdx in 1:C) {
hits <- 0
p.vals <- rep(NA, B)
for (bdx in 1:B) {
  # First generate our data using the time constant effects
  constantdata <- data.frame(id = 1:nperson, fe = rnorm(nperson, sd = fe.sd))
  
  # We expand our data by nobs
  fulldata <- constantdata[rep(1:nperson, each = nobs), ]
  
  # Using the generalized apply function coded above
  fulldata$t <- gapply(rep(1, length(fulldata$id)), 
    group = fulldata$id, 
    fun = cumsum)

  # Now we are ready to calculate the time variant xs
  # And our unobservable error
  fulldata$x <- fulldata$fe + rnorm(nobs * nperson, sd = x.sd)
  fulldata$x2 <- rnorm(nobs * nperson, mean = 2, sd = 2)
  fulldata$u <- rnorm(nobs * nperson, mean = 0, sd = u.sd)
  
  # Finally we are ready to simulate our y variables
  fulldata$y <- beta * fulldata$x +  .5 * fulldata$fe + fulldata$u

  fulldata.plm <- pdata.frame(fulldata, index = c("id", "t"))
  cwe.sim <- consistency.within.est(panel.set = fulldata.plm, 
    dep.var = c("y"), indep.vars = c("x"), 
    group.index = c("id"), time.index = "t")
  
  if (cwe.sim$wald < cwe.sim$chisq.bound) { 
    hits <- hits + 1 
  }
  p.vals[bdx] <- cwe.sim$p.value
}
p.vals.total[, cdx] <- p.vals
}



cwe.sim$wald
cwe.sim$chisq.bound
cwe.sim$gen.wald

model.plm <- plm(y ~ x, data = fulldata.plm, model = "within")
summary(model.plm)
model.lm <- lm(y ~ x + factor(id) - 1, data = fulldata.plm)
summary(model.lm)



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
