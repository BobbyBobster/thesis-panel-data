source("consistencytest.R")

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


mu <- c(1, 0)
Sigma <- diag(c(2, 1))
C <- 20
B <- 100
hits.total <- rep(0, C)
p.vals.total <- matrix(NA, nrow = B, ncol = C)
for (cdx in 1:C) {
hits <- 0
p.vals <- rep(NA, B)
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

# WHY FIXED EFFECTS ALSO IN X??
fulldata$x <- fulldata$fe + rnorm(nobs * nperson, sd = x.sd)
fulldata$x2 <- rnorm(nobs * nperson, mean = 2, sd = 2)
fulldata$u <- rnorm(nobs * nperson, mean = 0, sd = u.sd)

# Finally we are ready to simulate our y variables
fulldata$y <- beta * fulldata$x +  .5 * fulldata$fe + fulldata$u


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
p.vals.total[, cdx] <- p.vals
hits.total[cdx] <- hits
}