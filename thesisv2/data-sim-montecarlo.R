montecarlo <- function (data.gen.fun, nobs, rho.lag, nperson, num.sims) {
  #cl <- makeCluster(3, type="FORK")
  #registerDoParallel(cl)
  #clusterSetRNGStream(cl, iseed = NULL)
  
  booty <- foreach (1:num.sims) %dopar% {
    accepts <- gen.accepts <- agreements <- 0
    
    fulldata <- data.gen.fun(nobs, rho.lag, nperson) 
    fulldata <- pdata.frame(fulldata, index = c("id", "t"))
    cwe.sim <- consistency.within.est(
      panel.set = fulldata, 
      dep.var = c("y"), indep.vars = c("x2"), 
      group.index = c("id"), time.index = c("t"),
      alpha = 0.05, no.gen = FALSE)
    
    if (cwe.sim$wald < cwe.sim$chisq.bound) { 
      accepts <- 1
    }
    if (cwe.sim$gen.wald < cwe.sim$gen.chisq.bound) {
      gen.accepts <- 1
    }
    if (accepts == gen.accepts) {
      agreements <- 1
    }
    
    c(accepts, cwe.sim$p.value, 
      gen.accepts, cwe.sim$gen.p.value,
      cwe.sim$within.est$coefficients,
      plm(y ~ x2, data = fulldata, model = "pooling")$coefficients[2],
      agreements) 
  }
  booty <- matrix(unlist(booty), ncol = lengths(booty)[1], byrow = TRUE)
  colnames(booty) <- 
    c("wald.accepts", "p.vals", "gen.wald.accepts", "gen.p.vals", 
      "within.ests", "pooled.ests",
      "agreements")
  
  c(sum(booty[, 1]), mean(booty[, 2]),
    sum(booty[, 3]), mean(booty[, 4]), 
    mean(booty[, 5]), mean(booty[, 6]),
    sum(booty[, 7]))
  
  #stopCluster(cl)
  #rm(cl)
}


#{
#withins <- hist(booty[, 3], breaks = 10, plot = FALSE)
#pooleds <- hist(booty[, 4], breaks = 10, plot = FALSE)
#plot(withins, col=rgb(0,0,1,1/4), xlim = c(0,5), main = "", xlab = NA)
#plot(pooleds, col=rgb(1,0,0,1/4), xlim = c(0,5), add = TRUE)
#abline(v = beta2, col = "green", lty = 1)
##legend(
##  x = 2.5, y = 20,
##  legend = c("withins", "pooleds", "true beta2 == 2.0"),
##  col = c("blue", "red", "green"),
##  lty = c(1, 1))
#title(
#  main = paste(
#    "model y ~ x2 \n",
#    "test accepts:", total.accepts / B, " gen test accepts:", total.gen.accepts / B, "\n"
#    #"DGP: y = b1 + b2 * x2 + u, \n",
#    #"DGP: y = b1 + b2 * x2 + b3 * x3 + fe + u, \n",
#    #"where: x3 = d1 + d2 * x2 + e \n",
#  ),
#  sub = paste(
#    "bootstrap samples:", B, " sample size:", B.size, "\n",
#    "time periods:", nobs))
#}
