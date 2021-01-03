montecarlo <- function (data.gen.fun, nobs, rho.lag, nperson, num.sims, make.cluster = FALSE) {
  if (make.cluster) {
    cl <- makeCluster(3, type="FORK")
    registerDoParallel(cl)
    clusterSetRNGStream(cl, iseed = NULL)
  }
  
  
  booty <- foreach (1:num.sims) %dopar% {
    accepts <- gen.accepts <- agreements <- 0
    
    fulldata <- data.gen.fun(nobs, rho.lag, nperson) 
    fulldata <- pdata.frame(fulldata, index = c("id", "t"))
    cwe.sim <- consistency.within.est(
      panel.set = fulldata, 
      dep.var = c("y"), indep.vars = c("x2"), 
      group.index = c("id"), time.index = c("t"),
      alpha = 0.05, no.gen = FALSE)
    
    if (cwe.sim$p.value >= 0.05) { 
      accepts <- 1
    }
    if ((cwe.sim$p.value >= 0.05 && cwe.sim$wald < cwe.sim$chisq.bound) || 
        (cwe.sim$p.value < 0.05 && cwe.sim$wald >= cwe.sim$chisq.bound)) {
      print("wtf is wrong with the p-value?")
    }
      
      
    if (cwe.sim$gen.p.value >= 0.05) {
      gen.accepts <- 1
    }
    if (accepts == gen.accepts) {
      agreements <- 1
    }
    
    c(accepts, cwe.sim$p.value, 
      gen.accepts, cwe.sim$gen.p.value,
      cwe.sim$within.est$coefficients,
      plm(y ~ x2, data = fulldata, model = "pooling")$coefficients[2],
      agreements, cwe.sim$R.sq, cwe.sim$rss, cwe.sim$tss) 
  }
  booty <- matrix(unlist(booty), ncol = lengths(booty)[1], byrow = TRUE)
  colnames(booty) <- 
    c("wald.accepts", "p.vals", "gen.wald.accepts", "gen.p.vals", 
      "within.ests", "pooled.ests",
      "agreements", "R.sq", "rss", "tss")
  
  if (make.cluster) {
    stopCluster(cl)
    rm(cl)
  }
  
  rsq <- 1 - sum(booty[, 9])/sum(booty[, 10])
  c(sum(booty[, 1]), mean(booty[, 2]),
    sum(booty[, 3]), mean(booty[, 4]), 
    mean(booty[, 5]), mean(booty[, 6]),
    sum(booty[, 7]), rsq)
}
