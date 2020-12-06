library(doParallel)

sim.type <- c("OVB size MC") # change function
file.name <- c("me-size-output-mc.txt")
B <- 500

print.and.cat <- function (output.str, file = file.name) {
  print(output.str)
  cat(c(output.str, "\n"), file = file.name, append = TRUE)
}

ns <- c(100, 500, 1000)
rhos <- c(0.6, 0.9)
Ts <- c(5, 10)
all.montes <- matrix(NA, nrow = 3 * 2 * 2, ncol = 4 + 2 + 1 + 3)
colnames(all.montes) <- 
  c("wald.accepts", "p.vals", "gen.wald.accepts", "gen.p.vals", 
    "within.ests", "pooled.ests", 
    "agreements",
    "Tt", "rho.lag", "samplesize")

cl <- makeCluster(3, type="FORK")
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = NULL)

print.and.cat("Running Monte Carlo simulation")
print.and.cat(sim.type)
print.and.cat(paste("B:", B))

monte.counter <- 1
for (nobs in Ts) {
  print.and.cat(paste("nobs", nobs))
  for (rho.lag in rhos) {
    print.and.cat(paste("rho.lag", rho.lag))
    for (B.size in ns) {
      print.and.cat(paste("B.size", B.size))
      all.montes[monte.counter, 1:7] <- 
        montecarlo(data.gen.ovb.size, nobs = nobs, rho.lag = rho.lag, nperson = B.size, num.sims = B)
      all.montes[monte.counter, 8:10] <- c(nobs, rho.lag, B.size)
      print.and.cat(all.montes[monte.counter, ])
      monte.counter <- monte.counter + 1
      print.and.cat(monte.counter)
    }
  }
}

stopCluster(cl)
rm(cl)
