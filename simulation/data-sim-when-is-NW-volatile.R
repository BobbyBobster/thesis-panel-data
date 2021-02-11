library(doParallel)

sim.type <- c("When is NW volatile?") # change function
file.name <- c("when-is-NW-volatile-ouput.txt")
B <- 50

print.and.cat <- function (output.str, file = file.name) {
  print(output.str)
  cat(c(output.str, "\n"), file = file.name, append = TRUE)
}

ns <- c(100, 50, 10)
rho.lag <- 0.6
Ts <- 5
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
  print.and.cat(paste("rho.lag", rho.lag))
  print.and.cat(paste("Tt", Tt))
      all.montes[monte.counter, 1:7] <- 
        montecarlo(data.gen.ovb.size, nobs = nobs, rho.lag = rho.lag, nperson = Tt, num.sims = B)
      all.montes[monte.counter, 8:10] <- c(nobs, rho.lag, B.size)
      print.and.cat(all.montes[monte.counter, ])
      monte.counter <- monte.counter + 1
      print.and.cat(monte.counter)
    }

stopCluster(cl)
rm(cl)
