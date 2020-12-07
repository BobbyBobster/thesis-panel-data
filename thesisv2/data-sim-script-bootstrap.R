library(doParallel)



sim.type <- c("ME size") # change function
file.name <- c("me-size-output.txt")
nperson <- 5000
B <- 500
replace <- FALSE

print.and.cat <- function (output.str, file = file.name) {
  print(output.str)
  cat(c(output.str, "\n"), file = file.name, append = TRUE)
}

ns <- c(100, 500, 1000)
rhos <- c(0.6, 0.9)
Ts <- c(5, 10)
all.bootys <- matrix(NA, nrow = 3 * 2 * 2, ncol = 4 + 2 + 1 + 3)
colnames(all.bootys) <- 
  c("wald.accepts", "p.vals", "gen.wald.accepts", "gen.p.vals", 
    "within.ests", "pooled.ests", 
    "agreements",
    "Tt", "rho.lag", "samplesize")

cl <- makeCluster(3, type="FORK")
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = NULL)

print.and.cat("Running simulation")
print.and.cat(sim.type)
print.and.cat(paste("nperson:", nperson))
print.and.cat(paste("B:", B))
print.and.cat(paste("replace:", replace))

booty.counter <- 1
for (nobs in Ts) {
  print.and.cat(paste("nobs", nobs))
  for (rho.lag in rhos) {
    print.and.cat(paste("rho.lag", rho.lag))
    fulldata <- data.gen.me.size(nobs, rho.lag, nperson = nperson) # power vs. size
    for (B.size in ns) {
      print.and.cat(paste("B.size", B.size))
      all.bootys[booty.counter, 1:7] <- 
        bootstrapper(fulldata, B.size, B = B, nperson = nperson, replace = replace)
      all.bootys[booty.counter, 8:10] <- c(nobs, rho.lag, B.size)
      print.and.cat(all.bootys[booty.counter, ])
      booty.counter <- booty.counter + 1
      print.and.cat(booty.counter)
    }
  }
}

stopCluster(cl)
rm(cl)



#plot(
#  x = ns,
#  y = sim.results[1:3, 4], 
#  type ="l", col = 1)
#for (idx in c(4, 7, 10)) {
#  lines(
#    x = ns,
#    sim.results[idx:(idx+2), 4], col = idx)
#}
#title("wald")
#legend(
#  x = 800,
#  y = 0.25,
#  legend = c("5, 0.6", "5, 0.9", "10, 0.6", "10, 0.9"),
#  col = c(1, 4, 7, 10), 
#  lty = c(1, 1, 1, 1))





