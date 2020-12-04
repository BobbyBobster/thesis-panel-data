tab <- rbind(
c(0.780,             0.838,  5,     0.6,        100),
c(0.200,             0.330,  5,     0.6,        500),
c(0.010,             0.038,  5,     0.6,       1000),
c(0.478,             0.540,  5,     0.9,        100),
c(0.000,             0.000,  5,     0.9,        500),
c(0.000,             0.000,  5,     0.9,       1000),
c(0.636,             0.906,  10,     0.6,        100),
c(0.014,             0.444,  10,     0.6,        500),
c(0.000,             0.030,  10,     0.6,       1000),
c(0.050,             0.128,  10,     0.9,        100),
c(0.000,             0.000,  10,     0.9,        500),
c(0.000,             0.000,  10,     0.9,       1000)
)

truetab <- cbind(1 - tab[, 1], 1 - tab[, 2], tab[, 3:5])
colnames(truetab) <- 
  c("wald.power", "gen.wald.power", 
    "Tt", "rho.lag", "samplesize")
xtable::print.xtable(xtable::xtable(truetab, type = "latex"), filename = "powertab.tex")

f.d <- data.creation(5, 0.6, 1000)
bootys.fiddle <- bootstrapper(fulldata = f.d, B.size = 500, B = 100, nperson = 1000, replace = FALSE)


tab.sim <- `data-sim-script_results`$all.bootys


  
