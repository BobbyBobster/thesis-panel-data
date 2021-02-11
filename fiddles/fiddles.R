N <- 100
pRsq <- rep(NA, N)
sRsq <- rep(NA, N)
prss <- rep(NA, N)
srss <- rep(NA, N)
ptss <- rep(NA, N)
stss <- rep(NA, N)
for (idx in 1:N) {
  fdp <- pdata.frame(data.gen.ovb.power(nobs = 5, rho.lag = 0.6, nperson = 100), index = c("id", "t"))
  fds <- pdata.frame(data.gen.ovb.size(nobs = 5, rho.lag = 0.6, nperson = 100), index = c("id", "t"))
  
  cwep <- consistency.within.est(fdp, "y", "x2", "id", "t")
  cwes <- consistency.within.est(fds, "y", "x2", "id", "t")
  
  pRsq[idx] <- cwep$R.sq 
  sRsq[idx] <- cwes$R.sq 
  prss[idx] <- cwep$rss 
  srss[idx] <- cwes$rss 
  ptss[idx] <- cwep$tss 
  stss[idx] <- cwes$tss 
}

mean(pRsq) 
mean(sRsq) 
psumRsq <- 1 - sum(prss)/sum(ptss)
ssumRsq <- 1 - sum(srss)/sum(stss)
mean(psumRsq)
mean(ssumRsq)


ht <- plm(lwage ~ wks + south + smsa + married + exp + I(exp ^ 2) + 
    bluecol + ind + union + sex + black + ed |
    bluecol + south + smsa + ind + sex + black |
    wks + married + union + exp + I(exp ^ 2), 
  data = Wages, index = 595,
  model = "within", inst.method = "baltagi")
summary(ht)



eigen(t(cwe.wages$R) %*% cwe.wages$Sigma.hat %*% cwe.wages$R)$values


pbltest(y ~ x2, data = fds)

pbltest(lwage ~ wks + south + smsa + married + exp + I(exp ^ 2) + 
    bluecol + ind + union + sex + black + ed,
  data = Wages, index = 595)

N <- 10
ps <- rep(NA, N)
gp <- rep(NA, N)
for (idx in 1:N) {
  fulldata <- data.gen.me.power(5, 0.9, 100)
  cwe <- consistency.within.est(pdata.frame(fulldata, index = c("id","t")), "y","x2","id","t")
  ps[idx] <- cwe$p.value
  gp[idx] <- cwe$gen.p.value
}








{
object <- cwe.wages
ylabs <- 
  c("experience", "experience^2", "weeks worked", 
    "bluecollar dummy", "industry dummy", "southern dummy",
    "std metro area dummy", "married dummy", "union dummy")
  
  # Create nice rectangle of plots 
  plt.row <- round(sqrt(object$k))
  plt.col <- ceiling(sqrt(object$k))
  par(
    mfrow = c(plt.row, plt.col),
    mar = c(5.1, 4.1, 0.1, 2.1),
    oma = c(0, 0, 4, 0)
  )
          
  for (row in 1:object$k) {
    plot(object$beta.hat.unstacked[row, ], 
      type = "b", 
      xlab = "j",
      #ylab = rownames(object$beta.hat.unstacked)[row], 
      ylab = ylabs[row]
      )
    abline(h = object$within.est$coefficients[row], lty = 3)
    abline(h = 0, lty = 2)
    arrows(
      x0 = 1:(object$Tt - 1), 
      y0 = object$beta.hat.unstacked[row, ] - 1.96 * object$beta.std.err.unstacked[row, ],
      x1 = 1:(object$Tt - 1),
      y1 = object$beta.hat.unstacked[row, ] + 1.96 * object$beta.std.err.unstacked[row, ],
      code = 3,
      angle = 90,
      length = 0.05)
  }
}













sws <- 0
gens <- 0
for (idx in 1:20) {
  cwe <- consistency.within.est(
    pdata.frame(data.gen.si.power(5, 0.6, 100), index = c("id", "t")), "y", "x2", "id", "t", try.normal.inv = FALSE)
  if (cwe$p.value <- 0.05) { sws <- sws + 1 }
  if (cwe$gen.p.value <- 0.05) { gens <- gens + 1 }
}
sws
gens

me.power.b50 <- rbind(
c(       0.936,            0.950,        5,     0.6,        100),    
c(       0.956,            0.952,        5,     0.6,        500),
c(       0.944,            0.930,        5,     0.6,       1000),
c(       0.926,            0.944,        5,     0.9,        100),
c(       0.962,            0.968,        5,     0.9,        500),
c(       0.924,            0.936,        5,     0.9,       1000),
c(       0.886,            0.954,       10,     0.6,        100),
c(       0.940,            0.942,       10,     0.6,        500),
c(       0.934,            0.932,       10,     0.6,       1000),
c(       0.904,            0.942,       10,     0.9,        100),
c(       0.946,            0.966,       10,     0.9,        500),
c(       0.934,            0.956,       10,     0.9,       1000)
       )
cbind(1 - me.power.b50[, c(1,2)], me.power.b50[, 3:5])


par(mfrow = c(1, 2))
fd.si.power <- pdata.frame(data.gen.si.power(5, 0.6, 100), index = c("id", "t"))
cwe.si <- consistency.within.est(fd.si.power, "y", "x2", "id", "t", try.normal.inv = FALSE)

N <- 100
pvls <- rep(NA, N)
gpvls <- rep(NA, N)
for (idx in 1:N) {
  fd.me.power <- pdata.frame(data.gen.me.power(5, 0.9, 500), index = c("id", "t"))
  cwe.me <- consistency.within.est(fd.me.power, "y", "x2", "id", "t", try.normal.inv = FALSE)
  pvls[idx] <- cwe.me$p.value
  gpvls[idx] <- cwe.me$gen.p.value
}
length(pvls[pvls < 0.05])
length(gpvls[gpvls < 0.05])


sigma.xs <- 1.5
sigma.eps <- 0.5
xs  <- rnorm(5000, mean = 2, sd = sigma.xs)
eps <- rnorm(5000, mean = 0, sd = sigma.eps)
beta <- 2
ys <- beta * xs + eps
linmod <- lm(ys ~ xs)
plot(ys ~ xs); abline(linmod)
summary(linmod)
1 - var(eps) / var(ys) # sample R^2
1 - (sigma.eps)**2 / ((beta * sigma.xs)**2 + (sigma.eps)**2) # simulation R^2





