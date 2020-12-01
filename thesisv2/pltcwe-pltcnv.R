plt.cwe <- function (cwe, step.size = 1, step.end = 100) {
  df <- cwe$k * (cwe$Tt - 2)
  
  xs.chisq <- seq(0, step.end - step.size, by = step.size)
  pts.chisq <- dchisq(xs.chisq, df = df)
  
  plot(
    x = xs.chisq,
    y = pts.chisq,
    type = "l", 
    lty = 1)
  abline(
    v = cwe$chisq.bound, 
    col = "red", 
    lty = 2)
  abline(
    v = cwe$wald, 
    col = "green",
    lty = 1)
  title(
    paste("cwe wald ", cwe$wald, cwe$within.est$formula[2]))
  #legend(
  #  x = 75,
  #  y = 0.04,
  #  legend = c("chisq", "Stat. value", "Crit. value single 95%", "Crit. value double 95%"),
  #  col = c("black", "red", "blue", "green"), 
  #  lty = c(1, 1, 2, 2))
}

plt.cnv <- function (cnv, step.size = 0.5, step.end = 40) {
  Sigma.hat <- cnv$Sigma.hat
  R <- cnv$R
  Tt <- cnv$Tt
    
  Sigma.R <- t(R) %*% Sigma.hat %*% R
  evs <- eigen(Sigma.R)$values
  
  # package CompQuadForm for the "true" distribution
  pts.imhof <- rep(NA, step.end / step.size)
  xs.imhof <- seq(0, step.end - step.size, by = step.size)
  for (point in xs.imhof) {
    pts.imhof[point * (1 / step.size)] <- 
      (1 - CompQuadForm::imhof(q = point, lambda = evs)$Qq)
  }
  
  if (step.size >= 1) {
    crit.val <- length(pts.imhof) - length(pts.imhof[pts.imhof > 1 - cnv$alpha])
  } else if (step.size < 1 & step.size > 0) {
    crit.val <- (length(pts.imhof) - length(pts.imhof[pts.imhof > 1 - cnv$alpha])) * step.size
  }
  
  plot(
    x = xs.imhof,
    y = pracma::gradient(pts.imhof, h1 = step.size), 
    col = "black",
    type = "l",
    lty = 1)
  abline(
    v = crit.val,
    col = "red",
    lty = 2)
  abline(
    v = cnv$gen.wald, 
    col = "green",
    lty = 1)
  title(
    paste("cnv wald ", cnv$wald, cnv$within.est$formula[2]))
  #legend(
  #  x = 150,
  #  y = 0.025,
  #  legend = c("Gen. Chisq dist", "Crit. value", "Stat. value"),
  #  col = c("blue", "red", "green"), 
  #  lty = c(1, 2, 1))
}
