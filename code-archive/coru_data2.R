library(MASS)
library(nlme)
library(plm)
library(extRC)
library(tidyverse)

data("Wages", package = "plm")

Wages$bluecol <- ifelse(Wages$bluecol == "no", 0, 1)
Wages$south <- ifelse(Wages$south == "no", 0, 1)
Wages$smsa <- ifelse(Wages$smsa == "no", 0, 1)
Wages$married <- ifelse(Wages$married == "no", 0, 1)
Wages$sex <- ifelse(Wages$sex == "male", 0, 1)
Wages$union <- ifelse(Wages$union == "no", 0, 1)
Wages$black <- ifelse(Wages$black == "no", 0, 1)

wages <- pdata.frame(Wages, index = 595)

form <- lwage ~ 
  exp + I(exp**2) + wks + bluecol + ind + south +
  smsa + married + union + ed + sex + black

fit.fe <- plm(
  form,
  data = wages,
  model = "within")
summary(fit.fe)

fit.gls <- update(fit.fe, 
  model = "random")
summary(fit.gls)

phtest(fit.fe, fit.gls)

fit.fd <- update(fit.fe, 
  model = "fd")
summary(fit.fd)


# CoRu 
y <- wages %>% filter(id == 1) %>% select(lwage) %>% data.matrix()
x <- wages %>% filter(id == 1) %>% select(exp) %>% data.matrix()
V <- wages %>% filter(id == 1) %>%  select(-lwage, -exp, -wks, -id, -time) %>%
  select(ind) %>% 
  data.matrix()

proj.V <- V %*% solve(t(V) %*% V) %*% t(V)
Q.V <- diag(dim(V)[1]) - proj.V
b <- solve(t(x) %*% Q.V %*% x) %*% t(x) %*% Q.V %*% y



## Table 8.1 from Greene, Econometric Analysis 7th ed
greene.ols <- plm(
  wks ~ lwage + ed + union + sex,
  data = wages,
  model = "pooling")
summary(greene.ols)
greene.1 <- plm(
  wks ~ lwage + ed + union + sex |
    ind + ed + union + sex,
  data = wages,
  model = "pooling")
summary(greene.1)
greene.2 <- plm(
  wks ~ lwage + ed + union + sex |
    ind + ed + union + sex + smsa,
  data = wages,
  model = "pooling")
summary(greene.2)

lms <- NA
for (idx in 1:595) {
  lms[idx] <- lm(lwage ~ exp, data = wages, subset = id == idx)$coefficients
}
mean(lms)


fixed.dum <- lm(lwage ~ exp + factor(id) - 1, data = wages)
summary(fixed.dum)
yhat <- fixed.dum$fitted
sampling.ids <- sample(1:595, size = 10)
wages.sample <- wages %>% filter(id %in% sampling.ids)
car::scatterplot(yhat ~ wages$exp | wages$id, 
  boxplots = FALSE, 
  xlab = "exp",
  ylab = "yhat",
  smooth = FALSE)

# Wald tests ====
pwaldtest(fit.fe,
  test = "Chisq")

wald.b <- fit.fe$coefficients
wald.var <- fit.fe$vcov
wald.test.stat <- t(wald.b) %*% solve(wald.var) %*% wald.b


# Old Wald test ====
restr.mat <- diffMatrix(6, 1)
X <- wages %>% 
  filter(time == 1) %>% 
  dplyr::select(-id, -time, -lwage) %>% 
  data.matrix()
W <-t(wald.b) %*% restr.mat %*% solve(wald.var) %*% t(restr.mat) %*% wald.b


### OLD ESTIMATES ====
# Timespans estimator beta hat
estimates <- matrix(nrow = k, ncol = Tt - 1)
rownames(estimates) <- c("exp", "expsq", "wks", "bluecol", "ind", "south", "smsa", "married", "union")
for (tdx in 1:(Tt - 1)) {
  D <- diffMatrix(Tt, tdx)
  Delta <- D %*% t(D)
  
  first.term <- 0
  second.term <- 0
  for (idx in 1:n) {
    x <- wages %>% filter(id == idx) %>% select(!!!dep.vars) %>%  data.matrix()
    y <- wages %>% filter(id == idx) %>% select(!!indep.var) %>% data.matrix()
    first.term <- first.term + t(x) %*% Delta %*% x
    second.term <- second.term + t(x) %*% Delta %*% y
  }
  estimates[, tdx] <- solve(first.term) %*% second.term
}
estimates

# u.i's
#beta.hat <- as.vector(t(estimates))
u.i <- matrix(nrow = k, ncol = (Tt - 1))
for (tdx in 1:(Tt - 1)) {
  D <- diffMatrix(Tt, tdx)
  Delta <- D %*% t(D)
  
  first.term.u <- matrix(0, nrow = k, ncol = k)
  for (idx in 1:n) {
    x <- wages %>% filter(id == idx) %>% select(!!!dep.vars) %>% data.matrix()
    y <- wages %>% filter(id == idx) %>% select(!!indep.var) %>% data.matrix()
    
    first.term.u <- first.term.u + t(x) %*% Delta %*% x
  }
  
  second.term.u <- matrix(0, nrow = k, ncol = 1)
  for (idx in 1:n) {
    x <- wages %>% filter(id == idx) %>% select(!!!dep.vars) %>% data.matrix()
    y <- wages %>% filter(id == idx) %>% select(!!indep.var) %>% data.matrix()
    
    epsilon <- y - x %*% estimates[, tdx]
    
    second.term.u <- second.term + t(x) %*% Delta %*% epsilon
  }
  
  u.i[, tdx] <- first.term.u %*% second.term.u
}

# vcov matrix
B <- rbind(diag(Tt - 2), 0) - rbind(0, diag(Tt - 2))
R <- kronecker(B, diag(k))
outer.sum.vcov <- matrix(0, nrow = Tt, ncol = k)
for (idx in 1:n) {
  x <- wages %>% filter(id == idx) %>% select(!!!dep.vars) %>% data.matrix()
  y <- wages %>% filter(id == idx) %>% select(!!indep.var) %>% data.matrix()
  #u <- y - x %*% estimates[1] ################ Wat moet hier eigenlijk komen ipv estimates[1]
}

# Sigma <- solve(outer.vcov) %*% inner.vcov %*% solve(outer.vcov) 





u.i.stacked <- as.vector(t(u.i))
middle.term <- t(R) %*% u.i.stacked %*% t(u.i.stacked) %*% R

# test
est.stacked <- as.vector(t(estimates))
q.w <- t(est.stacked) %*% R %*% MASS::ginv(middle.term) %*% t(R) %*% est.stacked
q.w

plot(x = NA, xlim = c(1, 6), ylim = c(-0.06, 0.15))
for (row in 1:9) {
  lines(estimates[row, ], col = row)
  abline(h = mean(estimates[row, ]), col = row, lty = 2)
}

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(3, 3))
for (row in 1:9) {
  plot(estimates[row, ], 
    type = "l", 
    xlab = "timespan",
    ylab = rownames(estimates)[row])
  abline(h = mean(estimates[row, ]), lty = 3)
  abline(h = 0, lty = 2)
}
par(oldpar)

