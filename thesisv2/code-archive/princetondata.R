library(tidyverse)
library(plm)

Panel101Princeton <- haven::read_dta("Panel101Princeton.dta")
Panel101Princeton <- transform(Panel101Princeton, country = as.numeric(country))

prince <- pdata.frame(Panel101Princeton,
  index = c("country", "year"))

form.pr <- y ~ x1
  
plot(y ~ x1, data = prince, col = country)
pr.ols <- plm(form.pr,
  data = prince,
  model = "pooling")
abline(pr.ols)

pr.re <- update(pr.ols, model = "random")
pr.fe <- update(pr.ols, model = "within")
pr.effects <- fixef(pr.fe, type = "level")

for (idx in 1:7) {
  current <- prince %>% filter(country == idx)
  clip(
    x1 = min(current$x1),
    x2 = max(current$x1),
    y1 = min(current$y),
    y2 = max(current$y))
  abline(a = pr.effects[idx], b = pr.fe$coefficients["x1"], col = idx)
}



## Setup ====
panel.set <- prince

Tt <- 10
n <- 7

indep.var <- quo(y)
dep.vars <- quos(x1, x2, x3)
group.index <- quo(country)
time.index <- quo(year)

k <- length(dep.vars)
#B <- rbind(diag(Tt - 1), 0) - rbind(0, diag(Tt - 1))
B <- diffMatrix(Tt, 1)
R <- kronecker(B, diag(k))

