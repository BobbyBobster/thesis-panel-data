library(tidyverse)
library(nlme)
library(plm)

form <- lwage ~ 
  wks + south + smsa + married + exp + I(exp ^ 2) + 
  bluecol + ind + union + sex + black + ed

fit.gls <- plm(form,
  data = wages,
  model = "random")
summary(fit.gls)

fit.fe <- plm(form,
  data = wages,
  model = "within")
summary(fit.fe)


# hausman test
phtest(fit.fe, fit.gls)

# HT - Hausman Taylor (1981)
ht <- plm(
  lwage ~ wks + south + smsa + married + exp + I(exp ^ 2) + 
    bluecol + ind + union + sex + black + ed |
    bluecol + south + smsa + ind + sex + black |
    wks + married + union + exp + I(exp ^ 2), 
  data = wages,
  random.method = "ht", model = "random", inst.method = "baltagi")
summary(ht)

# AM - Amemiya and MaCurdy (1986)
fit.am <- plm(
  lwage ~ wks + south + smsa + married + exp + I(exp ^ 2) + 
    bluecol + ind + union + sex + black + ed |
    bluecol + south + smsa + ind + sex + black |
    wks + married + union + exp + I(exp ^ 2), 
  data = wages,
  random.method = "ht",
  model = "random",
  inst.method = "am")
summary(fit.am)
