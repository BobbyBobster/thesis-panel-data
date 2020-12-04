library(plm)
library(lmtest)

data("Fatalities", package = "AER")
Fatalities$frate <- with(Fatalities, fatal / pop * 10000)
fm <- frate ~ beertax

poolmod <- plm(fm, 
  Fatalities, model = "pooling")

dmod <- plm(diff(frate, 5) ~ diff(beertax, 5), 
  Fatalities, model = "pooling")

femod <- plm(fm,
  Fatalities, model = "within")
