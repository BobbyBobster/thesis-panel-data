
#### Wages ####
{
  data("Wages", package = "plm")
  Wages <- Wages %>% mutate(expsq = exp**2)
  
  wages <- pdata.frame(Wages, index = 595)
  
  dv.wages <- c("lwage")
  iv.wages <- c("exp", "expsq", "wks", "bluecol", "ind", "south", "smsa", "married", "union")
  gid.wages <- c("id")
  tid.wages <- c("time")
  
  cwe.wages <- 
    consistency.within.est(wages, dv.wages, iv.wages, gid.wages, tid.wages)
  #plot(cwe.wages)
  cnv.wages <- 
    consistency.within.est(wages, dv.wages, iv.wages, gid.wages, tid.wages, 
      no.ginv = TRUE)
  
  # within estimator
  wages.form <- as.formula(paste(dv.wages, paste(iv.wages, collapse=" + "), sep=" ~ "))
  wages.est <- plm(
    wages.form,
    data = wages,
    model = "within")
}



#### Practical Test for Strict Exogeneity ... by Liangjun Su, Yonghui Zhang et al. ####
{
  output.full <- 
    readr::read_csv(file = "../datasets/faostat/Macro-Statistics_Key_Indicators_E_All_Data_(Normalized).csv") %>% 
    filter(Item == "Value Added (Agriculture, Forestry and Fishing)") %>% 
    filter(Element == "Value US$") %>% 
    select(AreaCode = "Area Code", Area, Year, Output = Value) 
  area.full <-
    readr::read_csv(file = "../datasets/faostat/Inputs_LandUse_E_All_Data_(Normalized).csv") %>% 
    filter(Item == "Agricultural land") %>% 
    select(AreaCode = "Area Code", Area, Year, AgrArea = Value) 
  pop.full <- 
    readr::read_csv(file = "../datasets/faostat/Employment_Indicators_E_All_Data_(Normalized).csv") %>% 
    filter(Indicator == "Employment in agriculture") %>% 
    select(AreaCode = "Area Code", Area, Year, Pop = Value) 
  govt.full <- 
    readr::read_csv(file = "../datasets/faostat/Investment_GovernmentExpenditure_E_All_Data_(Normalized).csv") %>% 
    filter(Item == "Agriculture, forestry, fishing (Central Government)") %>% 
    filter(Element == "Value US$") %>% 
    select(AreaCode = "Area Code", Area, Year, Govt = Value) 
  
  from.Year <- 2002
  till.Year <- 2013
  output <- output.full %>% filter(Year >= from.Year & Year <= till.Year)
  area <- area.full %>% filter(Year >= from.Year & Year <= till.Year)
  pop <- pop.full %>% filter(Year >= from.Year & Year <= till.Year)
  govt <- govt.full %>% filter(Year >= from.Year & Year <= till.Year)
  faostat <- output %>% left_join(area) %>% left_join(pop) %>% left_join(govt) %>% 
    group_by(Area) %>% filter(all(!is.na(AgrArea)) & all(!is.na(Pop)) & all(!is.na(Govt))) %>% 
    ungroup()
  
  faostat <- faostat %>% mutate(Id = as.numeric(factor(AreaCode)))
  
  faostat <- faostat %>% 
    mutate(loutput = log(Output), lagrarea = log(AgrArea), lpop = log(Pop), lgovt = log(Govt))
    
  faostat <- pdata.frame(faostat, 
    index = c("Id", "Year"), 
    drop.index = FALSE)
  dv.faostat <- c("loutput")
  iv.faostat <- c("lagrarea", "lpop", "lgovt")
  #iv.faostat <- c( "lgovt")
  gid.faostat <- c("Id")
  tid.faostat <- c("Year")
  
  cwe.faostat <-
    consistency.within.est(faostat, dv.faostat, iv.faostat, gid.faostat, tid.faostat)
  cnv.faostat <-
    consistency.within.est(faostat, dv.faostat, iv.faostat, gid.faostat, tid.faostat, no.ginv = TRUE)
  #plot(cwe.faostat)
  
  
  faostat.form <- loutput ~ lagrarea + lpop + lgovt
  faostat.fe <- plm(faostat.form, data = faostat, model = "within")
  faostat.re <- update(faostat.fe, model = "random")
  faostat.fd <- update(faostat.fe, model = "fd")
  summary(faostat.fe)
  summary(faostat.fd)
}






#### Grunfeld ####
{
  data("Grunfeld", package = "plm")
  grunfeld <- pdata.frame(Grunfeld, 
    index = c("firm", "year"), 
    drop.index = FALSE)
  
  dv.grunfeld <- c("inv")
  iv.grunfeld <- c("value", "capital")
  gid.grunfeld <- c("firm")
  tid.grunfeld <- c("year")
  
  cwe.grunfeld <-
    consistency.within.est(grunfeld, dv.grunfeld, iv.grunfeld, gid.grunfeld, tid.grunfeld)
  cnv.grunfeld <-
    consistency.within.est(grunfeld, dv.grunfeld, iv.grunfeld, gid.grunfeld, tid.grunfeld, no.ginv = TRUE)
  #plot(cwe.grunfeld)
}

#### Princeton ####
{
  Panel101Princeton <- haven::read_dta("../datasets/Panel101Princeton.dta")
  Panel101Princeton <- transform(Panel101Princeton, country = as.numeric(country))
  
  prince <- pdata.frame(Panel101Princeton,
    index = c("country", "year"))
  
  dv.prince <- c("y")
  iv.prince <- c("x1", "x2", "x3")
  gid.prince <- c("country")
  tid.prince <- c("year")
  
  cwe.prince <-
    consistency.within.est(prince, dv.prince, iv.prince, gid.prince, tid.prince)
  cnv.prince <-
    consistency.within.est(prince, dv.prince, iv.prince, gid.prince, tid.prince, no.ginv = TRUE)
  #plot(cwe.prince)
}


#### Produc ####
{
  data("Produc", package="plm")
  
  Produc <- Produc %>% 
    mutate(
      id = as.integer(state),
      lgsp = log(gsp),
      lpcap = log(pcap),
      lpc = log(pc),
      lemp = log(emp),
      lunemp = log(unemp))
  
  produc <- pdata.frame(Produc, index = c("id", "year"))
  
  dv.produc <- c("lgsp")
  iv.produc <- c("lpcap", "lpc", "lemp", "lunemp")
  gid.produc <- c("id")
  tid.produc <- c("year")
  
  cwe.produc <- 
    consistency.within.est(produc, dv.produc, iv.produc, gid.produc, tid.produc)
  cnv.produc <- 
    consistency.within.est(produc, dv.produc, iv.produc, gid.produc, tid.produc, no.ginv = TRUE)
  #plot(cwe.produc)
}

