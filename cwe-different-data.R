


{
  plot(cwe.wages)
 plt.cwe(cwe.wages)
 plt.cnv(cnv.wages, 0.0002, 0.2)
  plot(cwe.faostat)
 plt.cwe(cwe.faostat)
 plt.cnv(cnv.faostat, 0.03, 3)
  plot(cwe.grunfeld)
 plt.cwe(cwe.grunfeld)
 plt.cnv(cnv.grunfeld, 0.01, 0.5)
 #plot(cwe.prince)
 #plt.cwe(cwe.prince)
 #plt.cnv(cnv.prince)
  plot(cwe.produc)
 plt.cwe(cwe.produc)
 plt.cnv(cnv.produc, 0.01, 2)
  
  # save all open plots
  {
    plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
    plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
    file.copy(from=plots.png.paths, to="plot-pakket")
    
    plots.png.detials <- file.info(plots.png.paths)
    plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
    sorted.png.names <- gsub(plots.dir.path, "plot-pakket", row.names(plots.png.detials), fixed=TRUE)
    numbered.png.names <- paste0("plot-pakket/", 1:length(sorted.png.names), ".png")
    
    # Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
    file.rename(from=sorted.png.names, to=numbered.png.names)
  }
}



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
  
  # within estimator
  wages.form <- lwage ~ 
    exp + I(exp ** 2) + wks + bluecol + ind + south + smsa + married + union +
    sex + black + ed
  #wages.form.old <- as.formula(paste(dv.wages, paste(iv.wages, collapse=" + "), sep=" ~ "))
  wages.po <- plm(wages.form, data = wages, model = "pooling")
  wages.fe <- plm(wages.form, data = wages, model = "within")
  wages.re <- plm(wages.form, data = wages, model = "random")
  wages.fd <- plm(wages.form, data = wages, model = "fd")
  # HT - Hausman Taylor (1981)
  ht <- plm(
    lwage ~ 
      exp + I(exp ** 2) + wks + bluecol + ind + south + smsa + married + union +
      sex + black + ed |
      bluecol + south + smsa + ind + sex + black |
      wks + married + union + exp + I(exp ^ 2), 
    data = wages,
    random.method = "ht", model = "random", inst.method = "baltagi")
  summary(ht)
  summary(wages.po)
  summary(wages.fe)
  summary(wages.re)
  summary(wages.fd)
}



#### Practical Test for Strict Exogeneity ... by Liangjun Su, Yonghui Zhang et al. ####
{
  output.full <- 
    readr::read_csv(file = "../datasets/faostat/Macro-Statistics_Key_Indicators_E_All_Data_(Normalized).csv") %>% 
    dplyr::filter(Item == "Value Added (Agriculture, Forestry and Fishing)") %>% 
    dplyr::filter(Element == "Value US$") %>% 
    select(AreaCode = "Area Code", Area, Year, Output = Value)
      
  area.full <-
    readr::read_csv(file = "../datasets/faostat/Inputs_LandUse_E_All_Data_(Normalized).csv") %>% 
    dplyr::filter(Item == "Agricultural land") %>% 
    select(AreaCode = "Area Code", Area, Year, AgrArea = Value) 
      
  pop.full <- 
    readr::read_csv(file = "../datasets/faostat/Employment_Indicators_E_All_Data_(Normalized).csv") %>% 
    dplyr::filter(Indicator == "Employment in agriculture") %>% 
    select(AreaCode = "Area Code", Area, Year, Pop = Value) 
      
  govt.full <- 
    readr::read_csv(file = "../datasets/faostat/Investment_GovernmentExpenditure_E_All_Data_(Normalized).csv") %>% 
    dplyr::filter(Element == "Value US$") %>% 
    dplyr::filter(Item == "Agriculture, forestry, fishing (General Government)") %>% 
    select(AreaCode = "Area Code", Area, Year, Govt = Value)
    
  
  from.Year <- 2002
  till.Year <- 2007
  output <- output.full %>% dplyr::filter(Year >= from.Year & Year <= till.Year)
  area <- area.full     %>% dplyr::filter(Year >= from.Year & Year <= till.Year)
  pop <- pop.full       %>% dplyr::filter(Year >= from.Year & Year <= till.Year)
  govt <- govt.full     %>% dplyr::filter(Year >= from.Year & Year <= till.Year)
  faostat <- output %>% left_join(area) %>% left_join(pop) %>% left_join(govt) %>% 
    group_by(Area) %>% dplyr::filter(all(!is.na(AgrArea)) & all(!is.na(Pop)) & all(!is.na(Govt))) %>% 
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
  plot(cwe.faostat)
  
  
  faostat.form <- loutput ~ lagrarea + lpop + lgovt
  faostat.fe <- plm(faostat.form, data = faostat, model = "within")
  faostat.po <- update(faostat.fe, model = "pooling")
  faostat.re <- update(faostat.fe, model = "random")
  faostat.fd <- update(faostat.fe, model = "fd")
  summary(faostat.po)
  summary(faostat.re)
  summary(faostat.fe)
  summary(faostat.fd)
}



#### Labor Supply ####
{
  data(LaborSupply, package = "Ecdat")
  LaborSupply$age.sq <- LaborSupply$age ** 2
  labor <- pdata.frame(LaborSupply,
    index = c("id", "year"),
    drop.index = FALSE)
  
  dv.labor <- c("lnhr")
  iv.labor <- c("lnwg", "age", "age.sq", "kids", "disab")
  gid.labor <- c("id")
  tid.labor <- c("year")
  
  cwe.labor <-
    consistency.within.est(labor, dv.labor, iv.labor, gid.labor, tid.labor, try.normal.inv = TRUE)
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
  #plot(cwe.grunfeld)
  
  grunfeld.form <- inv ~ value + capital
  grunfeld.fe <- plm(grunfeld.form, data = grunfeld, model = "within")
  grunfeld.po <- update(grunfeld.fe, model = "pooling")
  grunfeld.re <- update(grunfeld.fe, model = "random")
  grunfeld.fd <- update(grunfeld.fe, model = "fd")
  summary(grunfeld.po)
  summary(grunfeld.re)
  summary(grunfeld.fe)
  summary(grunfeld.fd)
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
    consistency.within.est(prince, dv.prince, iv.prince, gid.prince, tid.prince, try.normal.inv = TRUE)
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
  iv.produc <- c("lpcap", "lpc", "lemp", "unemp")
  gid.produc <- c("id")
  tid.produc <- c("year")
  
  cwe.produc <- 
    consistency.within.est(produc, dv.produc, iv.produc, gid.produc, tid.produc)
  #plot(cwe.produc)
}

