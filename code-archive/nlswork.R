data("nlswork", package = "sampleSelection")
nlswork <- nlswork %>% 
  mutate(ttl_exp_sq = ttl_exp**2) %>% 
  mutate(id = idcode) %>% 
  select(!idcode)

# drop non-full time participants
ids <- nlswork %>% 
  group_by(id) %>% 
  summarise(yrs = n()) %>% 
  filter(yrs == 15) %>% 
  select(id)

full.yrs <- nlswork %>% 
  filter(id %in% ids$id) %>% 
  select(year)
  

  
nlswork <- pdata.frame(nlswork, 
  c("id", "year")) 

dv.nlswork <- c("ln_wage")
iv.nlswork <- c("ttl_exp", "ttl_exp_sq")
gid.nlswork <- c("id")
tid.nlswork <- c("year")

cwe.nlswork <- 
  consistency.within.est(nlswork, dv.nlswork, iv.nlswork, gid.nlswork, tid.nlswork)
plot(cwe.nlswork)
cnv.nlswork <- 
  consistency.novar(nlswork, dv.nlswork, iv.nlswork, gid.nlswork, tid.nlswork)