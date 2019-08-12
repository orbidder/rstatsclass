library(lubridate)
library(raster)
library(amt)

data("deer")
deer

summarize_sampling_rate(deer)

data("sh_forest")
sh_forest

deer %>% 
  steps_by_burst() %>% random_steps(n = 15) %>% 
  extract_covariates(sh_forest) %>% 
  mutate(forest = factor(sh.forest, levels = 1:2, labels = c("forest", "non-forest")), 
         cos_ta = cos(ta_), 
         log_sl = log(sl_)) -> ssf1

m0 <- ssf1 %>% fit_clogit(case_ ~ forest + strata(step_id_))
m1 <- ssf1 %>% fit_clogit(case_ ~ forest + forest:cos_ta + forest:log_sl + log_sl * cos_ta + strata(step_id_))
m2 <- ssf1 %>% fit_clogit(case_ ~ forest + forest:cos_ta + forest:log_sl + log_sl + cos_ta + strata(step_id_))

summary(m0)
summary(m1)
summary(m2)

QIC.coxph(m0, details=TRUE)
QIC.coxph(m1, details=TRUE)
QIC.coxph(m2, details=TRUE)

