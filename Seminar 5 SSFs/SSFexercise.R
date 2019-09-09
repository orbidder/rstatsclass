library(lubridate)
library(raster)
library(amt)
library(sf)
library(dplyr)
library(tidyverse)

envtrasters <- stack("/Users/justinesmith/Documents/UCB/Data/Rasters/SGNP_all_stack_south_crop.tif") 
names(envtrasters) <- c("aspect", "dem", "rough",  "slope",  "tri", "max_ndvi")

read_csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma_data 3hr2.csv") %>% 
  # First, select just the columns you need for the analysis
  dplyr::select(ID = animals_id, 5:7) %>% 
  # Then, make sure your date/times are in the correct time zone
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan")) %>% 
  # Make a "track", which is used by package amt for movement analysis
  mk_track(., longitude, latitude, timestamp, id=ID, crs = CRS("+init=epsg:4326")) %>% 
  # Very fast way to determine day or night from movement data!
  time_of_day() %>%
  # Transform our latlongs into UTMs
  transform_coords(., CRS("+init=epsg:32719")) -> puma_tr1

trk.class<-class(puma_tr1)

puma_tr1 %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_=map(data, nsd)) %>% unnest() -> puma_tr1

puma_tr1 %>%
  mutate(dir_abs = as_degree(dir_abs),
         dir_rel = as_degree(dir_rel)) -> puma_tr1

puma_tr1%>% 
  mutate(
    week=week(t_),
    month = month(t_, label=TRUE), 
    year=year(t_),
    hour = hour(t_)) -> puma_tr1

class(puma_tr1)

class(puma_tr1)<-trk.class

summarize_sampling_rate_many(puma_tr1, id)

ggplot(puma_tr1, aes(x = dir_rel, y=..density..)) + geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Angles Direct") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20), 
                     labels = seq(-180, 180, by=20))+
  facet_wrap(~id)

ggplot(puma_tr1, aes(x = dir_rel)) +  geom_histogram(breaks = seq(-180,180, by=20))+
  theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Angles Relative") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20),
                     labels = seq(-180, 180, by=20))+facet_wrap(~id, scales="free")

ggplot(puma_tr1, aes(x = tod_, y = log(sl))) + 
  geom_boxplot()+geom_smooth()+facet_wrap(~id)



puma_tr1 %>% group_by(id) %>% nest() %>% 
  mutate(data = map(
    data, ~ .x %>% 
      track_resample(rate = hours(3), tolerance = minutes(15)) %>%
      filter_min_n_burst(min_n = 2) %>% 
      steps_by_burst() %>% 
      random_steps(n = 10) %>% 
      extract_covariates(envtrasters))) %>% unnest() -> ssfdat

# A note on extracting covariates:
#https://rdrr.io/cran/amt/man/extract_covariates.html

class(ssfdat)<-trk.class


ssfdat %>% mutate(elev.s = scale(dem), 
                  slope.s = scale(slope),
                  tri.s = scale(tri),
                  ndvi.s = scale(max_ndvi),
                  step.id = with(ssfdat, paste0(id, step_id_, sep = "")),
                  case_ = as.numeric(case_),
                  cos_ta = cos(ta_), 
                  log_sl = log(sl_)) -> ssfdat


m0 <- clogit(case_ ~ elev.s + tri.s + ndvi.s + strata(step_id_),method = "efron", robust = TRUE, data = ssfdat)
m1 <- clogit(case_ ~ elev.s + + tri.s + ndvi.s + elev.s:cos_ta + elev.s:log_sl + log_sl * cos_ta + strata(step_id_), method = "efron", robust = TRUE, data = ssfdat)
m2 <- clogit(case_ ~ elev.s + + tri.s + ndvi.s + elev.s:cos_ta + elev.s:log_sl + log_sl + cos_ta + strata(step_id_), method = "efron", robust = TRUE, data = ssfdat)

summary(m0)
summary(m1)
summary(m2)

QIC.coxph <- function(mod, details = FALSE) {
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  return(-2*quasi + 2*trace)
}

QIC.coxph(m0)
QIC.coxph(m1)
QIC.coxph(m2)








fitted_ssf <- function(data){
  fit_issf(case_ ~ elev.s + tri.s + ndvi.s + strata(step_id_),method = "efron", robust = TRUE, data=data)
}
m1 <-ssfdat %>%  group_by(id) %>% nest() %>% 
  mutate(mod = map(data, fitted_ssf)) 

m1 <- m1 %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow))

m1$tidy

#' Now, create data frame w/ the coefficients, etc
ssf_coefs <- m1 %>%
  unnest(tidy) 

ssf_coefs %>% spread(term, estimate)

#' Plot coefficients
#+ fig.width=12, fig.height= 8
ssf_coefs %>% 
  ggplot(., aes(x=term, y=estimate,group = id, col = id, pch = sex)) + 
  geom_dotplot(binaxis="y", stackdir="center")+geom_hline(yintercept=0)+
  facet_wrap(~term, scales="free")

p1 <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
  select(id, sex, coef) %>% unnest %>%
  mutate(id = factor(id)) %>%
  ggplot(., aes(x = term, y = estimate, group = id, col = id, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin,
                          ymax = ymax), data = d2, inherit.aes = FALSE,
            fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,
                             y = mean, yend = mean), data = d2, inherit.aes = FALSE,
               size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 0.7), size = 0.8)  +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Habitat", y = "Relative Selection Strength") +
  theme_light() +
  scale_x_discrete(labels = c("Developed (open)", "Developed (other)", "Natural", "Crops"))
p1