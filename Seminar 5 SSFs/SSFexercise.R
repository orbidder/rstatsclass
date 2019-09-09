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

# Save the class of data for later!
trk.class<-class(puma_tr1)

# With our track, we can calculate a number of movement parameters.
# Here we are calculating the absolute direction of the step, the relative direction
#   of the step, the step length, and the net squared displacement
puma_tr1 %>% nest(-id) %>% 
  mutate(dir_abs = map(data, direction_abs,full_circle=TRUE, zero="N"), 
         dir_rel = map(data, direction_rel), 
         sl = map(data, step_lengths),
         nsd_ = map(data, nsd)) %>% unnest() -> puma_tr1

# You might have noticed that the bearing covariates were in radians
#   For ease of use, we'll turn them into degrees
puma_tr1 %>%
  mutate(dir_abs = as_degree(dir_abs),
         dir_rel = as_degree(dir_rel)) -> puma_tr1

# If we want to subset the data by year or month, or use time as an interacting
#   term with a covariate in our model, we can add some columns for time into our track
# When might you expect time to interact with a covariate?
puma_tr1%>% 
  mutate(
    week = week(t_),
    month = month(t_, label=TRUE), 
    year = year(t_),
    hour = hour(t_)) -> puma_tr1

# You'll notice that now our data are not a track anymore
class(puma_tr1)

# So we need to change our tibble back to a track
class(puma_tr1)<-trk.class

# Check out the sampling rate of each of our study animals
# Do the medians look about right? What do you notice about the max values?
# Take a close look at the 3rd quartile. Which animals do you think have 
#   lower fix success?
summarize_sampling_rate_many(puma_tr1, id)

# First, we can look at the distribution of our turning angles as a rose plot
ggplot(puma_tr1, aes(x = dir_rel, y=..density..)) + geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Angles Direct") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20), 
                     labels = seq(-180, 180, by=20))+
  facet_wrap(~id)

# Does the rose plot make your head spin? You can aslo look at turning angles as 
#   a traditional histogram
ggplot(puma_tr1, aes(x = dir_rel)) +  geom_histogram(breaks = seq(-180,180, by=20))+
  theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Angles Relative") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20),
                     labels = seq(-180, 180, by=20))+facet_wrap(~id, scales="free")

# We can also look at how movement parameters might differ between day and night
# Here, we're looking at the log of step length
# Why do we take the log of the step length? Modify the plot code to see!
ggplot(puma_tr1, aes(x = tod_, y = log(sl))) + 
  geom_boxplot()+geom_smooth()+facet_wrap(~id)


# Next we eliminate steps that are longer than 3 hours apart, filter out bursts
#   that only have one point, simulate steps from our distribution of step lengths
#   and turn angles, and extract the covariates for each step. 
#   All in one pipe!
puma_tr1 %>% group_by(id) %>% nest() %>% 
  mutate(data = map(
    data, ~ .x %>% 
      # Eliminate steps longer than 3 hours
      track_resample(rate = hours(3), tolerance = minutes(15)) %>%
      # Eliminate bursts with fewer than 2 points
      filter_min_n_burst(min_n = 2) %>% 
      steps_by_burst() %>% 
      # Randomly sample 10 steps per real step
      random_steps(n = 10) %>% 
      # Extract covariates from our raster stack
      extract_covariates(envtrasters))) %>% unnest() -> ssfdat

# A note on extracting covariates:
# Above, we've extracted the covariates at the end points of each step
# In some cases, it may make more sense to extract covariates ALONG a step
# So, rather than ask: "did the animal end up in forest habitat?" you can ask:
#   "what proportion of forest did they move through on their movement path?"
#https://rdrr.io/cran/amt/man/extract_covariates.html

# Scale the covariates, turn "case" into a binary 1/0 variable, 
#   and add some movement parameters that we can use as covariates
ssfdat %>% mutate(elev.s = scale(dem), 
                  slope.s = scale(slope),
                  tri.s = scale(tri),
                  ndvi.s = scale(max_ndvi),
                  case_ = as.numeric(case_),
                  cos_ta = cos(ta_), 
                  log_sl = log(sl_)) -> ssfdat

# Write out your candidate models
m0 <- clogit(case_ ~ elev.s + tri.s + ndvi.s + strata(step_id_),method = "efron", robust = TRUE, data = ssfdat)
m1 <- clogit(case_ ~ elev.s + + tri.s + ndvi.s + elev.s:cos_ta + elev.s:log_sl + log_sl * cos_ta + strata(step_id_), method = "efron", robust = TRUE, data = ssfdat)
m2 <- clogit(case_ ~ elev.s + + tri.s + ndvi.s + elev.s:cos_ta + elev.s:log_sl + log_sl + cos_ta + strata(step_id_), method = "efron", robust = TRUE, data = ssfdat)

summary(m0)
summary(m1)
summary(m2)

# Time for model comparison!
# Instead of AIC, for conditional logistic regression models we use QIC
# QIC is better suited to accommodate comparisons within strata while correcting
#   for autocorrelation among strata
# Here's our function for calculating QIC from an SSF model
QIC.coxph <- function(mod, details = FALSE) {
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  return(-2*quasi + 2*trace)
}

# Which model is best supported?
QIC.coxph(m0)
QIC.coxph(m1)
QIC.coxph(m2)


#__________********_________

# The global model is a good start, expecially if we think animals respond 
#   similarly to their environment while moving
# But, we might want to look at individual variation among animals in their movement behavior
# Why might we want to do this? Discussion break!
#...............
#...............
# So let's fit individual SSFs to each of our animals.
# To do this, we can apply the conditional logistic regression model to our nested data
fitted_ssf <- function(data){
  fit_issf(case_ ~ elev.s + tri.s + ndvi.s + strata(step_id_),method = "efron", robust = TRUE, data=data)
}
ssfdat %>%  group_by(id) %>% nest() %>% 
  mutate(mod = map(data, fitted_ssf)) -> m1

# Package broom can help us tidy up out models into a tibble
m1 %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow)) -> m1

# To vizualise our model output, we can reveal the estimates of all iSSFs
m1$tidy

# Next, create data frame with the coefficients from all the SSFs
ssf_coefs <- m1 %>%
  unnest(tidy) 

# Just to make it a little more interesting, we can add the sex of the animals to our vizualization
n.covs -> 3
mutate(ssf_coefs, sex = c(rep(c("f","f","m","m","m","f","m","m","f"),each = n.covs))) -> ssf_coefs

# Plot the coefficients!
ssf_coefs %>% 
  ggplot(., aes(x=term, y=estimate, group = id, col = factor(id), pch = sex)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 0.7), size = 0.8) +
  geom_hline(yintercept = 0, lty = 2)+
  facet_wrap(~term, scales="free") + 
  labs(x = "Covariate", y = "Relative Selection Strength") +
  theme_light()

# Now that you can see the individual variation, what modifications might you make to our model?
