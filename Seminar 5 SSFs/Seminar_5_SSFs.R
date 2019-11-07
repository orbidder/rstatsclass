# Today we'll be going over step-selection functions to model habitat selection
# Specifically, we'll prepare the data, sample available locations, fit three different kinds
#   of models, and assess model performance with cross-validation
# We will talk about SSFs vs. iSSFs, fixed effects models vs. GEEs vs. mixed effects models,
#   and how to think about autocorrelation and correlation within individuals

# Fundimentally, an SSF measures habitat selection along the movement path of an animal
# In comparison to an RSF, it better represents actual selection decisions made by animals 
#   in real time, as they choose to move to specific locations out of a limited possibilies 
#   available to them at each step

# First, install your packages

install.packages("lubridate")
install.packages("raster")
install.packages("amt")
install.packages("sf")
install.packages("tidyverse")
install.packages("survival")
install.packages("lmtest")
install.packages("MuMIn")

# Load your packages

# Syntax
library(tidyverse)
library(lubridate)
# Spatial
library(sf)
library(raster)
# New SSF package!
library(amt)
# Fitting conditional logistic regression models (in a Cox PH framework)
library(survival)
# Model selection using mixed-effects cox models
library(lmtest)
library(MuMIn)

# Load in your environmental covariates
# Unlike last week, when all our rasters were already in a stack, today we have 4 individual files
# So we list the files and then bring them in as a stack using paste0
raster_files <- list.files("Seminar 5 SSFs/",pattern = ".tif") 
envtrasters <- raster::stack(paste0("Seminar 5 SSFs/", raster_files))
names(envtrasters) <- c("dem", "slope",  "tri", "max_ndvi")

# Next, we'll bring in our GPS data, but we won't yet deal with making it spatial
puma_tr1 <- read_csv("Seminar 5 SSFs/puma_data_2015.csv") %>% 
  # First, select just the columns you need for the analysis
  dplyr::select(ID = animals_id, 5:7) %>% 
  # Then, make sure your date/times are in the correct time zone
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan")) 

# To treat each animal differently, we will nest the data by animal ID
puma_tr1 <- puma_tr1 %>%  nest(-ID) 
# Check out the structure of the data
puma_tr1

# Now we'll make a "track", which is used by package amt for movement analysis
# This format helps amt to manage with variable fix rates and fix success so you don't have to!
puma_tr1 <- puma_tr1 %>% 
  mutate(trk = map(data, function(d){
    mk_track(d, longitude, latitude, timestamp, crs = CRS("+init=epsg:4326")) %>% 
      # Transform our latlongs into UTMs
      transform_coords(CRS("+init=epsg:32719"))
  }))

# Are our fix rates normal? We should have steps 3 hours apart
puma_tr1 %>%
  mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  dplyr::select(ID, sr) %>%  unnest

# Nope! So let's eliminate steps that are longer than 3 hours apart and filter out bursts
#   that only have one point
puma_tr1 %>% 
  mutate(steps = map(trk, function(x){
    x %>% 
      # Eliminate steps longer than 3 hours
      track_resample(rate = hours(3), tolerance = minutes(15)) %>%
      # Eliminate bursts with fewer than 2 points
      filter_min_n_burst(min_n = 2)})) -> ssfdat

# Did we get rid of any locations?
# Compare the dimentions of our individual tracks in comparison to the steps we're keeping
ssfdat

# Check out the sampling rate of each of our study animals
# Do the medians look about right? What do you notice about the max values?
# Take a close look at the means. Which animals do you think have 
#   lower fix success?
ssfdat %>%
  mutate(sr = lapply(steps, summarize_sampling_rate)) %>% 
  dplyr::select(ID, sr) %>%  unnest

# Next we simulate steps from our distribution of step lengths and turn angles, 
#   and extract the covariates for each step. 

# A note on extracting covariates:
# Generally we extract covariates at the END points of each step
# However, in some cases, it may make sense to extract covariates ALONG a step
# So, rather than ask: "did the animal end up in forest habitat?" you can ask:
#   "what proportion of forest did they move through on their movement path?"
#https://rdrr.io/cran/amt/man/extract_covariates.html
# In other cases, it might make sense to extract covariates at the BEGINNING of the step
# WHY??? 
#...
#...
# As an interaction term - to see if the start location influences the end location or a 
#   movement parameter
#   e.g. a habitat that is hard to move through
#   e.g. due to group/herd effects
# In the extract_covariates command, use where = "start", "end", or "both" depending on 
#   your goals

# Here we make a function that will apply a bunch of commands to each of our individual animals
# In this one piped function, we can sample random points, make a day/night covariate, and 
#   extract covariates!
ssfdat %>%
  mutate(moddata = map(steps, function (x){
    x %>% 
      steps_by_burst() %>% 
      # Randomly sample 10 steps per real step
      random_steps(n = 10) %>% 
      # Determine day or night from movement track data
      time_of_day(where = "start") %>% 
      # Extract covariates from our raster stack
      amt::extract_covariates(envtrasters,where="both")})) -> ssfdat

# We need to scale our covariates, but right now our data are in different tibbles by individual
# If we want to scale the covariates from all the data, we need to make a single dataframe
# To make one dataframe, we need a column for ID to tell the animals apart
# To do that, we'll pull the IDs and the number of rows from each component of our data list
ssfdat.all <- do.call(rbind,ssfdat$moddata)
ID<-c()
for (i in 1:length(ssfdat[[1]])) {
   id <- rep(ssfdat[[i,1]],dim(ssfdat$moddata[[i]])[1])
   ID<-c(ID,id)
}
ssfdat.all$ID <- ID

ssfdat.all

# Scale the covariates, turn "case" into a binary 1/0 variable, 
#   and add some movement parameters that we can use as covariates
ssfdat.all %>% mutate(elev_s_start = scale(dem_start), 
                  slope_s_start = scale(slope_start),
                  tri_s_start = scale(tri_start),
                  ndvi_s_start = scale(max_ndvi_start),
                  elev_s_end = scale(dem_end), 
                  slope_s_end = scale(slope_end),
                  tri_s_end = scale(tri_end),
                  ndvi_s_end = scale(max_ndvi_end),
                  case_ = as.numeric(case_),
                  cos_ta = cos(ta_), 
                  log_sl = log(sl_)) -> ssfdat.all

#----------***------------
# Before we run our models, let's visualize our step lengths and turn angles

# Our turn angles are in radians
#   For ease of vizualization, we'll turn them into degrees
ssfdat.all %>%
  mutate(ta_degree = as_degree(ta_)) -> ssfdat.all

# First, we can look at the distribution of our turning angles as a rose plot
ssfdat.all %>% 
  filter(case_ == 1) %>% 
  ggplot(., aes(x = ta_degree, y=..density..)) + geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Angles Direct") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=20), 
                     labels = seq(-180, 180, by=20))+
  facet_wrap(~ID)

# Also check out the distribution of step lengths
ssfdat.all %>% 
  filter(case_ == 1) %>% 
  ggplot(., aes(x = sl_)) +  geom_histogram() + theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Step Lengths") +
  facet_wrap(~ID, scales="free")

# We can also look at how movement parameters might differ between day and night
# Here, we're looking at the log of step length
# Why do we take the log of the step length? Modify the plot code to see!
ssfdat.all %>% 
  filter(case_ == 1) %>% 
  ggplot(., aes(x = tod_start_, y = log(sl_))) + 
  geom_boxplot() + geom_smooth() + 
  facet_wrap(~ID)

# Let's also look at what used and available steps look like
# Here, the black point is our starting location and the colored points are our end locations
ssfdat.all %>% 
  filter(step_id_ == 1, ID == 6) %>% 
  ggplot(.) + geom_point(aes(x = x2_,y = y2_,color = as.factor(case_))) + 
  geom_point(aes(x = x1_, y = y1_))

# Check out some other strata to see how they look by modifying the ID and step_id_

#----------***------------
# A note on behaviors
# You might imagine that animals select habitat really differently if they are in 
#   different behavioral states
#   (i.e. resting, feeding, meandering, directed travel)
# Therefore, it is often wise to seperate your data into states before running SSF analyses
# You can determine the state of your animal in a few different ways
# With just GPS data, you can fit a hidden markov model (HMM) to determine state with the moveHMM package
# If you have fine-scale accelerometer data, you can also use that to determine state
# We won't seperate by states today, but there are some great resources that do this in 
#   the github folder for today
#   (Abrahms et al. 2017, Suraci et al. 2019)

#----------***------------

# Now to the models!
# As with RSFs, start by writing out your candidate models
# Here are just a few examples of ways you can structure your models based on different hypotheses

# m0 just has habitat covariates at the end
# note: amt has a wrapper for clogit called fit_issf. These commands should produce identical outcomes
m0 <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               strata(step_id_),method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)
m0.b <- fit_issf(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               strata(step_id_),method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)

summary(m0)
summary(m0.b)

# Do animals select for ndvi differently between day and night? 
# m1 has habitat covariates at the end with an interaction term between NDVI and time of day
m1 <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               ndvi_s_end:tod_start_ +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)

# Do animals select rugged habitats depending on their speed of movement?
# m2 has habitat covariates the end with an interaction term between TRI at the start and movement rate
m2 <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               tri_s_end:log_sl +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)

# Are animals more likely to stay in the same NDVI as where they started?
m3 <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               ndvi_s_end:ndvi_s_start +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)

# Let's compare our models
summary(m0)
summary(m1)
summary(m2)
summary(m3)

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
QIC.coxph(m3)


# From hab package (in development)
# https://github.com/basille/hab/blob/master/R/kfold.r
library(devtools)
install_github("basille/hab")
library(hab)

# Normally we would want at least 100 repetitions, but in the interest of time we'll use a 
#   smaller number in class today
kfold.CV <- kfold(m2, k = 5, nrepet = 10, jitter = FALSE,
            reproducible = TRUE, details = TRUE)

# The correct validation value is just that of the observed points
kfold.CV %>%
  group_by(type) %>%
  summarize(mean_cor = mean(cor))


#__________********_________

# The global model is a good start, expecially if we think animals respond 
#   similarly to their environment while moving
# But, we might want to look at individual variation among animals in their movement behavior
# Why might we want to do this? Discussion break!
#...............
#...............

# So let's fit individual SSFs to each of our animals.
# First, let's make a function for individual SSF models
fitted_ssf <- function(issf_model){
  fit_issf(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + strata(step_id_),method = "efron", robust = TRUE, data=issf_model)
}

# Next, we can apply the conditional logistic regression model to our nested data
ssfdat.all %>% 
  nest(-ID) %>% 
  mutate(mod = map(data, fitted_ssf)) -> m.ind
    
# Package broom can help us tidy up out models into a tibble
m.ind %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow)) -> m.ind

# To vizualise our model output, we can reveal the estimates of all iSSFs
m.ind$tidy

# To vizualize or report the results, we may want to create a data frame with the 
#   coefficients from all the SSFs
ssf_coefs <- m.ind %>%
  unnest(tidy) 

# Just to make it a little more interesting, we can add the sex of the animals for our vizualization
n.covs <- 3
unique(m.ind$ID)
mutate(ssf_coefs, sex = c(rep(c("f","f","m","m","m","f","m","m","f"),each = n.covs))) -> ssf_coefs

# Plot the coefficients!
ssf_coefs %>% 
  ggplot(., aes(x=term, y=estimate, group = ID, col = factor(ID), pch = sex)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 0.7), size = 0.8) +
  geom_hline(yintercept = 0, lty = 2)+
  facet_wrap(~term, scales="free") + 
  labs(x = "Covariate", y = "Relative Selection Strength") +
  theme_light()

# Now that you can see the individual variation, what modifications might you make to our model?

#-----------***--------------
# Above, we made models for different individuals, but we might want a single population model
# There are multiple ways to integrate correlated data into a single model
# The first way is to build in a correlation structure in your model
# We can do this by adding a cluster() parameter that seperates your data into clusters 
#   of correlated data
# The robust standard errors are now fit using generalized estimating equations (GEE), 
#   which correct for autocorrelated data
# Note that using GGEs are only appropriate if your data are temporally autocorrelated (see Prima et al. 2017)
# Important!!: For GGEs to do an adequate job fitting robust SEs, there need to be 
#   enough clusters (>30) (see Prima et al. 2017)
# Therefore, when we use cluster(), we might want to seperate individual animals into multiple clusters
#   by breaking individual animal data into groups 
# However, becasue there is  temporal autocorrelation WITHIN individual animals,
#   we need to use destructive sampling to make sure the groups are not temporally autocorrelated, 
#   so we remove data between groups for the time period in which temp AC exists based on the 
#     ACF (autocorrelation function)
# This is called "destructive sampling" because you are actually throwing away data
# The GEE approach DOES NOT impact the fit of the coefficent values, but only helps to accurately 
#   calculate the robust standard errors 
#     (by accounting for heteroskedasticity, i.e. differential variance among sub-populations)

# Let's revisit our fixed effects models and add a cluster term

# Although we should run an ACF analysis to make eliminate autocorrelation between clusters, 
#   for now (just to see how the model works) we're just going to make clusters by animal ID and year
ssfdat.all %>% 
  mutate(year = year(t1_)) %>% 
  unite("clust_id_yr",c(ID,year),remove = FALSE) -> ssfdat.all
table(ssfdat.all$clust_id_yr)

m0_c <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + cluster(clust_id_yr) +
               strata(step_id_),method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)
m1_c <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + ndvi_s_end:tod_start_ + cluster(clust_id_yr) +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)
m2_c <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               tri_s_end:log_sl + cluster(clust_id_yr) +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)
m3_c <- clogit(case_ ~ elev_s_end + tri_s_end + ndvi_s_end + 
               ndvi_s_end:ndvi_s_start + cluster(clust_id_yr) +
               strata(step_id_), method = "efron", robust = TRUE, data = ssfdat.all, model = TRUE)


# Compare our cluster models to our original fixed effects conditional logistic regression models
# What is different?
summary(m0)
summary(m0_c)

summary(m1)
summary(m1_c)


# The second way to deal with correlated data is using random effects (like we did with RSFs)
# The assumptions here are different! We don't need to have temporal autocorrelation within 
#   individuals but we expect that habitat selection varies among individuals
# Unlike the GEE models that just correct your SEs, mixed-effects models will also change 
#   your coefficient estimates

# Random effects SSF with the coxme package

m0_me<-coxme(Surv(rep(1, length(ssfdat.all$ID)), case_) ~ elev_s_end + tri_s_end + ndvi_s_end +
                 strata(step_id_) + (1|ID),
               data = ssfdat.all, na.action = na.fail)
m1_me<-coxme(Surv(rep(1, length(ssfdat.all$ID)), case_) ~ elev_s_end + tri_s_end + ndvi_s_end + 
                ndvi_s_end:tod_start_ +
                strata(step_id_) + (1|ID),
                data = ssfdat.all, na.action = na.fail)
m2_me<-coxme(Surv(rep(1, length(ssfdat.all$ID)), case_) ~ elev_s_end + tri_s_end + ndvi_s_end + 
                tri_s_end:log_sl +
                strata(step_id_) + (1|ID),
                data = ssfdat.all, na.action = na.fail)
m3_me<-coxme(Surv(rep(1, length(ssfdat.all$ID)), case_) ~ elev_s_end + tri_s_end + ndvi_s_end + 
               ndvi_s_end:ndvi_s_start +
               strata(step_id_) + (1|ID),
             data = ssfdat.all, na.action = na.fail)

summary(m0_me)

lrtest(m0_me,m1_me,m2_me,m3_me)
model.sel(m0_me,m1_me,m2_me,m3_me,rank=AIC)

# You might want to know if a mixed effects model is better than your fixed effects one
# Use the below function to test to see which one is more supported
TestRanef <- function(coxphModel,coxmeModel) {
  if (!class(coxphModel) == "coxph" | !class(coxmeModel) == "coxme") {
    stop("Wrong models")
  }
  ## Degrees of freedom
  coxphDf <- sum(!is.na(coef(coxphModel)))
  coxmeDf <- coxmeModel$df
  names(coxmeDf) <- c("Integrated","Penalized")
  ## DF differnces
  dfDiff <- coxmeDf - coxphDf
  ## Log likelihodds
  coxphLogLik <- coxphModel$loglik[2]
  coxmeLogLik <- coxmeModel$loglik + c(0, 0, coxmeModel$penalty)
  coxmeLogLik <- coxmeLogLik[2:3]
  ## -2 logLik difference
  logLikNeg2Diff <- c(-2 * (coxphLogLik - coxmeLogLik["Integrated"]),
                      -2 * (coxphLogLik - coxmeLogLik["Penalized"]))
  ## p-values
  pVals <- pchisq(q = logLikNeg2Diff, df = dfDiff, lower.tail = FALSE)
  ## Combine
  outDf <- data.frame(dfDiff, logLikNeg2Diff, pVals)
  colnames(outDf) <- c("df diff", "-2logLik diff", "p-values")
  outDf
}

# For the function to work, we have to make our model only class "coxph"
class(m0)
class(m0) <- "coxph"
TestRanef(m0,m0_me)
