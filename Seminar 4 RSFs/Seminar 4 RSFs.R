# Today we'll be doing some basic habitat selection analyses by fitting resource selection functions (RSFs)
# Because animal movement data only gives us used (and non unused) locatations, we generally can't
#   fit resource selection probability functions, or RSPFs. 
# However, by randomly sampling available locations, we can approximate habitat selection with RSFs.
# We will go through how to sample available locations and model habitat selection today using 
#   GPS data from vicuñas.

install.packages("tidyverse")
install.packages("lubridate")
install.packages("maptools")
install.packages("adehabitatHR")
install.packages("rgeos")
install.packages("rgdal")
install.packages("sf")
install.packages("raster")
install.packages("lme4")
install.packages("AICcmodavg")
install.packages("AUC")
install.packages("ROCR")
install.packages("mapview")
install.packages("RColorBrewer")
install.packages("sp")
install.packages("gridExtra")

# First, load the packages that you will be using

# Syntax and formatting
library(tidyverse)
library(lubridate)
library(maptools)
# Home ranges
library(adehabitatHR)
# Spatial data
library(rgeos)
library(rgdal)
library(raster)
library(sf)
# Model fitting and selection
library(lme4)
library(AICcmodavg)
# ROC curves
library(AUC)
library(ROCR)
# Visualizations
library(mapview)
library(RColorBrewer)
library(sp)
library(gridExtra)

# Read in, clean, and project all your data - in a single pipe!
# use read_csv to bring in data
vicuna <- read_csv("Seminar 4 RSFs/vicuna_data_2015.csv") %>% 
  # We'll make a timestamp column with the correct time zone
  # in this file, I've already corrected the time to be in Argentina hours, so we just need to assign the tz
  # we will also create a day/night column to differentiate day and night
  mutate(timestamp = force_tz(acquisition_time,tz="America/Argentina/San_Juan"),
         sunrise = sunriset(SpatialPoints(cbind(longitude,latitude), proj4string=CRS("+init=epsg:4326")), timestamp, direction="sunrise", POSIXct.out=TRUE)[,2],
         sunset = sunriset(SpatialPoints(cbind(longitude,latitude), proj4string=CRS("+init=epsg:4326")), timestamp, direction="sunset", POSIXct.out=TRUE)[,2],
         daynight = ifelse(timestamp>sunrise&timestamp<sunset,1,0)) %>% 
  # then select only the columns we need for analysis - animal ID, datetime, coordinates, and daynight
  dplyr::select(ID = animals_id, 12, 6:7,15) %>% 
  # then create an sf object by calling the coordinates columns
  # it's very important that you know the coordinate reference system that your data are in!
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>% 
  # next, let's transform our latitude and longitude columns to UTMs
  st_transform("+init=epsg:32719")
  
# Always look at your data!
vicuna

# Now we can pull in our environmental data, which is in the form of a raster stack with 4 raster layers
# Our environmental data include layer for digital elevation model (or dem, which is just elevation),
#   slope, terrain ruggedness index (mean difference between a central pixel and its surrounding cells),
#   and max NDVI from 2015 (normalized difference vegetation index, or a measure of greenness)
envtrasters <- stack("Seminar 4 RSFs/vicuna_envt_layers.tif")
names(envtrasters) <- c("dem","slope","tri","max_ndvi")

# Let's look at the structure of the data
envtrasters

# Visualize environmental layers to make sure they make sense
NDVI.raster <- as.data.frame(envtrasters[[4]], xy = TRUE)
ggplot() +
  geom_raster(data = NDVI.raster, aes(x = x, y = y, fill = max_ndvi)) +
  scale_fill_viridis_c() +
  coord_sf()

# You can also visualize all the layers at once in base R
plot(envtrasters)

# Choosing an available range
# A conservative option is to use a 99% KUD (kernel utilization distribution)
# We can also use a 95% KUD if we want to make sure to focus on areas used more frequently by individuals
# To get unique home ranges, we'll use the kernel density estimation that you learned in Seminar 2
vicuna_ud <- adehabitatHR::kernelUD(as(vicuna, "Spatial")[1],  
                                   grid = 450)
vicuna_hr <- adehabitatHR::getverticeshr(vicuna_ud, 95)

# Check out how the home ranges stack onto the dem layer
mapview(envtrasters[[2]]) + mapview(vicuna_hr[1])
# Or visualize individual home ranges on the NDVI layer
mapview(envtrasters[[4]]) + mapview(vicuna_hr[18,1])

# Keep in mind temporal effects!
# RSF may be a spatial analysis but if we think selection or availability may change over time, 
#   we should consider that in our modeling approach
# We can control for time either via appropriate covariate inclusion, or simply by limiting our scope.
# Let's say we only want to know how vicuñas select habitat in winter (June-Sept) during the day

vicuna %>% 
  # filter by months June-Sept and daytime
  filter(lubridate::month(timestamp) %in% 6:9, daynight == 1) %>% 
  # make a new column to show that these are real GPS data from collared vicuñas
  mutate(Used = 1) -> vicuna


#--------
# Now that we have created our home ranges, we can sample available points within each home range. 

# First, determine how many points to sample in each home range. 
# Here we will calculate a 1:1 or 10:1 ratio of available:used points and randomly sample the available points from the individual UDs

# Extract a list of vicuña IDs
print(vicuna_ids <- unique(vicuna$ID))
# Make an empty list to store your new locations by individual vicuña
availables <- list()

# Randomly sample points from within the home range of each individual according to your chosen ratio
for(i in 1:length(vicuna_ids)){
  st_sample(st_as_sf(vicuna_hr)[i,], nrow(filter(vicuna, ID == vicuna_ids[i]))) %>% # e.g. 1:1 ratio
    #st_sample(st_as_sf(vicuna_hr)[i,], 10*nrow(filter(vicuna, ID == vicuna_ids[i])))  %>%  # e.g. 1:10 ratio
    st_sf(geometry = .) %>%
    # We'll save the list with columns that match our available data, including a "Used" column populated with zeroes
    mutate(ID = vicuna_ids[i], 
           timestamp = NA,
           daynight = NA,
           Used = 0) -> availables[[i]] 
}

# Then combine individual lists into one data.frame and combine it with the real vicuña GPS points
availables %>% 
  do.call(rbind,.) %>%
  rbind(vicuna, .) -> vicuna_all

# Do we have a 1:1 ratio of used:available?
table(vicuna_all$ID,vicuna_all$Used)

# Let's take a look at some of our used and available data
plot(crop(envtrasters[[4]],extent(vicuna_hr[3,])))
plot(vicuna_hr[3,],add=T)
plot(st_geometry(filter(vicuna_all,ID == 16,Used == 0)),col="red",add=T)
plot(st_geometry(filter(vicuna_all,ID == 16,Used == 1)),add=T)


#--------

# Next we have to extract environmental covariates and process them for analysis
# We will use raster::extract to get environmental values for each used and available location
# We will also scale and center environmental covariates
# There are two main reasons to scale covariates:
#   1. it allows to to compare the relative strength of each variable by looking at the coefficient values alone
#   2. sometimes if your variables are on really different scales, the model won't converge
vicuna_all %>% mutate(NDVI = raster::extract(envtrasters[[4]], as(., "Spatial")),
                     elev = raster::extract(envtrasters[[1]], as(., "Spatial")),
                     slope = raster::extract(envtrasters[[2]], as(., "Spatial")),
                     tri = raster::extract(envtrasters[[3]], as(., "Spatial")),
                     NDVI.scaled = scale(NDVI, center = TRUE, scale = TRUE),
                     elev.scaled = scale(elev, center = TRUE, scale = TRUE),
                     slope.scaled = scale(slope, center = TRUE, scale = TRUE),
                     tri.scaled = scale(tri, center = TRUE, scale = TRUE)) -> vicuna_full

# Save your scale parameters - you'll need them for plotting later!
print(elevscalelist <- list(scale = attr(vicuna_full$elev.scaled, "scaled:scale"),center = attr(vicuna_full$elev.scaled, "scaled:center")))
print(slopescalelist <- list(scale = attr(vicuna_full$slope.scaled, "scaled:scale"),center = attr(vicuna_full$slope.scaled, "scaled:center")))
print(triscalelist <- list(scale = attr(vicuna_full$tri.scaled, "scaled:scale"),center = attr(vicuna_full$tri.scaled, "scaled:center")))
print(NDVIscalelist <- list(scale = attr(vicuna_full$NDVI.scaled, "scaled:scale"),center = attr(vicuna_full$NDVI.scaled, "scaled:center")))


#-----------
# Fitting the RSF

# Always check the structure of the data
# For example, look at the center parameter for NDVI.scaled. Are the values positive when NDVI is above that value?
vicuna_full

# When we run our model, we want to make sure that the animal ID is treated as a factor
vicuna_full$ID <- as_factor(vicuna_full$ID)

# Before we run any models, we need to check for colinearity among habitat covariates 
# If two covariates are correlated, including them in the same model can overinflate their standard errors 
#   and make it appear that they are not good predictors (when they acutally might be!)  
# General practice is to not include two covariates for which r > 0.7 (or, conservatively, 0.5)
vicuna_cor <- vicuna_full
st_geometry(vicuna_cor) <- NULL
cor(vicuna_cor[,c(5:8)])
rm(vicuna_cor)

# Are any of our covariates correlated?


# As a first pass, let's fit a population-level RSF (using fixed effects only)
# For example, we can compare the fit of slope and ruggedness to decide which one we want to include
summary(fixed.global <- glm(Used ~ NDVI.scaled + elev.scaled + slope.scaled, data=vicuna_full, family=binomial(link="logit")))
summary(fixed.global <- glm(Used ~ NDVI.scaled + elev.scaled + tri.scaled, data=vicuna_full, family=binomial(link="logit")))

# Does slope or ruggedness improve model fit more?

# Next, we'll fit a mixed-effects binomial logistic regression
# i.e. an RSF that accounts for variation among individuals

#---------
# BREAK - random effects exercise! 

# First: what is a random effects model?
#   A type of hierarchical linear model, whereby the data being analyzed are drawn from a hierarchy of populations
#   Or, more simply, a model that contains random variables, which fits the coefficient values with random
#     effects independently for each sub-populations

# Let's explore how random effects alter model shape
# To do this, we'll vary the slope and intercept of a logistic regression model and visualize the differences

# First make a function for a binomial logistic regression model
fake.data<-function(x,b,m){
  1/(1+exp(-(b + x*m)))
}

# Set up fake data
fakem0.5<-rep(0.5,times=101)
fakem1.5<-rep(1.5,times=101)
fakem1<-rep(1,times=101)
fakex<-seq(0,20,by=0.2)
fakeb1<-rep(-1,times=101)
fakeb10<-rep(-10,times=101)
fakeb20<-rep(-20,times=101)

# Random intercept plot
fakeyb1<-fake.data(fakex,fakeb1,fakem1)
fakeyb10<-fake.data(fakex,fakeb10,fakem1)
fakeyb20<-fake.data(fakex,fakeb20,fakem1)
plot(fakeyb1~fakex,ylim=c(0,1))
points(fakeyb10~fakex,ylim=c(0,1),col="blue")
points(fakeyb20~fakex,ylim=c(0,1),col="red")

# Random slope plot
fakeym1<-fake.data(fakex,fakeb10,fakem1)
fakeym0.5<-fake.data(fakex,fakeb10,fakem0.5)
fakeym1.5<-fake.data(fakex,fakeb10,fakem1.5)
plot(fakeym1~fakex,ylim=c(0,1))
points(fakeym0.5~fakex,ylim=c(0,1),col="blue")
points(fakeym1.5~fakex,ylim=c(0,1),col="red")

# Discussion
# When might you want to use a random intercept? What about a random slope?

# OK - back to the models!

#---------

# Random intercept model
#   Most common way to include a random effect
#   Addresses variation in sample sizes, moves logistic curve left and right
# Example hypothesis: general strength of selection varies by individual
random.int.global <- glmer(Used ~ NDVI.scaled + elev.scaled + slope.scaled + (1|ID), 
                           data=vicuna_full, family=binomial(link="logit"))
# Look at the AIC and strength of individual covariates
summary(random.int.global)
# Look at coefficient values for all covariates for each individual
coef(random.int.global)
# Isolate just the random effects. How is "ranef" related to the first column of "coef"?
ranef(random.int.global)

# Random slope model
#   Use a random slope if you think individuals respond differently to environmental covariates
#   Risk in over-fitting, but the good news is that you can test for that in the cross-validation process
# Example hypothesis: some individuals select for vegetation, whereas others do not, 
#   because of variation in internal state and risk-foraging tradeoffs
random.slp_ndvi.global <- glmer(Used ~ NDVI.scaled + elev.scaled + slope.scaled + (0+NDVI.scaled|ID), 
                                data=vicuna_full, family=binomial(link="logit"))
summary(random.slp_ndvi.global)
coef(random.slp_ndvi.global)
ranef(random.slp_ndvi.global)

# Both random slope and intercept
# Example hypothesis: some individuals select for ruggedness, whereas others avoid it, 
#   because of variation in boldness, and general strength of selection varies by individual
random.int.slp_tri.global <- glmer(Used ~ NDVI.scaled + elev.scaled + slope.scaled + (1+tri.scaled|ID), 
                                  data=vicuna_full, family=binomial(link="logit"))
summary(random.int.slp_tri.global)
coef(random.int.slp_tri.global)
ranef(random.int.slp_tri.global)


#-----------
#Model selection

# Your turn! Come up with some alternative hypotheses
# List your candidate models below based on your alternative hypotheses
# Run your candidate models
model1<-glmer(Used ~ NDVI.scaled + elev.scaled + slope.scaled + (1|ID), 
              data=vicuna_full, family=binomial(link="logit"))
model2<-
model3<-
model4<-

# Make a list of your models and compare the deltaAIC, log likelihook, and akaike weights
cand.models<-c(model1,model2,model3,model4)
print(df.kill.models <- data.frame(aictab(cand.models)))

#----------
# How do we decide which points are predicted to be used or available?
# Do determine this, we calculate the optimal cutoff between used and available
#   that maximizes sensitivity and specificity
# Sensitivity: true positive rate 
# Specificity: true negative rate 

# Determine model cutoff between "used" and "available". This will help you interpret the visualizations.
# Enter your top model name in the "fitted" command below
fittedmod<-fitted(model1)
target_pred<-as.matrix(fittedmod)
target_class<-as.matrix(as.factor(vicuna_full$Used))
pred<-prediction(target_pred,target_class)
perf<-performance(pred,"tpr","fpr")
fpr = perf@x.values[[1]] 
tpr = perf@y.values[[1]] 
sum = tpr + (1-fpr) 
index = which.max(sum) 
print(cutoff <- perf@alpha.values[[1]][[index]])

# To get an initial feel for how good your model performed, you can look at the sensitivity and specificity
#   in an ROC curve
# A strong model has an area-under-the-curve (AUC) of over 0.8.
# An AUC of 0.5 means the model cannot discriminate between a true and false positive
plot(perf,col="black",lty=3, lwd=3)
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
minauc<-min(round(auc, digits = 2))
minauct <- paste(c("AUC  = "),minauc,sep="")
legend(0.4,0.4,c(minauct,"\n"),border="white",cex=1.7,box.col = "white")

# How does our model look?

# IMPORTANT NOTE ON CROSS-VALIDATION
# This just shows us a first pass. We still need to do cross-validation to check our model fit
# Cross-validation is a very important part of the modeling process
# While selecting based on AIC tells you which of your candidate models is the best
#   it doesn't tell you if your "best" model is actually a GOOD model!
# Cross-validation allows you to see how well your model performs using your validated data
# We will go over cross-validation in Seminar 8, but remember that you should always conduct a
#   CV analysis on your final model to examine model fit


#------------
# Visualizing your model in space:
#   predict your model output to a predictive surface

# Ideally, we would be able to create a spatial predictive surface of habitat selection using our 
#   environmental layers
# To do this, we could try to fit our model using the unscaled variables to match the original raster scales
random.int.raw <- glmer(Used ~ NDVI + elev + slope + (1|ID), 
                          data=vicuna_full, family=binomial(link="logit"))

# You'll note that the model doesn't converge because of the very different scales of the covariates
# When this happens, you have to instead convert the rasters to the scale 
#   of the centered and normalized covariates in the model for mapping
# To do this, we'll revisit our scale parameters for our covariates

# Make a new raster stack with just the covariates in the model
env.rasters <- stack(envtrasters[[4]], envtrasters[[1]], envtrasters[[2]])

# Scale the covariates to be centered and normalized, based on the scale parameters from our model covariates
env.rasters[[1]]<-(env.rasters[[1]]-NDVIscalelist$center)/NDVIscalelist$scale
env.rasters[[2]]<-(env.rasters[[2]]-elevscalelist$center)/elevscalelist$scale
env.rasters[[3]]<-(env.rasters[[3]]-slopescalelist$center)/slopescalelist$scale

# Name the raster layers to match the model covariate names
names(env.rasters) <- c("NDVI.scaled", "elev.scaled", "slope.scaled")

# Time to make a predicted raster layer based on your model!
predictionmap<-raster::predict(object=env.rasters,model=model1,
                               re.form = NA,type="response")

# To visualize the differences among individuals, let's also map the individuals with the 
#   greatest deviations from the model intercept
ranef(model1)
predictionmap.2<-raster::predict(object=env.rasters,model=model1,
                               const=(data.frame(ID="16")),type="response")
predictionmap.3<-raster::predict(object=env.rasters,model=model1,
                               const=(data.frame(ID="20")),type="response")

# Three types of visualizations!

# Visualization 1 in base::plot
plot(predictionmap)

# Visualization 2 in sp::spplot
breaks.plot<-c(seq(0,1,by=0.1))
my.palette <- rev(brewer.pal(n = 10, name = "RdYlBu"))
spplot(stack(predictionmap.2,predictionmap,predictionmap.3),
       col.regions=my.palette,at=breaks.plot)

# Visualization 3 in ggplot2::ggplot
prediction.raster <- as.data.frame(predictionmap, xy = TRUE)
prediction.raster.2 <- as.data.frame(predictionmap.2, xy = TRUE)
prediction.raster.3 <- as.data.frame(predictionmap.3, xy = TRUE)
gg.0<-ggplot() +
  geom_raster(data = prediction.raster, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() +theme(legend.position = "none") +
  coord_sf()
gg.low<-ggplot() +
  geom_raster(data = prediction.raster.2, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() +theme(legend.position = "none") +
  coord_sf()
gg.high<-ggplot() +
  geom_raster(data = prediction.raster.3, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() + theme(legend.position = "none") +
  coord_sf()
grid.arrange(gg.low,gg.0,gg.high,ncol=3)
