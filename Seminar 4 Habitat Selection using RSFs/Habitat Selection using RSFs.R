#this is a test line!

#modified Dana Seidel RSF code

#load packages
library(raster)
library(lme4)
library(sf)
library(tidyverse)
library(adehabitatHR)
library(mapview)
library(lubridate)
library(rgeos)
library(rgdal)

# Read in, clean, and project all your data
vicuna <- read_csv("vicuna_data.csv") %>% 
  dplyr::select(ID = animals_id, 4:6) %>% 
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>% 
  st_transform("+init=epsg:32719") #%>% 
  #mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan"))

envtrasters <- stack("/Users/justinesmith/Documents/UCB/Data/Rasters/SGNP_all_stack_south_crop.tif") %>% 
#  projectRaster(crs="+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
names(envtrasters) <- c("aspect", "dem", "rough",  "slope",  "tri", "max_ndvi")

# Visualize environmental layers to make sure they make sense
NDVI.raster <- as.data.frame(envtrasters[[6]], xy = TRUE)
ggplot() +
  geom_raster(data = NDVI.raster, aes(x = x, y = y, fill = max_ndvi)) +
  scale_fill_viridis_c() +
#  geom_sf(data = vicuna) +
  coord_sf()

# Keep in mind temporal effects!
# RSF may be a spatial analysis but if we think selection or availability may change over time, we must control for that in our models
# We can control for time either via appropriate covariate inclusion, or simply by limiting our scope.
# Let's say we only want to know how vicuñas select habitat in the winter (June - September)

vicuna %>% 
  filter(lubridate::month(timestamp) %in% 6:9) %>% 
  dplyr::select(-acquisition_time) %>%
  mutate(Used = 1) -> vicuna

# Choosing an available range. 
#conservative option - 99% UD
#Use below command to draw from range of all zebra locations
vicuna_ud <- adehabitatHR::kernelUD(as(vicuna, "Spatial"), 
                                   grid = 450)
vicuna_hr <- adehabitatHR::getverticeshr(vicuna_ud, 99)
mapview(vicuna_hr)

#Use below command instead to get unique home ranges by adding the id column
vicuna_ud <- adehabitatHR::kernelUD(as(vicuna, "Spatial")[1],  
                                   grid = 450)
vicuna_hr <- adehabitatHR::getverticeshr(vicuna_ud, 99)

# Visualize individual home ranges
homerange <- mapview(vicuna_hr)
homerange

#--------
# Sample available points within each home range. 

# First, determine how many points to sample in each home range. 
# Here we will calculate a 1:1 or 10:1 ratio of available:used points and randomly sample the available points from the population or individual UDs

vicuna_ids <- unique(vicuna$ID)
availables <- list()
for(i in 1:length(vicuna_ids)){
  st_sample(st_as_sf(vicuna_hr)[i,], nrow(filter(vicuna, ID == vicuna_ids[i]))) %>% # e.g. 1:1 ratio; sample individual UDs
    #st_sample(st_as_sf(vicuna_hr), 10*nrow(filter(vicuna, ID == vicuna_ids[i])))  %>%  # e.g. 1:10 ratio; sample population UD
    st_sf(geometry = .) %>%
    mutate(ID = vicuna_ids[i], 
           timestamp = NA, 
           Used = 0) -> availables[[i]] 
}

availables %>% 
  do.call(rbind,.) %>%
  rbind(vicuna, .) -> vicuna_all

#--------

# Scale and center environmental covariates. 
# Scaling covariates allows to to compare the relative strength of each variable by looking at the coefficient values alone
vicuna_all %>% mutate(NDVI = raster::extract(envtrasters2[[6]], as(., "Spatial")),
                     elev = raster::extract(envtrasters2[[2]], as(., "Spatial")),
                     slope = raster::extract(envtrasters2[[4]], as(., "Spatial")),
                     tri = raster::extract(envtrasters2[[5]], as(., "Spatial")),
                     NDVI.scaled = scale(NDVI, center = TRUE, scale = TRUE),
                     elev.scaled = scale(elev, center = TRUE, scale = TRUE),
                     slope.scaled = scale(slope, center = TRUE, scale = TRUE),
                     tri.scaled = scale(tri, center = TRUE, scale = TRUE)) %>% 
  group_by(ID) -> vicuna_full

# Save your scale parameters - you'll need them later!
print(elevscalelist <- list(scale = attr(vicuna_full$elev.scaled, "scaled:scale"),center = attr(vicuna_full$elev.scaled, "scaled:center")))
print(slopescalelist <- list(scale = attr(vicuna_full$slope.scaled, "scaled:scale"),center = attr(vicuna_full$slope.scaled, "scaled:center")))
print(triscalelist <- list(scale = attr(vicuna_full$tri.scaled, "scaled:scale"),center = attr(vicuna_full$tri.scaled, "scaled:center")))
print(NDVIscalelist <- list(scale = attr(vicuna_full$NDVI.scaled, "scaled:scale"),center = attr(vicuna_full$NDVI.scaled, "scaled:center")))


#-----------
# Fitting the RSF

# Always check the structure of the data
head(vicuna_full)
table(vicuna_full$Used) #check for 1:1 ratio
table(vicuna_full$Used,vicuna_full$ID) #check for 1:1 ratio by individual

# First we need to check for colinearity. 
# If two covariates are correlated, including them in the same model can overinflate their standard errors and make it appear that they are not good predictors (when they acutally might be!)  
# General practice is to not include two covariates for which r > 0.7 (or, conservatively, 0.5)
vicuna_cor <- vicuna_full
st_geometry(vicuna_cor) <- NULL
cor(vicuna_cor[,c(4:8)])

vicuna_full$ID <- as_factor(vicuna_full$ID)

# Population-level RSF (fixed effects only)
# Let's compare the fit of slope and ruggedness to decide which one we want to include
fixed.global <- glm(Used ~ NDVI.scaled + elev.scaled + slope.scaled, data=vicuna_full, family=binomial(link="logit"))
summary(fixed.global)
fixed.global <- glm(Used ~ NDVI.scaled + elev.scaled + tri.scaled, data=vicuna_full, family=binomial(link="logit"))
summary(fixed.global)

# Does slope or ruggedness improve model fit more?

# Next, we'll fit a mixed-effects binomial logistic regression
# i.e. an RSF that accounts for variation among individuals

#---------
# BREAK - random effects exercise!

# Let's explore how random effects alter model shape
# To do this, we'll vary the slope and intercept of a logistic regression model and visualize the differences

# First make a function for a binomial logistic regression model
fake.data<-function(x,b,m){
  1/(1+exp(-(b + x*m)))
}

# Set up data
fakem0.5<-rep(0.5,times=101)
fakem1.5<-rep(1.5,times=101)
fakem1<-rep(1,times=101)
fakex<-seq(0,20,by=0.2)
fakebn1<-rep(-1,times=101)
fakebn10<-rep(-10,times=101)
fakebn20<-rep(-20,times=101)

# Random intercept plot
fakeyn1<-fake.data(fakex,fakebn1,fakem1)
fakey0<-fake.data(fakex,fakebn10,fakem1)
fakey05<-fake.data(fakex,fakebn20,fakem1)
plot(fakey0~fakex,ylim=c(0,1))
points(fakey05~fakex,ylim=c(0,1),col="blue")
points(fakeyn1~fakex,ylim=c(0,1),col="red")

# Random slope plot
fakeym1<-fake.data(fakex,fakebn10,fakem1)
fakeym0.5<-fake.data(fakex,fakebn10,fakem0.5)
fakeym1.5<-fake.data(fakex,fakebn10,fakem1.5)
plot(fakeym1~fakex,ylim=c(0,1))
points(fakeym0.5~fakex,ylim=c(0,1),col="blue")
points(fakeym1.5~fakex,ylim=c(0,1),col="red")

# OK - back to the models!
#---------

# Random intercept
# Most common way to include a random effect
# Addresses variation in sample sizes, moves logistic curve left and right
random.int.global <- glmer(Used ~ NDVI.scaled + elev.scaled + tri.scaled + (1|ID), 
                           data=vicuna_full, family=binomial(link="logit"))
summary(random.int.global)
coef(random.int.global)
ranef(random.int.global)

# Random slope
# If you think individuals respond differently to environmental covariates
# Risk in over-fitting
# Example hypothesis: some individuals select for vegetation, whereas others do not, because of variation in internal state and risk-foraging tradeoffs
random.slp_ndvi.global <- glmer(Used ~ NDVI.scaled + elev.scaled + tri.scaled + (0+NDVI.scaled|ID), 
                                data=vicuna_full, family=binomial(link="logit"))
summary(random.slp_ndvi.global)
coef(random.slp_ndvi.global)
ranef(random.slp_ndvi.global)

# Both random slope and intercept
# Example hypothesis: some individuals select for ruggedness, whereas others avoid it, because of variation in boldness
random.int.slp_tri.global <- glmer(Used ~ NDVI.scaled + elev.scaled + tri.scaled + (1+tri.scaled|ID), 
                                  data=vicuna_full, family=binomial(link="logit"))
summary(random.int.slp_tri.global)
coef(random.int.slp_tri.global)
ranef(random.int.slp_tri.global)


#-----------
#Model selection

# Your turn! List your candidate models below based on your alternative hypotheses
# Run your candidate models
model1<-glmer(Used ~ NDVI.scaled + elev.scaled + tri.scaled + (1|ID), 
              data=vicuna_cor, family=binomial(link="logit"))
model2<-
model3<-
model4<-

# Make a list of your models and compare the deltaAIC, relative log likelihook, and akaike weights
cand.models<-c(model1,model2,model3,model4)
aic.func<-function(x){
  AIC(logLik(x))
}
aics<-unlist(map(cand.models,aic.func),use.names=FALSE)

library(qpcR)
akaike.weights(aics)

#----------
# Find optimal cutoff between used and available
library(AUC)
library(ROCR)

# Determine model cutoff between "used" and "available". This will help you interpret the visualizations.
# (We'll revisit ROC/AUC to assess model fit using cross-validation later)
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

#------------
# Predict to map

# Ideally, we would be able to create a spatial predictive surface of habitat selection. 
# To do this, we have to fit our model using the unscaled variables to match the original raster scales
random.int.raw <- glmer(Used ~ NDVI + elev + tri + (1|ID), 
                          data=vicuna_full, family=binomial(link="logit"))

# You'll note that the model doesn't converge because of the very different scales of the covariates
# When this happens, you have to instead convert the rasters to the scale of the centered and normalized covariates in the model for mapping
# To do this, we'll revisit our scale parameters for our covariates

# Make a new raster stack with just the covariates in the model
env.rasters <- stack(envtrasters[[6]], envtrasters[[2]], envtrasters[[5]])

# Scale the covariates to be centered and normalized, based on the scale parameters from our model covariates
env.rasters[[1]]<-(env.rasters[[1]]-NDVIscalelist$center)/NDVIscalelist$scale
env.rasters[[2]]<-(env.rasters[[2]]-elevscalelist$center)/elevscalelist$scale
env.rasters[[3]]<-(env.rasters[[3]]-triscalelist$center)/triscalelist$scale

# Name the raster layers to match the model covariate names
names(env.rasters) <- c("NDVI.scaled", "elev.scaled", "tri.scaled")

# Time to make a predicted raster layer based on your model!

# The predict function requires that all the parameters in the model are specified. That means we need to give it a value for "ID"
# We want to pick the individual most representative of the population. Therefore, we want the individual that deviates the least from the global intercept coefficient estimate
# To select this individual, look at the random effects and identify the individual that has a deviation from the intercept closest to zero
ranef(random.int)
predictionmap<-raster::predict(object=env.rasters,model=random.int,
                               const=(data.frame(ID="13")),type="response")

# To visualize the differences among individuals, let's also map the individuals with the greatest deviations from the model intercept
predictionmap.2<-raster::predict(object=env.rasters,model=random.int,
                               const=(data.frame(ID="32")),type="response")
predictionmap.3<-raster::predict(object=env.rasters,model=random.int,
                               const=(data.frame(ID="25")),type="response")

# Three types of visualizations!

# Visualization 1 in base::plot
plot(predictionmap)

# Visualization 2 in sp::spplot
library(RColorBrewer)
library(sp)
breaks.plot<-c(seq(0,1,by=0.1))
my.palette <- rev(brewer.pal(n = 10, name = "RdYlBu"))
spplot(stack(predictionmap.2,predictionmap,predictionmap.3),col.regions=my.palette,cuts=9,at=breaks.plot)

# Visualization 3 in ggplot2::ggplot
library(gridExtra)
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
