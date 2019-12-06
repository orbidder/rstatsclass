# Today we will be working with cluster data. Animal create "clusters" of GPS locations
#   in places that they spend a significant amount of time (relevant to the fix rate
#   of their tracking device). We care about GPS clusters because they can help us track 
#   animals, and in particular, to find out what animals were doing at specific, important
#   locations. We use GPS clusters to find feeding locations, den sites, sleeping sites,
#   and sometimes other fun things like important communication locations.

# Our class today will go over what a cluster is, how to specify it, and how to extract
#   clusters from GPS data. 
# Next, we will predict a certain kind of cluster (kill sites) from all GPS data using 
#   a model that we will fit with only field-investigated clusters
# Finally, we will apply what we learned in lesson 4 (RSFs) to see which habitat 
#   characteristics are selected for at a specific cluster type (the kill sites modeled
#   in the above step)

# Syntax and datetimes
install.packages("tidyverse")
install.packages("lubridate")
# For day/night column and making sure data have a standardized fix rate
install.packages("amt")
# Spatial packages
install.packages("rgeos")
install.packages("raster")
install.packages("rgdal")
install.packages("maptools")
install.packages("sf")
install.packages("mapview")
# Multicolinearity
install.packages("usdm")
install.packages("corrr")
# Modeling
install.packages("lme4")
# Cross-validation (createDataPartition)
install.packages("caret")

# Load your libraries
library(maptools)
library(rgdal)
library(tidyverse)
library(raster)
library(lme4)
library(sf)
library(mapview)
library(lubridate)
library(rgeos)
library(caret)
library(amt)
library(usdm)
library(corrr)

###*** Portion 1 of today's exercise is on how to generate clusters from GPS data
# This is really useful to do BEFORE or DURING fieldwork so you can investigate
#   the clusters on the ground
# (Break to discuss how to design field sampling protocos with clusters)

# Read in the data from file "puma_data_2015.csv" in the "Seminar 5 SSFs" folder. 
# Set up a pipe with the following steps:
# 1. bring in the csv
# 2. select only columns called animals_id, acquisition_time, longitude, and latitude
#   2b. when you select column animals_id, rename it as "ID"
# 3. make a new column called "timestamp" that is in the correct tz of "America/Argentina/San_Juan"
#   IMPORTANT: the acquisition_time column is in UTC, so you don't need to use force_tz
#     to correct it
# 4. sort the data by ID, then timestamp
# 5. In the interest of time, we're only going to look at data collected in the winter,
#     so subset the data to include only June and July dates.
# 6. name the final tibble "puma.data"

# Let's look at your tibbles!
puma.data

# ____________***____________
# Concept break!

# First let's bring in our environmental data to use as a basemap
raster_files <- list.files("Seminar 5 SSFs/",pattern = ".tif") 
envtrasters <- raster::stack(paste0("Seminar 5 SSFs/", raster_files))
names(envtrasters) <- c("dem", "max_ndvi", "slope",  "tri")

# So what is a cluster?
# Let's look at a puma track to get a feel for it

# We'll start by just grabbing about a week of data from one puma for visualization
puma.data %>% 
  filter(ID == 6) %>% 
  slice(., 1:60) -> puma_points
puma.data %>% 
  filter(ID == 6) %>% 
  slice(., 1:60) %>% 
  summarise(do_union = FALSE) %>%
  st_cast(., "LINESTRING") -> puma_lines

# Vizualize your tracks!
plot(crop(envtrasters[[2]],extent(puma_lines)))
plot(st_geometry(puma_lines),add=T)
plot(st_geometry(puma_points),add=T)

# Or you can zoom in on a track in Mapview!
mapview(puma_lines)

# We need to make a track to make sure we eliminate extraneous fixes 
#   if there are multiple fix rates in the data
# We can also figure out if each location is day or night within the track
puma.data %>% 
  # make a track
  mk_track(., longitude, latitude, id = ID, timestamp, crs = CRS("+init=epsg:4326")) %>% 
  # day or night?
  time_of_day() %>%
  # make sure the track is arranged correctly by ID and time
  arrange(id,t_)-> puma_track

# Next we need to make sure the fix rate is consistent. 
# This step is essential if you have a variable fix rate in your data!
# In our full dataset, we have 3 hr and 1 hr fixes
# To compare, let's first see how many data rows we have
nrow(puma_track)
# Then, we resample our data based on the 3 hr fix rate
puma_track %>% group_by(id) %>% nest() %>% 
  mutate(data = map(
    data, ~ .x %>% 
      track_resample(rate = hours(3), tolerance = minutes(15)))) %>% unnest() -> puma_track
# Are there fewer data points now?
nrow(puma_track)

# If you know that the fix rate is consistent then making a track is unnecesary
# In that case, you can assign points to day/night using the maptools package
# e.g.
#puma.data %>%
#  mutate(sunrise = sunriset(SpatialPoints(cbind(longitude,latitude), proj4string=CRS("+init=epsg:4326")), timestamp, direction="sunrise", POSIXct.out=TRUE)[,2],
#         sunset = sunriset(SpatialPoints(cbind(longitude,latitude), proj4string=CRS("+init=epsg:4326")), timestamp, direction="sunset", POSIXct.out=TRUE)[,2],
#         daynight = ifelse(timestamp>sunrise&timestamp<sunset,1,0)) -> puma.data

# But, if we had a variable fix rate and identified day/night 
#   from our track, we could now attach the day/night data to the original dataframe
#   using a right join
# Why do we use a right join??
puma.data %>%
  # We attach the day/night covariate by joining by time and puma id
  right_join(.,puma_track,by=c("timestamp" = "t_", "ID" = "id")) %>%
  # Then we make a new variable for day/night as a binary variable: 1 for day and 0 for night
  mutate(daynight = ifelse(tod_=="night",0,1)) %>%
  # Because the join will add all columns from each tibble (many of which are duplicates)
  #   we will select only the columns we want to keep
  dplyr::select(1:5,10) %>%
  # Finally, we can transform our data to UTMs from latlongs
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>%
  st_transform("+init=epsg:32719") -> puma.data


# Set the time (days) and distance (meters) windows, as well as the fix rate
s_time <- 1.33
s_dist <- 20
fix_r <- 3
npuma <- length(unique(puma.data$ID))

# Specify the time zone the data are in (timezone_1) 
#   and the time zone they should be in (timezone_2)
# Because we've already corrected the timezone in the timestamp column, we use the same
#   timezone for both
timezone_1 <- "America/Argentina/San_Juan"
timezone_2 <- "America/Argentina/San_Juan"

source("Seminar 8 Clusters/cluster_script.R")

# Now run the cluster function!
# This will give you a centroids csv file and a points csv file
# If you only have one animal, you can just run the cluster and points functions outside of the loop
centroids_df <- list()
points_df <- list()
for(j in 1:npuma){
  puma.data %>% 
    filter(ID == unique(puma.data$ID)[j]) -> subsetpuma
  centroids_df[[j]] <- cluster_centroids_func(df = subsetpuma, stime = s_time, sdist = s_dist, 
                                              fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)
  centroids_df[[j]]$puma = unique(puma.data$ID)[j]
  points_df[[j]] <- cluster_points_func(df = subsetpuma, stime = s_time, sdist = s_dist, 
                                        fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)
  points_df[[j]]$puma = unique(puma.data$ID)[j]
}
centroids_df<-data.frame(do.call(rbind,centroids_df))
# Remove rows left over from merging overlapping clusters by subsetting by clusters with duration > 0
centroids_df<-subset(centroids_df,time>0)
points_df<-data.frame(do.call(rbind,points_df))

# Always look at the data to make sure they make sense
head(centroids_df)
subset(points_df,ID==6&puma==1)
puma.data[52:62,]

# Which cluster characteristics would you add to the cluster function?

## Save the data files
#write.csv(centroids_df, "clustercentroids.csv")
#write.csv(points_df, "clusterpoints.csv")

## Save the shape files
#coordinates(centroids_df)<-~x+y
#crs(centroids_df) <- ("+init=epsg:32719")
#writeOGR(centroids_df,".","centroids_df",driver="ESRI Shapefile")

#coordinates(points_df)<-~x+y
#crs(points_df) <- ("+init=epsg:32719")
#writeOGR(points_df,".","points_df",driver="ESRI Shapefile")


#__________********_________

###*** Portion 2 of today's exercise is on how to use data from field-investigated
#         clusters to model the probability of a certain cluster type from 
#         uninvestigated clusters
# For this next exercise, we will predict kills by modeling investigated clusters
#   and applying the kill prediction model to uninvestigated clusters

# Let's bring in our field-collected data
centroids_df <- read_csv("Seminar 8 Clusters/investigatedclusters.csv")%>% 
  st_as_sf(., coords = 3:4, crs = "+init=epsg:32719") %>% 
  mutate(first = lubridate::ymd_hms(first,tz="America/Argentina/San_Juan"),
         last = lubridate::ymd_hms(last,tz="America/Argentina/San_Juan"),
         killYN = as.factor(killYN))

# Check out the spatial distribution of the data to get a better feel for it
plot(centroids_df["killYN"])
plot(centroids_df["nightratio"])
plot(centroids_df["actratio"])
plot(centroids_df["time"],breaks = c(0,0.1,0.3,1,3,9))

# How many of our clusters were investigated, and how many were kills?
nrow(subset(centroids_df,is.na(killYN))) #number of uninvestigated clusters
table(centroids_df$killYN) #number of kills and non-kills among investigated clusters

# In order to fit our model, we only want to use investigated clusters, so let's
#   get rid of rows with "NA" under the killYN column
centroids_df %>%
  drop_na() -> centroids_df_inv


#__________********_________

# ACTIVITY: Now, use what you learned from fitting an RSF model in Seminar 4 to fit a
#   kill prediction model (also a mixed-effects binomial logistic regression model)
#   with cluster characteristics as predictor variables
# Steps:
#   1. assess how to include correlated variables
#   2. write your candidate models
#   3. perform model selection by comparing AIC values

# Let's start by looking at the distribution and correlation among the covariates
#   for the non-kill clusters
centroids_df_inv %>% 
  filter(killYN == 0) %>% 
  dplyr::select(time:points,bin) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales = "free")
# Now edit the above code to look at the distribution of kill clusters 
# What do these distributions tell us? 
# What might we change if we were to do the study again?

# Use the "cor" command to look at correlations between our 6 covariates
#...
#...
#...


# To vizualize these correlations, we can also use package corrr
library(corrr)
centroids_df_inv %>% 
  dplyr::select(time:points,bin) %>% 
  correlate() %>% 
  network_plot()

# You can also use variance inflation factor (VIF) rather than correlation 
#   to look at mutlicollinearity
# If there are variables with values higher than 4, that suite of variables 
#   shouldn't be in the same model
vif(as.data.frame(centroids_df_inv[,c(3:7,10)]))
# What happens when we get rid of one of the variables with a >4 VIF score?
vif(as.data.frame(centroids_df_inv[,c(4:7,10)]))
vif(as.data.frame(centroids_df_inv[,c(3,5:7,10)]))
# What if we only keep one of the variables with a >4 VIF score?
vif(as.data.frame(centroids_df_inv[,c(3,5:6,10)]))
vif(as.data.frame(centroids_df_inv[,c(4:6,10)]))
vif(as.data.frame(centroids_df_inv[,c(5:7,10)]))

# Compare the correlation value of the "bin" covariate with others, and then look
#   at its VIF scores with the suites of covariates
# Becasue "bin" is a binomial covariate, it responds a little differently between approaches
# How do you think you should proceed with this covariate?

# Scale and center your variables
#....
#....
#....


# Now write out your candidate models and select the best model

pumamodel1<-
  
pumamodel2<-
  
pumamodel3<-
  
  
#__________********_________ 

# Time to see how accurate our best model is! 
# For this, we will use bootstrapping cross-validation
  
# First, we set a seed so our results are reproducible
set.seed(6)

# Then, we decide on a number of simulations we want to perform
n.sim = 100

# Next, make a new dataframe for testing the model and 
#   create an empty dataframe to keep track of your results
traintest<-centroids_df_inv
kfold.kills<-data.frame(matrix(vector(), n.sim, 1,
                               dimnames=list(c(), c("accuracy"))),
                        stringsAsFactors=F)

# Now for the cross validation simulation. We will fit the model to some of the data and see how 
#   well the model predicts the remaining data

for (i in 1:n.sim){
  
  #First, split your data into two dataframes, one with 80% of the data for training
  #   and one with 20% for testing
  inTrain <- createDataPartition(y=traintest$killYN, p=0.8, list=FALSE)
  training <- traintest[inTrain,]
  testing <- traintest[-inTrain,]
  
  # ***Write out your best model below to fit using the training data***
  test.model<-glmer()
  
  # Next, predict the testing data with the model fit from the training data
  # To predict binomial data, make sure to include type="response" in your command
  predict1<-predict(test.model,testing,type="response")
  
  # We also need to examine the fit of the training data to the model output 
  PM<-predict(test.model,training,type="response")
  target_pred<-as.matrix(PM)
  target_class<-as.matrix(as.factor(training$killYN))
  
  # To determine the optimal cutoff to distinguish a kill between 0 and 1, 
  #   we calculate the false positive rate (1 - specificity) and 
  #   the true positive rate (sensitivity)
  pred<-prediction(target_pred,target_class)
  perf<-performance(pred,"tpr","fpr")
  fpr = perf@x.values[[1]] 
  tpr = perf@y.values[[1]] 
  sum = tpr + (1-fpr) 
  index = which.max(sum) 
  
  # We use the FPR and TPR to find the optimal cutoff for the test model
  cutoff = perf@alpha.values[[1]][[index]] 
  
  # Then we can assign our model output values to be kills (1) or not (0)
  predict2<-ifelse(predict1>cutoff,1,0)
  
  
  # Finally, we summarize the results by assessing the accuracy of the 
  #   model predictions in the testing dataset
  confM<-confusionMatrix(as.factor(predict2), as.factor(testing$killYN))
  kfold.kills$accuracy[i]<-confM$overall[1]
}

# Calculate the mean overall accuracy from our cross-validation simulations
mean(kfold.kills$accuracy)

# Now let's look at a confusion matrix from our last simulation. 
# What do the Sensitivity (TPR) and Specificity (TNR) tell us about our data?
# Are we more likely to falsy call a real empty cluster a kill, 
#   or a real kill an empty cluster?
confM


#__________********_________

# ACTIVITY: Let's now identify the best cutoff for our full investigated dataset 
#   to apply to our uninvestigated clusters
#   Modify the cutoff code in the cross validation loop to find the best cutoff
#     for all of the investigated cluster data

#...
#...
#...


# Predict model output for uninvestigated clusters
centroids_df$predictkill<-predict(pumamodel1,newdata=centroids_df,
                                  type="response",allow.new.levels=T)

# And apply the optimal cutoff to the predicted values
centroids_df$modkillYN<-ifelse(centroids_df$predictkill>cutoff,1,0)

# Look at the distribution of your modeled data from all generated and 
#   investigated clusters
table(centroids_df$modkillYN)
table(centroids_df$killYN)

# Compare your model output to your investigated cluster data
confusionMatrix(as.factor(centroids_df$killYN), as.factor(centroids_df$modkillYN))

# How might a different cutoff value change your results?



#__________********_________
# ACTIVITY: Use the environmental covariates in the Seminar 5 folder to fit an 
#   RSF with your new predicted kill dataset
# Hint: You can use the code in the Seminar 4 exercise!

