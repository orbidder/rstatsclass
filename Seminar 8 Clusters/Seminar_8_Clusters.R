library(maptools)
library(rgdal)
library(tidyverse)
library(raster)
library(lme4)
library(sf)
library(lubridate)
library(rgeos)
library(caret)

# Load up the GPS location data from one individual animal
puma.data <- read_csv("Seminar 8 Clusters/puma_data_P4.csv") %>% 
  dplyr::select(ID = animals_id, 5:7) %>% 
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan")) %>%
  arrange(., timestamp)

# figure out if each location is day or night
puma.data %>% 
  mk_track(., longitude, latitude, timestamp, crs = CRS("+init=epsg:4326")) %>% 
  time_of_day() -> puma_track

# Make sure the fix rate is consistent. 
# This is essential if you have a variable fix rate. 
# In this dataset, we have 3 hr, 1 hr, and 30 minute fixes
# To compare, let's first see how many data rows we have
nrow(puma_track)
# Then, we resample our data based on the 3 hr fix rate
puma_track %>% 
  track_resample(rate = hours(3), tolerance = minutes(15)) -> puma_track
# Are there fewer data points now?
nrow(puma_track)

# attach the day/night data to the original dataframe
puma.data %>%
  right_join(.,puma_track,by=c("timestamp" = "t_")) %>%
  mutate(daynight = ifelse(tod_=="night",0,1)) %>%
  dplyr::select(1:5,10) %>%
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>%
  st_transform("+init=epsg:32719") -> puma.data

# In the interest of time, let's only look at data collected in the winter of 2014
puma.data %>%
  filter(lubridate::month(timestamp) %in% 6:9,lubridate::year(timestamp)==2014) -> puma.data

# Set the time (days) and distance (meters) windows, as well as the fix rate
s_time <- 1.33
s_dist <- 20
fix_r <- 3

# Specify the time zone the data are in (timezone_1) 
#   and the time zone they should be in (timezone_2)
timezone_1 <- "America/Los_Angeles"
timezone_2 <- "America/Argentina/San_Juan"

# Now run the cluster function!
centroids_df <- cluster_centroids_func(df = puma.data, stime = s_time, sdist = s_dist, 
                                        fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)

# For field work, you will likely also want the points associated with each cluster
points_df <- cluster_points_func(df = puma.data, stime = s_time, sdist = s_dist, 
                                  fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)

# Always look at the data to make sure they make sense
head(centroids_df)
puma.data[13:25,]
head(points_df,11)

# Which cluster characteristics would you add to the cluster function?

#__________********_________

# For this next exercise, we will predict kills by modeling investigated clusters
#   and applying the model to uninvestigated clusters

centroids_df <- read_csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma_clusters_halfsubset.csv")%>% 
  st_as_sf(., coords = 2:3, crs = "+init=epsg:32719") %>% 
  mutate(first = lubridate::ymd_hms(first,tz="America/Argentina/San_Juan"),
         last = lubridate::ymd_hms(last,tz="America/Argentina/San_Juan"),
         killYN = as.factor(killYN))

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

# ACTIVITY: Now, use what you learned from fitting an RSF model to fit a
#   kill prediction model (also a mixed-effects binomial logistic regression model)
#   with cluster characteristics as predictor variables
# Assess how to include correlated variables, write your candidate models, 
#   and perform model selection by comparing AIC values

# Let's start by looking at the distribution and correlation among the covariates
centroids_df_inv %>% 
  dplyr::select(time:points,bin) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales = "free")

cor(centroids_df_inv[,c(2:6,9)])

# To vizualize these correlations in space, we can use package corrr
library(corrr)
centroids_df_inv %>% 
  dplyr::select(time:points,bin) %>% 
  correlate() %>% 
  network_plot()


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

# Now for the simulation. We will fit the model to some of the data and see how 
#   well the model predicts the remaining data

for (i in 1:n.sim){
  
  #First, split your data into two dataframes, one with 80% of the data for training
  #   and one with 20% for testing
  inTrain <- createDataPartition(y=traintest$killYN, p=0.8, list=FALSE)
  training <- traintest[inTrain,]
  testing <- traintest[-inTrain,]
  
  # ***Write out your best model below to fit using the training data
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




# Predict model output for uninvestigated clusters
centroids_df$predictkill<-predict(pumamodel1,newdata=centroids_df,
                        type="response",allow.new.levels=T)

# And apply the optimal cutoff to the predicted values
centroids_df$modkillYN<-ifelse(centroids_df$predictkill>cutoff,1,0)

# Compare your model output to your investigated cluster data
table(centroids_df$killYN,centroids_df$modkillYN)

# And look at the distribution of your modeled data from all generated clusters
table(centroids_df$modkillYN)


#__________********_________
# ACTIVITY: Use the environmental covariates in the Seminar 8 folder to fit an 
#   RSF and and RSPF with your kill / no kill data



