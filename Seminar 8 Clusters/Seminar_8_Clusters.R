library(maptools)
library(rgdal)
library(tidyverse)
library(raster)
library(lme4)
library(sf)
library(lubridate)
library(rgeos)
library(StreamMetabolism)
library(caret)

#We are selecting only Puma ID = 4 today to examine one puma's clusters
puma4 <- read_csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma_data 3hr2.csv") %>% 
  dplyr::select(ID = animals_id, 5:7) %>% 
  filter(ID == 4) %>%
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>% 
  st_transform("+init=epsg:32719") %>% 
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan"))

#In the interest of time, let's only look at data collected in the winter of 2014
puma4 %>%
  filter(lubridate::month(timestamp) %in% 6:9,lubridate::year(timestamp)==2014) -> puma4T

#Determine sunrise and sunset times by day. 
#The lat/longs should be chosen based on a representative location in your study area
#Determine if each location is between sunrise and sunset (day = 1) or not (night = 0)
puma4T %>%
  rowwise() %>% 
  mutate(sunrise = as.POSIXct(sunrise.set(-29.18351, -69.34985, timestamp)[1,1],origin="1970-01-01",format="%H:%M:%S",force_tz("America/Argentina/San_Juan")),
         sunset = as.POSIXct(sunrise.set(-29.18351, -69.34985, timestamp)[1,2],origin="1970-01-01",format="%H:%M:%S",force_tz("America/Argentina/San_Juan")),
         daynight = ifelse(timestamp<sunrise|timestamp>sunset,0,1)) -> puma4T

#Set the time (days) and distance (meters) windows, as well as the fix rate
s_time <- 1.33
s_dist <- 20
fix_r <- 3

#Specify the time zone the data are in (timezone_1) 
# and the time zone they should be in (timezone_2)
timezone_1 <- "America/Los_Angeles"
timezone_2 <- "America/Argentina/San_Juan"


centroids_df <- cluster_centroids_func(df = puma4T, stime = s_time, sdist = s_dist, 
                                        fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)
points_df <- cluster_points_func(df = puma4T, stime = s_time, sdist = s_dist, 
                                  fixr = fix_r, timezone1 = timezone_1, timezone2 = timezone_2)

#Always look at the data to make sure they make sense
head(centroids_df)
puma4T[14:24,]
head(points_df)



#-----------------
# Predict kills from investigated clusters and GPS data

centroids_df <- read.csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma1_clusters.csv")

centroids_df_inv<-centroids_df[complete.cases(centroids_df), ]
cor(as.data.frame(centroids_df[,c(3:7,10)]))

# Now, use what you learned from fitting an RSF model to fit a
#   kill prediction model (also a mixed-effects binomial logistic regression model)
# Assess how to include correlated variables, write your candidate models, 
#   and perform model selection by comparing AIC values

pumamodel1<-

pumamodel2<-
  
pumamodel3<-

  
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
  # Write out your best model below to fit using the training data
  test.model<-glmer(killYN~points+actratio+(1|puma),family=binomial,data=training)
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

# Let's now identify the best cutoff for our full investigated dataset 
#   to apply to our uninvestigated clusters
PM <- fitted(pumamodel1)
target_pred <- as.matrix(PM)
target_class <- as.matrix(as.factor(centroids_df_inv$killYN))
pred <- prediction(target_pred,target_class)
perf <- performance(pred,"tpr","fpr")
fpr <- perf@x.values[[1]] 
tpr <- perf@y.values[[1]] 
sum <- tpr + (1-fpr) 
index <- which.max(sum) 
cutoff <- perf@alpha.values[[1]][[index]] 
cutoff

# Predict model output for uninvestigated clusters
centroids_df$predictkill<-predict(pumamodel1,newdata=centroids_df,
                        type="response",allow.new.levels=T)

# And apply our optimal cutoff to the predicted values
centroids_df$modkillYN<-ifelse(centroids_df$predictkill>cutoff,1,0)
