library(maptools)
library(rgdal)
library(tidyverse)
library(raster)
library(lme4)
library(sf)
library(lubridate)
library(rgeos)
library(StreamMetabolism)

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

centroids_df <- read.csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma1_clusters.csv")

centroids_df_inv<-centroids_df[complete.cases(centroids_df), ]
cor(as.data.frame(centroids_df[,c(3:7,10)]))

#try with points/time/night/bin interchanged
pumamodel1<-glmer(killYN~points+actratio+(1|puma),
                  family=binomial(link=logit),data=centroids_df)
