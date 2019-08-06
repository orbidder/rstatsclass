#Spatial Data in R#
#O. R. Bidder#
#August 2019#
#orbidder@berkeley.edu#

library(sp)
library(tidyverse)
library(rgdal)
library(raster)
library(adehabitatLT)
library(rgeos)

#R has a rich ecosystem of packages to perform spatial analysis. In this session, we will run quickly over
#the different spatial data types, how they are handled in R, and some simple functions that we use
#often in movement ecology analysis. This is not an exhaustive list, so I encourage you to look at the 
#documentation for these packages, and check out the following free ebooks...#

#Geocomputation with R; Lovelace, Nowosad and Muenchow. 
#https://geocompr.robinlovelace.net#

#Applied Spatial Data Analysis with R; Bivand, Pebesma, Gomez-Rubio. 
#http://gis.humboldt.edu/OLM/r/Spatial%20Analysis%20With%20R.pdf#

#There are currently two major spatial R packages that you should become accustomed with, sp and sf.
#sp is older, but its data structures are still the primary input for a lot of packages used in movement
#ecology (e.g. adehabitat). sf is the tidyverse update on sp, includes some great functionality, and works
#nicely with other tidyverse tools (e.g. ggplot2). We will start with sp...

#SP#
#sp is a package that contains methods and classes for dealing with spatial data, primarily vector data
#Of these vector data, there are three main types we will concentrate on today, points, lines and polygons.
#Those of you that have used GIS should be quite familiar with these. Let's look at these individually;

#Points#
#Points are features with just a single set of coordinates (x,y). An example of points we often run in to
#in movement ecology are GPS telemetry data. These can be loaded in to R using the following code...

#load the data as a dataframe from csv
elk_data <- read.csv("Spatial_Wiggins_Elk.csv", header = T, stringsAsFactors = F)
head(elk_data) #take a quick look at whats in the df

#set time to POSIXct
elk_data$Date_Time <- as.POSIXct(elk_data$Date_Time, "%Y-%m-%d %H:%M:%S", tz = "MST")

#set the coordinates
coordinates(elk_data) <- ~ Longitude + Latitude #always in the format x + y

#Once we set the coordinates, notice that R has changed the object type to sp's 'SpatialPointsDataFrame'

#set the CRS
proj4string(elk_data) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #set to WGS84
#proj4string(elk) <- CRS('+init=epsg:4326') #sometimes you will see folks use the epsg code instead

#We can then use the plot function as normal, and sp will act as a wrapper function
plot(elk_data, axes = T) #axes show lon/lat

#Notice that the elk data contain an outlier (usually produced because of a collar error)
#Luckily, we can manipulate sp's dataframes the same way we would any other data frame in R

#Exercise 1: Write some code that would remove the outlier, but leave valid data intact (5 mins)
#Answer 1:
elk_data <- elk_data[elk_data$Latitude > 20,] #keep only points above 20 degrees Latitude
#or
elk_data <- subset(elk_data, elk_data$Latitude > 20)

#Take a look at the data now it has been corrected...
plot(elk_data, axes = T)

#We can also look at the bounding box that contains all of our points
elk_data@bbox

#Many functions require our coordinates in m, so we can reproject using spTransform()
elk_utm <- spTransform(elk_data, CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(elk_utm, axes = T) #Notice now axes show metres!

#Lines#
#lines are similar to points, but now each feature is represented by a vector of coordinates (c(x,x,x),c(y,y,y))
#common examples we might encounter in movement ecology are roads and rivers.
#Lines can be loaded and set up manually as above, but given a lot of this information will be provided to you
#in ESRI shapefiles, we will load an example using the readOGR function in rgdal

rivers <- readOGR('wiggins_rivers.shp')
proj4string(rivers)

#We can access the attribute table for this layer by looking in the data slot
rivers@data

#We can plot the layers with the elk locations, by using the add=T argument 
plot(elk_utm, col = 'black', axes=T) ; plot(rivers, col = 'blue', add=T) #plot our elk data with the rivers

#Polygons#
#polygons are the next increment in vector data. Like lines, they have a vector of coordinates, but now
#their first and last coordinates are always the same, so they end up where they started, and discribe a 
#discrete region or area. Common examples of features you might see as polygons are country borders, 
#administrative boundaries etc. 

pub_lands <- readOGR('wiggins_lands.shp') #load in an example polygon data frame, parcels of public land
proj4string(pub_lands) #check the projection

#we can check that all of our projections are the same using
proj4string(elk_utm) == proj4string(pub_lands)

plot(pub_lands)
plot(pub_lands, axes = T, lwd = 2) ; plot(rivers, col = 'blue', axes=T, add=T) #plot our elk data with the lands

#Common spatial functions#

#When modelling animal movement, you might want to calculate displacement distances and turn angles
#between successive relocations. This can be done quickly using as.ltraj from adehabitatLT

elk_ltraj <- as.ltraj(xy = data.frame(elk_utm$Longitude,elk_utm$Latitude),
                       id = elk_utm$Elk,
                       date = elk_utm$Date_Time, burst = elk_utm$id, 
                       infolocs = data.frame(pkey = paste(elk_utm$id, elk_utm$Date_Time, sep = "."), 
                                             row.names = row.names(data.frame(elk_utm$Longitude,elk_utm$Latitude))),
                       proj4string = CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
elk_df <- ld(elk_ltraj) #creates a dataframe from an ltraj object

colnames(elk_df) #column "dist" holds the displacement distances and "rel.angle" the turn angles in radians

#The rgeos package contains a number of useful spatial functions for dealing with sp objects...
#gDistance measures the distance between two features

gDistance(elk_utm[1,], pub_lands[pub_lands$SURFACE == "State",]) #finds the nearest parcel of state land to first location

#gBuffer makes a buffer around an sp object, output is a polygon
elk_500m <- gBuffer(elk_utm[1:5,], width = 500)
plot(elk_500m); plot(elk_utm[1:5,], col = 'blue', add = T)

#gIntersects returns true or false whether two features touch, which is useful for checking if elk are
#on public land
gIntersects(elk_utm[1,], pub_lands, byid = T) #include argument byid=T to check each polygon
pub_lands$SURFACE[gIntersects(elk_utm[1,], pub_lands, byid = T)] #returns the value of SURFACE at elk_utm[1,]

#Another common task is to measure the distance of each GPS location to a linear feature like a river, this can be
#done using the gDistance function in rgeos...

dechrau <- Sys.time()
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (i in 1:length(elk_utm)) {
  elk_utm$dist_river[i] <- gDistance(elk_utm[i,], rivers)
  setTxtProgressBar(pb, i)
}
close(pb)
diwedd <- Sys.time()
period <- round(diwedd-dechrau,2)
print(paste("Distance measurement took",round(diwedd-dechrau,2),attributes(period)$units,sep = " "))



#SF#
#Mixed Objects#

#Rasters#
#Loading rasters, raster stacks#
#Plotting#
#Extracting values#
#Raster calculations#