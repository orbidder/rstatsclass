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
library(sf)
library(mapview)
library(ggspatial)
library(plyr)

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
############################################################################################################
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
#----------------------------------------------------------------------------------------------------------
#Answer 1:
elk_data <- elk_data[elk_data$Latitude > 20,] #keep only points above 20 degrees Latitude
#or
elk_data <- subset(elk_data, elk_data$Latitude > 20)
#----------------------------------------------------------------------------------------------------------

#Take a look at the data now it has been corrected...
plot(elk_data, axes = T)

#We can also look at the bounding box that contains all of our points
elk_data@bbox

#Many functions require our coordinates in meters, so we can reproject using spTransform()
elk_utm <- spTransform(elk_data, CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(elk_utm, axes = T) #Notice now axes show metres!

#Lines#
#Lines are similar to points, but now each feature is represented by a vector of coordinates (c(x,x,x),c(y,y,y))
#common examples we might encounter in movement ecology are roads and rivers.

#Lines can be loaded and set up manually as above with sp, but given a lot of this information will be provided to you
#in ESRI shapefiles, we will load an example using the readOGR function in rgdal. You should note that a
#lot of spatial reading and writing is done using rgdal, which is the r wrapper package for the GDAL
#library.

rivers <- readOGR('wiggins_rivers.shp') #note this creates the SpatialLinesDataFrame, an sp class object
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
plot(pub_lands, axes = T, lwd = 2) ; plot(rivers, col = 'blue', axes=T, add=T) #plot our rivers and public lands

#A quick note on projections...
#Sometimes you'll be given a CSV file by someone and you won't know what projection it's in. In cases like
#these it helps to know some common projections for your study site. You can check if you have the projection
#right by plotting data over a world map in the same projection.

world <- readOGR('Countries_WGS84.shp')

plot(world); plot(elk_data, add = T, col = 'red') #data should appear as a red cross in WY

proj4string(world); proj4string(elk_data) #are these the same?

#if data don't appear in the right place, the projection might be incorrect...

#Exercise 2: Load the csv file 'elk_other_proj.csv', take a look at the coordinates in columns X and Y, 
#what do you notice about them? Can you make them a spatial object, which projection do you need to assign
#in order to get them to appear in the right place on the world map? Hint, you can either assign the projection
#and convert to WGS84 to match the map, or reproject the map using spTransform(). What projections are
#appropriate for this part of North America? (10 mins)

#----------------------------------------------------------------------------------------------------------
#Answer 2:
elk_unk <- read.csv('elk_other_proj.csv')
coordinates(elk_unk) <- ~ X + Y
proj4string(elk_unk) <- CRS('+init=epsg:5070')
world <- spTransform(world, CRS(proj4string(elk_unk)))
plot(world, axes = T); plot(elk_unk, col = 'blue', add=T)
#----------------------------------------------------------------------------------------------------------

############################################################################################################
#Common spatial functions#
#So far we've seen how rgdal helps get spatial data in to R and sp provides methods and classes to handle 
#the data once its there...

#A third package, rgeos, provides some basic functions to manipulate the data. These three packages
#form the core of the sp "R-Spatial" package ecosystem.

#The rgeos package contains a number of useful spatial functions for dealing with sp objects...

#gDistance measures the distance between two features

gDistance(elk_utm[1,], pub_lands[pub_lands$SURFACE == "State",]) #finds distance between point #1 and nearest state land

#gBuffer makes a buffer around an sp object, output is a polygon
elk_500m <- gBuffer(elk_utm[1:5,], width = 500) #just do the first 5 points because it's quick
plot(elk_500m, lwd = 2, col = 'green'); plot(elk_utm[1:5,], col = 'blue', add = T)

#gIntersects returns true or false whether two features touch, which is useful for checking if elk are
#on public land
gIntersects(elk_utm[1,], pub_lands, byid = T) #include argument byid=T to check each polygon
pub_lands$SURFACE[gIntersects(elk_utm[1,], pub_lands, byid = T)] #returns the value of SURFACE at elk_utm[1,]

#Exercise 3:
#A common task is to measure the distance of each GPS location to a linear feature like a river. Running
#gDistance(elk_utm, rivers, byid=T) would return an nrow(x) by nrow(y) matrix, and we aren't too bothered
#about comparing each location to each point on every river (takes a long time on large layers!). 

#Write a for loop that fills a new column, 'dist_river', in the elk_utm attribute table, 
#with the distance between each location and the nearest point on the rivers object (10 minutes)

#I've included a template for loop for you to fill in below, which will monitor progress using a
#progress bar, and also time the calculations for us to compare our solutions...

#For Loop Template#
dechrau <- Sys.time()
total <- #enter an appropriate max number here
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (i in 1:total) {
 
  #Enter your code here#
  
  setTxtProgressBar(pb, i)
}
close(pb)
diwedd <- Sys.time()
amser <- diwedd-dechrau
print(paste("Measurements took",round(diwedd-dechrau,2),attributes(amser)$units,sep = " "))

#----------------------------------------------------------------------------------------------------------

#Answer 3:
dechrau <- Sys.time()
total <- length(elk_utm )
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (i in 1:total) {
  
  #some code to speed up calcs when dealing with very big lines and polygons!
  #buffer <- gBuffer(elk_utm[i,], width = 10000)
  #plot(buffer)
  #cut_river <- rivers[c(gIntersects(rivers, buffer, byid=T)),]
  #plot(buffer);plot(rivers, add=T)
  #plot(buffer);plot(cut_river, add=T)
  
  elk_utm$dist_river[i] <- gDistance(rivers, elk_utm[i,])
  setTxtProgressBar(pb, i)
}
close(pb)
diwedd <- Sys.time()
amser <- diwedd-dechrau
print(paste("Distance measurement took",round(diwedd-dechrau,2),attributes(amser)$units,sep = " "))

head(elk_utm$dist_river)

#----------------------------------------------------------------------------------------------------------

############################################################################################################
#SF#
#As mentioned above, the sf package has many of sp's functionality, but has a couple of important advantages,
#namely, that it is faster, works with tidyverse packages, and can combine different geometry types in 
#the same object

#reading in csv files with sf is similar to sp
elk_data <- read.csv("Spatial_Wiggins_Elk.csv", header = T, stringsAsFactors = F)
elk_data$Date_Time <- as.POSIXct(elk_data$Date_Time, "%Y-%m-%d %H:%M:%S", tz = "MST")
elk_data <- elk_data[elk_data$Latitude > 20,] #don't forget we corrected that outlier!
elk_sf <- st_as_sf(elk_data, 
                     coords = c('Longitude', 'Latitude'), 
                     crs = '+init=epsg:4326')

#Converting an sp object to an sf object is easy!
elk_sf <- st_as_sf(elk_utm)

#As is loading in ESRI Shapefiles
rivers_sf <- st_read('wiggins_rivers.shp')
pub_lands_sf <- st_read('wiggins_lands.shp')


#Plotting with sf shows a summary for each layer
plot(elk_sf)

#...and we can plot specific layers by name
plot(elk_sf['dist_river'])

#and use ggplot2 to make nice maps quite quickly...
ggplot() + geom_sf(data = elk_sf, aes(color = dist_river)) + 
  geom_sf(data = rivers_sf) +
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  theme_light()

#...and even fancier interactive plots with mapview!
map <- mapview(pub_lands_sf['SURFACE'],
               layer.name = 'Land Ownership')
map

#note: you can add multiple layers to mapview by entering them as a list(), and including a vector of layer names

#and of course sf has some useful spatial functions. For instance joins can take attributes from one
#object and add them to another.

pub_lands_sf$id <- as.factor(1:312) #make a unique index for each polygon
elk_sf <- st_join(elk_sf, pub_lands_sf, join = st_intersects, left = T) #adds all pub_lands_sf attributes
plot(elk_sf['SURFACE'])

#and we can count how many points are in each polygon
point_in_poly <- st_join(elk_sf, pub_lands_sf, join = st_within)
point_poly_count <- count(point_in_poly$id)
tract_point_sf <- left_join(pub_lands_sf, point_poly_count, by = c("id" = "x"))
map <- mapview(tract_point_sf['freq'],
               layer.name = '# locations')
map


#sf also has its own distance function, st_distance(), you may see some tutorials with the following method
#Dont run...
#dist <- st_distance(y=elk_sf, x=rivers_sf) # a matrix of all pair-wise distances (all pnts to all lines)
#elk_sf$dist_river_sf <- apply(dist, 2, min) #run through, finding the minimum for each one

#...but I find the for loop much faster. Let's compare with the earlier sp version by timing it, and we'll
#make a new column so we can check the values between the methods
dechrau <- Sys.time() #start the timer
for (i in 1:length(elk_sf)) {
  elk_sf$dist_river_sf[i] <- min(st_distance(y=elk_sf[i,], x=rivers_sf, by_element = T))
}
diwedd <- Sys.time() #stop the timer
amser <- diwedd-dechrau #calc the time elapsed
print(paste("Measurements took",round(diwedd-dechrau,2),attributes(amser)$units,sep = " "))

head(elk_sf) #values are the same, but its much faster (I hope!)

############################################################################################################
#Rasters#
#So far we have concentrated on vector data; points, lines and polygons. A fourth and very important data
#type is the raster. These are essentially 2D matrices that represent a grid of the earth's surface. For
#satellite imagery, each cell represents a pixel from the image. A common use of the raster datatype is
#elevation data...

#Loading Digital Elevation Model (DEM) for our area...
elev <- raster('wiggins_dem.tif')
plot(elev) ; plot(elk_sf[2], add=T) #plot elevation and overlay the elk locations

#hist() will give you a summary histogram of values in the raster
hist(elev)

#we can get some info about our raster easily...
extent(elev)
res(elev)
crs(elev)
elev #calling elev will give us a lot of summary info at once

#The raster package is one of the most comprehensive packages available for R. The documentation is excellent
#and I encourage you to read it to see the breadth of things you can do.

#For instance, raster contains useful terrain functions for calculating aspect, slope etc.
asp <- terrain(elev, opt = "aspect")
slo <- terrain(elev, opt = "slope")

plot(asp)
plot(slo)

#These different layers can be combined in to a 'Raster Stack'
dem <- stack(elev,asp,slo)
plot(dem)

#Individual layers in a stack can be called using the $ operator
dem$wiggins_dem

#think of a stack as a 3D matrix, where each layer or slice represents a new dimension of the same
#area in space (in our example, elevation, aspect and slope). For multi-spectral imagery, the layers
#often correspond with the different spectral bands (e.g. IR, visible light bands). When measuring things 
#like NDVI, we might stack revisits by the satellite and index them by time, allowing us to look at trends 
#over the seasons (i.e. phenology) or to look at conditions on the ground in time and space for use in our 
#movement models (see session on SSFs).

#Some restrictions on stacks are; they must have the same extent, resolution and projection.

#Notice that our rivers and public lands objects are smaller than the dem raster...
st_bbox(pub_lands_sf)
st_bbox(elev)
plot(elev) ; plot(pub_lands_sf[1], add=T) #plot elevation and overlay the elk locations

#But we can use the crop function to cut the object to the extent of another
sm_dem <- crop(dem, pub_lands_sf)
plot(sm_dem$wiggins_dem, axes=T) ; plot(pub_lands_sf[1], add=T) #plot elevation and overlay the elk locations

#or even mask the raster to the outline of the polygon
mask_dem <- mask(dem, pub_lands_sf)
plot(mask_dem$wiggins_dem, axes=T) #There are no polygons in the middle (private land), so the data are NA


#For movement modelling, you might want to sample the raster values at each location, which can be easily done
#using raster's extract() function

elk_sf$id = seq_len(nrow(elk_sf)) #extract needs sf objects to have a unique index called id
elk_sf$elevation <- extract(sm_dem$wiggins_dem, elk_sf)

plot(elk_sf['elevation'])

#and remember the data in sf objects can be treated just like regular data frames...
plot(elk_sf$Date_Time, elk_sf$elevation, axes = T) #annual cycle of elk elevation


#Exercise 4:
#Using the functions in sf and raster, sample the mean elevation in each polygon and add it to the 
#pub_lands_sf attribute table as 'mean_elev (10 mins)

#HINT: This should be done in two steps. Assign the output of the first step to a temp object, 
#what is the output type and structure? Look at the sapply function docs, how can we use it to get 
#the mean for each polyon id back in to our polygon sf object?

#----------------------------------------------------------------------------------------------------------
#Answer 4:
lands_dem <- extract(sm_dem$wiggins_dem, pub_lands_sf)
pub_lands_sf$mean_elev <- sapply(lands_dem, mean)
#----------------------------------------------------------------------------------------------------------

#This final bit of code makes an aesthetically pleasing map with contour lines and hillshade, feel free to
#use it in your own work...

# create hillshade
hs = hillShade(slope = dem$slope, aspect = dem$aspect)
plot(hs, col = gray(0:100 / 100), legend = FALSE)
# overlay with DEM
plot(dem$wiggins_dem, col = terrain.colors(25), alpha = 0.5, legend = FALSE, add = TRUE)
# add contour lines
contour(dem$wiggins_dem, col = "white", add = TRUE, nlevels = 5)

#or if you ever need to make contours to plot in ggplot...
cl <- rasterToContour(dem$wiggins_dem, nlevels = 5)
cl <- st_as_sf(cl)

#Note rasters have to be turned in to data frames to work with ggplot
dem_df <- as.data.frame(dem$wiggins_dem, xy=T)
hs <- as.data.frame(hs, xy=T)
str(dem_df)

ggplot() +
  geom_raster(data = dem_df , aes(x = x, y = y, fill = wiggins_dem)) +
  geom_raster(data = hs , aes(x = x, y = y, alpha = layer)) +
  geom_sf(data = cl, color = 'white') +
  scale_alpha(range = c(0.15, 0.65), guide = "none") +  
  scale_fill_viridis_c() +
  coord_sf()

#Finally, our maps can be saved using the function ggsave()
ggsave(filename = 'wiggins_dem_with_hill.png', plot = last_plot())

#Exercise 5:
#Using everything you've learned in this session, pick a new terrain metric (look at the docs) and 
#sample it for each elk location. Then produce a map, including a dem basemap, hillshade, rivers and
#scale bar and north arrows. Make it as aesthetically pleasing as possible using scale_fill_viridis_c() or
#terrain.colors(). (15 mins)

#----------------------------------------------------------------------------------------------------------
#Answer 5:
tri <- terrain(dem$wiggins_dem, opt = "TRI")
elk_sf$tri <- extract(tri, elk_sf)

ggplot() +
  #geom_sf(aes(color = tri)) +
  geom_raster(data = dem_df , aes(x = x, y = y, fill = wiggins_dem)) +
  geom_raster(data = hs , aes(x = x, y = y, alpha = layer)) +
  scale_alpha(range = c(0.15, 0.65), guide = "none") +  
  scale_fill_viridis_c(direction = -1) +
  #scale_fill_distiller(palette = "Spectral") + 
  #scale_fill_brewer(palette = "Greens") +
  coord_sf() + 
  #annotation_scale(location = "br", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       style = north_arrow_fancy_orienteering) +
  theme_light()
#----------------------------------------------------------------------------------------------------------