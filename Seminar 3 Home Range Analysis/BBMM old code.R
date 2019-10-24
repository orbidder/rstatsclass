##Making a SpatialPolygonsDataFrame from a list of kernelbb objects
library(adehabitatHR)
library(sf)
library(plyr)
library(move)
library(ggplot2)
library(imager)
library(tidyverse)

#This code is a carry over from an earlier version of the seminar, but may be useful if you ever find
# yourself having to deal with a list of bbmm objects.


#Load the data
puma <- read.csv("puma_data1.csv", header = T, stringsAsFactors = F)
puma$acquisition_time <- as.POSIXct(puma$acquisition_time, "%Y-%m-%d %H:%S:%S", tz = "America/Buenos_Aires")

puma <- as.data.frame(puma)

puma <- puma[complete.cases(puma$acquisition_time),] #have to remove NAs because Daylight Savings sometimes produces them when clocks roll forward

puma_ltraj <- as.ltraj(xy = data.frame("x" = puma$x, "y"=puma$y),
                       date = as.POSIXct(puma$acquisition_time), id = puma$animals_id,
                       burst = puma$animals_id, 
                       proj4string = CRS('+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))


lik <- liker(puma_ltraj, sig2 = 60, rangesig1 = c(1, 40))

bbmm_puma <- list()
for (i in 1:3) {
  print(paste('Calculating BBMM for animal: ',ld(puma_ltraj[i])$id[1], sep = ""))
  bbmm_puma[[i]] <- kernelbb(puma_ltraj[i], sig1 = lik[[i]]$sig1, sig2 = 60, grid = 50)#calculate the bbmm
  bbmm_puma[[i]]$id <- ld(puma_ltraj[i])$id[1]
}

bbmm_contours <- list()
for (i in 1:3) {
  #we can get the contours just as we did with the KUD, but only one at a time...
  bbmm_contours[[i]] <- getverticeshr(bbmm_puma[[i]], 95) 
  #To get the IDs back, we have to paste them in to the appropriate slots...
  slot(slot(bbmm_contours[[i]], "polygons")[[1]], "ID") <- as.character(bbmm_puma[[i]]$id[1])
}

#This code will combine a list of spatialpolygons, in to a spatialpolygonsdataframe
#Getting polygon IDs
IDs <- sapply(bbmm_contours, function(x)
  slot(slot(x, "polygons")[[1]], "ID"))

#Checking
length(unique(IDs)) == length(bbmm_contours)

#Making SpatialPolygons from list of polygons
bbmm_contours <- SpatialPolygons(lapply(bbmm_contours,
                                        function(x) slot(x, "polygons")[[1]]))

plot(bbmm_contours, lwd=2)