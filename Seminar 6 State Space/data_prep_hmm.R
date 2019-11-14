##Data Preg for HMM Seminar
#O. R. Bidder, Sept 2019#
library(adehabitatLT)
library(raster)

URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/",
              "10255/move.373/Elliptical%20Time-Density%20Model%20%28Wall%",
              "20et%20al.%202014%29%20African%20Elephant%20Dataset%20%",
              "28Source-Save%20the%20Elephants%29.csv")
rawData <- read.csv(url(URL))

write.csv(rawData, "elephants.csv")

wigs <- read.csv('C:/Data/Targhee Wiggins Fork Elk/wiggins_fork/Wiggins Fork Apr _19/All_Wiggins.csv')
unique(wigs$Elk)
wigs <- wigs[wigs$Elk == 'Wiggins1',]
wigs <- wigs[wigs$Year == 2016,]
wigs <- wigs[wigs$Latitude > 30,]
wigs <- wigs[,-1]
head(wigs)

coordinates(wigs) <- c("Longitude","Latitude")
proj4string(wigs) <- CRS('+init=epsg:4326')
wigs <- spTransform(wigs, CRSobj = CRS('+init=epsg:32612'))
wigs <- as.data.frame(wigs)
colnames(wigs)[10] <- "Longitude"
colnames(wigs)[11] <- "Latitude"

colnames(wigs)[13] <- "UTM_X"
colnames(wigs)[14] <- "UTM_Y"


write.csv(wigs, file = 'Wiggins1.csv')

#Prep for the grey seal example
seal <- read.csv('mcconnell_seal.csv')

head(seal)
unique(seal$individual.local.identifier)

coordinates(seal) <- c("long","lat")
proj4string(seal) <- CRS('+init=epsg:4326')
seal <- spTransform(seal, CRSobj = CRS('+init=epsg:32630'))

seal <- as.data.frame(seal)

seal <- as.ltraj(xy = data.frame("x" = seal$long, "y"= seal$lat), 
                   date = as.POSIXct(seal$timestamp), 
                   id = seal$individual.local.identifier, 
                   burst = seal$individual.local.identifier, 
                   proj4string = CRS('+init=epsg:32630'))

plot(seal)
seal <- ld(seal)
hist(seal$dist, breaks = 100)
hist(seal$dist[seal$dist < 1000], breaks = 100)
head(seal)

bedrock <- raster('etopo1_bedrock.tif')
bedrock <- projectRaster(bedrock, crs = '+init=epsg:32630')
proj4string(bedrock)

head(seal)
coordinates(seal) <- c("x","y")
proj4string(seal) <- CRS('+init=epsg:32630')

plot(bedrock); plot(seal, add=T)

seal$depth <- raster::extract(bedrock, seal)

plot(seal$depth, seal$dist)
seal <- as.data.frame(seal)

orig_seal <- read.csv('mcconnell_seal.csv')
head(seal)
seal <- seal[,c('id','date','x','y','dt','dist','rel.angle','depth')]

seal <- cbind(seal, orig_seal[,c('long','lat')])
head(seal)
length(seal$id) #2535
seal <- seal[seal$depth < 0,]

plot(seal$depth, seal$dist, ylim = c(0,8000))

chloro <- raster('chloroa.nc')
plot(chloro)
chloro <- projectRaster(chloro, crs = '+init=epsg:32630')

head(seal)
coordinates(seal) <- c("x","y")
proj4string(seal) <- CRS('+init=epsg:32630')

plot(chloro, xlim = c(487921, 962241.9), ylim = c(6092524, 6529550)) ; plot(seal, add=T)
max(seal$y)
min(seal$y)
seal$chl <- raster::extract(chloro, seal)

plot(seal$depth, seal$chl)
plot(seal$dist, (seal$depth*seal$chl))

seal$sc_depth <- scale(seal$depth)
seal$sc_chl <- scale(seal$chl)

plot(seal$dist, (seal$sc_depth*seal$sc_chl))
seal$dep_chl <- (seal$sc_depth*seal$sc_chl)

plot(seal$x, seal$y, col = seal$dep_chl)

library(ggplot2)
library(sf)

seal_sf <- st_as_sf(seal)

ggplot(data = seal_sf) +
  geom_sf(aes(color = (sc_depth + sc_chl)))

ggplot(data = seal_sf) +
  geom_sf(aes(color = (sc_depth)))

ggplot(data = seal_sf) +
  geom_sf(aes(color = (sc_chl)))

ggplot(data = seal_sf) +
  geom_sf(aes(color = (norm_depth)))

ggplot(data = seal_sf) +
  geom_sf(aes(color = (norm_chl)))

seal_sf$norm_depth <- (seal$sc_depth - min(seal$sc_depth))/(max(seal$sc_depth) - min(seal$sc_depth))
seal_sf$norm_chl <- (seal$sc_chl - min(seal$sc_chl))/(max(seal$sc_chl) - min(seal$sc_chl))

seal_sf$dt <- seal_sf$dt/3600

#This isn't working because I'm not accounting for the step length!!! Add it in here!!
seal_sf$speed <- seal_sf$dist/seal_sf$dt

seal_sf$dive_mu <- seal_sf$speed*-25 + seal_sf$norm_depth*50
ggplot(data = seal_sf) +
  geom_sf(aes(color = (dive_mu)))

seal_sf$n_dives <- rpois(length(seal_sf$depth), seal_sf$dive_mu)
seal_sf$n_dives[is.na(seal_sf$n_dives)] <- 0

ggplot(data = seal_sf) +
  geom_sf(aes(color = (n_dives)))

seal_out <- as.data.frame(seal_sf)

seal <- as.data.frame(seal)
colnames(seal)
colnames(seal_out)

seal_data <- seal[,c('id', 'date', 'x', 'y', 'long', 'lat', 'depth', 'chl')]
seal_data <- cbind(seal_data,seal_out[,'n_dives'])
colnames(seal_data)[9] <- 'n_dives'
colnames(seal_data)[10] <- 'dive'

seal_data$depth <- abs(seal_data$depth)
seal_data <- seal_data[,-1]
write.csv(seal_data, 'seal_data.csv')


obsData <- miExample$obsData

# error ellipse model
err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)

# Fit crwMLE models to obsData and predict locations 
# at default intervals for both individuals
crwOut <- crawlWrap(obsData=seal_ltraj, timeStep="hour")

# create data frame with fake data stream
data <- data.frame(ID=rep(factor(c(1,2)),times=c(753,652)),
                   time=c(1:753,1:652),
                   fake=rpois(753+652,5))

# merge fake data stream with crwOut
crwOut <- crawlMerge(crwOut,data,"time")

seal_ltraj <- as.ltraj(xy = data.frame('x' = seal_data$x[1:1800], 'y' = seal_data$y[1:1800]),
                       date = seal_data$time[1:1800], id = seal_data$ID[1:1800], 
                       burst = seal_data$ID[1:1800], 
                       infolocs = data.frame('ID' = seal_data$ID[1:1800],
                                             'time' = seal_data$time[1:1800],
                                             'dive' = seal_data$dive[1:1800],
                                             'depth' = seal_data$depth[1:1800]))
plot(seal_ltraj)
seal_ltraj <- ld(seal_ltraj)
plot(seal_ltraj$dist)

colnames(seal_ltraj)[11] <- "ID"
colnames(seal_ltraj)[3] <- "time"

help('as.ltraj')

##Prep for exercise 3, Morales et al 2004#
trackData <- read.table("http://www.esapubs.org/archive/ecol/E085/072/elk_data.txt",
                        sep="\t",header=TRUE)
colnames(trackData)
elk_data <- trackData[,c("Individual","Easting", "Northing", "Habitat", "dist_water..meters.",
                         "dist_swamp..meters.", "dist_otw..meters.", "dist_openfor..meters.",
                         "dist_ntw..meters.", "dist_devel...meters.", "dist_ddf..meters.",
                         "dist_conifer..meters.")]
head(elk_data)
colnames(elk_data) <- gsub("..meters.", "", colnames(elk_data))
colnames(elk_data) <- gsub("l.", "l", colnames(elk_data))
colnames(elk_data)[1] <- "ID"
colnames(elk_data)[2] <- "X"
colnames(elk_data)[3] <- "Y"
elk_data$time <- 1:length(elk_data$ID)
elk_data <- elk_data[,c(1,13,2,3,4:12)]
elk_data[is.na(elk_data$dist_water),]
elk_data <- elk_data[1:735,]

elk_data$X <- elk_data$X/1000
elk_data$Y <- elk_data$Y/1000

elk_data[,6:13] <- scale(elk_data[,6:13])

write.csv(elk_data, 'elk_data.csv', row.names = F)
