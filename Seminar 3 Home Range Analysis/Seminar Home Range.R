#Seminar 3: Home Range Analysis
#O. R. Bidder, July 2019#

#Please install the following package while you read the introduction!#
library(adehabitatHR)
library(sf)
library(plyr)
library(move)
library(ggplot2)
library(imager)
library(tidyverse)

#Intro#
#Burt (1943) defines a home range as "that area traversed by the animal during its normal activities of food
#gathering, mating and caring for young". In ecology, its often useful to describe an animals home range 
#in space, so that we can measure it, model it, and better understand what dictates where animals establish
#home ranges. In this session, we're going to look at three methods used to delineate an animal's home-range,
#Minimum Convex Polygon (MCP), Kernel Utilization Distribution (KUD) and the Brownian Bridge Movement Model.
#We will do this using the adehabitat and move packages.

#MCP#
#Load the data
puma <- read.csv("puma_data1.csv", header = T, stringsAsFactors = F)
puma$acquisition_time <- as.POSIXct(puma$acquisition_time, "%Y-%m-%d %H:%S:%S", tz = "America/Buenos_Aires")


#Exercise 1: 
# Check the data frame for the following;
#1. What columns are present? Where are our coordinates held?
#2. How many animals are present? What are their ids?
#3. Convert puma from a data frame in to a SpatialPointsDataFrame ready for use in adehabitat. Use the head()
#   function to look for which columns hold the coordinates. Do they look like they are in WGS84 or UTM 19S 
#   projection? If WGS84, use the proj4string '+init=epsg:4326', if UTM 19S, use
#   '+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
#(5 mins)

#Ex1 Answer:
colnames(puma) #what columns have we got?
unique(puma$animals_id) #which animals have we got?
head(puma)
coordinates(puma) <- c("x","y")
proj4string(puma) <- CRS('+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#The adehabitatHR package gives us the mcp() function for calculating the MCP. It takes a spatialpointsdataframe
#with coordinates in UTM meters, the percentage of points we want to use for computation of the MCP,
#unin, which is the units of input, and unout, the units for output.

#mcp() will output an sp object class, spatialpolygonsdataframe, which is a polygon that encompasses
#95% of all the puma locations.

#Minimum Convex Polygons
mcp_all <- mcp(puma, percent = 95, unin = "m", unout = "km2") #an MCP calculated for all animals
mcp_ind <- mcp(puma[,4], percent = 95, unin = "m", unout = "km2") #MCP for each animal in puma
plot(mcp_ind); points(puma, col=puma$animals_id)

#We can also work out the area of each MCP using mcp.area
area <- mcp.area(puma[,4], percent = c(15:95), unin = "m", unout = "km2") 
#ideally you want the slope to level off to show you have a good estimation of the animal's home range

mcp_ind$area #the areas are also stored as a vector in the mcp object

#KUD#
#Let's plot our MCPs again...
plot(mcp_ind); points(puma, col=puma$animals_id) 

#Throwing a rubber band around all of our points is a useful but primitive form of home-range estimation,
#we see in our plot a lot of white space, meaning the area is in our MCP, but the animal never visited
#these sites! Hardly what Burt meant by "area traversed during normal activities"!

#The next method, Kernel Utilization Distribution, is a more formal description of a home range, by using
#a bivariate probability density function, the UD, which describes the probability that we will find
#an animal in any given area. A neat visualization that is useful when conceptualizing the UD can be found
#here; https://mathisonian.github.io/kde/

#A seminal paper on using KUD for animal movements is Worton (1989), you can find it in the folder
#for this seminar.

#The KUD can be calculated using adehabitat's kernelUD() function, which takes a spatialpointsdataframe,
#the smoothing parameter h (see bandwidth in that animation above), as well as other arguements like
#'kern' that we can use to specify which kernel we want to use ("bivorm" is default).

#The 'href' value for h calculates the 'reference bandwidth' and uses it for the KUD estimation, the
#reference bandwidth is calculated as h = sigma * n^(-1/6), 
#with sigma = 0.5 * (sigma_x + sigma_y)

#This reference bandwidth is a quick shorthand for finding h, but assumes the true distribution is
#bivariate normal, which is almost always not the case. The alternatives are h = "LSCV" to robustly find
#h by least-square cross validation (slow) or by picking a value for h by subjective trial and error
#(actually can give sensible results). For now, we will stick with 'href'...

kud_all <- kernelUD(puma, h = 'href')
kud_ind <- kernelUD(puma[,4], h = 'href')

#If you ever want to try and use LSCV to pick h, I suggest working out from the 'href' value using the 'hlim'
#arguement, which specifies the search space for h as hlim[1]*href to hlim[2]*href 

#kud_ind <- kernelUD(puma[,4], h = 'LSCV', hlim = c(0.1,2)) #Example, don't run!
#plotLSCV(kudl) #plot the LSCV results, look for convergence (flattening of the line).

#Other issues you may run in to are errors concerning your grid and extent parameters. Grid is essentially
#the resolution of your output matrix, the higher the number, the greater the number of values in the 
#grid x grid matrix. Extent is the margin around the locations that R allocates before starting the calculations.
#If you pick too small an extent, or have a lot of locations, the UD can run over the edge of the matrix.
#extent = 0.2 will add a 20% margin around, extent = 2 will add a 200% one etc.

#Let's visualize our KUDs
image(kud_all) #notice the output is a matrix, with dimensions for x and y, as well as the KUD values z
image(kud_ind) #for multiple animals, kernelUD outputs adehabitat's estUDm class, which is a list of
#the matrices, one element for each animal. Because it's a list, we need to use the double square brackets
# e.g. kud_ind[[1]], to call the individual elements

#KUD also has its own area function...
area <- kernel.area(kud_ind)
plot(area)

kud_ind$`1`@h #you can always see what 'h' was used by running the following, replacing `1` with your animal_id

#Sometimes, it's easier to show a home range as a polygon. For this we take the contour of the UD, or the 
#area in xy space that represent a proportion of the underlying distribution (i.e. probability of finding
#the animal there). For instance, the 95% contour represents the area in which we are likely to observe the 
#animal in 95% of cases.
kud_contour <- getverticeshr(kud_ind, percent = 95)
plot(kud_contour); points(puma, col=puma$animals_id)

#Exercise 2: It might be useful to plot each home range seperately, along with the points that produced them,
#in order to assess each one. Write a for loop that runs through all the animals and plots each one.
#Hint: I would do this in the following steps, 1) plot the first contour only by calling kud_contour[1,]
#2) use the points function to plot only the first animal's locations, you will need to slice puma by
#the animals_id column 3) replace the number 1 with i, and index the loop for the number of pumas
#we have. You will get bonus points if you can also write an appropriate title for each plot that tells us
#which animal is plotted (try using the paste() function) (15 mins)

#Ex2 Answer:
for (i in 1:3) {
  plot(kud_contour[i,], main = paste("KUD for Puma", unique(puma$animals_id)[i], sep = " "))
  points(puma[puma$animals_id == unique(puma$animals_id)[i],]) #plot one KUD at a time...
}

#Notice this polygon contains a lot less white space than the MCP, and provides a more realistic description
#of the home range. If you ever need to export your home range (or sp object in R), you can use the
#writeOGR() function from the rgdal package.
writeOGR(kud_contour, "Puma_95_KUD.shp", driver = "ESRI Shapefile", layer = "puma_kud")


#Dynamic Brownian Bridge Movement Model#
#The KUD gave a nice estimation of the home range, but that probabilistic method doesn't account for time
#or the sequence of locations in any way. The next development in home ranges comes from the use of brownian
#bridges, which place the bivariate kernel, not just over the locations, as with KUD, but also over the 
#space between two sequential locations, to give the bridge. This gives an indication of the latent position of the
#animal between positional fixes.

#Let's pull in an example taken from the adehabitat manual...
im<-load.image("BBMM_demo.PNG")
plot(im, axes = F)

#So the classic BBMM had two parameters of note;

#sigma_1 controls the width of the "bridge" that connects successive locations. In effect, sigma_1 is 
#related to the speed of the animal, the faster it moves, the more uncertainty we have about the 
#latent position between locations, the wider the bridge must be to capture all possible positions 
#between the recorded locations.

#sigma_2 controls the width of the "bumps" over the locations themselves. It is similar to the smoothing 
#parameter 'h' in KUD, but conceptually is related to the error when recording locations.

#In classic BBMM, sigma_1 is fixed for the duration of the track, but sigma_2 can vary, because lots 
#of GPS collars return a estimate of error at each recording, derived from factors like number of satellites
#etc. Adehabitat has functions to estimate these parameters and produce the classic BBMM.
puma <- as.data.frame(puma)

puma <- puma[complete.cases(puma$acquisition_time),] #have to remove NAs because Daylight Savings sometimes produces them when clocks roll forward

puma_ltraj <- as.ltraj(xy = data.frame("x" = puma$x, "y"=puma$y),
                       date = as.POSIXct(puma$acquisition_time), id = puma$animals_id,
                       burst = puma$animals_id, 
                       proj4string = CRS('+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

#Let's search for the optimum sig1, assuming our GPS error is 60 m
#Here's a for loop, going through each puma in puma_ltraj, finding the optimum sig1, which is stored in
#an object called lik, and then calculating the kernelbb with that sig1 value
lik <- liker(puma_ltraj, sig2 = 60, rangesig1 = c(1, 40)) #find the optimum sig_1
bbmm_puma <- kernelbb(puma_ltraj, sig1 = lik[[1]]$sig1, sig2 = 60, grid = 50)#calculate the bbmm
image(bbmm_puma)

# Notice the sig1 value is different for each animal. So usually we'd have to pick a value to share
#among all animals, OR we can write a loop that calculates the sig1 for each animal, calculates a 
#bbmm, and combines the outputs in a list...

#Exercise 3: Write another for loop that repeats these steps for the other 2 pumas, and places the
#results in to the bbmm_puma object *HINT*: assign using the list notation bbmm_puma[[i]] <- kernelbb()
#Answer 3:

bbmm_puma <- list()
for (i in 1:3) {
  print(paste('Calculating BBMM for animal: ',ld(puma_ltraj[i])$id[1], sep = ""))
  lik <- liker(puma_ltraj[i], sig2 = 60, rangesig1 = c(1, 40)) #find the optimum sig_1
  bbmm_puma[[i]] <- kernelbb(puma_ltraj[i], sig1 = lik[[1]]$sig1, sig2 = 60, grid = 50)#calculate the bbmm
  bbmm_puma[[i]]$id <- ld(puma_ltraj[i])$id[1]
}

#Now kernelbb isn't optimized for doing multiple animals at once like kernelUD is (not to my knowledge!),
#but we can use this iterative approach to make a list of results, and then combine at the end to give
#us a spatialpolygonsdataframe, like kernelUD does...

#First, change the class of bbmm_puma to adehabitats matrix class...
class(bbmm_puma) <- "estUDm"
image(bbmm_puma) #Notice the ids are gone! We can get them back later

#Let's get the contours for each animal and get the polygons for output
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
#Now we can output using the writeOGR() method

#So where does this 'dynamic' part of BBMM come in? Well with the calculations above, we assume a constant
#value for sig1 throughout the entire track. Of course, animals change their movement rates as they change
#behaviour over time. The Dynamic BBMM, or dBBMM, uses likelihood statistics to find when the animal performs
#a switch in behaviour. The process works a little like above, but by breaking the track down in to multiple
#segments and finding the optimum sig1 at each segment (for more details, see Kranstauber et al 2012, in folder)
#Of course, this all adds to the computation time, but luckily the move package has this process optimised
#for faster calculation.

#first we turn the puma dataframe in to a move object
puma_move <- ld(puma_ltraj)
puma_move <- move(x=puma_move$x, y=puma_move$y, 
               animal=puma_move$id, time=puma_move$date, data=puma_move, 
               proj='+proj=utm +zone=19 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#we can plot the object
plot(puma_move)

#I find the dBBMM works best when you set up the raster on which the calculations are made first (YMMV!)
r<-raster(puma_move, resolution=100) #an empty raster, with extent the size of puma_move, and 100 m resolution
r<-extend(r,200) #extend out the raster so the dBBMM doesn't run over the edge!

#brownian.bridge.dyn() has some optional parameters, like the margin and window.size for the change point
#analysis. I will leave these at default for now, but feel free to play around with these later.

#We will also only do this for the first puma, because computation takes a while...
plot(puma_move[[1]])
dbbmm_puma <- brownian.bridge.dyn(puma_move[[1]], 
                                 raster = r,
                                 #dimSize=150, 
                                 location.error=60  #can also provide vector with changing positional errors
                                 #margin=11,
                                 #ext=2,
                                 #window.size=21, 
                                 #time.step=12
)

poly_dbbmm <- raster2contour(dbbmm_puma, levels=.95) #Make in to a polygon if needed
plot(poly_dbbmm, add=T)

#As you can see, the home range is much tigher than the MCP and KUD methods, and provides a nice 
#characterization of the latent positions between fixes.

#Optional if time#
#Sometime we may need to know how much two or more home ranges overlap, for example, to see how much 
#two animals interact/exclude each other, or to look at site fidelity over multiple periods. Luckily,
#adehabitatHR has a functions for doing this...

bbmm_overlap <- kerneloverlaphr(bbmm_puma, method = "HR")

#There are multiple indices of overlap, the simplest being "HR", the proportion of the home range
#covered by another. Adehabitat has a few, type help('kerneloverlaphr') for a list, and check Fieberg
#and Kochanny (2005, in folder) for a review.

#Exercise 3#
#For the final exercise in this class, perform a home range analysis from start to finish. Please do the 
#following;
#1. Load in the file 'puma_data2.csv" and convert it to a spatialpointsdataframe, note the data are
#   already in the UTM19S projection.
#2. Calculate either a KUD or BBMM, either for all animals in the data, or for one animal across multiple years.
#   For multiple years, either split the data frame by year, or assign a burst id using the paste() function.
#3. Use the kerneloverlaphr() function above, what does it's output tell us. Can you deduce anything about
#   the animals from the results?


