#Seminar 6: Hidden Markov Models for Animal Movement Ecology
#O. R. Bidder, September 2019#

#Please install the following packages while you read the introduction!#
library(ggplot2)
library(imager)
library(sp)
library(plyr)
library(move)
library(tidyverse)
library(rgdal)

#Install momentuHMM from GitHub
#library(devtools)
#install_github("bmcclintock/momentuHMM")

#Note, momentuHMM is now on CRAN, so you can try installing as you would any other package!
library(momentuHMM)

#Intro


#In this session we will learn how to use the momentuHMM package. Some of the data and code has been
#adapted from this package's vignette to ensure this session runs as smoothly as possible and 
#that we can talk about some interesting results.

##HMM Example - Elephants##

#In this example we will make use of some publicly available data (collected by Wall et al. 2014, 
#see paper in folder) from bull African elephants. This data has previously been downloaded from
#movebank.org and is ready to be loaded in to R.

rawData <- read.csv('elephants.csv')
rawData <- rawData[,-1]
  
#There are two tracks in this data set and for this example we will consider only one

# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])

#let's look at the data...
head(rawData)

# convert times from factors to POSIX
rawData$time <- as.POSIXct(rawData$time,tz="GMT")

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(rawData[,3:4],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))

utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))

# add UTM locations to data frame
rawData$x <- attr(utmcoord,"coords")[,1]
rawData$y <- attr(utmcoord,"coords")[,2]

#We've now appended the coordinates to the rawData df...
head(rawData)

#A note about crawl...

# fit crawl model
crwOut <- crawlWrap(obsData=rawData, timeStep="hour",
                    theta=c(6.855, -0.007), fixPar=c(NA,NA))

# create momentuHMMData object from crwData object
elephantData <- prepData(data=crwOut, covNames="temp")

# add cosinor covariate based on hour of day
elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))

acf(elephantData$step[!is.na(elephantData$step)],lag.max=300)

#Here we're aiming to fit a 2 state HMM model to the elephant track. The first state we expect to be
#a 'exploratory' state, with longer step length and a constrained distribution of turn angles biased toward 0.
#The second state would be an 'encamped' state, with smaller step lengths and an unconstrained distribution
#of turn angles.

#Let's start by fitting a model with no covariates. This should give us some good starting parameters to
#fit our more complex model later on#

# label states
stateNames <- c("encamped","exploratory")

# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")

# initial parameters
Par0_m1 <- list(step=c(100,500,100,200),angle=c(0.3,0.7))

# fit model
m1 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames)     

#Note, if the model won't converge, you can use the argument retryFits to up the number of attempts to
#minimize the negative log-likelihood

#We can plot the model distributions for step length and turn angle using the plot function...
plot(m1)

#So we have a nice model for two elephant behavioral 'states'. Obviously, animals don't cycle between
#states at random, they are usually responding to their environment in some way.

#To illustrate this, let's incorporate our first covariates; time of day and temperature.

#First, let's think about the biological justification for each.

#time of day; elephants need to rest at some point, at which point their movement rates should be reduced.
#Without assuming anything about this animal's activity pattern (nocturnal, crepuscular, diurnal), we can
#guess that time of day may play some kind of role in determining if the animal is moving a lot, or not 
#much at all.

#temp; animals may have very different needs in hot weather (water, shady areas), or cool temperatures
#may be good for finding food etc.

#We may expect that there is some kind of interaction between temp and time of day. If so, the following model
#formula might be appropriate for modelling the transition probabilities...

# formula for transition probabilities
formula <- ~ temp * cosinor(hour, period = 24)

#Notes: two things to look at here, the cosinor function and our astericks (*) in the model formula;

#Cosinor
#Here we model time of day as the cos and sin of hour using the following formulae
#cos(2 * pi * hour/period)
#sin(2 * pi * hour(period)

#we do this because time is a circular variable, and 1 hour either side of 0 needs to have the same
#value. If we modelled time as a linear variable, we would run in to problems. Other temporal variables
#like day of the year need to be handled in the same way.

#Astericks
#the * signifies an interaction, as well as adding first order variables. Another way to write this would be
#formula <- ~ temp + cosinor(hour, period = 24) + temp:cosinor(hour, period = 24)


#As stated earlier, we can use our first model to get some parameter estimates from which to start from
#when fitting our second model

# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)

#Exercise 1:
# Using the code above, fit our second model using the fitHMM function, assigning our model object with the name
#"m2". You will need to find a way to specify the formula, as well as using the parameters obtained from 
#the first model (5 mins)

#############################################################################################
#Answer 1:
# fit model
m2 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)

#############################################################################################

#In model 'm2', we are modelling the transition probabilities as a function of the covariates. We can
#look at how this works by plotting the model object.
plot(m2)

#First, notice that for our plot of temp, hour is being help at its mean (Trans prob: hour = 11.5)
#The animal seems to have the highest probability of switching 1 -> 2 at low temperatures.
#The higher the temp, the more likely we are to observe a 2 -> 1 state transition.

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
                       sd = ~ temp * cosinor(hour, period = 24)),
           angle = list(concentration = ~ temp))

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)

# fit model
m3 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
             formula = formula)
