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
library(lubridate)
library(adehabitatLT)
library(setRNG)
library(momentuHMM)

#Note, momentuHMM is now on CRAN, so you can try installing as you would any other package.
#However, in case that doesn't work we can install from GitHub.
#Install momentuHMM from GitHub
#library(devtools)
#install_github("bmcclintock/momentuHMM")

#Intro

#In this session we will learn how to use the momentuHMM package to fit Hidden Markov Models to animals
#movement data. Some of the data and code has been adapted from this package's vignette to ensure this 
#session runs as smoothly as possible and that we can talk about some interesting results. With that in mind
#a detailed vignette for this package is available here;

# https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf

#There's also a paper that describes the package, along with its applications in the folder of this
#seminar (McClintock et al. 2018).

#In this first section I want to describe the motivation and justification behind HMMs. I think it would
#be better to SHOW you rather than TELL you.

#Section 1: high-level behaviors can be characterized by their step lengths and turning angles.

#Let's load in some GPS location data from an elk of the Wiggins Fork herd in Wyoming.
wigs <- read.csv('Wiggins1.csv')
wigs <- wigs[,-1] #Get rid of the extra index column
wigs$Date_Time <- as.POSIXct(wigs$Date_Time, "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

#Exercise 1: I would like you to do two things;

# 1) subset wigs in to two specific time periods; 
#     one starting on January 1st 2016 and ending on February 1st 2016, 
#     and a second that starts on the 11th of May 2016 and ends on the 18th of May 2016. 

#To do this, define a lubridate::interval for these two periods, and then subset wigs using these 
# interval object using the boolean operator %within%. Call the resulting subsets 'winter' and 'spring'
# respectively. 
# Remember an interval is defined with int_obj <- interval(start, end), with dates in YYYY-MM-DD format
#
# 2) convert both 'winter' and 'spring' subsets in to an adehabitat:ltraj object. To do this, you can use the
# following template...
#ltraj_obj <- as.ltraj(xy = data.frame("x" = <X coordinates>, "y"= <Y Coordinates>), 
#                       date = as.POSIXct(<datetime column>), 
#                       id = <animal ID>, 
#                       burst = <animal ID>, 
#                       proj4string = CRS('+init=epsg:32612'))
#
# 3) plot the ltraj objects 'winter' and 'spring'. What do you notice about the track. What can you infer
# about this animal's behavior from the shape of the track. Which features inform your answers?

#Answer 1:
#####################################################################################################
winter_int <- interval("2016-01-01", "2016-01-07")
spring_int <- interval("2016-05-12", "2016-05-14")

winter <- subset(wigs, Date_Time %within% winter_int) #One way using subset
spring <- subset(wigs, Date_Time %within% spring_int)


winter <- wigs[wigs$Date_Time %within% winter_int,] #another by slicing the DF
spring <- wigs[wigs$Date_Time %within% spring_int,]

winter <- as.ltraj(xy = data.frame("x" = winter$UTM_X, "y"= winter$UTM_Y), 
                       date = as.POSIXct(winter$Date_Time), 
                       id = winter$Elk, 
                       burst = winter$Elk, 
                       proj4string = CRS('+init=epsg:32612'))

spring <- as.ltraj(xy = data.frame("x" = spring$UTM_X, "y"= spring$UTM_Y), 
                   date = as.POSIXct(spring$Date_Time), 
                   id = spring$Elk, 
                   burst = spring$Elk, 
                   proj4string = CRS('+init=epsg:32612'))

plot(winter) #Elk has a mixture of steps, some of them small, going back on itself.
plot(spring) #track is much more linear, animal goes in one direction. Looks like a migration
#####################################################################################################

#Now if I use the ld function from adehabitat to turn the ltraj objects back in to data frames, I can see
# adehabitat has added a few new columns...
winter <- ld(winter)
spring <- ld(spring)

colnames(winter)
#Of these new columns two are of particular interest;
# dist: the distance in m between that location and the next one
# rel.angle: the relative change in direction needed to hit the next location

#If we compare the distributions of dist and rel.angle between our winter and spring tracks...
hist(spring$dist, col= 'red', breaks = 50, freq = F); hist(winter$dist, col = 'blue', add=T, breaks = 50, freq = F)

#We see that although spring (red) has quite a lot of small steps, it also has some very large ones, >2 km
#In winter (blue), these large steps distances are less common, with the majority being 0-500 m

#Let's look at relative angle
hist(spring$rel.angle, col= 'red', breaks = 50); hist(winter$rel.angle, col = 'blue', add=T, breaks = 50)

#Spring seems concentrated around 0, and winter movements have a more even distribution of turning angles

#Essentially what we are looking at here is the difference between directional movement conducted during
# a migration, and the more 'brownian motion' type movement of foraging during winter. So I hope this
# demonstrates that we can characterize behaviors or 'movement states' in terms of dist (or step lengths)
# and rel.angle (or turning angles). This principle forms the basis for the set of models we will look at today
# 'state-space' models or Hidden Markov Models, HMM.

#Now in this session, I don't think it would be helpful to dive too deeply in to the math behind HMMs, 
#but there are a couple of key concepts I want you to come to grips with before we learn how to apply HMMs 
#to animal movement data.

#Let's load a figure that demonstrates the process by which we assume the data have arisen, 
# the Markov Chain of states according to time...
im<-load.image("hmm_diagram.jpg")
plot(im, axes = F)

# The unobservable 'state':
# In a HMM, we assume that the data (y in the diagram, our locations and the step lengths and turn angles 
# we derive from them) are produced according to whatever state (x in diag) our animal is in at the time. 
# Example states could be foraging or directional displacement (like our migration above). When the animal is 
# in each state, we will observe differemt patterns in step lengths and turn angles, but we have to use 
# these patterns to infer the current state because we can't measure that directly. Thus, depending on the step 
# length and turn angle, we can calculate a likelihood that the animal is in any given state at that 
# point in time (observation probability).

# States are autocorrelated:
# If an animal is in the 'foraging' state at time 1, then one minute later they are still likely to be foraging.
# However, the longer time goes on, the less and less likely the state at time 1 will persist. Go one hour
# forward and the animal might still be foraging, one day though? One week? Animals can perform transitions 
# between all of their possible states (this matrix of state transition probabilities is the 'state space'),
# and these transition probabilities may differ according to their environment e.g. an animal that is in a 
# displacement state may be more likely to transition to a foraging state when they enter a region that has
# a lot of food available.

# The take home here being that we can use HMMs to estimate the number of behavioral states that could
# produce our observed data, the most likely sequence of states that would produce our data, and how
# changes in the environment might affect the state transition probabilities. Powerful stuff.

# The package we'll we using to fit these models has a simple workflow;
im<-load.image("momentuHMM_workflow.png")
plot(im, axes = F)

#Let's learn some more through an example...

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

# If we look at the locations obtained between March 22nd 2008 and September 30th 2010, we see we have some
# hourly locations missing...
int_missing <- interval("2008-03-22", "2010-09-30") #define an interval for this period
miss_period <- rawData[rawData$time %within% int_missing,] #subset rawData to times in this period
no_days <- as.numeric(as.Date("2010-09-30")-as.Date("2008-03-22")) #get the number of days

length(miss_period$time)/(no_days*24) # % of hourly samples take for this time period, around 1% missing!

# This 1% of missing locations isn't massive compared to some data sets, but HMM models do need continuous
# location data at regular intervals. Traditionally we would discount this period or use a different method.

# However, momentuHMM incorporates the functions of the crawl package to impute the missing data and make
# it a regularly sampled track. It does this using the CTCRW model described by Johnson et al. (2008), which
# is in this seminar's folder for more detail.

# fit crawl model
crwOut <- crawlWrap(obsData=rawData, timeStep="hour",
                    theta=c(6.855, -0.007), fixPar=c(NA,NA))
plot(crwOut)

# create momentuHMMData object from crwData object, taking one of the calculated tracks
elephantData <- prepData(data=crwOut, covNames="temp")

# add cosinor covariate based on hour of day
elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))

# Let's look for any periodicity in the data, in this next step, we will calculate the Autocorrelation Function
acf(elephantData$step[!is.na(elephantData$step)],lag.max=300)

# So here we are plotting the correlation between the step lengths (y) as a function of the distance in time
# between them (the Lag on x). What do you notice about the pattern? We'll need to take this pattern in to
# consideration in our analyses.

#Here we're aiming to fit a 2 state HMM model to the elephant track. The first state we expect to be
#a 'exploratory' state, with longer step length and a constrained distribution of turn angles biased toward 0.
#The second state would be an 'encamped' state, with smaller step lengths and an unconstrained distribution
#of turn angles.

#Let's start by fitting a model with no covariates. This should give us some good starting parameters to
#fit our more complex model later on#

# label states
stateNames <- c("encamped","exploratory")

# We are going to fit our step lengths to a gamma distribution, with parameters mean and standard deviation
# and a wrapped cauchy distribution to our turn angles, with the parameter concentration

# Distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")


#To get some good starting values, let's look at our data
plot(elephantData) 

#Note that the gamma has two parameters, shape and scale, is a continous distribution, many low values
#The wrapped cuachy describes turning angles from -pi to +pi, it has one parameter, concentration.

# initial parameters
Par0_m1 <- list(step=c(100, #shape encamped
                       500, #shape exploratory
                       100, #scale encamped
                       200 #scale exploratory
                       ),
                angle=c(0.3, #concentration 1
                        0.7 #concentration 2
                        ))

# We use fitHMM to fit the model
m1 <- fitHMM(data = elephantData, #the data
             nbStates = 2, #number of states
             dist = dist, #the distributions we are using for the data streams
             Par0 = Par0_m1, #the initial values for the distributions
             estAngleMean = list(angle=FALSE), #whether or not we want to calculate the angle mean
             stateNames = stateNames) #labels for the states    

#Note, if the model won't converge, you can use the argument retryFits to up the number of attempts to
#minimize the negative log-likelihood

#We can plot the model distributions for step length and turn angle using the plot function...
plot(m1)

# First plot, we see the modeled distributions for step length, with encamped having lots of short steps
# and exploratory having more at 500+ m.
# Second plot, we have the distributions for turning angle, encamped with an even distribution around 0,
# exploratory with a lot more bias towards 0 turning angle, more directed movement,
# Third plot is our most likely state sequence according to our HMM model.

#So we have a nice model for two elephant behavioral 'states'. Obviously, animals don't cycle between
#states at random, they are usually responding to their environment in some way.

#To illustrate this, let's incorporate our first covariates; time of day and temperature.

#First, let's think about the biological justification for each.

#time of day; elephants need to rest at some point, at which point their movement rates should be reduced.
#Without assuming anything about this animal's activity pattern (nocturnal, crepuscular, diurnal), we can
#guess that time of day may play some kind of role in determining if the animal is moving a lot, or not 
#much at all. Our investigation of ACF yielded a pattern that seemed to cycle every 24 hrs.

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

#The parameters themselves are accessed within this object
Par0_m2$Par

#Exercise 2:
# Adapting the code above, fit our second model using the fitHMM function, assigning our model object with 
# the name "m2". You will need to find a way to include the formula, as well as using the parameters obtained 
# from the first model (5 mins)

#Hint: use help('fitHMM') to learn more about the arguments that can be passed to this function

#############################################################################################
#Answer 2:
# fit model
m2 <- fitHMM(data = elephantData, 
             nbStates = 2, 
             dist = dist, 
             Par0 = Par0_m2$Par,
             #beta0=Par0_m2$beta, 
             stateNames = stateNames, 
             formula=formula)

#############################################################################################

#In model 'm2', we are modelling the transition probabilities as a function of the covariates. We can
#look at how this works by plotting the model object.
plot(m2)

#First, notice that for our plot of temp, hour is being help at its mean (Trans prob: hour = 11.5)
#The animal seems to have the highest probability of switching 1 -> 2 at low temperatures.
#The higher the temp, the more likely we are to observe a 2 -> 1 state transition.

#Second, what can we see about how the state probabilities change with respect to time? 
# 1 -> 2 transitions are most likely to occur between 1500 and 2000 hrs
# 2 -> 1 transitions are most likely to occur between 0500 and 1000 hrs

#We can extract the most likely state sequence from our model using the viterbi function
help('viterbi')
elephantData$State <- viterbi(m2)

#or if we want the names themselves
elephantData$State <- stateNames[viterbi(m2)]

# Now there is one last version of the model that I'd like to demonstrate on the elephant data, but its
# a little more complex. So far we have allowed the state transition probabilities to differ according
# to the hour and temp covariates. In addition to this, we can allow the parameters of the distributions
# for step length and turn angle (step length; gamma with mean and sd, turn angle; concentration) themselves
# to differ according to the covariates...

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24), #step shape a function of temp and hour
                       sd = ~ temp * cosinor(hour, period = 24)), #step scale a function of temp and hour
           angle = list(concentration = ~ temp)) #angle concentration a function of temp only

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)

# fit model
load('m3_model.RData') #this will be faster than fitting the model again on your machine, but 
#the code below will re-run the model... (Don't run it right now though!)

m3 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
             formula = formula)

#save(m3, file = 'm3_model.RData') #in case you didn't know how to save R environment objects!

# decode most likely state sequence
elephantData$State_m3 <- stateNames[viterbi(m3)]

#Are they the same?
table(elephantData$State == elephantData$State_m3)
(1686/20448) *100 #about 8% disagreement

# derive percentage of time spent in each state
table(elephantData$State_m3)/nrow(elephantData)

# Compare models using AIC
AIC(m1, m2, m3) #lowest AIC indicates best fit, m3 far lower

# visualize the model
plot(m3, plotCI = TRUE, covs = data.frame(hour=12))

# Interestingly, this model suggests step lengths (plot 1) and directional persistence (plot 10) for the
#'encamped' state decreased as temperature increased. Step lengths for both states
#tended to decrease in the late evening and early morning (plot 2 and 6), and transition probabilities
#from the 'encamped' to 'exploratory' state decreased as temperature increased (plot 13).

# compute pseudo-residuals for the steps and the angles
pr <- pseudoRes(m3)

# plot the ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

#These pseudo residuals show that we've managed to explain a lot of the periodicity in step with our model, 
# but there's still some periodicity (at both 12 and 24 hr lags) that isn't explained (producing residuals
# between our model prediction and the observed data). This might be because of missing covariates that
# operate on this temporal scale.

#Section 3: Additional data streams

#So far we fitted our models to two data streams, step lengths and turning angles. There are of course other
#data that biologgers can record from animals. Can you think of a few?

#One of the great things about momentuHMM is that it can incorporate additonal data streams, we'll now run 
# a quick example, using northern fur seal data from McClintock et al. (2014, see folder). These data are 
# interesting because they are collected by Fastloc GPS (sea water attenuates GPS signals, so locations 
# are obtained with error and at irregular intervals), and also have corresponding Time Depth Recorder 
# data, which are represented by the number of foraging dives detected between locations.


#So let's go back and look at the workflow we could follow
im<-load.image("momentuHMM_workflow.png")
plot(im, axes = F)

#We have irregular data, so will need to use crawlWrap again, but because we have very irregular data
# (more than just the 1% missing in the elephant example), we'll need to run crawlWrap many times
# (i.e. Multiple Imputation), then fit our HMM model to all of these tracks, averaging over all the models
# to get our parameter estimates that take in to account our uncertainty in the true track. Exciting!

# First let's set the seed so we all get the same results...
oldRNG<-setRNG::setRNG()
setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=30)

# load the data from the folder...
load('nfsData.RData')

#Note here I'm loading in an .RData file, which contains 3 R environment objects;
# 1) nfsData - the telemetry location data, projected to UTM 2N (EPSG:32602)
# 2) foragedives - the TDR data, with time period and number of dives recorded
# 3) Par - the starting parameters for the distributions that we'll use

#nfsFits and crwOut objects are just results from later in this script, which I'm loading to speed
#up this demonstration.

nSims <- 30 # number of imputatons
retryFits <- 30 # number attempt to re-fit based on random perturbation
ncores <- 7 # number of CPU cores

# So because of the irregular location data, we'll need to use the crawl method again. Because we have hourly
# dive records, we need to make sure we are predicting the locations at times that match that data.
# Luckily, we can specify the times we'd like location estimates for by passing them to crawlWrap.

# set time steps based on overlapping time period of location and dive data
predTimes <- seq(max(foragedives$time[1],nfsData$time[1]),
                 min(foragedives$time[nrow(foragedives)],
                     nfsData$time[nrow(nfsData)])+60*60,"hour")

head(predTimes, 10) #a vector of times to match TDR times

# use default starting values (theta) and explore likelihood surface using retryFits
crwOut<-crawlWrap(nfsData, #data as dataframe with time and ID columns
                  ncores=ncores, #number of CPUs to use
                  retryFits=30, #number of times to try and converge, using different starting values
                  err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr), #error models to use
                  fixPar=c(1,1,NA,NA), #fixing first two parameters, estimate the NAs
                  attempts=100, #number of times to try and converge, using each parameter combination
                  predTime=predTimes) #the times we want to extract locations for
plot(crwOut)

crwOut$crwFits$NFSF0907$nms #names of estimated params

#tau.x and tau.y are scaling factors for the error in x and y planes, but are data are UTM so we can 
#fix those to 1
#sigma is the velocity variation
#beta is the autocorrelation in velocity between samples

crwOut$crwFits$NFSF0907$par #fixed first to params

# We can merge the crwData object with forage dive data using crawlMerge
crwOut <- crawlMerge(crwOut,foragedives,"time")

nbStates <- 3
stateNames <- c("resting", "foraging", "transit")

#We are using similar distributions used in elephant example for step and angle, and adding another, the
# poisson distribution for the dive counts. A quick recap on the poisson distribution.
hist(rpois(1000,10), col = 'red', breaks = 50) #Poisson has a parameter 'rate'
hist(rpois(1000,5), col = 'blue', breaks = 50) #It is a discrete distribution, all values are integers
hist(rpois(1000,1), col = 'green', breaks = 50) #All values are above 0, so it's good for counts!

#Ok so we specify which distributions to use in a list as before.
dist <- list(step = "gamma", angle = "wrpcauchy", dive = "pois") 

#So the next few lines of code are set up just to put boundaries on the values for the parameters for each
# state. Otherwise our 1st state might drift and end up being our 'foraging' or 'transit' state when we
# intended it to be our 'resting' state. The values here were decided upon by McClintock et al., so I'm
# just reproducing them so that our model behaves itself.

### construct pseudo-design matrix constraining parameters (to avoid label switching across imputations)
# constrain step length mean parameters: transit > resting
stepDM<-matrix(c(1,1,0,0,0,0,
                 0,0,1,0,0,0,
                 1,0,0,0,0,0,
                 0,0,0,1,0,0,
                 0,0,0,0,1,0,
                 0,0,0,0,0,1),2*nbStates,6,byrow=TRUE,
               dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),
                             c("mean_13:(Intercept)","mean_1","mean_2:(Intercept)",
                               paste0("sd_",1:nbStates,":(Intercept)"))))

stepworkBounds <- matrix(c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,Inf,0,rep(Inf,4)),6,2,
                         dimnames=list(colnames(stepDM),c("lower","upper")))

# constrain turning angle concentration parameters: transit > resting
angleDM<-matrix(c(1,1,0,
                  0,0,1,
                  1,0,0),nbStates,3,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_13:(Intercept)","concentration_1","concentration_2:(Intercept)")))

angleworkBounds <- matrix(c(-Inf,-Inf,-Inf,Inf,0,Inf),3,2,dimnames=list(colnames(angleDM),c("lower","upper")))

# constrain dive lambda parameters: foraging > transit
diveDM<-matrix(c(1,0,0,
                 0,1,0,
                 0,1,1),nbStates,3,byrow=TRUE,
               dimnames=list(paste0("lambda_",1:nbStates),
                             c("lambda_1:(Intercept)","lambda_23:(Intercept)","lambda_3")))

diveworkBounds <- matrix(c(-Inf,-Inf,-Inf,rep(Inf,2),0),3,2,
                         dimnames=list(colnames(diveDM),c("lower","upper")))

DM<-list(step=stepDM,angle=angleDM,dive=diveDM)

workBounds<-list(step=stepworkBounds,angle=angleworkBounds,dive=diveworkBounds)

Par0 <- getParDM(nbStates = nbStates, dist = dist,
                 Par = Par, DM = DM, workBounds = workBounds,
                 estAngleMean = list(angle = FALSE))

fixPar <- list(dive = c(-100, NA, NA))

# set prior to help prevent working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#So this time we will pass a crwOut object to the function MIfitHMM, which is a wrapper for fitHMM that;
#    crwSimulator - draw parameters from the posteriors of the previously fitted crawl
#    crwPostID - use those parameters to predict a new track
#    prepData - create a new momentuHMM object with the relevant data
#    fitHMM - feed that object to fitHMM and run the space state model

#This process is repeated for the number of times specified in nSims. We do this to account for the 
# uncertainty in the track estimates. It's a pragmatic approach to deal with the fact that HMMs need
# discrete time data obtained at regular intervals.

#Now we fit the model again, but this time we pass some different args...
nfsFits <- MIfitHMM(crwOut, #the crawl object itself (not prepData object as before)
                    nSims = nSims, #the number of simulated tracks to use
                    ncores = ncores, #the number of processors we'd like to use for parallel processing
                    nbStates = nbStates, #number of states
                    dist = dist, #the distributions we want to use for each data stream
                    Par0 = Par0, #the Par0 object that contains the initial values for the parameters
                    DM = DM, #Design matrix for the prob distribution parameters for each stream
                    workBounds = workBounds, #bounds for the prob distribution, transition probs and initial parameters
                    estAngleMean = list(angle = FALSE), #whether or not to estimate the angle means
                    fixPar = fixPar, #list of parameters for which we assume a value
                    retryFits = retryFits, #how many times to iteratively refit the model
                    prior = prior, #stop the parameter estimates from getting stuck on the boundary
                    stateNames=stateNames) #our statenames

print(nfsFits)
plot(nfsFits,legend.pos="topright")

#In the plots we can see that seals in 'transit' take bigger steps between locations, followed by foraging,
# with resting seals taking overall smaller steps.

#The turning angles of 'transit' are also far more concentrated around 0, indicating much more directional
# movement. Foraging has some bias towards 0 but also includes a fair amount of 90 degree turns (1.57 radians).
# Resting seals have an almost uniform distribution of turning angles.

#Our third data streams, the number of dives, clearly distinguishes foraging states from the others.

#So, we see that adding a third data stream provides extra means to differentiate states. Now if we had
# additional environmental covariates, we might be able to investigate which of those causes seals to 
# transition to feeding states, and thus infer which resources or landscape features are important to 
# seal feeding behavior. So I hope you can see, combining data streams gives us the ability to do some
# pretty cool science.

#########################################################################################################
#Exercise 3#
# Providing there is time, I'd like us all to work with the folks near us to run a simple space-state model
# from start to finish. We will use a new data set from elk.

elk_data <- read.csv('elk_data.csv')

#elk_data <- subset(elk_data, ID != "elk-115")
head(elk_data)

#Steps to follow;
#1) Use the prepData function from momentuHMM to convert the data to the appropriate format.
# Note: To speed the model fitting up, UTM coordinates have been divided by 1000 to convert to Km
#      Covariates have already been scaled and centered using scale()
#      You'll need to specify where the coordinates are using the coordNames = c("X","Y") argument 
#       in prepData

#2) Fit a 2-state model with fitHMM. Take a look at the plots and parameter estimates.

# Notes: Use a gamma ("gamma") distribution for step and Von Mises ("vm") distribution for angle
#        Use plot() on the prepData output to guess at some good initial values for the parameters 
#        There are many 0's in the step lengths, so you will need to add a pair of zeromass parameters,
#        which should be between 0 and 1 (i.e. the proportion of values that are 0s)
#        The vm will need initial values to be in the format mean1, mean2, concentration 1, concentration 2
#        Because we are including the mean parameters, we will need to set estAngleMean = list(angle=TRUE)
#         when running fitHMM

#3) Fit a final model that includes dist_water, dist_openfor and dist_devel as covariates in the transition 
# probabilities part of the model

#4) Plot and discuss the results.

#Answers 3#
########################################################################################################
#1
colnames(elk_data)
elkData <- momentuHMM::prepData(data=elk_data, 
                                covNames= c('dist_water',
                                            'dist_openfor',
                                            'dist_devel'), 
                                coordNames = c("X","Y"))

#2
# label states
stateNames <- c("encamped","exploratory")

# Distributions for observation processes
dist = list(step = "gamma", angle = "vm")

#To get some good starting values, let's look at our data
plot(elkData)

# initial parameters
Par0_m1 <- list(step=c(0.1, #shape encamped
                       1, #shape exploratory
                       0.1, #scale encamped
                       1, #scale exploratory
                       0.05, #encamped zeromass
                       0.01 #exploratory zeromass
),
angle=c(pi,
        0,
        1, #concentration 1
        1 #concentration 2
))

# We use fitHMM to fit the model
m1 <- momentuHMM::fitHMM(data = elkData, #the data
             nbStates = 2, #number of states
             dist = dist, #the distributions we are using for the data streams
             Par0 = Par0_m1, #the initial values for the distributions
             estAngleMean = list(angle=TRUE), #whether or not we want to calculate the angle mean
             stateNames = stateNames,
             formula = ~1) #labels for the states    
m1
plot(m1)

m2 <- momentuHMM::fitHMM(data = elkData, #the data
                         nbStates = 2, #number of states
                         dist = dist, #the distributions we are using for the data streams
                         Par0 = Par0_m1, #the initial values for the distributions
                         estAngleMean = list(angle=TRUE), #whether or not we want to calculate the angle mean
                         stateNames = stateNames,
                         formula = ~ dist_water + dist_openfor + dist_devel) 
plot(m2)
AIC(m1, m2)

#Interpretation:
#state1 has small steps (~350 m), high concentration around pi
#state2 has larger steps (~3.5 km), high concentration around 0


#Looking at state transitions;
#Animals more likely to transition to exploratory when further from water (1->2)
#Unlikely to stay in exploratory (2->2) while still close to water

#Animals are likely to stay encamped near open forests (1->1)
#But those in exploratory more likely to slow when near open forests (2->1) and more likely
#to stay in exploratory state if far from an open forest.

#Animals more likely to stay in exploratory when near development (2->2)
#Switch to encamped once distance to development increases
#However, animals that are already encamped unlikley to perform transition as a result of dist_devel

#############################################################################################################
##Extra Code##
#In case you want to fit a three state model during the session...
stateNames2 <- c("encamped","exploratory", "directed")

# initial parameters
Par0_3state <- list(step=c(0.1, #shape encamped
                       0.5, #shape exploratory
                       3, #shape directed
                       0.05,#scale encamped
                       0.5, #scale exploratory
                       1, #scale directed
                       0.05, #encamped zeromass
                       0.01, #exploratory zeromass
                       0.01 #directed zeromass
),
angle=c(pi, #mean encamped
        pi, #mean exploratory
        0, #mean migration
        1, #concentration encamped
        1, #concentration exploratory
        1 #concentration migration
))

# We use fitHMM to fit the model
m_3state <- momentuHMM::fitHMM(data = elkData, #the data
                         nbStates = 3, #number of states
                         dist = dist, #the distributions we are using for the data streams
                         Par0 = Par0_3state, #the initial values for the distributions
                         estAngleMean = list(angle=TRUE), #whether or not we want to calculate the angle mean
                         #stateNames = stateNames2,
                         formula = ~1) #labels for the states    
plot(m_3state)

AIC(m_3state) #the three state model has the lowest AIC
##########################################################################################################