#Seminar 7: Net Squared Displacement Modelling
#O. R. Bidder, June 2019#

#Please install the following package while you read the introduction!#
install.packages("devtools") #if you don't have devtools installed, otherwise skip this step
library(devtools)
install_github("dbspitz/migrateR/migrateR", build_vignettes = T)

library(migrateR)
library(plyr)

#Introduction#
#In this class we're going to look at how we can characterize patterns in animal movements using 
#Net Squared Displacement (NSD) modelling. NSD is the squared straight line distance between two points.
#When studying animal movements, we calculate the NSD between an arbitrary starting point and each
#successive location in our data. By tracking NSD over time, and looking at the patterns, we can see 
#if our animal exhibits migratory movement patterns, is nomadic, resident etc. Knowing this about
#our animal not only provides an insight to their ecology, but by fitting curves to the NSD data, we can
#pull objective estimates of parameters such as migration distance, timing etc.

#In this session we will simulate our own NSD data, see how different equations produce patterns in the data,
#the fit non-linear models to the data to determine which model best describes our NSD data.

#After we have tried fitting models to simulated data, we will shift to the migrateR package, and fit
#these same models to real animal collar data. These exercises should then equip you with the code needed
#to perform these analyses on your own data in future.

#Below are the 5 common equations used to characterize patterns in NSD data, which were originally proposed
#by Bunnefeld et al. (2011).

# Migrant: NSD* = delta/(1+exp((theta-t)/phi)) - delta/(1+exp((theta+2*phi+2*phi2+rho-t)/phi2))
# 
# Mixed Migrant: NSD* = delta/(1+exp((theta-t)/phi))- (delta*zeta)/(1+exp((theta+2*phi+2*phi2+rho-t)/phi2))
# 
# Disperser: NSD* = delta/(1+exp((theta-t)/phi))
# 
# Nomad: NSD* = beta * t
# 
# Resident: NSD* = gamma*(1-exp(kappa*t)

#Part 1#
#Let's take a look at each parameter and generate some data from these equations

t <- 1:365 #time in julian days, 1-365 (notice this is a vector of numbers)

delta <- 5000 #distance between ranges, in NSD (km^2)

theta <- 60 #mid-point in time of outward migration

phi <- 10 #duration between 50% and 75% of first outward migration
phi2 <- 10 #duration between 50% and 75% of second outward migration

rho <- 160 #duration of summer range occupancy in julian days

zeta <- 0.8 #proportion of return leg completed by a mixed migrant

gamma <- 200 #mean NSD around a resident animal's range
beta <- 0.2 # a constant, roughly the movement rate
kappa <- -0.01 #a constant

#Now lets punch these values in to the equations above and look at the NSD data

NSD <- (delta/(1+exp((theta-t)/phi)) - delta/(1+exp((theta+2*phi+2*phi2+rho-t)/phi2)))  #migrant
plot(t, NSD) 

#See how the animal increases the NSD until the delta value, then after a period of time equivalent
#to rho, returns back to the initial point? This elegant pattern characterizes the important process of animal
#migration, in which millions of tons of biomass move across our planet in response to seasonal changes
#in resouce availability.

#Exercise 1: Play with the values for delta, theta, phi and rho, and observe how the NSD data change (5 mins)
#Exercise 2: Edit the equation to make it produce a mixed-migration pattern (hint: look at the zeta parameter),
# how does this curve differ to the migrant one? Which values of zeta might make it difficult to discern between
# the two model types?
#Ex1 Solution:
###################################################################################################

###################################################################################################

#Now we will produce a more realistic NSD curve, with some noise to represent the movement an animal
#might make to seek out resources, or avoid predation

noise <- rnorm(365, mean = 0, sd = 200) #vector of 365 numbers, around zero so NSD might move up or down
NSD <- (delta/(1+exp((theta-t)/phi)) - delta/(1+exp((theta+2*phi+2*phi2+rho-t)/phi2))) + noise  #migrant with noise
plot(t, NSD) 

#Now let's fit a non-linear model to this data, using nls() function and the equation above

#nls() takes a formula and some starting values. This helps the function converge quickly.
#Notice here, the only data we are passing are the NSD values, and time (t). We've substituted in letters
#to represent the different parameters (delta, theta etc.), because if we were doing this with real data we 
#wouldn't know their values, and would need to estimate them with our model!

mig_mod <-nls(NSD ~ d/(1+exp((th-t)/p1)) - d/(1+exp((th+2*p1+2*p2+rh1-t)/p2)), 
              start = list(d=5000,
                            th=60,
                            p1=10,
                            p2=10,
                            rh1=150))
summary(mig_mod) #are the estimates for the parameters close to the ones above? The added noise has muddied
#our estimates a bit (we don't get the exact values above). 
#The more imprecise our measurements (i.e. the more noise), the worse they'll be.
#Luckily, GPS data are pretty accurate, but often migratory processes aren't the only thing dictating where
#animals travel on a fine scale. If other influences are strong, the more noise there will be.

#Let's plot the fit
plot(t, NSD)
lines(t, predict(mig_mod), col = "red", lty=2, lwd=3) #looks good!

#Exercise 3: Now use nls to fit the other models to our NSD data (mix-mig, disperser, nomadic and resident).
# Hint: you'll need to swap in the equations, replace the parameters, and give some appropriate initial values.
# Some of the models might not fit the data well. If you get an error saying the model has reached the max
# iterations, add the following code to your call to nls(), control = list(maxiter = 500)
# If you get a 'singular gradient' error, you need to choose better initial values.
# For good intial values, try plotting and trying to gauge a likely estimate from the plot.
# (15 mins).
#Ex3 Solution:
####################################################################################################


#####################################################################################################

#Now let's look at the fit of each model to our data
plot(t, NSD)
lines(t, predict(mig_mod), col = "red", lty=2, lwd=3)
lines(t, predict(mix_mod), col = "green", lty=2, lwd=3)
lines(t, predict(dis_mod), col = "blue", lty=2, lwd=3)
lines(t, predict(nom_mod), col = "orange", lty=2, lwd=3)
lines(t, predict(res_mod), col = "purple", lty=2, lwd=3)

#clearly the migrant model is the best, although sometimes depending on the level of noise we add, 
#mixed-migration can also be a good candidate model. At zeta > 0.8, the two models start to look very
#similar.

#To compare these models more rigoroursly, we can calculate the AIC for each model and look for the best
#(i.e. lowest) score
AIC(mig_mod)
AIC(mix_mod)
AIC(dis_mod)
AIC(nom_mod)
AIC(res_mod)

#So that's the process going on behind NSD modelling. Of course, writing scripts to calculate the NSD
#for every point for every animal and then fitting the models and performing model selection for each animal
#individually is arduous. Luckily, Derek Spitz and colleagues have written a great package that automates
#a lot of these tasks, which we will now learn how to use.

#First we'll load in some elk GPS location data, collected for the Wiggins Fork herd in the South-Eastern
#part of the Greater Yellowstone Ecosystem
df_elk1 <- read.csv("elk1_data.csv", header = T, stringsAsFactors = F)
df_elk1$Date_Time <- as.POSIXct(df_elk1$Date_Time, "%Y-%m-%d %H:%M:%S", tz = "MST")

plot(df_elk1$UTM_X, df_elk1$UTM_Y) #make a quick plot of the data
unique(df_elk1$id) #which IDs in file?
unique(df_elk1$Year) #which years in file?

#First we will convert the data in to adehabitat's ltraj object type
xy1 <- as.data.frame(cbind("x" = df_elk1$UTM_X, "y" = df_elk1$UTM_Y)) #needs to be a dataframe
date1 <- df_elk1$Date_Time #needs to be POSIXct class, which we made it in to above
id1 <- paste(df_elk1$id,df_elk1$Year,sep = "_") #combine to make a new id, which is the animal_year format

#Note here we can split our data by year, because our elk migrates in spring and fall, if either
#migration was split across Dec 31st, we might have to use another method.

#We can use as.ltraj to build the ltraj object.
tracks1 <- as.ltraj(xy1, date1, id1, burst = id1, typeII = TRUE,
                    slsp = c("remove", "missing"),
                    infolocs = data.frame(pkey = paste(id1, date1, sep="."),
                                          row.names=row.names(xy1)),
                    proj4string = CRS('+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

plot(tracks1) 

#Notice the data is now split in to seperate 'bursts'. The NSD models are fitted for each
#burst seperately. Do you notice anything else about these three bursts?

#To fit the NSD models, we can now pass the ltraj object as is to the mvmtClass() function, which is a part of
#migrateR. The function will fit the models to all the bursts in the data, without us having to split anything.

#In this instance, the function should be able to fit the models without any additional arguments.
#However, sometimes the models have trouble converging (just as they did when we fit the models manually).
#We can suggest lower, starting and upper values for every parameter using the p.est arguement. e.g.
#elk.nsd1 <- mvmtClass(tracks1, p.est=pEst(l.d = 50, s.d=500, u.d=10000))
#Here I'm setting the lowest delta as 50, the starting value as 500, and the upper delta to 10,000.
#Every parameter can be specified in this way, so check the help file for this function for further details.

#Another set of arguments you might need to use on your own data are stdt and ecut. These set the date
#from which to begin fitting the curve, and are specified in MM-DD format, e.g. stdt "3-1".

#Let's continue fitting the NSD models to our elk data...
elk.nsd1 <- mvmtClass(tracks1)

#migrateR is having trouble fitting the migrant model for 2017, let's take a look at the output to see why...
plot(elk.nsd1)

#Exercise 4: using the output plot of elk.nsd1, what change the call to mvmtClass() 
#so that all models converge (10 minutes).
#Hint: look at the additional arguements above!

#Ex4 Solution:
###################################################################################################

###################################################################################################
#Now that all models have converged, let's take a look at the plots again
plot(elk.nsd1)

#You can see that for years 2015 and 2016, we see a migratory pattern, with outward and return
#migrations. Elk track the green-up of forage across the landscape in spring, so these oscillations in their
#movements are typical.

#For 2017, we see something else. The animal moves to the summer range, and then continues on. The
#AIC ouput in the figure legend indicates that disperser was the model that fit the data best.
#We can also plot these data in 2D with the spatmig() function.

spatmig(tracks1, elk.nsd1)

#Now we see that in 2017, our elk has returned to a different winter range. This winter range switching is
#quite rare for elk in WY, only about 1-3% of migrations result in winter range switching.

#We might also be interested in the parameter estimates, or the estimated migration timings. These results
#could be used in other analyses, such as time-to-event modelling (to see what environmental factors trigger
#migration).

#Exctract migration times
time_df <- mvmt2dt(elk.nsd1, mod = "mixmig") #extracts the timing estimates from our mixed migration models

#Here's some code to tidy the output in to a dataframe for saving
time_df <- ldply(time_df, data.frame) #convert to DF for export, dday is julian day etc...
time_df$dday <- yday(time_df$date) #I often prefer julian days to decimal days
time_df$attr <- rep(c("str1", "end1", "str2", "end2"), each = 1) #label the times with their respective attributes
time_df <- time_df[,c(1,4,2,3)] #reorganize the DF
write.csv(time_df, paste("nsd_times_wiggins",df_elk1$id[1],".csv",sep = "")) #save for later

#Extract parameter estimates
top_mod <- topmvmt(elk.nsd1) #gets the top model (in terms of AIC) for each burst
params_df <- mvmt2df(top_mod) #exctract the parameters in to a list

#Some code to tidy the output for saving
for (i in 1:length(params_df)){
  params_df[[i]]$burst <- rownames(params_df[[i]])
  params_df[[i]] <- params_df[[i]][,order(colnames(params_df[[i]]))]
}
params_df <- do.call(rbind.fill, params_df) #here rbind.fill() is a function from the plyr package
#it performs an rbind and fills missing columns with NA. It's useful here because some parameters
#are specific to the top model for that burst. e.g. nomadic contains only a beta parameter!

write.csv(params_df, paste("nsd_params_wiggins",df_elk1$id[1],".csv",sep = "")) #save for later

#Exersize 5: Using everything you have learned, load in the file XXX, fit some NSD models and observe
#the results. How do they differ to the first animal?
#Is this animal performing a migration every year? How does it compare to the firse animal we analysed? 
#Think about whether this animal has non-overlapping seasonal ranges, distance travelled (e.g. sqrt(delta)).
#(20 mins). 

#Try your best to get all models to converge, but ask if you need a hand.

df_elk2 <- read.csv("elk2_data.csv", header = T, stringsAsFactors = F)
df_elk2$Date_Time <- as.POSIXct(df_elk2$Date_Time, "%Y-%m-%d %H:%M:%S", tz = "MST")

#...the rest is up to you!

#Ex5 Solution:
#############################################################################################################

##############################################################################################################