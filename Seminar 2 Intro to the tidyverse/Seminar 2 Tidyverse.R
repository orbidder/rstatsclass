# Loading library "tidyverse" will load the core packages in the tidyverse
# These include:
#   readr: the package that imports data into tibbles (enhanced dataframes)
#     key commands: read_csv
#   ggplot2: the data vizualization package
#     key commands: ggplot
#   tidyr: helps to create tidy data (whereby each variable is in its own column and each observation is in its own row)
#     key commands: gather, spread
#   dplyr: assists in data manipulation
#     key commands: mutate, filter, summarise, arrange, select, *joins* 
#   purrr: apply functions to lists or vectors
#     key commands: map, nest, unnest

library(tidyverse)
library(lubridate)

# base::read.csv brings in your data as a dataframe
vicuna <- read.csv("Seminar 4 RSFs/vicuna_data_2015.csv")

class(vicuna)
vicuna

# readr::read_csv creates a tibble, which works well with piping other tidyverse commands
vicuna <- read_csv("Seminar 4 RSFs/vicuna_data_2015.csv")

class(vicuna)
vicuna

# Is the base dataframe is harder to review than the tibble? Why?
# Tibbles are enhanced data frames that force your tables to be cleaner
# Review the data types of each column. Do they make sense?

#-----------
#Some basic functions to get to know your data

# How many observations do you have? This command is applied to an entire tibble or dataframe
nrow(vicuna)

# How many variables do you have? This command is also applied to an entire tibble or dataframe
ncol(vicuna)

# How many values are there? This command can be applied to any vector
length(vicuna$X1)
length(vicuna[1,])

# How many of each unique value are there? This command can be applied to one or more vectors
# "table" is a shortcut that can be great to check data for errors
# A little later, we'll also go over how to summarize data in dplyr
table(vicuna$animals_id)
# To vizualize how this can be applied to multiple vectors, let's simulate some binomial data
mutate(vicuna, random = rbinom(nrow(vicuna),1,0.5)) -> vicuna
table(vicuna$animals_id,vicuna$random)
vicuna$random = NULL

# What is the distribution of my data?
# To vizualize this, let's simulate some normally distributed data
mutate(vicuna, normal = rnorm(nrow(vicuna),0,1)) -> vicuna
summary(vicuna$normal)
hist(vicuna$normal)
vicuna$normal = NULL

# What is the structure of my data? This can be applied to vectors or dataframes/tibbles
str(vicuna)
str(vicuna$animals_id)
str(vicuna$acquisition_time)

# It is always a good idea to look at your data and make sure it is formatted correctly before proceeding with analyses!
# Look at histograms of your variables. Do they make sense? Are there weird outliers?
# 

# -----------
# Using tidyr - gather and spread
# We don't often use "gather" and "spread" in animal movement data because we want every observation to count independently
# But, just for good measure, let's go over these commands
# First, we'll create a tibble with three people and their heart rates before and after aerobic exercise
data_hr <- tibble(
  name = c("Wilbur", "Petunia", "Gregory"),
  before = c(67, 80, 64),
  after = c(100, 115, 107)
)
data_hr

# gather allows us to consolidate values into one column and distinguish them by a factor
data_hr %>%
  gather(drug, heartrate, before:after) -> data_hr_tidy

#spread allows us to split a column by a factor
data_hr_tidy %>%
  spread(drug,heartrate)

#-------------
# Time to get to know the dplyr commands! dplyr helps with data manipulation
# Let's go back to our vicuna GPS tibble

# What if we want to only look at data from individual #16?
filter(vicuna, animals_id==16)

# What if we only want to see certain columns from the tibble?
select(vicuna,4:7)
# Or exlude some columns specifically?
select(vicuna,-(1:3))
# We can also rename columns within the "select" command
select(vicuna,4:7, vicunaID = animals_id)
# Or, you can use "rename" to rename without selecting columns
rename(vicuna, vicunaID = animals_id)

# What if we want to sort our data by a variable?
# We can use a "-" to sort in decreasing order
arrange(vicuna,-animals_id)
# We can also sort by multiple variables
arrange(vicuna,animals_id,-gps_data_animals_id)

# Often we want to add new columns based on functions applied to other columns
# For this, we use "mutate". mutate allows us to make multiple new columns in one command
library(geosphere)
mutate(vicuna,hav.distance.to.first = distHaversine(cbind(longitude,latitude),(cbind(longitude[1],latitude[1]))))

# Whelp, we can't see our new column. Time to pipe!
# Piping allows you to do mulitple commands in one...pipe
vicuna %>%
  mutate(hav.distance.to.first = distHaversine(cbind(longitude,latitude),(cbind(longitude[1],latitude[1])))) %>%
  select(c(5:7,12))

# We can also select the columns we want to see using a base approach
mutate(vicuna,hav.distance.to.first = distHaversine(cbind(longitude,latitude),(cbind(longitude[1],latitude[1]))))[,c(5:7,12)]
  
# Or we can make a column entirely using a base approach
vicuna$hav.distance.to.first <- distHaversine(cbind(vicuna$longitude,vicuna$latitude),(cbind(vicuna$longitude[1],vicuna$latitude[1])))
vicuna[,c(5:7,12)]

# Making one new column is pretty much just as easy in base as in dplyr
# The value of dplyr is that you can do a lot at once without making a bunch of new objects
#   and you can visualize columns before you add them
# Here we'll add two columns, distance and a random number from a normal distribution
vicuna %>%
  mutate(hav.distance.to.first = distHaversine(cbind(longitude,latitude),(cbind(longitude[1],latitude[1]))),
         random.number = rnorm(nrow(vicuna),0,1))

# We can choose to save these columns by reassigning the tibble name
vicuna %>%
  mutate(hav.distance.to.first = distHaversine(cbind(longitude,latitude),(cbind(longitude[1],latitude[1]))),
         random.number = rnorm(nrow(vicuna),0,1)) -> vicuna

# Or, we can save a new tibble with our new columns and other specificiations
# We can link up multiple commands all in one pipe!
# For example:
vicuna %>%
  mutate(random.bin = rbinom(nrow(vicuna),1,0.5)) %>%
  filter(animals_id > 14)-> vicuna2

vicuna2

#__________********_________
# Activity!!

# Link up mutate, filter, select, and arrange all in one pipe and save a new tibble



#__________********_________
# Often times we have data from two different sources that we want to combine into one tbl
# This is called joining, whereby we link data from two sources based on a common vector
# If we want to join x to y, we have a number of options
#   inner_join returns all rows from x that are also in y, and all columns from both
#   left_join returns all rows from x, and all columns from both
#     If there is an x with no y match, those rows from y will be filled in with "NA"
#   right_join returns all rows from y, and all columns from both
#     If there is a y with no x match, those rows from x will be filled in with "NA"
#   full_join returns all rows and columns from x and y
#     Any missing values are filled in with NAs

# Let's make up some data:
tibble(animalID = c(seq(1,5),seq(7,10)),
  sex = c(rbinom(9,1,0.5))) -> animal.sex

tibble(animalID = c(seq(1,6),seq(8,10)),
  site = c(rbinom(9,1,0.5))) %>%
  mutate(site.factor = ifelse(site == 1, "A", "B")) %>%
  select(1,3) -> animal.site

inner_join(animal.sex,animal.site)

right_join(animal.sex,animal.site)

left_join(animal.sex,animal.site)

full_join(animal.sex,animal.site)

# But what if our columns are named different things? Then we need to use "by" in our command
tibble(name.of.animal = c(seq(1,5),seq(7,10)),
       sex = c(rbinom(9,1,0.5))) -> animal.sex

inner_join(animal.sex,animal.site, by = c("name.of.animal" = "animalID"))


# How might we integrate a join into our vicuna data?
#__________********_________
# Activity!!
# Generate a tibble with two columns: the first with each vicuna ID, 
#   and the second with some hypothetical value of your choosing
# Then join the data from your new tibble to our vicuna movement data
# Hint: you can use rbinom, runif, or rnorm to generate your fake data
# Hint2: use "unique" to get a vector of vicuna IDs



#__________********_________
# purrr helps to replace for loops by applying functions to subsets of data within a tbl
# purrr relies on "map" commands

vicuna %>%
  split(.$animals_id) %>% # from base R
  map(~ lm(hav.distance.to.first ~ random.number, data = .)) %>%
  map(summary) %>%
  map_dbl("r.squared")

# In comparison to a base R approach:
df<- data.frame(vicID = unique(vicuna$animals_id),
                rsq = 0)
for (i in 1:nrow(df)){
  animalid <- df$vicID[i]
  data <- subset(vicuna,animals_id == animalid)
  model <- lm(hav.distance.to.first ~ random.number, data = data)
  df$rsq[i] <- summary(model)$r.squared
}
df

# In this course, we will mostly use map to deal with splitting our data into individual animals
#   and applying some model or function to animals seperately


#__________********_________
# Time for TIME
# Dealing with times in R can be magical, or it can be a pain in the a**
# It's all about knowing what format and time zone your data are in
#   and knowing how to transform them into the format and time zone you want
# We will be very dependent on the package lubridate for all of our time needs
library(lubridate)

# First, just for fun...what time is it?
now()

# Let's revisit our vicuna data, where "acquisition_time" is the time of the GPS location
vicuna %>%
  select(4:7) -> vicuna

# One of the truly mind-blowing things about lubridate is that it can handle multiple data formats at once
#   As long as the parts of the date are in the right order!
x <- c(20090101, "2009-01-02", "2009 01 03", "2009-1-4",
       "2009-1, 5", "Created on 2009 1 6", "200901 !!! 07")
ymd(x)

# However, when we use read_csv, the tidy gods will automatically read our datetimes in as POSIX objects
#   as long as they are formatted in a way the package can understand
str(vicuna$acquisition_time)

# Before we do ANYTHING with GPS data, we need to confirm the time zone our GPS data are in
tz(vicuna$acquisition_time)

# If the time is right and the time zone is right, leave it!
# If the time is right and the time zone is wrong, use "force_tz" to assign the correct time zone
# If the time is wrong, use "with_tz" to change the time and time zone

# If we assume that our data show the correct time in Argentina, we need to correct the time zone
# To do that, we can use "force_tz"
vicuna %>% 
  mutate(timestamp = force_tz(acquisition_time,tz="America/Argentina/San_Juan")) -> vicuna
vicuna

# The times look the same right? But let's check out the time zones
tz(vicuna$acquisition_time)
tz(vicuna$timestamp)

# And if we look at the vectors alone, you can see the difference in the data
head(vicuna$acquisition_time)
head(vicuna$timestamp)

# What if the collars were programmed in California, but then brought down to Argentina?
# Then we have to first assign the correct time zone to the datetimes (California),
#   and then transform the times to the correct time zone (Argentina)
vicuna %>%
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan")) -> vicunaCA2Arg
vicunaCA2Arg

# Now the times look different, but the time zone for the new timestamp is still accurate
tz(vicunaCA2Arg$acquisition_time)
tz(vicunaCA2Arg$timestamp)

# Now that we feel good about our timezones, we can start manipulating our datetimes
# With datetimes, we can extract a lot of data that is important to us
#   including year, month, week of year, day of year, hour, etc.

# First, let's go through ways of formatting datetimes
# Generally, our movement data will be formatted as ymd_hms objects - year-month-day hour:minute:second
# Once we have our ymd_hms object, we can piece it apart
tail(vicuna)

tail(date(vicuna$timestamp))

tail(month(vicuna$timestamp))

tail(week(vicuna$timestamp))

tail(wday(vicuna$timestamp))

tail(day(vicuna$timestamp))

tail(hour(vicuna$timestamp))

tail(minute(vicuna$timestamp))

tail(am(vicuna$timestamp))
tail(pm(vicuna$timestamp))

# If we want to round our datetimes, we have a couple of options
# Let's say seconds are just confusing us and we want them out
# We can use round_date to get rid of them by rounding to the nearest minute, while keeping the same format
vicuna %>%
  mutate(timestamp = round_date(timestamp, unit = "minute"))
# OR we can get rid of the seconds altogether by reformatting the datetime
vicuna %>%
  mutate(timestamp = format(timestamp,format='%Y-%m-%d %H:%M'))

