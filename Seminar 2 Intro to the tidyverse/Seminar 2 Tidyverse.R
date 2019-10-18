# Loading library "tidyverse" will load the core packages in the tidyverse
# These include:
#   readr: the package that imports data into tibbles (enhanced dataframes)
#     key commands: read_csv
#   tidyr: helps to create tidy data (whereby each variable is in its own column and each observation is in its own row)
#     key commands: gather, spread
#   dplyr: assists in data manipulation
#     key commands: mutate, filter, summarise, arrange, select, *joins* 
#   purrr: apply functions to lists or vectors
#     key commands: map, nest, unnest
#   ggplot2: the data vizualization package
#     key commands: ggplot
# We will also go over times and dates today with lubridate, which works well with tidy commands

install.packages("tidyverse",dependencies=T)
install.packages("lubridate",dependencies=T)
install.packages("data.table",dependencies=T)
install.packages("geosphere",dependencies=T)
install.packages("RColorBrewer",dependencies=T)
install.packages("viridis",dependencies=T)
install.packages("wesanderson",dependencies=T)

library(tidyverse)
library(lubridate)
library(data.table)

# base::read.csv brings in your data as a dataframe
vicuna <- read.csv("Seminar 4 RSFs/vicuna_data_2015.csv")
class(vicuna)
vicuna

# data.table::fread creates a data.table
# data.table is faster at processing very large datasets, and can also be used in tidy pipes
vicuna <- fread("Seminar 4 RSFs/vicuna_data_2015.csv",sep=",")
class(vicuna)
vicuna
# however, if you have to read in a large file but still want to work with a tibble, 
#   you can just transform the data structure
vicuna <- as_tibble(vicuna)

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

# FILTER
# What if we want to only look at data from individual #16?
filter(vicuna, animals_id==16)

# SELECT
# What if we only want to see certain columns from the tibble?
select(vicuna,4:7)
# Or exlude some columns specifically?
select(vicuna,-(1:3))
# We can also rename columns within the "select" command
select(vicuna,4:7, vicunaID = animals_id)
# Or, you can use "rename" to rename without selecting columns
rename(vicuna, vicunaID = animals_id)

# ARRANGE
# What if we want to sort our data by a variable?
# We can use a "-" to sort in decreasing order
arrange(vicuna,-animals_id)
# We can also sort by multiple variables
arrange(vicuna,animals_id,-gps_data_animals_id)

# MUTATE
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

# SUMMARIZE
# What if we want to know summary information about our data?
summarize(vicuna, mean_dist = mean(hav.distance.to.first))
# We can also use summarize to look at summary info by a factor
summarize(group_by(vicuna, animals_id), mean_dist = mean(hav.distance.to.first))
# Or we can make a new tibble with multiple summary columns
summarize(group_by(vicuna,animals_id), 
          count = n(),
          mean_dist = mean(hav.distance.to.first),
          sd_dist = sd(hav.distance.to.first),
          se_dist = sd(hav.distance.to.first)/sqrt(n()),
          upperse = mean_dist + se_dist,
          lowerse = mean_dist - se_dist) -> summary_dist
summary_dist

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

# First, just for fun...what time is it?
now()

# Let's revisit our vicuna data, where "acquisition_time" is the time of the GPS location
vicuna %>%
  select(4:7) -> vicuna
# Becasue we're only doing one command, this can actually be done without a pipe
# e.g.
# select(vicuna,4:7) -> vicuna
# Here, I'm showing it as a pipe to help get used to the piping syntax

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
mutate(vicuna,timestamp = round_date(timestamp, unit = "minute"))
# OR we can get rid of the seconds altogether by reformatting the datetime
mutate(vicuna,timestamp = format(timestamp,format='%Y-%m-%d %H:%M'))

#__________********_________
# ACTIVITY 1!
# Make a tibble of GPS location data from only vicuna #23 for the month of February at 3 AM
# See if you can do it in one line of code!
# How many rows are there in your new tibble?

# ACTIVITY 2!
# Save a new tibble with added columns for year, week of year, and hour that only has April datetimes

#__________********_________
# Data visualization with ggplot

# ggplot is your best friend for pretty visualizations!
# For our ggplot exploration, we'll use existing datasets available in R
data(mpg)
# ggplot allows you to specify themes and color palettes to improve your aesthetics
# You can also group by a factor to easily plot data from multiple groups on the same plot
# Let's get started!

# First, we need to get familiar with general ggplot syntax
#  "+" is your key to adding lines to your plot code!
# In ggplot we call the data first
# We can either include the x and y in the initial command or as we add visualization types (e.g. points or lines)
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy))
ggplot(data = mpg, aes(x = displ, y = hwy)) + 
  geom_point()

# What if we want to plot by a factor? Again, we can add our variables up top or in geom_point
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class))
ggplot(data = mpg,aes(x = displ, y = hwy, color = class)) + 
  geom_point()

# We can also distinguish a factor by size...
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, size = class))
# ...or shape...
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, shape = class))
# ...or transparency. This makes the most sense to use if you have an ordinal factor
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, alpha = class))

# To change the color off all points rather than by factor, take the color command out of the aes statement
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy), color = "red")

# Once we start using shapes other than the default, we have other commands we can consider
# For example, if we use a shape with a border, we can use "stroke" to specify the border width
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point(shape = 21, color = "darkorange4", fill = "#EEE8CD", size = 3, stroke = 2)
# *** Try modifying the border color, fill, size, and stroke

# You can also plot based on a factor with facet wrap
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) + 
  facet_wrap(~ class, nrow = 2)

# **COLOR BREAK***
# A note on color!
# In r, you can use color names (e.g. "red") or hexadecimal RGB triplets (e.g. "#FFCC00")
# There are also a number of color palettes available to you, through base, RColorBrewer,
#   Viridis, or, my favorite, wesanderson
# palettes can be one of three types:
#   1. sequential: gradients that indicate high-to-low variables (e.g. temperature)
#   2. diverging: gradients that diverge from a middle neutral point (e.g. temperature anomaly)
#   3. qualitative: unrelated colors for nominal or categorical data (e.g. dog breed)

# Let's try some palettes!
library(RColorBrewer)
display.brewer.all()
library(viridis)
library(wesanderson)

#First, categorical/qualitative data
# Base
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class))
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class)) +
  scale_color_manual(values = c("blue", "red", "green","purple","black","chocolate1","coral3"))
# R color brewer
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class)) +
  scale_color_brewer(palette = "Accent")
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class)) +
  scale_color_brewer(palette = "Pastel2")
# wesanderson
#   these palettes are so fun, but most only have 4-5 colors in them
#   so to illustrate how to use wesanderson, we'll use a factor with fewer options, cyl
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = factor(cyl))) +
  scale_color_manual(values = wes_palette("Moonrise2", n = 4))
# Viridis
#   these palettes are inherently continuous, but can be used by adding discrete = T
#   generally, it can be misleading to use a continuous color palette for discrete data
#   for categorial data, base/RColorBrewer/wesanderson will be better bets
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = class)) +
  scale_color_viridis(discrete = T, option = "inferno")

#Now, continuous data
# Base
#   in base, we have a number of palettes available including heat.colors, terrain.colors, and rainbow
#   and custom gradients
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty))
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_gradientn(colors = terrain.colors(10))
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_gradient(low = "lightblue", high = "darkgreen")
# RColorBrewer
#   diverging
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_gradientn(colors = brewer.pal(n = 3, name = "RdBu"))
#   sequential
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_gradientn(colors = brewer.pal(n = 3, name = "OrRd"))
#wesanderson
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous"))
#Viridis
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  scale_color_viridis(discrete = F, option = "magma", direction = -1)

# Lastly, let's talk themes
# There are a lot of manual modifications you can do in ggplot in regard to axes, labels, background, etc.
# But themes make it a lot more efficient to make the kind of plot you want
# Here are a few to check out: theme_bw, theme_classic, theme_dark, theme_minimal
#   Replace the theme command below with these to see how they look!
ggplot(data = mpg) + 
  geom_point(aes(x = displ, y = hwy, color = cty)) +
  theme_classic()

# That's it for colors and themes for now! It's a bottomless hole you may never escape, 
#   but it's also a lot of fun
# Things we didn't cover today: labels, legends, axes/scales 
# If you are interested in learning more about customization, you can find a ton of information
#   at https://ggplot2.tidyverse.org/index.html
#-------***----------

# We've been working with point data so far, but we might also want to plot summarized data

# One thing you might want to represent is how many records you have of different kinds of observations
# This can be done by only setting an x and not a y with the geom_bar command
ggplot(data = mpg, aes(x = class)) + 
  geom_bar()
# We can also look at how this breaks down by a second factor
ggplot(data = mpg) + 
  geom_bar(aes(x = class, fill = factor(year)), position = "dodge")

# How are our data distributed within factor values? There are a few easy ways to look at this
# Boxplots
ggplot(data = mpg) + 
  geom_boxplot(aes(x = class, y = hwy))
# Violin plots
ggplot(data = mpg) + 
  geom_violin(aes(x = class, y = hwy))
# Violin PLUS data mean with standard error
ggplot(data = mpg,aes(x = class, y = hwy)) + 
  geom_violin() + 
  stat_summary()
# We can also use geom_bar to illustrate means and standard errors, 
#   but we have to calculate those values first

# *** ACTIVITY ***
# Use "summarize" to calculate the mean and standard errors of hwy ~ class from the mpg dataset
# Make a new tibble with these values called hwyplot

# Now we can make a barplot with means and error bars
# To avoid geom_bar counting our data, we have to tell it to use the actual y value using stat = "identity"
# Then we can add error bars using geom_errorbar
ggplot(data = hwyplot,aes(x = class, y = meanhwy)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = meanhwy - sehwy, ymax = meanhwy + sehwy), width=0.2)

# Bonus! Change the color palette of your plot to an RColorBrewer palette
#   Change the theme to your theme of choice!
# Hint: because we are using shapes (bars) rather than points/lines, use "fill" instead of "color" in the aes
#   and use "scale_fill_brewer" rather than "scale_color_brewer"
