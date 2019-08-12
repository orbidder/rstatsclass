library(maptools)
library(rgdal)
library(tidyverse)
library(raster)
library(lme4)
library(sf)
library(lubridate)
library(rgeos)
library(StreamMetabolism)

#We are selecting only Puma ID = 4 today to examine one puma's clusters
puma4 <- read_csv("/Users/justinesmith/Documents/UCB/Data/Puma_data/puma_data 3hr2.csv") %>% 
  dplyr::select(ID = animals_id, 5:7) %>% 
  filter(ID == 4) %>%
  st_as_sf(., coords = 3:4, crs = "+init=epsg:4326") %>% 
  st_transform("+init=epsg:32719") %>% 
  mutate(timestamp = lubridate::with_tz(ymd_hms(acquisition_time,tz="America/Los_Angeles"),"America/Argentina/San_Juan"))

#Determine sunrise and sunset times by day. The lat/longs should be chosen based on a representative location in your study area
#Determine if each location is between sunrise and sunset (day = 1) or not (night = 0)
puma4 %>%
  rowwise() %>% 
  mutate(sunrise = as.POSIXct(sunrise.set(-29.18351, -69.34985, timestamp)[1,1],origin="1970-01-01",format="%H:%M:%S",force_tz("America/Argentina/San_Juan")),
         sunset = as.POSIXct(sunrise.set(-29.18351, -69.34985, timestamp)[1,2],origin="1970-01-01",format="%H:%M:%S",force_tz("America/Argentina/San_Juan")),
         daynight = ifelse(timestamp<sunrise|timestamp>sunset,0,1)) -> puma4

#set our time (days) and distance (meters) windows
stime <- 1.33
sdist <- 20
fixr <- 3

# find all the points that are fewer than stime days apart
timearray <- array(0,c(nrow(puma4),nrow(puma4)))
timedif<- array(0,c(nrow(puma4),nrow(puma4)))
for (i in 1:nrow(puma4)){
  for (j in 1:nrow(puma4)){
    timearray[i,j] = as.numeric(difftime(puma4$timestamp[j],puma4$timestamp[i],units="days"))
    timedif[i,j] = ifelse(timearray[i,j]<=stime,1,0)                        
  }
}

# create distance matrix of all points
puma4 <- st_as_sf(puma4)
distance <- st_distance(puma4, puma4, by_element = FALSE, which = "Euclidean")
puma4$x<-st_coordinates(puma4[1])[,1]
puma4$y<-st_coordinates(puma4[1])[,2]

# empty distance array
D <- array(0,c(nrow(puma4),nrow(puma4)))

# empty cluster array
C <- array(0,c(nrow(puma4),nrow(puma4)))

# populate cluster array based on time and distance thresholds.
for(i in 1:nrow(puma4)){
    for(j in 1:nrow(puma4)){
        if(timedif[i,j]==1){
            D[i,j]=distance[i,j]}
        else{
            D[i,j]=0}
        if(D[i,j]>0 & D[i,j]<sdist){
            C[i,j]=1}
    }
}

# create centroids from clusters...
# holder array for centroids
centroidtemp <- array(0,c(nrow(puma4),2))

# empty centroid array
#centroid <- array(NA,c(nrow(puma4),12))
N <- nrow(puma4)
centroid <- data.frame(x=double(N),
                 y=double(N),
                 time=double(N),
                 night=integer(N),
                 actratio=double(N),
                 nightratio=double(N),
                 points=integer(N),
                 first=.POSIXct(double(N)),
                 last=.POSIXct(double(N)),
                 adj=double(N),
                 bin=integer(N),
                 fix.success=double(N),
                 clustID=integer(N),
                 stringsAsFactors=FALSE)

# Next we create an empty array of point ids associated with each centroid
# This array will have a row for each centroid with all the point IDs associated with it
# 400 is just an arbitrary large number well over the likely number of points in a cluster
ID=array(0,c(nrow(puma4),400))

# used to ignore points already associated with a centroids
excess=numeric(nrow(puma4))

# Populate centroid and id arrays
# For our centroid data frame, we calculate cluster characteristics like duration and night ratio
b=1 #counter
for(i in 1:nrow(puma4)){
    for(j in 1:nrow(puma4)){
        if((excess[i]!=1)&(excess[j]!=1)){
            if(C[i,j]==1){
                a=1
                excess[i]=1
                excess[j]=1
                centroidtemp[a,1]=(puma4$x[j]+puma4$x[i])/2
                centroidtemp[a,2]=(puma4$y[j]+puma4$y[i])/2
                ID[b,c(a,a+1)]=c(i,j) #fill in first two points of centroid
                centroid[b,1]=centroidtemp[a,1]
                centroid[b,2]=centroidtemp[a,2]
                for(k in 1:nrow(puma4)){
                    if(((timedif[i,k]==1)&(k!=i)&(k!=j))|((timedif[max(ID[b,]),k]==1)&(k!=j)&(k!=i))){
                        if(sqrt((centroidtemp[a,1]-puma4$x[k])^2+(centroidtemp[a,2]-puma4$y[k])^2)<sdist){
                            a=a+1
                            centroidtemp[a,1]=(a*centroidtemp[(a-1),1]+puma4$x[k])/(a+1)#weighted average
                            centroidtemp[a,2]=(a*centroidtemp[(a-1),2]+puma4$y[k])/(a+1)
                            centroid[b,1]=as.numeric(centroidtemp[a,1])
                            centroid[b,2]=as.numeric(centroidtemp[a,2])
                            ID[b,a+1]=k
                            excess[k]=1
                        }
                    }
                }
				cluster=nnzero(ID[b,])
				z=0
				for(l in 1:(nnzero(ID[b,]))){
				  if(puma4$daynight[ID[b,l]]==0){
					z=z+1
					centroid[b,4]=z
					} 
					else{
					centroid[b,4]=z
					} #total number of night points
				}
				centroid[b,5]=(cluster)/(max(ID[b,1:cluster])-min(ID[b,1:cluster])+1) #ratio of active to total cluster points
				centroid[b,6]=centroid[b,4]/(cluster) #ratio of night to total cluster points
				centroid[b,7]=cluster #total number of points in cluster
				centroid[b,8]=as.POSIXct(puma4$timestamp[min(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz="America/Los_Angeles") #first fix
				centroid[b,9]=as.POSIXct(puma4$timestamp[max(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz="America/Los_Angeles") #last fix
				centroid[b,3]=as.numeric(difftime(centroid[b,9],centroid[b,8],units="days")) #total time
				centroid[b,10]=centroid[b,7]/((max(ID[b,1:cluster])-min(ID[b,1:cluster])+1)/((round(centroid[b,3]*24)/fixr)+1)) #adjusted fix 
				if(centroid[b,3]>1){
					centroid[b,11]=1
				} 
				else{
					centroid[b,11]=0 #binary variable
				} 
				centroid[b,12]=round((centroid[b,7]/centroid[b,5])/(centroid[b,3]*8+1),digits=2) #fix success over cluster period
				centroid[b,13]=b #cluster ID number
				b=b+1      
			}
		}
	}
}
centroid$first<-with_tz(centroid$first,tzone="America/Argentina/San_Juan")
centroid$last<-with_tz(centroid$last,tzone="America/Argentina/San_Juan")

# find and combine clusters that have overlapping points
overlap <- 1
L <- nnzero(ID[,1])
for(i in 1:L){
	for(j in 1:L){
		if((centroid[j,1]!=centroid[i,1])&(centroid[j,2]!=centroid[i,2])){
			total=nnzero(ID[i,])
			if(total>0){
				intersect=intersect(ID[i,],ID[j,]) #identify any shared points
				union=sort(union(ID[i,],ID[j,]),decreasing=FALSE)
				union=union[-1]
				if((length(intersect)-1)>=overlap){#always has 0 as intersect)
					ID[i,1:length(union)]=union
					ID[j,]=0
					centroid[i,1]=mean(puma4$x[union])
					centroid[i,2]=mean(puma4$y[union])
					centroid[j,c(1:7,10:12)]=999 #used to identify and eliminate old cluster j that has been combined with cluster i
					cluster=nnzero(ID[i,])
					z=0
					for(l in 1:(nnzero(ID[i,]))){
						if(puma4$daynight[ID[i,l]]==0){
							z=z+1
			  				centroid[i,4]=z
							} 
						else{
							centroid[i,4]=z
							} #total number of night points
								}
					centroid[i,5]=(cluster)/(max(ID[i,1:cluster])-min(ID[i,1:cluster])+1) #ratio of active to total cluster points
					centroid[i,6]=centroid[i,4]/(cluster) #ratio of night to total cluster points
					centroid[i,7]=cluster #total number of points in cluster
					centroid[i,8]=as.POSIXct(puma4$timestamp[min(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz="America/Los_Angeles") #first fix
					centroid[i,9]=as.POSIXct(puma4$timestamp[max(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz="America/Los_Angeles") #last fix
					centroid[i,3]=as.numeric(difftime(centroid[i,9],centroid[i,8],units="days")) #total time
					centroid[i,10]=centroid[i,7]/((max(ID[i,1:cluster])-min(ID[i,1:cluster])+1)/((round(centroid[i,3]*24)/fixr)+1)) #adjusted fix 
					if(centroid[i,3]>1){
						centroid[i,11]=1}
					else{
						centroid[i,11]=0} #binary variable
					centroid[i,12]=round((centroid[i,7]/centroid[i,5])/(centroid[i,3]*8+1),digits=2)
					centroid[i,13]=i
				}
			}
		}
	}
}
centroid[centroid == 999] <- NA
centroid<-na.omit(centroid)

#create ID array with x and y data, cluster ID, day or night
N=nnzero(ID)
IDcentroid <- data.frame(ID=double(N),
                       x=double(N),
                       y=double(N),
                       daynight=integer(N),
                       datetime=.POSIXct(double(N)),
                       stringsAsFactors=FALSE)
p=0
for(i in 1:nrow(ID)){
	for(j in 1:ncol(ID)){
		if(ID[i,j]!=0){
			p=p+1      
			IDcentroid[p,1]=i #centroid ID
			IDcentroid[p,2]=puma4$x[(ID[i,j])] #X coordinate
			IDcentroid[p,3]=puma4$y[(ID[i,j])] #Y coordinate
			IDcentroid[p,4]=puma4$daynight[(ID[i,j])] #Daynight
			IDcentroid[p,5]=as.POSIXct(puma4$timestamp[(ID[i,j])],format="YYYY-MM-DD hh:mm:ss",tz="America/Los_Angeles") #time
		}
	}
}
IDcentroid$datetime<-with_tz(IDcentroid$datetime,tzone="America/Argentina/San_Juan")

write.csv(centroid, "centroids.csv", row.names=FALSE)
write.csv(IDcentroid, "centroidIds.csv", row.names=FALSE)
