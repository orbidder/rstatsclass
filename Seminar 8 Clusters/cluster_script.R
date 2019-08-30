cluster_centroids_func <- function(df, stime, sdist, fixr, timezone1, timezone2){
  
  # find all the points that are fewer than stime days apart
  timearray <- array(0,c(nrow(df),nrow(df)))
  timedif<- array(0,c(nrow(df),nrow(df)))
  for (i in 1:nrow(df)){
    for (j in 1:nrow(df)){
      timearray[i,j] = abs(as.numeric(difftime(df$timestamp[j],df$timestamp[i],units="days")))
      timedif[i,j] = ifelse(timearray[i,j]<=stime,1,0)                        
    }
  }
  
  # create distance matrix of all points
  df <- st_as_sf(df)
  distance <- st_distance(df, df, by_element = FALSE, which = "Euclidean")
  df$x<-st_coordinates(df[1])[,1]
  df$y<-st_coordinates(df[1])[,2]
  
  # empty distance array
  D <- array(0,c(nrow(df),nrow(df)))
  
  # empty cluster array
  C <- array(0,c(nrow(df),nrow(df)))
  
  # populate cluster array based on time and distance thresholds.
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
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
  centroidtemp <- array(0,c(nrow(df),2))
  
  # empty centroid array
  #centroid <- array(NA,c(nrow(df),12))
  N <- nrow(df)
  centroid <- data.frame(x=double(N),
                         y=double(N),
                         time=double(N),
                         night=integer(N),
                         actratio=double(N),
                         nightratio=double(N),
                         points=integer(N),
                         first=.POSIXct(double(N)),
                         last=.POSIXct(double(N)),
                         bin=integer(N),
                         clustID=integer(N),
                         stringsAsFactors=FALSE)
  
  # Next we create an empty array of point ids associated with each centroid
  # This array will have a row for each centroid with all the point IDs associated with it
  # 400 is just an arbitrary large number well over the likely number of points in a cluster
  ID=array(0,c(nrow(df),400))
  
  # used to ignore points already associated with a centroids
  excess=numeric(nrow(df))
  
  # Populate centroid and id arrays
  # For our centroid data frame, we calculate cluster characteristics like duration and night ratio
  b=1 #counter
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
      if((excess[i]!=1)&(excess[j]!=1)){
        if(C[i,j]==1){
          a=1
          excess[i]=1
          excess[j]=1
          centroidtemp[a,1]=(df$x[j]+df$x[i])/2
          centroidtemp[a,2]=(df$y[j]+df$y[i])/2
          ID[b,c(a,a+1)]=c(i,j) #fill in first two points of centroid
          centroid[b,1]=centroidtemp[a,1]
          centroid[b,2]=centroidtemp[a,2]
          for(k in 1:nrow(df)){
            if(((timedif[i,k]==1)&(k!=i)&(k!=j))|((timedif[max(ID[b,]),k]==1)&(k!=j)&(k!=i))){
              if(sqrt((centroidtemp[a,1]-df$x[k])^2+(centroidtemp[a,2]-df$y[k])^2)<sdist){
                a=a+1
                centroidtemp[a,1]=(a*centroidtemp[(a-1),1]+df$x[k])/(a+1)#weighted average
                centroidtemp[a,2]=(a*centroidtemp[(a-1),2]+df$y[k])/(a+1)
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
            if(df$daynight[ID[b,l]]==0){
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
          centroid[b,8]=as.POSIXct(df$timestamp[min(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #first fix
          centroid[b,9]=as.POSIXct(df$timestamp[max(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #last fix
          centroid[b,3]=as.numeric(difftime(centroid[b,9],centroid[b,8],units="days")) #total time
          if(centroid[b,3]>1){
            centroid[b,10]=1
          } 
          else{
            centroid[b,10]=0 #binary variable
          } 
          centroid[b,11]=b #cluster ID number
          b=b+1      
        }
      }
    }
  }
  centroid$first<-with_tz(centroid$first,tzone=timezone2)
  centroid$last<-with_tz(centroid$last,tzone=timezone2)
  
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
            centroid[i,1]=mean(df$x[union])
            centroid[i,2]=mean(df$y[union])
            centroid[j,c(1:7,10:11)]=999 #used to identify and eliminate old cluster j that has been combined with cluster i
            cluster=nnzero(ID[i,])
            z=0
            for(l in 1:(nnzero(ID[i,]))){
              if(df$daynight[ID[i,l]]==0){
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
            centroid[i,8]=as.POSIXct(df$timestamp[min(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #first fix
            centroid[i,9]=as.POSIXct(df$timestamp[max(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #last fix
            centroid[i,3]=as.numeric(difftime(centroid[i,9],centroid[i,8],units="days")) #total time
            if(centroid[i,3]>1){
              centroid[i,10]=1}
            else{
              centroid[i,10]=0} #binary variable
            centroid[i,11]=i
          }
        }
      }
    }
  }
  centroid[centroid == 999] <- NA
  centroid<-na.omit(centroid)
  return(centroid)
}



cluster_points_func <- function(df, stime, sdist, fixr, timezone1, timezone2){

  # find all the points that are fewer than stime days apart
  timearray <- array(0,c(nrow(df),nrow(df)))
  timedif<- array(0,c(nrow(df),nrow(df)))
  for (i in 1:nrow(df)){
    for (j in 1:nrow(df)){
      timearray[i,j] = abs(as.numeric(difftime(df$timestamp[j],df$timestamp[i],units="days")))
      timedif[i,j] = ifelse(timearray[i,j]<=stime,1,0)                        
    }
  }
  
  # create distance matrix of all points
  df <- st_as_sf(df)
  distance <- st_distance(df, df, by_element = FALSE, which = "Euclidean")
  df$x<-st_coordinates(df[1])[,1]
  df$y<-st_coordinates(df[1])[,2]
  
  # empty distance array
  D <- array(0,c(nrow(df),nrow(df)))
  
  # empty cluster array
  C <- array(0,c(nrow(df),nrow(df)))
  
  # populate cluster array based on time and distance thresholds.
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
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
  centroidtemp <- array(0,c(nrow(df),2))
  
  # empty centroid array
  #centroid <- array(NA,c(nrow(df),12))
  N <- nrow(df)
  centroid <- data.frame(x=double(N),
                         y=double(N),
                         time=double(N),
                         night=integer(N),
                         actratio=double(N),
                         nightratio=double(N),
                         points=integer(N),
                         first=.POSIXct(double(N)),
                         last=.POSIXct(double(N)),
                         bin=integer(N),
                         clustID=integer(N),
                         stringsAsFactors=FALSE)
  
  # Next we create an empty array of point ids associated with each centroid
  # This array will have a row for each centroid with all the point IDs associated with it
  # 400 is just an arbitrary large number well over the likely number of points in a cluster
  ID=array(0,c(nrow(df),400))
  
  # used to ignore points already associated with a centroids
  excess=numeric(nrow(df))
  
  # Populate centroid and id arrays
  # For our centroid data frame, we calculate cluster characteristics like duration and night ratio
  b=1 #counter
  for(i in 1:nrow(df)){
    for(j in 1:nrow(df)){
      if((excess[i]!=1)&(excess[j]!=1)){
        if(C[i,j]==1){
          a=1
          excess[i]=1
          excess[j]=1
          centroidtemp[a,1]=(df$x[j]+df$x[i])/2
          centroidtemp[a,2]=(df$y[j]+df$y[i])/2
          ID[b,c(a,a+1)]=c(i,j) #fill in first two points of centroid
          centroid[b,1]=centroidtemp[a,1]
          centroid[b,2]=centroidtemp[a,2]
          for(k in 1:nrow(df)){
            if(((C[i,k]==1)&(k!=i)&(k!=j))|((C[max(ID[b,]),k]==1)&(k!=j)&(k!=i))){
                a=a+1
                centroidtemp[a,1]=(a*centroidtemp[(a-1),1]+df$x[k])/(a+1)#weighted average
                centroidtemp[a,2]=(a*centroidtemp[(a-1),2]+df$y[k])/(a+1)
                centroid[b,1]=as.numeric(centroidtemp[a,1])
                centroid[b,2]=as.numeric(centroidtemp[a,2])
                ID[b,a+1]=k
                excess[k]=1
            }
          }
          cluster=nnzero(ID[b,])
          z=0
          for(l in 1:(nnzero(ID[b,]))){
            if(df$daynight[ID[b,l]]==0){
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
          centroid[b,8]=as.POSIXct(df$timestamp[min(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #first fix
          centroid[b,9]=as.POSIXct(df$timestamp[max(ID[b,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #last fix
          centroid[b,3]=as.numeric(difftime(centroid[b,9],centroid[b,8],units="days")) #total time
          if(centroid[b,3]>1){
            centroid[b,10]=1
          } 
          else{
            centroid[b,10]=0 #binary variable
          } 
          centroid[b,11]=b #cluster ID number
          b=b+1      
        }
      }
    }
  }
  centroid$first<-with_tz(centroid$first,tzone=timezone2)
  centroid$last<-with_tz(centroid$last,tzone=timezone2)
  
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
            centroid[i,1]=mean(df$x[union])
            centroid[i,2]=mean(df$y[union])
            centroid[j,c(1:7,10:11)]=999 #used to identify and eliminate old cluster j that has been combined with cluster i
            cluster=nnzero(ID[i,])
            z=0
            for(l in 1:(nnzero(ID[i,]))){
              if(df$daynight[ID[i,l]]==0){
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
            centroid[i,8]=as.POSIXct(df$timestamp[min(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #first fix
            centroid[i,9]=as.POSIXct(df$timestamp[max(ID[i,1:cluster])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #last fix
            centroid[i,3]=as.numeric(difftime(centroid[i,9],centroid[i,8],units="days")) #total time
            if(centroid[i,3]>1){
              centroid[i,10]=1}
            else{
              centroid[i,10]=0} #binary variable
            centroid[i,11]=i
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
        IDcentroid[p,2]=df$x[(ID[i,j])] #X coordinate
        IDcentroid[p,3]=df$y[(ID[i,j])] #Y coordinate
        IDcentroid[p,4]=df$daynight[(ID[i,j])] #Daynight
        IDcentroid[p,5]=as.POSIXct(df$timestamp[(ID[i,j])],format="YYYY-MM-DD hh:mm:ss",tz=timezone1) #time
      }
    }
  }
  IDcentroid$datetime<-with_tz(IDcentroid$datetime,tzone=timezone2)
  return(IDcentroid)
}
