## tripSplit    #####################################################################################################################

## MAIN UPDATE: tidyverse, simple features, CROSSING DATELINE workaround
## included main new wrap-around feature to work across all IDs in a dataset

## Steffen oppel, 5 March 2019; based on work by Phil Taylor and Mark Miller in 2011
## update 6 March - plotting is weird in RStudio, so removed plot and included a panel figure of all trips

## this script splits central place foraging animal movement data
## into individual trips away from the colony based on distance and time
## away from a defined colony. A distance buffer is set, under which data is
## assumed to be either roosting or device error and is ignored.

## tracks must be a DataFrame with Latitude, Longitude, ID and DateTime fields
## Colony must be a DataFrame with Latitudes and Longitudes
## InnerBuff is a number indicating the distance in km that must be travelled for the
## movement to be considered a trip
## ReturnBuff is a number indicating the proximity in km that is required for a trip
## to be considered as returning.
## Duration is the length of time, in hours, that the birds must be at large for for the
## movement to be considered a trip.
## the calculations will be done projected on the data's mean latitude and longitude
## plotit=T will plot the trips of 20 individuals


#### MAIN WRAPPER FUNCTION THAT INCLUDES DATA PREP AND LOOP OVER EACH ID

tripSplit <- function(tracks, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, nests=FALSE, plotit=T)
  {
  
  ## load required packages ##
  require(sp)
  require(maps)
  require(mapdata)
  require(maptools)
  require(rgdal)
  require(geosphere)
  require(ggplot2)
  require(tidyverse)
  
  ## provide error messages ##
  if(!"Latitude" %in% names(tracks)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(tracks)) stop("Longitude field does not exist")
  if(!"ID" %in% names(tracks)) stop("ID field does not exist")
  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")
  if(!(is.double(InnerBuff) & is.double(ReturnBuff))) stop ("InnerBuff and ReturnBuff should be numbers")
  
  ## set required fields
  cleantracks <- tracks %>%
    mutate(DateTime = ymd_hms(DateTime)) %>%   ### needs some clever trick to convert to POSIXct if it isn't already POSIXct
    mutate(TrackTime = as.double(DateTime)) %>%
    mutate(trip_id = ID) %>%
    dplyr::select(ID, trip_id, Latitude, Longitude,DateTime, TrackTime) %>%
    arrange(ID, TrackTime)
  range(cleantracks$Longitude)
  
  ### CREATE PROJECTED DATAFRAME ###
  DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(cleantracks$Longitude, cleantracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = cleantracks, match.ID=F)
  mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
  
  ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
  if (min(cleantracks$Longitude) < -170 &  max(cleantracks$Longitude) > 170) {
    longs=ifelse(cleantracks$Longitude<0,cleantracks$Longitude+360,cleantracks$Longitude)
    mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))}
  
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
  DataGroup.Projected <- spTransform(DataGroup, CRS=proj.UTM)

  
### LOOP OVER EVERY SINGLE ID ###
for(nid in 1:length(unique(tracks$ID))){
  TrackIn <- subset(DataGroup.Projected, ID == unique(DataGroup.Projected$ID)[nid])
  TrackOut<-splitSingleID(Track=TrackIn,Colony=Colony,InnerBuff = InnerBuff, ReturnBuff = ReturnBuff, Duration = Duration, nests=nests, proj.UTM=proj.UTM)
  if(nid == 1) {Trips <- TrackOut} else {Trips <- spRbind(Trips,TrackOut)}
}
  
  
### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE
  if(plotit == TRUE)
    {  
  
    if(length(unique(Trips@data$ID))>25){
      selectIDs<-unique(Trips@data$ID)[1:25]
      plotdat<-  Trips@data %>% filter(ID %in% selectIDs)
      warning("Too many individuals to plot. Only the first 25 ID's will be shown")
      }else{plotdat<-Trips@data}
  
    TRACKPLOT<-plotdat %>% mutate(complete=ifelse(Returns=="N","no","yes")) %>% 
      arrange(ID,TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
      ggplot(aes(x=Longitude, y=Latitude, col=complete)) +
      geom_path() +
      geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
      facet_wrap(~ID) +
      theme(panel.background=element_rect(fill="white", colour="black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="black", fill="white"),
          panel.border = element_blank())
    
    
    ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ###
    if (min(cleantracks$Longitude) < -170 &  max(cleantracks$Longitude) > 170) {
      plotdat<-  Trips@data %>% 
        mutate(Longitude=ifelse(Longitude<0,Longitude+360,Longitude))
      Colony$Longitude<-ifelse(Colony$Longitude<0,Colony$Longitude+360,Colony$Longitude)
      longlimits<-c(min(plotdat$Longitude)-2, max((plotdat$Longitude)+2))
      longbreaks<-round(seq(longlimits[1],longlimits[2],by=10)/10,0)*10
      longlabels<-ifelse(longbreaks>180,longbreaks-360,longbreaks)
      
      TRACKPLOT<-plotdat %>% mutate(complete=ifelse(Returns=="N","no","yes")) %>% 
        arrange(ID,TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
        ggplot(aes(x=Longitude, y=Latitude, col=complete)) +
        geom_path() +
        geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
        scale_x_continuous(limits = longlimits,
                           breaks = longbreaks,
                           labels = longlabels) +
        facet_wrap(~ID) +
        theme(panel.background=element_rect(fill="white", colour="black"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              panel.border = element_blank())}
    
    
    print(TRACKPLOT)
  } ## end plotit=T loop
  
  
return(Trips)
}



#### ACTUAL TRIP SPLIT FUNCTION THAT WORKS ON SINGLE ID ONLY, this reflects the original tripSpllit function as
## originally written by Mark Miller and Phil Taylor
## wrapped in wrapper function above for convenience


splitSingleID <- function(Track, Colony,InnerBuff = 15, ReturnBuff = 45, Duration = 12, nests=FALSE,proj.UTM){

  
  ### facilitate nest-specific distance calculations ###
  if(nests == TRUE)
  {  if(!"ID" %in% names(Colony)) stop("Colony missing ID field")
    nest<- Colony[match(unique(tracks$ID), Colony$ID),]
    Colony.Wgs <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } else{
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } 		## ends the else loop for nests=FALSE
  Colony.Projected <- spTransform(Colony.Wgs, CRS=proj.UTM)
  
  ### set up data to include in output ###
  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]
  
  Track$Returns <- ""
  Track$trip_id <- 0
  Track$ColDist <- spDistsN1(Track, Colony.Projected)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  ReturnBuff <- ReturnBuff * 1000   ### convert from km into UTM units (m)
  InnerBuff <- InnerBuff * 1000   ### convert from km into UTM units (m)
  
  
  # ### plot data (DEPRECATED) ###
  # if(plotit == TRUE)
  # {
  #   plot(Track, pch=1, cex=0.5)
  #   legend("topleft", paste(Track$ID[1]))
  #   points(Colony.Projected, pch=18, cex=1.5, col=2)
  # }
  
  
  ### SPLIT THE DATA INTO DISCRETE TRIPS ###
  i <- 0
  while(i < nrow(Track))
  {
    i <- i + 1
    if(Track$ColDist[i] < InnerBuff) {Track$trip_id[i] <- -1} else {
      k <- i
      if(i == nrow(Track)) {Track$trip_id[i] <- -1; break}      ### need to look at how these breaks affect the DataGroup loop
      Dist <- Track$ColDist[i]
      while(Dist >= InnerBuff)
      {
        if(k == nrow(Track) & Dist < ReturnBuff) {break} else {
          if(k == nrow(Track))
          {
            print(paste("track ", Track$ID[1], Trip.Sequence + 1, " does not return to the colony", sep=""))
            Track$Returns[i:k] <- "N" ; break
          }
        }
        k <- k + 1
        #if(plotit == TRUE){points(Track[k,], col=2, pch=16, cex=0.5)}
        Dist <- Track$ColDist[k]
      }
      Time.Diff <- (Track$TrackTime[k] - Track$TrackTime[i]) / 3600
      Max.Dist <- max(Track$ColDist[i:k])
      if(Time.Diff < Duration |  Max.Dist < InnerBuff)
      {
        Track$trip_id[i:k] <- -1;
        i <- k;
        print(paste("trip ", Track$ID[1], Trip.Sequence + 1, " is too small a trip"))
        next
      }
      Trip.Sequence <- Trip.Sequence + 1
      Track$trip_id[i:k] <- paste(Track$ID[1], Trip.Sequence, sep="")
      i <- k
      print(paste(Track$ID[1], Trip.Sequence, sep=""))
    }
  }
  #if(plotit == TRUE){points(Track, pch=16, cex=0.75, col=as.factor(Track$trip_id))}
  return(Track)
}






