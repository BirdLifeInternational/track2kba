## tripSplit    #####################################################################################################################

## MAIN UPDATE: tidyverse, simple features, CROSSING DATELINE workaround
## included main new wrap-around feature to work across all IDs in a dataset

## Steffen oppel, 5 March 2019; based on work by Phil Taylor and Mark Miller in 2011

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
## if plotit = TRUE a map will be drawn.
## the calculations will be done projected on the data's mean latitude and longitude


#### MAIN WRAPPER FUNCTION THAT INCLUDES DATA PREP AND LOOP OVER EACH ID

tripSplit <- function(tracks, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, nests=FALSE,plotit = FALSE)
  {
  
  ## load required packages ##
  require(sp)
  require(maps)
  require(mapdata)
  require(maptools)
  require(rgdal)
  require(geosphere)
  
  ## provide error messages ##
  if(!"Latitude" %in% names(tracks)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(tracks)) stop("Longitude field does not exist")
  if(!"ID" %in% names(tracks)) stop("ID field does not exist")
  if(!"TrackTime" %in% names(tracks)) stop ("TrackTime field does not exist")
  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")
  if(!(is.double(InnerBuff) & is.double(ReturnBuff))) stop ("InnerBuff and ReturnBuff should be numbers")
  
  ## set required fields
  tracks <- tracks %>%
    #mutate(DateTime = dmy_hms(paste(DateGMT,TimeGMT, sep = " "))) %>%   ### needs some clever trick to convert to POSIXct if it isn't already POSIXct
    mutate(TrackTime = as.double(DateTime)) %>%
    mutate(trip_id = ID) %>%
    dplyr::select(ID, trip_id, Latitude, Longitude,DateTime, TrackTime)
  
  
  ### CREATE PROJECTED DATAFRAME ###
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mean(tracks$Longitude), " +lat_0=", mean(tracks$Latitude), sep=""))
  DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = tracks, match.ID=F)
  DataGroup.Projected <- spTransform(DataGroup, CRS=proj.UTM)


  
  
### LOOP OVER EVERY SINGLE ID ###
for(nid in 1:length(unique(tracks$ID))){
  TrackIn <- subset(DataGroup.Projected, ID == unique(DataGroup.Projected$ID)[nid])
  TrackOut<-splitSingleID(Track=TrackIn,InnerBuff = InnerBuff, ReturnBuff = ReturnBuff, Duration = Duration, nests=nests,plotit = plotit)
  if(nid == 1) {Trips <- TrackOut} else {Trips <- spRbind(Trips,TrackOut)}
  }
return(Trips)
}



#### ACTUAL TRIP SPLIT FUNCTION THAT WORKS ON SINGLE ID ONLY, this reflects the original tripSpllit function as
## originally written by Mark Miller and Phil Taylor
## wrapped in wrapper function above for convenience


splitSingleID <- function(Track, InnerBuff = 15, ReturnBuff = 45, Duration = 12, nests=FALSE,plotit = FALSE){

  
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
  Track$ColDist <- spDists(Track, Colony.Projected)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  ReturnBuff <- ReturnBuff * 1000   ### convert from km into UTM units (m)
  InnerBuff <- InnerBuff * 1000   ### convert from km into UTM units (m)
  
  
  ### plot data (OPTIONAL) ###
  if(plotit == TRUE)
  {
    plot(Track, pch=1, cex=0.5)
    legend("topleft", paste(Track$ID[1]))
    points(Colony.Projected, pch=18, cex=1.5, col=2)
  }
  
  
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
        if(plotit == TRUE){points(Track[k,], col=2, pch=16, cex=0.5)}
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
  if(plotit == TRUE){points(Track, pch=16, cex=0.75, col=as.factor(Track$trip_id))}
  return(Track)
}

