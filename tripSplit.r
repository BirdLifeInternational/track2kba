## tripSplit    #####################################################################################################################

## MAIN UPDATE: tidyverse, simple features, CROSSING DATELINE workaround

## Phil Taylor & Mark Miller, 2011

## this script splits central place foraging animal movement data
## into individual trips away from the colony based on distance and time
## away from a defined colony. A distance buffer is set, under which data is
## assumed to be either roosting or device error and is ignored.

## Track must be either a DataFrame or SpatialPointsDataFrame with Latitude, Longitude,
## ID and TrackTime fields
## Colony must be a DataFrame with Latitudes and Longitudes
## InnerBuff is a number indicating the distance in km that must be travelled for the
## movement to be considered a trip
## ReturnBuff is a number indicating the proximity in km that is required for a trip
## to be considered as returning.
## Duration is the length of time, in hours, that the birds must be atlarge for for the
## movement to be considered a trip.
## if plotit = TRUE a map will be drawn.
## if MidPoint = TRUE the calculations will be done projected on the data's centroid.


tripSplit <- function(Track, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, plotit = FALSE, MidPoint = FALSE, nests=FALSE)
  {

  if(!"Latitude" %in% names(Track)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(Track)) stop("Longitude field does not exist")
  if(!"ID" %in% names(Track)) stop("ID field does not exist")
  if(!"TrackTime" %in% names(Track)) stop ("TrackTime field does not exist")

  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")

  if(!(is.double(InnerBuff) & is.double(ReturnBuff))) stop ("InnerBuff and ReturnBuff should be numbers")

  require(sp)
  require(maps)
  require(mapdata)
  require(maptools)
  require(rgdal)
    require(geosphere)

  if(class(Track) != "SpatialPointsDataFrame")
    {
    Track.Wgs <- SpatialPoints(data.frame(Track$Longitude, Track$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    Track.Projected <- spTransform(Track.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", Colony$Longitude[1], " +lat_0=", Colony$Latitude[1], sep="")))
    Track <- SpatialPointsDataFrame(Track.Projected, data = Track)
    }

    ### added section by Steffen Oppel to facilitate nest-specific distance calculations######
if(nests == TRUE)
    {  if(!"ID" %in% names(Colony)) stop("Colony missing ID field")
    nest<- Colony[match(unique(Track$ID), Colony$ID),]
    Colony.Wgs <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    usedCRS<-CRS(proj4string(Track))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=usedCRS)
} else{

if(MidPoint == FALSE)
    {
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    #Colony.Projected <- spTransform(Colony.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", Colony$Longitude, " +lat_0=", Colony$Latitude, sep="")))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=CRS(proj4string(Track)))    ### added 5 June 2017 because midpoint caused problems
        } else
    {
   mid_point<-data.frame(centroid(cbind(Track$Longitude, Track$Latitude)))
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")))
    Track <- spTransform(Track, CRS=CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")))### Transform DataGroup too.
        }
    } 		## ends the else loop for nests=FALSE


  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]

  Track$Returns <- ""
  Track$trip_id <- 0
  Track$ColDist <- spDists(Track, Colony.Projected)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  ReturnBuff <- ReturnBuff * 1000
  InnerBuff <- InnerBuff * 1000

  if(plotit == TRUE)
    {
    plot(Track, pch=1, cex=0.5)
    legend("topleft", paste(Track$ID[1]))
    points(Colony.Projected, pch=18, cex=1.5, col=2)
    }

  i <- 0
  while(i < nrow(Track))
    {
    i <- i + 1
  if(Track$ColDist[i] < InnerBuff) {Track$trip_id[i] <- -1} else
    {
    k <- i
    if(i == nrow(Track)) {Track$trip_id[i] <- -1; break}      ### need to look at how these breaks affect the DataGroup loop
    Dist <- Track$ColDist[i]
    while(Dist >= InnerBuff)
      {
      if(k == nrow(Track) & Dist < ReturnBuff) {break} else
       {
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
