## tripSplit    #####################################################################################################################

#' Split tracking data into trips
#'
#' \code{tripSplit} employs \code{splitSingleID} to split data from multiple individuals' into discrete trips made from centrally-located places.
#'
#' This function splits central place foraging animal movement data into individual trips away from a central location based on distance and time.
#'
#' \emph{nests}=T may be used if it is desired, for example, to use specific nest locations instead of one central location for all individuals/dataGroup.
#'
#' @param dataGroup data.frame. Must contain 'Latitude', 'Longitude', 'ID' and 'DateTime' columns (correct format may be assured using \code{\link{formatFields}} function).
#' @param colony data.frame. Containing 'Latitude' and 'Longitude' fields specifying the central location(s) from which trips begin. If data are from MoveBank this information may be extracted using the \code{\link{move2KBA}} function. If \emph{nests}=T, each row should correspond to an appropriate location (Lat/Lon) for each ID value in \emph{dataGroup}.
#' @param innerBuff numeric (in kilometers). Indicate the distance that an animal must travel for the movement to be considered a trip. Note that this is also the metric that determines whether two subsequent trips are split - if your animal records locations > \code{innerBuff} (km) from its place of origin and no locations at the place of origin (e.g. for burrow-nesting species) then subsequent trips may be lumped into a single trip. Increase \code{innerBuff} to ensure correct splitting of trips. 
#' @param returnBuff numeric (in kilometers). Indicate the proximity required for a trip to be considered as returning. This is useful for identifying incomplete trips (i.e. where data storage/transmission failed during the trip).
#' @param duration numeric (in hours). The period of time that the animals must be at large for the movement to be considered a trip.
#' @param nests logical scalar (TRUE/FALSE). Should the central place used in trip-splitting be specific to each ID? If so, each place must be matched with an 'ID' value in both \emph{dataGroup} and \emph{colony} objects.
#' @param rmNonTrip logical scalar (TRUE/FALSE). Should periods not associated with trips be filtered out? Note that this does not filter out the trip starting and ending points which fall within innerBuff, to allow more accurate calculations of duration and distance covered with \code{tripSummary}. Default is TRUE.
#' @return Returns an un-projected (WGS84) SpatialPointsDataFrame, with the field 'trip_id' added to identify each unique trip-ID combination. If rmNonTrip=TRUE, then output has been filtered of points deemed not associated with trip movements.
#'
#' @seealso \code{tripSummary}, \code{mapTrips}
#'
#' @examples
#' \dontrun{Trips <- tripSplit(dataGroup, 
#' colony=colony, 
#' innerBuff=2, 
#' returnBuff=20, 
#' duration=1, 
#' nests = F, 
#' rmNonTrip = T)}
#' @export
#' @import sp
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats median

#### MAIN WRAPPER FUNCTION THAT INCLUDES DATA PREP AND LOOP OVER EACH ID

tripSplit <- function(dataGroup, colony, innerBuff = NULL, returnBuff = NULL, duration = NULL, nests=FALSE, rmNonTrip=TRUE)
{
  dataGroup <- dataGroup %>%
      mutate(DateTime = lubridate::ymd_hms(.data$DateTime)) %>% 
      mutate(TrackTime = as.double(.data$DateTime)) %>%
      mutate(trip_id = .data$ID) %>%
      arrange(.data$ID, .data$TrackTime)

  ### CREATE PROJECTED DATAFRAME ###
  DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(dataGroup$Longitude, dataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = dataGroup, match.ID=FALSE)

  DataGroup.Projected <- DataGroup

  ### LOOP OVER EVERY SINGLE ID ###
  for(nid in seq_len(length(unique(dataGroup$ID)))){
    TrackIn <- base::subset(DataGroup.Projected, DataGroup.Projected$ID == unique(DataGroup.Projected$ID)[nid])
    TrackOut <- splitSingleID(Track=TrackIn, colony=colony, innerBuff = innerBuff, returnBuff = returnBuff, duration = duration, nests=nests)
    if(nid == 1) {Trips <- TrackOut} else {Trips <- maptools::spRbind(Trips, TrackOut)}
  }

  if(rmNonTrip==TRUE) {
    Trips <- Trips[Trips$trip_id != "-1",]
    # Trips <- Trips[Trips$ColDist < innerBuff] # removes start and end points of trips 
  }
  return(Trips)
}


#' @rdname tripSplit
#' @param Track dataFrame. If using singleSplitID() directly, 'ID' field not needed.


splitSingleID <- function(Track, colony, innerBuff = 15, returnBuff = 45, duration = 12, nests=FALSE){

  ### facilitate nest-specific distance calculations ###
  if(nests == TRUE)
  {  if(!"ID" %in% names(colony)) stop("colony missing ID field")
    nest<- colony[match(unique(Track$ID), colony$ID),]
    colonyWGS <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } else{
    colonyWGS <- SpatialPoints(data.frame(colony$Longitude, colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } 		## ends the else loop for nests=FALSE

  ### set up data to include in output ###
  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]

  Track$Returns <- ""
  Track$StartsOut <- ""
  Track$trip_id <- 0
  Track$ColDist <- spDistsN1(Track, colonyWGS, longlat = T) * 1000 # distance calculated on great circle (WGS84)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  returnBuff <- returnBuff * 1000   ### convert from km into m
  innerBuff <- innerBuff * 1000   ### convert from km into m
  if(is.null(duration)) { duration <- 0.0001 }

  ### SPLIT THE DATA INTO DISCRETE TRIPS ###
  i <- 0
  while(i < base::nrow(Track))
  {
    i <- i + 1
    if(Track$ColDist[i] < innerBuff) {Track$trip_id[i] <- -1} else {
      k <- i
      Dist <- Track$ColDist[i]
      if(i == nrow(Track)) {Track$trip_id[i] <- -1
      break
      } 
      if(i>1) {i <- i-1}
      while(Dist >= innerBuff)
      {
        if(k == nrow(Track) & Dist < returnBuff) {break} else {
          if(k == nrow(Track))
          {
            message(paste("track ", Track$ID[1], Trip.Sequence + 1, " does not return to the colony", sep=""))
            Track$Returns[i:k] <- "No" 
            break
          }
        }
        k <- k + 1
        Dist <- Track$ColDist[k]
      }
      # print(paste(Track$ID[1], i))
      Time.Diff <- (Track$TrackTime[k] - Track$TrackTime[i]) / 3600
      Max.Dist <- max(Track$ColDist[i:k])
      if(Time.Diff < duration |  Max.Dist < innerBuff)
      {
        Track$trip_id[i:k] <- -1
        i <- k
        next
      }
      Trip.Sequence <- Trip.Sequence + 1
      if(i==1) { # if track starts outside innerBuff, print message
        message(paste0("track ", Track$ID[1], Trip.Sequence, " starts out on trip", sep=""))
        Track$StartsOut[i:k] <- "Yes" 
        Track$trip_id[i:k] <- paste(Track$ID[1], Trip.Sequence, sep="")
      } else {
        Track$trip_id[i:k] <- paste(Track$ID[1], Trip.Sequence, sep="")
      }
      i <- k
    }
  }
  
  Track$Returns <- ifelse(Track$Returns != "No" & Track$trip_id != "-1", "Yes", Track$Returns)
  
  return(Track)
}
