## tripSplit    ################################################################
#' Split tracking data into trips
#'
#' \code{tripSplit} employs \code{splitSingleID} to split data from multiple 
#' individuals' into discrete trips made from centrally-located places.
#'
#' This function splits central place foraging animal movement data into 
#' individual trips away from a central location based on distance and time.
#'
#' \code{nests=TRUE} may be used if it is desired, for example, to use specific 
#' nest locations instead of one central location for all individuals/dataGroup.
#'
#' @param dataGroup data.frame. Must contain 'Latitude', 'Longitude', 'ID' and 
#' 'DateTime' columns (correct format may be assured using 
#' \code{\link{formatFields}} function).
#' @param colony data.frame. Containing 'Latitude' and 'Longitude' fields 
#' specifying the central location(s) from which trips begin. If data are from 
#' MoveBank this information may be extracted using the \code{\link{move2KBA}} 
#' function. If \emph{nests}=TRUE, each row should correspond to an appropriate 
#' location (Lat/Lon) for each ID value in \emph{dataGroup}.
#' @param innerBuff numeric (in kilometers). Indicate the distance that an 
#' animal must travel for the movement to be considered a trip. Note that this 
#' is also the metric that determines whether two subsequent trips are split - 
#' if your animal records locations > \code{innerBuff} (km) from its place of 
#' origin and no locations at the place of origin (e.g. for burrow-nesting 
#' species) then subsequent trips may be lumped into a single trip. Increase 
#' \code{innerBuff} to ensure correct splitting of trips. 
#' @param returnBuff numeric (in kilometers). Indicate the proximity required 
#' for a trip to be considered as returning. This is useful for identifying 
#' incomplete trips (i.e. where storage/transmission failed during the trip).
#' @param duration numeric (in hours). The period of time that the animals must 
#' be at large for the movement to be considered a trip.
#' @param gapLimit numeric (in days). The period of time between points to be
#' considered too large to be a contiguous tracking event. Can be used to ensure
#' that deployments on the same bird in different years do not get combined into
#' extra long trips. Defaults to one year.
#' @param nests logical scalar (TRUE/FALSE). Should the central place used in 
#' trip-splitting be specific to each ID? If so, each place must be matched with
#'  an 'ID' value in both \emph{dataGroup} and \emph{colony} objects.
#' @param rmNonTrip logical scalar (TRUE/FALSE). Should periods not associated 
#' with trips be filtered out? Note that this does not filter out the trip 
#' starting and ending points which fall within innerBuff, to allow more 
#' accurate calculations of duration and distance covered with 
#' @param verbose logical scalar (TRUE/FALSE). Should the function print 
#' messages when trips start outside the innerBuffer or doesn't return to the 
#' 'colony'? Default is TRUE.
#' \code{tripSummary}. Default is TRUE.
#' @return Returns an un-projected (WGS84) SpatialPointsDataFrame, with the 
#' field 'tripID' added to identify each unique trip-ID combination. 
#' If rmNonTrip=TRUE, then output has been filtered of points deemed not 
#' associated with trip movements.
#'
#' @seealso \code{tripSummary}, \code{mapTrips}
#'
#' @examples
#' \dontrun{Trips <- tripSplit(dataGroup, 
#' colony=colony, 
#' innerBuff=2, 
#' returnBuff=20, 
#' duration=1, 
#' nests = FALSE, 
#' rmNonTrip = TRUE)}
#' @export
#' @import sp
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats median

tripSplit <- function(
  dataGroup, colony, innerBuff = NULL, returnBuff = NULL, 
  duration = NULL, gapLimit = NULL, nests=FALSE, rmNonTrip=TRUE, verbose=TRUE) {
  
  if(is.null(gapLimit)){gapLimit <- 365*24}
  
  dataGroup <- dataGroup %>%
      mutate(DateTime = lubridate::ymd_hms(.data$DateTime)) %>% 
      mutate(tripID = .data$ID) %>%
      arrange(.data$ID, .data$DateTime)

  # check for duplicated data
  dup_check <- dataGroup %>% group_by(.data$ID) %>% 
    mutate(duplicates = duplicated(.data$DateTime)) %>% ungroup() %>% 
    summarise(duplicates = sum(.data$duplicates))
  if(dup_check$duplicates > 0){message(
  "WARNING:dataset may contain duplicated data, this will affect trip-splitting"
  )}
  
  ### CREATE PROJECTED DATAFRAME ----------------------------------------------
  DataGroup <- SpatialPointsDataFrame(
    SpatialPoints(
      data.frame(dataGroup$Longitude, dataGroup$Latitude), 
      proj4string=sp::CRS(SRS_string = "EPSG:4326")
      ), 
    data = dataGroup, match.ID=FALSE)

  DataGroup.Projected <- DataGroup

  ### LOOP OVER EVERY SINGLE ID -----------------------------------------------
  for(nid in seq_len(length(unique(dataGroup$ID)))){
    TrackIn <- base::subset(
      DataGroup.Projected, 
      DataGroup.Projected$ID == unique(DataGroup.Projected$ID)[nid])
    
    TrackOut <- splitSingleID(
      Track=TrackIn, colony=colony, 
      innerBuff = innerBuff, returnBuff = returnBuff, 
      duration = duration, gapLimit = gapLimit, nests=nests, verbose = verbose)
    
    if(nid == 1) {Trips <- TrackOut} else {
      Trips <- maptools::spRbind(Trips, TrackOut)
      }
  }

  if(rmNonTrip==TRUE) {
    Trips <- Trips[Trips$tripID != "-1",]
  }
  return(Trips)
}


#' @rdname tripSplit
#' @param Track dataFrame.


splitSingleID <- function(
  Track, colony, innerBuff = 15, returnBuff = 45, duration = 12, gapLimit=gapLimit, nests=FALSE, verbose=verbose){
  
  ### facilitate nest-specific distance calculations ###
  if(nests == TRUE)
  {  if(!"ID" %in% names(colony)) stop("colony missing ID field")
    nest <- colony[match(unique(Track$ID), colony$ID),]
    colonyWGS <- SpatialPoints(
      data.frame(nest$Longitude, nest$Latitude), 
      #proj4string=CRS("+proj=longlat + datum=wgs84") ## causes error in latest rgdal release
      proj4string=CRS(SRS_string = "EPSG:4326")
      )
  } else{
    colonyWGS <- SpatialPoints(
      data.frame(colony$Longitude, colony$Latitude), 
      #proj4string=CRS("+proj=longlat + datum=wgs84") ## causes error in latest rgdal release
      proj4string=CRS(SRS_string = "EPSG:4326")
    )
  }

  ### set up data to include in output-----------------------------------------
  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]
  Track$Returns <- ""
  Track$StartsOut <- ""
  Track$tripID <- 0
  ### distance calculated on great circle (WGS84) -----------------------------
  Track$ColDist <- spDistsN1(Track, colonyWGS, longlat = TRUE) * 1000
  steptime <- c(abs(as.numeric(
    difftime(Track$DateTime[1:nrow(Track)-1],
             Track$DateTime[2:nrow(Track)], units="hours"))), NA)
  Trip.Sequence <- 0
  trip_dur <- 0
  Max.Dist <- 0
  returnBuff <- returnBuff * 1000 ### convert from km into m ------------------
  innerBuff <- innerBuff * 1000   
  if(is.null(duration)) { duration <- 0.0001 }

  ### SPLIT THE DATA INTO DISCRETE TRIPS ###
  i <- 0
  while(i < base::nrow(Track))
  {
    i <- i + 1
    if(Track$ColDist[i] < innerBuff) {Track$tripID[i] <- "-1"} else {
      k <- i
      Dist <- Track$ColDist[i]
      if(i == nrow(Track)) {Track$tripID[i] <- "-1"
      break
      } 
      if(i>1 & Track$tripID[i] == "-1") {i <- i-1}
      while( (Dist >= innerBuff) )
      {
        if(k == nrow(Track) & Dist < returnBuff) {break} else {
          if(k == nrow(Track)) {
            if(verbose == TRUE) {
              message(
                paste("track ", Track$ID[1], Trip.Sequence + 1, 
                      " does not return to the colony", sep="")
              )
            }
            Track$Returns[i:k] <- "No" 
            break
          } else if(steptime[k] > gapLimit) {
            if(Dist > returnBuff){
              Track$Returns[i:k] <- "No" 
            } else {Track$Returns[i:k] <- "Yes" }
            break
          }
        }
        k <- k + 1
        Dist <- Track$ColDist[k]
      }
      trip_dur <- as.numeric(
        difftime(Track$DateTime[k], Track$DateTime[i], units="hours")
        )
      Max.Dist <- max(Track$ColDist[i:k])
      if(trip_dur < duration |  Max.Dist < innerBuff)
      {
        Track$tripID[i:k] <- "-1"
        i <- k
        next
      }
      Trip.Sequence <- Trip.Sequence + 1
      if(i==1) { 
        if(verbose == TRUE){
          # if track starts outside innerBuff, print message --------------------
        message(
          paste0("track ", Track$ID[1], sprintf("%02d", Trip.Sequence), 
            " starts out on trip", sep="")
          )
        }
        Track$StartsOut[i:k] <- "Yes" 
        Track$tripID[i:k] <- paste(Track$ID[1], 
                                   sprintf("%02d", Trip.Sequence), sep="_")
      } else {
        Track$tripID[i:k] <- paste(Track$ID[1], 
                                   sprintf("%02d", Trip.Sequence), sep="_")
      }
      i <- k
    }
  }
  
  Track$Returns <- ifelse(
    Track$Returns != "No" & Track$tripID != "-1",
    "Yes", Track$Returns
    )
  
  return(Track)
}
