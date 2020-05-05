## tripSplit    #####################################################################################################################

#' Split tracking data into trips
#'
#' \code{tripSplit} employs \code{splitSingleID} to split data from multiple individuals' into discrete trips made from centrally-located places.
#'
#' This function splits central place foraging animal movement data into individual trips away from a central location based on distance and time.
#'
#' \emph{Nests}=T may be used if it is desired, for example, to use specific nest locations instead of one central location for all individuals/tracks.
#'
#' @param tracks data.frame. Must contain 'Latitude', 'Longitude', 'ID' and 'DateTime' columns (correct format may be assured using \code{\link{formatFields}} function).
#' @param Colony data.frame. Containing 'Latitude' and 'Longitude' fields specifying the central location(s) from which trips begin. If data are from MoveBank this information may be extracted using the \code{\link{move2KBA}} function. If \emph{Nests}=T, each row should correspond to an appropriate location (Lat/Lon) for each ID value in \emph{tracks}.
#' @param InnerBuff numeric (in kilometers). Indicate the distance that an animal must travel for the movement to be considered a trip. Note that this is also the metric that determines whether two subsequent trips are split - if your animal records locations > \code{InnerBuff} (km) from its place of origin and no locations at the place of origin (e.g. for burrow-nesting species) then subsequent trips may be lumped into a single trip. Increase \code{InnerBuff} to ensure correct splitting of trips. 
#' @param ReturnBuff numeric (in kilometers). Indicate the proximity required for a trip to be considered as returning. This is useful for identifying incomplete trips (i.e. where data storage/transmission failed during the trip).
#' @param Duration numeric (in hours). The period of time that the animals must be at large for the movement to be considered a trip.
#' @param Nests logical scalar (TRUE/FALSE). Should the central place used in trip-splitting be specific to each ID? If so, each place must be matched with an 'ID' value in both \emph{tracks} and \emph{Colony} objects.
#' @param plot logical scalar (TRUE/FALSE). Should trips be plotted? If so, the first 20 will be vizualized.
#' @param rmNonTrip logical scalar (TRUE/FALSE). Should periods not associated with trips be filtered out? Note that this does not filter out the trip starting and ending points which fall within InnerBuff, to allow more accurate calculations of duration and distance covered with \code{tripSummary}. Default is TRUE.
#' @return Returns an un-projected (WGS84) SpatialPointsDataFrame, with the field 'trip_id' added to identify each unique trip-ID combination. If rmNonTrip=TRUE, then output has been filtered of points deemed not associated with trip movements.
#'
#' @seealso 
#'
#' @examples
#' \dontrun{Trips <- tripSplit(tracks, 
#' Colony=Colony, 
#' InnerBuff=2, 
#' ReturnBuff=20, 
#' Duration=1, 
#' plot=T, 
#' Nests = F, 
#' rmNonTrip = T)}
#' ## needs to be set-up to use published example dataset
#' @export
#' @import ggplot2
#' @import sp
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats median

#### MAIN WRAPPER FUNCTION THAT INCLUDES DATA PREP AND LOOP OVER EACH ID

tripSplit <- function(tracks, Colony, InnerBuff = NULL, ReturnBuff = NULL, Duration = NULL, Nests=FALSE, plot=TRUE, rmNonTrip=TRUE)
{
  tracks <- tracks %>%
      mutate(DateTime = lubridate::ymd_hms(.data$DateTime)) %>% 
      mutate(TrackTime = as.double(.data$DateTime)) %>%
      mutate(trip_id = .data$ID) %>%
      arrange(.data$ID, .data$TrackTime)

  ### CREATE PROJECTED DATAFRAME ###
  DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = tracks, match.ID=FALSE)

  DataGroup.Projected <- DataGroup

  ### LOOP OVER EVERY SINGLE ID ###
  for(nid in seq_len(length(unique(tracks$ID)))){
    TrackIn <- base::subset(DataGroup.Projected, DataGroup.Projected$ID == unique(DataGroup.Projected$ID)[nid])
    TrackOut <- splitSingleID(Track=TrackIn, Colony=Colony, InnerBuff = InnerBuff, ReturnBuff = ReturnBuff, Duration = Duration, Nests=Nests)
    if(nid == 1) {Trips <- TrackOut} else {Trips <- maptools::spRbind(Trips, TrackOut)}
  }

  ### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE
  if(plot == TRUE)
  {
    if(length(unique(Trips@data$ID))>25){
      selectIDs <- unique(Trips@data$ID)[1:25]
      plotdat<-  Trips@data %>% dplyr::filter(.data$ID %in% selectIDs)
      warning("Too many individuals to plot. Only the first 25 ID's will be shown")
    }else{plotdat <- Trips@data}

    TRACKPLOT <- plotdat %>% mutate(complete=ifelse(.data$Returns=="No","No","Yes")) %>%
      arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
      ggplot(aes(.data$., x=.data$Longitude, y=.data$Latitude, col=.data$complete)) +
      geom_path() +
      geom_point(data=Colony, aes(x=.data$Longitude, y=.data$Latitude), col='red', shape=16, size=2) +
      facet_wrap(ggplot2::vars(.data$ID)) +
      theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_blank())

    ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ###
    if (min(tracks$Longitude) < -170 &  max(tracks$Longitude) > 170) {
      plotdat <-  Trips@data %>%
        mutate(Longitude=ifelse(.data$Longitude<0, .data$Longitude+360, .data$Longitude))

      Colony$Longitude  <- ifelse(Colony$Longitude<0, Colony$Longitude+360, Colony$Longitude)
      longlimits <- c(min(plotdat$Longitude)-2, max((plotdat$Longitude)+2))
      longbreaks <- round(seq(longlimits[1], longlimits[2], by=10)/10,0)*10
      longlabels <- ifelse(longbreaks>180, longbreaks-360, longbreaks)

      TRACKPLOT <- plotdat %>% mutate(complete=ifelse(.data$Returns=="No", "No", "Yes")) %>%
        arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
        ggplot(.data$., aes(x=.data$Longitude, y=.data$Latitude, col=.data$complete)) +
        geom_path() +
        geom_point(data=Colony, aes(x=.data$Longitude, y=.data$Latitude), col='red', shape=16, size=2) +
        scale_x_continuous(limits = longlimits,
          breaks = longbreaks,
          labels = longlabels) +
        facet_wrap(ggplot2::vars(.data$ID)) +
        theme(panel.background=element_rect(fill="white", colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="black", fill="white"),
          panel.border = element_blank())}

    base::print(TRACKPLOT)
  }

  if(rmNonTrip==TRUE) {
    Trips <- Trips[Trips$trip_id != "-1",]
    # Trips <- Trips[Trips$ColDist < InnerBuff] # removes start and end points of trips 
  }
  return(Trips)
}


#' @rdname tripSplit
#' @param Track dataFrame. If using singleSplitID() directly, 'ID' field not needed.


splitSingleID <- function(Track, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, Nests=FALSE){

  ### facilitate nest-specific distance calculations ###
  if(Nests == TRUE)
  {  if(!"ID" %in% names(Colony)) stop("Colony missing ID field")
    nest<- Colony[match(unique(Track$ID), Colony$ID),]
    ColonyWGS <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } else{
    ColonyWGS <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } 		## ends the else loop for Nests=FALSE

  ### set up data to include in output ###
  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]

  Track$Returns <- ""
  Track$StartsOut <- ""
  Track$trip_id <- 0
  Track$ColDist <- spDistsN1(Track, ColonyWGS, longlat = T) * 1000 # distance calculated on great circle (WGS84)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  ReturnBuff <- ReturnBuff * 1000   ### convert from km into m
  InnerBuff <- InnerBuff * 1000   ### convert from km into m
  if(is.null(Duration)) { Duration <- 0.0001 }

  ### SPLIT THE DATA INTO DISCRETE TRIPS ###
  i <- 0
  while(i < base::nrow(Track))
  {
    i <- i + 1
    if(Track$ColDist[i] < InnerBuff) {Track$trip_id[i] <- -1} else {
      k <- i
      Dist <- Track$ColDist[i]
      if(i == nrow(Track)) {Track$trip_id[i] <- -1
      break
      } 
      if(i>1) {i <- i-1}
      while(Dist >= InnerBuff)
      {
        if(k == nrow(Track) & Dist < ReturnBuff) {break} else {
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
      if(Time.Diff < Duration |  Max.Dist < InnerBuff)
      {
        Track$trip_id[i:k] <- -1
        i <- k
        next
      }
      Trip.Sequence <- Trip.Sequence + 1
      if(i==1) { # if track starts outside InnerBuff, print message
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
