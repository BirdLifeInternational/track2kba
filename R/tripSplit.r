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
#' @param plotit logical scalar (T/F). Should trips be plotted? If so, the first 20 will be vizualized.
#' @param rmColLocs logical scalar (T/F). Should points associated with the central location (e.g. colony) be filtered out of the output?
#' @param cleanDF logical scalar (T/F). Should columns which are non-essential for track2KBA analysis be removed from dataframe, or not? Removal will speed analysis up a bit. 
#' @return Returns a projected, SpatialPointsDataFrame, with the field 'trip_id' added to identify each unique trip-ID combination. If rmColLocs=T, then output has been filtered of points deemed not associated with trip movements.
#'
#' @seealso \code{\link{tripSummary}}
#'
#' @examples
#' \dontrun{Trips <- tripSplit(tracks, 
#' Colony=Colony, 
#' InnerBuff=2, 
#' ReturnBuff=20, 
#' Duration=1, 
#' plotit=T, 
#' Nests = F, 
#' rmColLocs = T)}
#' ## needs to be set-up to use published example dataset
#' @export
#' @import ggplot2
#' @import sp
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats median

#### MAIN WRAPPER FUNCTION THAT INCLUDES DATA PREP AND LOOP OVER EACH ID

tripSplit <- function(tracks, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, Nests=FALSE, plotit=T, rmColLocs=T, cleanDF=F)
{

  ## load required packages ##
  # require(sp)
  # #require(maps)
  # #require(mapdata)
  # require(maptools)
  # require(rgdal)
  # require(geosphere)
  # require(ggplot2)
  # require(tidyverse)
  # pkgs <-c('sp', 'tidyverse', 'geosphere', 'ggplot2', 'maptools', 'lubridate')
  # for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}


  ## provide error messages ##
  if(!"Latitude" %in% names(tracks)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(tracks)) stop("Longitude field does not exist")
  if(!"ID" %in% names(tracks)) stop("ID field does not exist")
  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")
  if(!(is.double(InnerBuff) & is.double(ReturnBuff))) stop ("InnerBuff and ReturnBuff should be numbers")

  if(cleanDF==TRUE){
    cleantracks <- tracks %>%
      mutate(DateTime = lubridate::ymd_hms(DateTime)) %>%   ### needs some clever trick to convert to POSIXct if it isn't already POSIXct
      mutate(TrackTime = as.double(DateTime)) %>%
      mutate(trip_id = ID) %>%
      dplyr::select(ID, trip_id, Latitude, Longitude,DateTime, TrackTime) %>%
      arrange(ID, TrackTime)
  } else {
    cleantracks <- tracks %>%
      mutate(DateTime = lubridate::ymd_hms(DateTime)) %>%   ### needs some clever trick to convert to POSIXct if it isn't already POSIXct
      mutate(TrackTime = as.double(DateTime)) %>%
      mutate(trip_id = ID) %>%
      arrange(ID, TrackTime)
  }
  

  ### CREATE PROJECTED DATAFRAME ###
  DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(cleantracks$Longitude, cleantracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = cleantracks, match.ID=F)
  mid_point<-data.frame(geosphere::centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))

  ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
  if (min(cleantracks$Longitude) < -170 &  max(cleantracks$Longitude) > 170) {
    longs=ifelse(cleantracks$Longitude<0, cleantracks$Longitude+360,cleantracks$Longitude)
    mid_point$lon <- ifelse(median(longs)>180, median(longs)-360, median(longs))}

  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
  DataGroup.Projected <- spTransform(DataGroup, CRS=proj.UTM)


  ### LOOP OVER EVERY SINGLE ID ###
  for(nid in 1:length(unique(tracks$ID))){
    TrackIn <- base::subset(DataGroup.Projected, ID == unique(DataGroup.Projected$ID)[nid])
    TrackOut <- splitSingleID(Track=TrackIn, Colony=Colony, InnerBuff = InnerBuff, ReturnBuff = ReturnBuff, Duration = Duration, Nests=Nests, proj.UTM=proj.UTM)
    if(nid == 1) {Trips <- TrackOut} else {Trips <- maptools::spRbind(Trips, TrackOut)}
  }


  ### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE
  if(plotit == TRUE)
  {

    if(length(unique(Trips@data$ID))>25){
      selectIDs <- unique(Trips@data$ID)[1:25]
      plotdat<-  Trips@data %>% dplyr::filter(ID %in% selectIDs)
      warning("Too many individuals to plot. Only the first 25 ID's will be shown")
    }else{plotdat <- Trips@data}

    TRACKPLOT <- plotdat %>% mutate(complete=ifelse(.data$Returns=="N","no","yes")) %>%
      arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
      ggplot(aes(., x=Longitude, y=Latitude, col=complete)) +
      geom_path() +
      geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
      facet_wrap(ggplot2::vars(ID)) +
      theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_blank())


    ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ###
    if (min(cleantracks$Longitude) < -170 &  max(cleantracks$Longitude) > 170) {
      plotdat <-  Trips@data %>%
        mutate(Longitude=ifelse(.data$Longitude<0, .data$Longitude+360, .data$Longitude))

      Colony$Longitude  <- ifelse(Colony$Longitude<0, Colony$Longitude+360, Colony$Longitude)
      longlimits <- c(min(plotdat$Longitude)-2, max((plotdat$Longitude)+2))
      longbreaks <- round(seq(longlimits[1], longlimits[2], by=10)/10,0)*10
      longlabels <- ifelse(longbreaks>180, longbreaks-360, longbreaks)

      TRACKPLOT <- plotdat %>% mutate(complete=ifelse(.data$Returns=="N", "no", "yes")) %>%
        arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
        ggplot(., aes(x=Longitude, y=Latitude, col=complete)) +
        geom_path() +
        geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
        scale_x_continuous(limits = longlimits,
          breaks = longbreaks,
          labels = longlabels) +
        facet_wrap(ggplot2::vars(ID)) +
        theme(panel.background=element_rect(fill="white", colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="black", fill="white"),
          panel.border = element_blank())}


    base::print(TRACKPLOT)
  } ## end plotit=T loop

  if(rmColLocs==T) { # optional argument to remove points not associated with trips (i.e colony and small trips)
    Trips <- Trips[Trips$trip_id != "-1",]
  }

  return(Trips)
}


#' @rdname tripSplit
#' @param Track dataFrame. If using singleSplitID() directly, 'ID' field not needed.
#' @param proj.UTM CRS (Coordindate Reference System) object. As generated by sp::CRS()


splitSingleID <- function(Track, Colony,InnerBuff = 15, ReturnBuff = 45, Duration = 12, Nests=FALSE,proj.UTM){


  ### facilitate nest-specific distance calculations ###
  if(Nests == TRUE)
  {  if(!"ID" %in% names(Colony)) stop("Colony missing ID field")
    nest<- Colony[match(unique(Track$ID), Colony$ID),]
    Colony.Wgs <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } else{
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
  } 		## ends the else loop for Nests=FALSE
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
  while(i < base::nrow(Track))
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
        # print(paste("trip ", Track$ID[1], Trip.Sequence + 1, " is too small a trip"))
        next
      }
      Trip.Sequence <- Trip.Sequence + 1
      Track$trip_id[i:k] <- paste(Track$ID[1], Trip.Sequence, sep="")
      i <- k
      # print(paste(Track$ID[1], Trip.Sequence, sep=""))
    }
  }

  #if(plotit == TRUE){points(Track, pch=16, cex=0.75, col=as.factor(Track$trip_id))}
  return(Track)
}
