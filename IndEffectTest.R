#### IndEffectTest ####

## Written by Virginia Morera (2018), adapted from varianceTest() written by Phil Taylor and Mark Miller (2012)

## Changes by Martin Beal (2019):
#   - set default UDLev value to 50
#   - Remove default Scale parameter of 186km (user must input)
#   - made it so the Matching Package doesn't print a message when the function is run (i.e. 'quietly)
#   - Allow for input of SPDF (projected or un-projected) instead of just dataframe
#     - added a snippet dealing with projecting data that crosses the dataline
#   - Changed Grouping_var name to GroupVar
#   - Added inGroupVar argument (i.e. within Grouping Variable, variable) this is the variable over which overlaps will be calculated
#     - Small changes to allow for subsetting by inGroupVar, rather than "ID" column. 

## ver3 (vm, 2019)
#   - added argument conditional from function kerneloverlap with defalut = T
#   - changed argument name from inGroupVar to tripID


### ARGUMENTS TO THE FUNCTION      
# Tracks: must be a dataframe or SpatialPointsDataFrame with at least following fields: 
#    - Latitude*
#    - Longitude*
#    - tripID: it can have any other name, the user will specify it in the tripID argument, but it must have a unique identifier for each trip
#    - GroupVar: it can have any other name, and there can be as many as necessary, the user will specify it in the GroupVar argument. Variable to make the within-group vs. between group comparison (Year, Bird, etc.)
# * If it is a dataframe it must be unprojected (lonlat). If it's a SPDF it can be in any projection but it must be specified in the proj4string slot
#
# method: character, one of the options from the adehabitatHR::kerneloverlap function to calculate overlap
# conditional: logical. if TRUE the function sets to 0 the pixels of the grid over which the UD is estimated, outside the home range of 
#              the animal estimated at a level of probability equal to UDLev. Practically, if TRUE the maximum overlap will be 
#              equal to ~UDLev, if FALSE the maximum overlap will be equal to 1. 
# UDLev: numeric, value containing the % of the UD at which to calculate home-range (irrelevant if conditional is set to FALSE)
# Scale is the smoothing factor to be used in the Kernel Density Estimation (in Km)
# grid is a number giving the size of the grid on which the UD should be estimated. 


IndEffectTest <- function(Tracks, tripID, GroupVar, plotit=T, #own arguments
                          method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), conditional = TRUE, UDLev=50, Scale, grid = 500, #from adehabitat::kerneloverlap
                          nboots = 1000) #from Matching::ks.boot
{
  
  # packages
  require(sp)
  require(adehabitatHR) 
  require(Matching, quietly = T)
  require(tidyverse)
  
  # initial chceks
  if (!"Latitude" %in% names(Tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(Tracks)) stop("Longitude field does not exist")
  if (!(tripID) %in% names(Tracks)) stop("Within-group field does not exist")
  if (!GroupVar %in% names(Tracks)) stop("Group field does not exist")
  
  
  # MB # Added this section which converts Tracks to spatial dataframe and projects it (and if already is SPDF it accepts this)
  if (class(Tracks) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanTracks <- Tracks %>%
      dplyr::select(GroupVar, tripID, Latitude, Longitude)
    mid_point <- data.frame(centroid(cbind(CleanTracks$Longitude, CleanTracks$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanTracks$Longitude) < -170 &  max(CleanTracks$Longitude) > 170) {
      longs = ifelse(CleanTracks$Longitude < 0,CleanTracks$Longitude + 360,CleanTracks$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180,median(longs) - 360,median(longs))}
    
    Tracks.Wgs <- SpatialPoints(data.frame(CleanTracks$Longitude, CleanTracks$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    Tracks.Projected <- spTransform(Tracks.Wgs, CRS = proj.UTM )
    TracksSpatial <- SpatialPointsDataFrame(Tracks.Projected, data = CleanTracks)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, tripID, Latitude, Longitude)
    Tracks.Wgs <- NULL
    Tracks.Projected <- NULL
    
  }else {## if data are already in a SpatialPointsDataFrame then check for projection
    if (is.projected(Tracks)) {
      if ("trip_id" %in% names(Tracks@data)) {
        TracksSpatial <- Tracks }
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, tripID, Latitude, Longitude)
    }else {## project data to UTM if not projected
      mid_point <- data.frame(centroid(cbind(Tracks@data$Longitude, Tracks@data$Latitude)))
      
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(Tracks@data$Longitude) < -170 &  max(Tracks@data$Longitude) > 170) {
        longs = ifelse(Tracks@data$Longitude < 0, Tracks@data$Longitude + 360, Tracks@data$Longitude)
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      TracksSpatial <- spTransform(Tracks, CRS = proj.UTM)
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, tripID, Latitude, Longitude)
    }
  }
  
  # remove tripID tracks with < 6 points as they can't be used to calculate kernel
  # MB edit # Changed this step to happen after SPDF set-up. Also added tripID column.
  UIDs <- names(which(table(TracksSpatial@data[, tripID]) > 5))             
  TracksSpatial <- TracksSpatial[TracksSpatial@data[, tripID] %in% UIDs, ]
  TracksSpatial@data[ ,tripID] <- droplevels(as.factor(TracksSpatial@data[ ,tripID]))
  
  # create vector with value of GroupVar for each trip
  gid <- TracksSpatial@data[!duplicated(TracksSpatial@data[, tripID]), ][[GroupVar]]
  
  # calculate overlap between tracks
  X <- kerneloverlap(xy = TracksSpatial[, tripID], method = method, percent = UDLev, conditional = conditional, h = Scale*1000, grid = grid)
  X[lower.tri(X, diag = T)] <- NA
  
  # assign value of GroupVar to rows and columns
  rownames(X) <- colnames(X) <- gid
  
  # separate within (WI) and between (BW) group overlaps
  WI <- NULL
  BW <- NULL
  for (i in seq_along(rownames(X))) {
    # i = 1
    x1 <- X[i,] 
    x2 <- x1[which(names(x1) == rownames(X)[i])]
    x3 <- x1[which(names(x1) != rownames(X)[i])]
    WI <- c(WI, x2)
    BW <- c(BW, x3)
  }
  BW <- BW[!is.na(BW)]
  WI <- WI[!is.na(WI)]
  
  ## VMP commented this out since the ks.boot function is robust to ties. Was leftover from when using stats::ks.test function
  # BW <- BW[BW != 0]
  # WI <- WI[WI != 0]
  
  # organize values in a dataframe for plotting
  Overlaps <- data.frame(Overlap = c(WI, BW), Type = c(rep("Within", length(WI)), rep("Between", length(BW))))
  
  if(plotit==T){
    print(ggplot(data = Overlaps, aes(x = Type, y = Overlap, fill = Type)) + geom_boxplot(notch = T) + theme_bw())
  }

  # ks <- ks.test(x = WI, y = BW)
  ks <- Matching::ks.boot(WI, BW, alternative = "two.sided", nboots = nboots) # more indicated when data don't come from continuous distr (ours have many 0s)
  
  # Organise output
  Result <- list()
  Result[1] <- list(X) # overlaps matrix
  Result[2] <- list(Overlaps) # df with overlap values (long format)
  Result[3] <- list(ks) # output from the ks.boot function
  names(Result) <- c("Overlap Matrix", "Overlaps", "Kolmogorov-Smirnov")
  return(Result)
}
