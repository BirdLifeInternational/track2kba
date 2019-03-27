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

# Tracks must be a dataframe with de following fields: 
#    Latitude: not projected (latlon)
#    Longitude: not projected (latlon)
#    ID: individual track ID 
#    Group: ID of the group we want to test the fidelity of (Bird, Year, etc.)
# UDLev is the quantile to be used for the Utilisation Distribution for the overlap
# method is one of the methods available from adehabitatHR::kerneloverlap
# conditional: if TRUE the function sets to 0 the pixels of the grid over which the UD is estimated, outside the home range of 
#              the animal estimated at a level of probability equal to UDLev. Practically, if TRUE the maximum overlap will be 
#              equal to UDLev, if FALSE the maximum overlap will be equal to 1. 
# Scale is the smoothing factor to be used in the Kernel Density Estimation (in Km)
# grid is a number giving the size of the grid on which the UD should be estimated. 


IndEffectTest <- function(Tracks, GroupVar, inGroupVar, UDLev=50, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), Scale, grid = 500, nboots = 1000) {
  
  # packages
  require(sp)
  require(adehabitatHR) 
  require(Matching, quietly = T)
  require(tidyverse)
  
  # initial chceks
  if (!"Latitude" %in% names(Tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(Tracks)) stop("Longitude field does not exist")
  if (!GroupVar %in% names(Tracks)) stop("Group field does not exist")
  if (!inGroupVar %in% names(Tracks)) stop("Within-group field does not exist")
  
  # MB # Added this section which converts Tracks to spatial dataframe and projects it (and if already is SPDF it accepts this)
  if(class(Tracks)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## filter DF to the minimum fields that are needed
    CleanTracks <- Tracks %>%
      dplyr::select(GroupVar, inGroupVar, Latitude, Longitude)
    mid_point <- data.frame(centroid(cbind(CleanTracks$Longitude, CleanTracks$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanTracks$Longitude) < -170 &  max(CleanTracks$Longitude) > 170) {
      longs=ifelse(CleanTracks$Longitude<0,CleanTracks$Longitude+360,CleanTracks$Longitude)
      mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))}
    
    Tracks.Wgs <- SpatialPoints(data.frame(CleanTracks$Longitude, CleanTracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    Tracks.Projected <- spTransform(Tracks.Wgs, CRS=proj.UTM )
    TracksSpatial <- SpatialPointsDataFrame(Tracks.Projected, data = CleanTracks)
    TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, inGroupVar, Latitude, Longitude)
    Tracks.Wgs<-NULL
    Tracks.Projected<-NULL

  }else{  ## if data are already in a SpatialPointsDataFrame then check for projection
    if(is.projected(Tracks)){
      if("trip_id" %in% names(Tracks@data)){
          TracksSpatial <- Tracks }
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, inGroupVar, Latitude, Longitude)
    }else{ ## project data to UTM if not projected
      mid_point <- data.frame(centroid(cbind(Tracks@data$Longitude, Tracks@data$Latitude)))
      
      ### MB  This part prevents projection problems around the DATELINE 
      if (min(Tracks@data$Longitude) < -170 &  max(Tracks@data$Longitude) > 170) {
        longs=ifelse(Tracks@data$Longitude<0,Tracks@data$Longitude+360,Tracks@data$Longitude)
        mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      TracksSpatial <- spTransform(Tracks, CRS=proj.UTM)
      TracksSpatial@data <- TracksSpatial@data %>% dplyr::select(GroupVar, inGroupVar, Latitude, Longitude)
    }
  }
  
  # remove inGroupVar tracks with < 6 points as they can't be used to calculate kernel
  # MB edit # Changed this step to happen after SPDF set-up. Also added inGroupVar column.
  UIDs <- names(which(table(TracksSpatial@data[, inGroupVar]) > 5))             
  TracksSpatial <- TracksSpatial[TracksSpatial@data[, inGroupVar] %in% UIDs, ]
  TracksSpatial@data[ ,inGroupVar] <- droplevels(as.factor(TracksSpatial@data[ ,inGroupVar]))

  ######
  gid <- TracksSpatial@data[!duplicated(TracksSpatial@data[, inGroupVar]), ][[GroupVar]]
  
  # calculate overlap between tracks
  X <- kerneloverlap(xy = TracksSpatial[, inGroupVar], method = method, percent = UDLev, conditional = F, h = Scale*1000, grid = grid)
  X[lower.tri(X, diag = T)] <- NA
  # assign Group ID to rows and columns
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
  BW <- BW[BW != 0]
  WI <- WI[WI != 0]
  
  # organize values in a dataframe for plotting
  Overlaps <- data.frame(Overlap = c(WI, BW), Type = c(rep("Within", length(WI)), rep("Between", length(BW))))
  # overlaps_key <- list(x = .97, y = .2, corner = c(1, 1),
  #                 text = list(c("Between", "Within")),
  #                 lines = list(type = c("l", "l"), col = c("blue", "magenta"),
  #                              lwd = 1, lty = 1))
  # print(ecdfplot(Overlaps$Overlap, groups = Overlaps$Type, ref = T, key = overlaps_key))
  # 
  # # plot boxplot
  # print(ggplot(data = Overlaps, aes(x = Type, y = Overlap, fill = Type)) + geom_boxplot(notch = T))
  # 
  # # plot densities
  # print(ggplot(data = Overlaps, aes(x = Overlap, fill = Type)) + geom_density(alpha = 0.5))
  # ks <- ks.test(x = WI, y = BW)
  ks <- Matching::ks.boot(WI, BW, alternative = "two.sided", nboots = nboots) # more indicated when data don't come from continuous distr (ours have many 0s)
  # # wt <- wilcox.test(x = WI, y = BW, paired = F, conf.int = T, conf.level = 0.95)
  Result <- list()
  Result[1] <- list(X)
  Result[2] <- list(Overlaps)
  Result[3] <- list(ks)
  names(Result) <- c("Overlap Matrix", "Overlaps", "Kolmogorov-Smirnov")
  return(Result)
}
