#### IndEffectTest ####

## Written by Virginia Morera (2018), adapted from varianceTest() written by Phil Taylor and Mark Miller (2012)

## Changes by Martin Beal (2019):
#   - set default UDLev value to 50
#   - Change grid argument to Res, to be consistent with other track2KBA fxns. 
#   - Allow for input of SPDF instead of just dataframe
#   - Tracks -> Trips? Since input ought to be colony-cleaned, trip-split data?
#   - seems to overlap whole-individual datasets with themselves, not between trips.
#   - Grouping_var --> GroupVar

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


IndEffectTest <- function(Tracks, GroupVar, UDLev=50, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), Scale = 186, grid = 500, nboots = 1000) {
  
  # packages
  require(sp)
  require(adehabitatHR) 
  require(Matching)
  
  # initial chceks
  if (!"Latitude" %in% names(Tracks)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(Tracks)) stop("Longitude field does not exist")
  if (!"ID" %in% names(Tracks)) stop("ID field does not exist")
  if (!Grouping_var %in% names(Tracks)) stop("Group field does not exist")
  
  # remove Tracks with < 6 in as they can't be used to calculate kernel
  UIDs <- names(which(table(Tracks$ID) > 5))             
  Tracks <- Tracks[Tracks$ID %in% UIDs,]
  Tracks$ID <- droplevels(as.factor(Tracks$ID))
  
  # convert to spatial and project
  TracksSpatial <- SpatialPointsDataFrame(coords = cbind(Tracks$Longitude, Tracks$Latitude), data = data.frame(ID = Tracks$trip_id), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  LAEAProj <- CRS(paste("+proj=laea +lon_0=", mean(Tracks$Longitude), " +lat_0=", mean(Tracks$Latitude), sep = ""))
  TracksSpatial <- spTransform(TracksSpatial, LAEAProj)
  gid <- Tracks[!duplicated(Tracks$ID),][[Grouping_var]]
  gid <- unique(Tracks$trip_id)
  
  # calculate overlap between tracks
  X <- kerneloverlap(xy = TracksSpatial, method = method, percent = UDLev, conditional = F, h = Scale*1000, grid = grid)
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