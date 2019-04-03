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


IndEffectTest <- function(Tracks, tripID, GroupVar, #own arguments
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
  X <- kerneloverlap(xy = TracksSpatial[, tripID], method = method, percent = UDLev, conditional = conditionalOverlap, h = Scale*1000, grid = grid)
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
  
  ## Several plotting options. Can be commented out/in at convenience
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
  
  # Organise output
  Result <- list()
  Result[1] <- list(X) # overlaps matrix
  Result[2] <- list(Overlaps) # df with overlap values (long format)
  Result[3] <- list(ks) # output from the ks.boot function
  names(Result) <- c("Overlap Matrix", "Overlaps", "Kolmogorov-Smirnov")
  return(Result)
}

#### simulate distribution ####

# ColonyGroup must be a data frame with at least Longitude and Latitude (not projected), Sp and Colony fields
# PopulationInfo must be a data frame with at least Colony, Sp and Pairs
# Sp is the species
# multi_factor is the number by which we want to multiply each population (i.e. simulate that number of positions from each animal). The higher the number, the larger the resulting point pattern

simulateDistribution <- function(ColonyGroup, PopulationInfo, multi_factor = 5){
  require(adehabitatHR)
  require(sp)
  require(spatstat)
  mapa <- rworldmap::getMap(resolution = "coarse")
  # ColonyGroup <- boot_list[[2]] #this can be used for testing if something goes wrong
  
  # convert to spatial points
  coordinates(ColonyGroup) <- ~ Longitude + Latitude 
  ColonyGroup@proj4string <- CRS(as.character(NA))
  
  # generate kernel
  Kernel.est <- kernelUD(ColonyGroup, h = 1.5, grid = 1000, extent = 0.2)
  
  # convert to pixel image
  r <- raster(as(Kernel.est, "SpatialPixelsDataFrame"))
  # raster.as.im function from Jeffrey Evans answer here: https://bit.ly/2TI0FXB
  raster.as.im <- function(im) {
    r <- raster::res(im)
    orig <- sp::bbox(im)[, 1] + 0.5 * r
    dm <- dim(im)[2:1]
    xx <- unname(orig[1] + cumsum(c(0, rep(r[1], dm[1] - 1))))
    yy <- unname(orig[2] + cumsum(c(0, rep(r[2], dm[2] - 1))))
    return(spatstat::im(matrix(raster::values(im), ncol = dm[1], 
                               nrow = dm[2], byrow = TRUE)[dm[2]:1, ], 
                        xcol = xx, yrow = yy))
  }
  kernel.im <- raster.as.im(r)
  
  # select colony size info
  SColony <- as.character(unique(ColonyGroup$Colony))
  SSp <- as.character(unique(ColonyGroup$Sp))
  Pop.size <- PopulationInfo[PopulationInfo$Sp == SSp & PopulationInfo$Colony == SColony,]$Pairs*2
  
  # we're going to simulate a nº of points equal to the pop size * multi_factor
  SimulateN <- Pop.size*multi_factor
  
  # this simulates the points as ppp
  SimPoints <- rpoint(SimulateN, kernel.im)
  
  # convert to dataframe, and from there to Spatial points
  SimPoints.df <- as.data.frame(SimPoints)
  SimPoints.sp <- SimPoints.df
  coordinates(SimPoints.sp) <- ~ x+y
  
  # plot to see everything has worked
  par(mfrow = c(1,2))
  plot(kernel.im, main = SColony)
  plot(mapa, add = T, border = "white")
  plot(kernel.im, main = SColony)
  plot(mapa, add = T, border = "white")
  plot(SimPoints.sp, add = T, pch = 20, col = "#ff000010", cex = 0.1)
  par(mfrow = c(1,1))
  
  # prepare output
  SimPoints.sp <- SpatialPointsDataFrame(coords = SimPoints.sp@coords, 
                                         data = data.frame(Colony = rep(SColony, length(SimPoints.sp)), 
                                                           Sp = rep(SSp, length(SimPoints.sp))),
                                         proj4string = CRS(projections$WGS84))
  return(SimPoints.sp)
}

### bootstrap colony effect ####
bootstrap_Colony <- function(DataGroup, SimulatedDistr, Scale=186, Iteration=50, UDLev)
{
  
  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  if (!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if (!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  if (!"Sp" %in% names(DataGroup)) stop("Sp field does not exist")
  
  
  if (class(DataGroup) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    mid_point <- data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS = DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
  }else{DgProj <- DataGroup@proj4string}
  
  if (class(SimulatedDistr) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    SimulatedDistr.Wgs <- SpatialPoints(data.frame(SimulatedDistr$Longitude, SimulatedDistr$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    SimulatedDistr.Projected <- spTransform(SimulatedDistr.Wgs, CRS = DgProj)
    SimulatedDistr <- SpatialPointsDataFrame(SimulatedDistr.Projected, data = SimulatedDistr)
  }else {stop("SimulatedDistr must be a data frame")}
  
  SimulatedDistr <- SimulatedDistr[SimulatedDistr$Sp == unique(DataGroup$Sp),]
  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  UIDs <- as.character(unique(DataGroup$ID))
  Ntrips <- length(UIDs)
  Nloop <- seq(1,(Ntrips - 1),ifelse(Ntrips > 100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each = Iteration), Iteration = rep(seq(1:Iteration),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- UDLev
  
  #setup parallel backend to use 4 processors
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  Result <- data.frame()
  
  Result <- foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp","adehabitatHR","geosphere","rgdal")) %dopar% {
    # output_list <- list()
    # for (j in seq_along(LoopNr))  {
    # j = 360
    # LoopN <- LoopNr[j]
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]
    # Coverage <- NULL
    # Inclusion <- NULL
    # History <- NULL
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration = i)
    # set.seed(123)
    RanNum <- sample(UIDs, N, replace = F)
    sink("selected_ids_colony_bootstrap.txt", append = TRUE)
    cat(RanNum, "\n", "\n")
    sink()
    SelectedCoords <- DataGroup[DataGroup$ID %in% RanNum,]
    # Ext <- (min(SelectedCoords@coords[,1]) + 3 * diff(range(SelectedCoords@coords[,1])))
    # if(Ext < (Scale * 1000 * 2)) {
    # BExt <- ceiling((Scale * 1000 * 3)/(diff(range(SelectedCoords@coords[,1])))) #} else {BExt <- 3}
    KDE.Surface <- kernelUD(SelectedCoords, h = Scale*1000, grid = 500,  same4all = FALSE)
    try(KDE.UD <- getverticeshr(KDE.Surface, percent = UDLev))
    if (isTRUE(class(KDE.UD) == "try-error")) {
      sink("errors.txt", append = T)
      cat(paste("Failed in iteration", i, "with sample", RanNum, sep = " "))
      sink()} else {
        KDE.UD@proj4string <- DgProj
        Overlain <- over(SimulatedDistr, KDE.UD)$area
        Output$InclusionMean <- length(Overlain[!is.na(Overlain)])/nrow(SimulatedDistr)
        return(Output)
        # output_list[j] <- list(Output)
      }
  }
  
  ## stop the cluster
  stopCluster(cl)
  closeAllConnections() 
  
  par(mfrow = c(1,1), mai = c(1,1,1,1))
  #Result <- Output[1:nrow(Output) - 1,]
  Result$Colony <- unique(DataGroup$Colony)
  try(M1 <- nls(Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize), data = Result, start = list(a = 1,b = 0.1)), silent = TRUE)
  if ('M1' %in% ls()) {       ### run this only if nls was successful
    Result$pred <- predict(M1)
    P2 <- aggregate(pred ~ SampleSize, Result, FUN = mean)
    P2$sd <- aggregate(InclusionMean ~ SampleSize, Result, FUN = sd)[,2]
    plot(InclusionMean ~ SampleSize, data = Result, 
         pch = 16, cex = 0.2, col = "darkgray", ylim = c(0, 1),  
         ylab = "Inclusion", xlab = "SampleSize", main = paste(unique(DataGroup$Colony), "UDLev", UDLev, sep = "_"))
    yTemp <- c((P2[,2] + P2[,3]), rev(P2[,2] - P2[,3]))
    xTemp <- c(P2[,1], rev(P2[,1]))
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data = Result, pch = 16, cex = 0.2, col = "darkgray")
    lines(P2, lty = 1,lwd = 2)
    Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
    RepresentativeValue <- max(P2$pred)/Asymptote*100
    Result$RepresentativeValue <- RepresentativeValue
    print(RepresentativeValue)
    text(x = 2.5, y = 0.9,paste(round(RepresentativeValue,2), "%", sep = ""), cex = 2, col = "gray45", adj = 0)
  } else{RepresentativeValue <- mean(Result$InclusionMean[Result$SampleSize == max(Result$SampleSize)])   ### if nls is unsuccessful then use mean output for largest sample size
  Result$RepresentativeValue <- (RepresentativeValue/(UDLev/100))*100}    ## added by Jono Handley to convert to same scale as nls output
  
  Result$Asymptote <- Asymptote
  write.table(Result, paste(unique(DataGroup$Colony), "UDLEv", UDLev, "bootout_temp.csv", sep = "_"), row.names = F, sep = ",")
  return(Result)
}
