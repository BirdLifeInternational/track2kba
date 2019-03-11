### findScale ###############################################################################
#########################################################################################

#### Description ###
# This script takes a tracking dataset, and outputs a one-row dataframe with smoothing parameter ('h') values calucluated in *three* different ways:
# 1. href: a simple, data-driven method which takes into account the number of points, and the variance in X and Y directions
# 2. Scale of Area-Restricted Search (ARS): this method tries to estimate the scale at which the animal interacts with the environment, using First-Passage Time Analysis (e.g.Fauchard & Taveraa; Pinuad & Weimerskirch).

# ## test ##
# findScale(tracks, 
#       ARSscale = T, 
#       max_TripDist = pull(trip_distances, "max_dist"), 
#       whichStage="Incubation")


findScale <- function(Trips, ARSscale=T, Colony) {
  
  ### Packages
  pkgs <- c('sp', 'tidyverse', 'geosphere', 'adehabitatHR')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}
  
  ### Warning
  if(!"Returns" %in% names(Trips)) warning("Your data may not have been split into trips AND filtered of non-trip periods. Many of these points may be stationary periods, and therefore skew the resulting ARSscale parameter. If unsure, run tracks through tripSplit() and specify that rmColLocs=T.")
  
  ##################################################################
  ### CREATE PROJECTED DATAFRAME ###  ***** NEED TO ADD CLEAN TRACKS BIT
  if(is(Trips, "SpatialPointsDataFrame") != TRUE ) {
    Trips.wgs <- SpatialPointsDataFrame(SpatialPoints(data.frame(Trips$Longitude, Trips$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = Trips, match.ID=F)
    mid_point<-data.frame(centroid(cbind(Trips.wgs$Longitude, Trips.wgs$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(Trips.wgs$Longitude) < -170 &  max(Trips.wgs$Longitude) > 170) {
      longs=ifelse(Trips.wgs$Longitude<0,Trips.wgs$Longitude+360,Trips.wgs$Longitude)
      mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))
    }
    
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    Trips.Projected <- spTransform(Trips.wgs, CRS=proj.UTM)
    
  } else { Trips.Projected <- Trips}
  
  ### If dataset are SPDF but user has NOT projected them,  do so!
  if(is.projected(Trips.Projected) != TRUE) {
    Trips.Projected <- spTransform(Trips.Projected, CRS=proj.UTM)
  }
  
  #### prep data frame to fill ####
  HVALS <- data.frame(
    href=0,
    ARSscale=0,
    stringsAsFactors=F
  )
  
  ##################################################################
  ##### Href calculation (code from adehabitat::kernelUD() ) ####
  ##################################################################
  
  xy <- Trips.Projected
  
  xy <- coordinates(xy)
  
  varx <- var(xy[, 1])
  vary <- var(xy[, 2])
  sdxy <- sqrt(0.5 * (varx + vary))
  n <- nrow(xy)
  ex <- (-1/6)
  href <- sdxy * (n^ex)
  
  ##################################################################
  ##### calculate scale of ARS ####
  ##################################################################
  
  if(ARSscale == T){
    
    if(!"ID" %in% names(Trips.Projected)) stop("ARSscale error: ID field does not exist")
    if(!"TrackTime" %in% names(Trips.Projected)) stop("ARSscale error: TrackTime field does not exist")
    if(!"Latitude" %in% names(Trips.Projected)) stop("ARSscale error: Latitude field does not exist")
    if(!"Longitude" %in% names(Trips.Projected)) stop("ARSscale error: Longitude field does not exist")
    
    Trips.Projected$X <- Trips.Projected@coords[,1]
    Trips.Projected$Y <- Trips.Projected@coords[,2]
    #Trips@data$ID <- as.numeric(as.character(Trips@data$ID))
    if(is.factor(Trips.Projected@data$ID)==T){Trips.Projected@data$ID <- droplevels(Trips.Projected@data$ID)} 		## avoids the error 'some id's are not present' in as.ltraj
    
    Tripslt <- as.ltraj(data.frame(Trips.Projected$X, Trips.Projected$Y), date=as.POSIXct(Trips.Projected$TrackTime, origin="1970/01/01", tz="GMT"), id=Trips.Projected$ID, typeII = TRUE)
    
    ##################################################
    ### Determination of scales ###
    ##################################################
    
    ## helper function to calculate distance unless no previous location
    poss_dist <- possibly(geosphere::distm, otherwise = NA)
    
    ## all summary in one pipe
    med_displace <- as.data.frame(Trips@data) %>% 
      nest(Longitude, Latitude, .key = "coords") %>%
      group_by(trip_id) %>% 
      mutate(prev_coords = lag(coords)) %>%
      mutate(Dist = map2_dbl(coords, prev_coords, poss_dist)) %>% 
      dplyr::summarise(value = round(median(na.omit(Dist)), 2) / 1000) ## convert to km
    
    #### NEED TO SET MINIMUM BASED ON RES/Hvalue!!!
    Scales <- c(seq(0.5, 20, 
      by = max(0.5, quantile(med_displace$value, 0.25))), 
      seq(20, 50, 
        by = max(1, quantile(med_displace$value, 0.5))), 
      seq(50, 100, 
        by = max(5, quantile(med_displace$value, 0.75))),
      seq(100, 200, 
        by = max(5, quantile(med_displace$value, 0.9)))
    )
    
    ## FPT analysis
    fpt.out <- fpt(Tripslt, radii = Scales, units = "seconds")
    fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
    Temp <- as.double(fpt.scales[1,])
    plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))
    
    ars.scales <- NULL
    UIDs <- unique(Trips.Projected$ID)
    for(i in 1:length(UIDs))
    {
      if(length(Scales) == length(which(is.na(fpt.scales[i,])))) {print(paste("Warning: ID", UIDs[i], "is smaller than smallest scale and will be ignored")); next}
      Temp <- as.double(fpt.scales[i,])
      #lines(Scales,Temp)
      plot(Scales, Temp, type="l")
      
      q <- which(!is.na(Temp))
      p <- 2
      while(!is.na(Temp[q[p]]) & Temp[q[p]] < Temp[q[p-1]] & q[p] != length(Temp)) {p <- p + 1}
      while(!is.na(Temp[q[p]]) & Temp[q[p]] > Temp[q[p-1]]) {p <- p + 1}
      
      rfpt <- Scales[q[p-1]]
      if(suppressWarnings(min(which(is.na(Temp))) == p)) {print(paste("ID", UIDs[i], "has no peak")); next}
      FirstPeak <- Scales[q[p-1]]
      MaxPeak <- Scales[which(Temp == max(Temp[q[p-1]:length(Temp)], na.rm=T))]
      if(Peak == "Flexible")
      {
        if(FirstPeak < MaxPeak[1])
        {
          MaxPeak <- MaxPeak[MaxPeak >= FirstPeak]
          ifelse(MaxPeak[1] < FirstPeak + (max(Scales)/3), ars.sc <- MaxPeak[1], ars.sc <- FirstPeak)
        }  else  {ars.sc <- FirstPeak}
      }
      if(Peak == "Max") {ars.sc <- MaxPeak}
      if(Peak == "First")  {ars.sc <- FirstPeak}
      if(Peak == "User")
      {
        print("Select Peak on Graph")
        N <- identify(Scales, Temp, n=1)
        ars.sc <- Scales[N]
      }
      abline(v=ars.sc, col="red", lty=2)
      ars.scales <- c(ars.scales, ars.sc)
      #print(ars.sc)
      #readline("proceed?")
    }
    
    AprScale <- round(median(ars.scales), 3)            ### changed from mean to median to make output less susceptible to choice of input scales
    plot((Scales), Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)), xlab="Scales (km)", ylab="")
    for(i in 1:length(UIDs))
    {
      Temp <- as.double(fpt.scales[i,])
      lines((Scales),Temp)
    }
    abline(v=ars.scales, col="red", lty=2)
    abline(v=AprScale, col="darkred", lty=1, lwd=3)
    #print(ars.scales)
    #print(AprScale)
    text(max(Scales)/2, 1, paste(AprScale, "km"), col="darkred", cex=3)
    
    HVALS$ARSscale <- AprScale ## add ARS scale to data frame
  }
  
  
  ##################################################################
  ##### calculate mean foraging range ####
  ##################################################################
  
  ### Use tripSummary
  trip_distances <- tripSummary(Trips, Colony = Colony, nests = F)
  
  ForRangeH <- trip_distances %>% 
    ungroup() %>% 
    summarise(med_max_dist = round(median(max_dist), 2), 
      mag = round(log(max(max_dist)), 2)) %>%
    #mutate(mag=ifelse(mag<1,1,mag)) %>%
    mutate(scaled_mag = round(med_max_dist/mag, 2)) %>%
    mutate(scaled_large = round(ifelse(scaled_mag > 15, scales::rescale(scaled_mag, to = c(15, 50)), scaled_mag),2)) %>%
    mutate(scaled_small = round(scales::rescale(mag, to = c(0.5, 50)), 2))
  ##################################################################
  ######### Compile dataframe
  HVALS$href <- round(href/1000, 2)
  HVALS <- cbind.data.frame(HVALS, ForRangeH) %>% 
    dplyr::select(med_max_dist, mag, scaled_mag, scaled_large, scaled_small, href, ARSscale)
  
  return(HVALS)
  
}
