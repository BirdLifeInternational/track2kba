### findScale ###############################################################################
#########################################################################################

#### Description ###
# This script takes a tracking dataset, and outputs a one-row dataframe with smoothing parameter ('h') values calucluated in *three* different ways:
# 1. href: a simple, data-driven method which takes into account the number of points, and the variance in X and Y directions
# 2. Scale of Area-Restricted Search (ARS): this method tries to estimate the scale at which the animal interacts with the environment, using First-Passage Time Analysis (e.g.Fauchard & Taveraa; Pinuad & Weimerskirch).

## test ##
findH(tracks, 
      ARSscale = T, 
      max_TripDist = pull(trip_distances, "max_dist"), 
      whichStage="Incubation")


findH <- function(DataGroup, ARSscale=T, max_TripDist, whichStage, Scales = c(seq(1, 25, 1), seq(30, 50, 5), 75, seq(100, 250, 50)), Peak = "Flexible")
{
  
  #### prep data ####
  
  S.df <- data.frame(
    stage=whichStage,
    href=0,
    ARSscale=0,
    foraging_range=0,
    stringsAsFactors=F
  )
  
  ##################################################################
  ### CREATE PROJECTED DATAFRAME ###  ***** NEED TO ADD CLEAN TRACKS BIT, AND RENAME DATAGROUP HERE
  DataGroup.wgs <- SpatialPointsDataFrame(SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = DataGroup, match.ID=F)
  mid_point<-data.frame(centroid(cbind(DataGroup.wgs$Longitude, DataGroup.wgs$Latitude)))
  
  ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
  if (min(DataGroup.wgs$Longitude) < -170 &  max(DataGroup.wgs$Longitude) > 170) {
    longs=ifelse(DataGroup.wgs$Longitude<0,DataGroup.wgs$Longitude+360,DataGroup.wgs$Longitude)
    mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))}
  
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
  DataGroup.Projected <- spTransform(DataGroup.wgs, CRS=proj.UTM)
  
  ##################################################################
  ##### Href calculation (code from adehabitat::kernelUD() ) ####
  ##################################################################
  
  xy <- DataGroup.Projected
  
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
    
    if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
    if(!"TrackTime" %in% names(DataGroup)) stop("TrackTime field does not exist")
    
    DataGroup.Projected$X <- DataGroup.Projected@coords[,1]
    DataGroup.Projected$Y <- DataGroup.Projected@coords[,2]
    #DataGroup@data$ID <- as.numeric(as.character(DataGroup@data$ID))
    if(is.factor(DataGroup.Projected@data$ID)==T){DataGroup.Projected@data$ID <- droplevels(DataGroup.Projected@data$ID)} 		## avoids the error 'some id's are not present' in as.ltraj
    
    DataGrouplt <- as.ltraj(data.frame(DataGroup.Projected$X, DataGroup.Projected$Y), date=as.POSIXct(DataGroup.Projected$TrackTime, origin="1970/01/01", tz="GMT"), id=DataGroup.Projected$ID, typeII = TRUE)
    
    Scales <- Scales * 1000
    
    fpt.out <- fpt(DataGrouplt, radii = Scales, units = "seconds")
    fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
    Temp <- as.double(fpt.scales[1,])
    plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))
    
    ars.scales <- NULL
    UIDs <- unique(DataGroup.Projected$ID)
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
    
    AprScale <- median(ars.scales)            ### changed from mean to median to make output less susceptible to choice of input scales
    AprScale <- round(AprScale/1000,3)
    plot((Scales/1000), Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)), xlab="Scales (km)", ylab="")
    for(i in 1:length(UIDs))
    {
      Temp <- as.double(fpt.scales[i,])
      lines((Scales/1000),Temp)
    }
    abline(v=ars.scales/1000, col="red", lty=2)
    abline(v=AprScale, col="darkred", lty=1, lwd=3)
    #print(ars.scales)
    #print(AprScale)
    text(max(Scales/1000)/2, 1, paste(AprScale, "km"), col="darkred", cex=3)

    S.df$ARSscale <- AprScale ## add ARS scale to data frame
  }
  
  
  ##################################################################
  ##### calculate mean foraging range ####
  ##################################################################
  
  forage_range <- mean(na.omit(max_TripDist))
  
  ##################################################################
  ######### Compile dataframe
  S.df$href[s] <- href/1000
  S.df$foraging_range <- forage_range
  
  return(S.df)
  
}
