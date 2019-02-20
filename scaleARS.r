## Jono H Test
## scaleARS     #####################################################################################################################

## MAIN UPDATE: tidyverse, simple features
## REVISION: make steps and sequence device and data dependent
## use device type to define seq (PTT: 5-250 km, GPS 0-250 km)
## use temporal data resolution to define intervals
## include error message if no ARS scale is sensible because no ARS

## REPLACE WITH DEFAULT FAMILY-SPECIFIC H-VALUE


## Phil Taylor & Mark Miller, 2012

## scaleARS undertakes First Passage Time (FPT) analysis on each trip in DataGroup (defined
## by the ID field) and identifies the scale at which each trip is interacting with the
## environment, based on the maximum variance in FPT value (Fauchard & Taveraa; Pinuad &
## Weimerskirch). The function relies on Adehabitat package for FPT calculation and returns a
## object of type numeric indicating the average scale across the datagroup.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields.
## Scales should be a vector with the scales (in kilometres) to be tested. Scales
## should be set according to the movements shown in the data (from 1:maximum distance travelled)
## and should not be too regular as this will begin to identify variances in sample rate rather
## than behaviour.
## Peak determines how the scale will be identified for each trip, this must be a character and from
## these options "Flexible", "First", "Max", "User". "Flexible" will make an automated decision,
## "First" will take the first peak in variance, "Max" will take the Maximum variance and "User"
## will allow the user to select the value on the graph.
    
## tested on 27 Dec 2016 with adehabitatLT, no issues encountered

scaleARS <- function(DataGroup, Scales = c(seq(1, 25, 1), seq(30, 50, 5), 75, seq(100, 250, 50)), Peak = "Flexible")
  {


  require(geosphere)
  require(sp)
  require(rgdal)
  require(rgeos)
  require(adehabitatLT)     ### updated to avoid loading deprecated adehabitat - tested and ok on 27 Dec 2016

  if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  if(!"TrackTime" %in% names(DataGroup)) stop("TrackTime field does not exist")

  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
    {
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    }else{DgProj<-DataGroup@proj4string}

  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  #DataGroup@data$ID <- as.numeric(as.character(DataGroup@data$ID))
  if(is.factor(DataGroup@data$ID)==T){DataGroup@data$ID <- droplevels(DataGroup@data$ID)} 		## avoids the error 'some id's are not present' in as.ltraj



  DataGrouplt <- as.ltraj(data.frame(DataGroup$X, DataGroup$Y), date=as.POSIXct(DataGroup$TrackTime, origin="1970/01/01", tz="GMT"), id=DataGroup$ID, typeII = TRUE)

  Scales <- Scales * 1000

  fpt.out <- fpt(DataGrouplt, radii = Scales, units = "seconds")
  fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
  Temp <- as.double(fpt.scales[1,])
  plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))

  ars.scales <- NULL
  UIDs <- unique(DataGroup$ID)
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
  return(AprScale)
  }
