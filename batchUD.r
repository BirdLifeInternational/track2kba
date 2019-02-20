## batchUD  #####################################################################################################################

## MAIN UPDATE: tidyverse, simple features (NOTE: kernelUD function and presumably all adehabitat functions DO NOT accept SF objects)
## REVISION: USE H-VALUE RELEVANT TO A SPECIES based on Oppel et al. (2018)
## INCLUDE WARNING FOR SPECIES NOT CONDUCIVE TO IBA



## Phil Taylor & Mark Miller, 2011

## batchUD calculates the Utilisation Distribution for groups of data in
## a larger dataset. The function relies on Adehabitat package for the
## calculations but returns the UD as SpatialPolygons so they can be exported
## as shapefiles and be more versatile in the script.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields. UD will be calculated for each unique ID value and
## each row should be a location. For each UD to be comparable, the data should be
## regularly sampled or interpolated.
## Scale should be the smoothing factor to be used in the Kernel Density Estimation
## and should be provided in Km.
## UDLev should be the quantile to be used for the Utilisation Distribution.
    
## UPDATED BY STEFFEN OPPEL ON 16 Dec 2016 to fix orphaned holes in output geometry
## UPDATED BY STEFFEN OPPEL ON 23 Dec 2016 to switch to adehabitatHR
## UPDATED BY ANA CARNEIRO AND STEFFEN OPPEL 31 JAN 2017 to same4all=T


batchUD <- function(DataGroup, Scale = 50, UDLev = 50)
    {
    require(sp)
    require(maptools)
    require(rgdal)
    require(adehabitatHR)   ### adehabitat is deprecated, switched to adehabitatHR on 23 Dec 2016
    require(geosphere)

    if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
    if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
    if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")

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
      
      
      ###### ENSURE THAT ONLY TRACKS WITH 5 OR MORE LOCATIONS ARE RETAINED (added 27 Feb 2017) ####
      UIDs <- names(which(table(DataGroup$ID)>5))             ### previous implementation restored
      #DataGroup@data$count<-1
      #nlocs<-aggregate(count~ID,DataGroup@data,sum)
      #retain<-as.character(nlocs$ID[nlocs$count>5])           ### kernelUD will fail for any ID with <5 locations
      DataGroup<-DataGroup[(DataGroup@data$ID %in% UIDs),]
      DataGroup@data$ID<-droplevels(as.factor(DataGroup@data$ID))        ### encountered weird error when unused levels were retained (27 Feb 2017)

    UIDs <- unique(DataGroup$ID)
    #note<-0          removed loop to check whether >5 data points exist per trip - considered unnecessary
    #KDE.Sp <- NULL
    TripCoords<-SpatialPointsDataFrame(DataGroup, data=data.frame(ID=DataGroup@data$ID,TrackTime=DataGroup@data$TrackTime))		
    TripCoords@data$TrackTime<-NULL 
    Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
    if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5} #changed from 3 to 5 on 23 Dec 2016 to avoid 'too small extent' error

    ### NEED TO EXPLORE: CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)  
      
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=1000, extent=BExt, same4all=TRUE)		## newer version needs SpatialPoints object and id no longer required in adehabitatHR, also removed 'extent' as it caused problems
KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")	## syntax differs from older version; THIS FUNCTION CAN FAIL WHEN same4all=FALSE

    
    KDE.Sp@proj4string <- DgProj
    KDE.Wgs <- spTransform(KDE.Sp, CRS=CRS("+proj=longlat +ellps=WGS84"))
    Tbl <- data.frame(Name_0 = rep(1, length(UIDs)), Name_1 = 1:length(UIDs), ID = UIDs)
    row.names(Tbl) <- UIDs
    KDE.Spdf <- SpatialPolygonsDataFrame(KDE.Sp, data=Tbl)

    plot(KDE.Spdf, border=factor(UIDs))
      
      ##### OVERLAYS IN THE polyCount FUNCTION WILL NOT WORK IF THE POLYGONS CONTAIN HOLES OR ARE ORPHANED
      ## simple fix to remove holes from polygon object
      va90a <- spChFIDs(KDE.Spdf, paste(KDE.Spdf$Name_0, KDE.Spdf$Name_1, KDE.Spdf$ID, sep = ""))
      va90a <- va90a[, -(1:4)]
      va90_pl <- slot(va90a, "polygons")
      va90_pla <- lapply(va90_pl, checkPolygonsHoles)
      p4sva <- CRS(proj4string(va90a))
      vaSP <- SpatialPolygons(va90_pla, proj4string = p4sva)
      va90b <- SpatialPolygonsDataFrame(vaSP, data = as(va90a, "data.frame"))   ### this returns an empty data frame
      va90b@data<-KDE.Spdf@data                                                   ### this adds the original data back into the data frame - may not work if entire polygons are removed
      
      
    return(va90b)     ## changed from KDE.Spdf to replace with cleaned version
    }
