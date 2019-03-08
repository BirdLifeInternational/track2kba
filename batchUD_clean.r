## batchUD  #####################################################################################################################

## MAIN UPDATE: tidyverse
## REVISION: USE H-VALUE RELEVANT TO A SPECIES based on Oppel et al. (2018)
## INCLUDE WARNING FOR SPECIES NOT CONDUCIVE TO IBA



## Phil Taylor & Mark Miller, 2011

## batchUD calculates the Utilisation Distribution for groups of data in
## a larger dataset. The function relies on adehabitatHR package for the
## calculations and returns the UD as estUDm object.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields. UD will be calculated for each unique ID value and
## each row should be a location. For each UD to be comparable, the data should be
## regularly sampled or interpolated.
## Scale should be the smoothing factor to be used in the Kernel Density Estimation
## and should be provided in Km.
## UDLev should be the quantile to be used for the Utilisation Distribution.
## polyOut = logical, gives user the option to return a simple feature with kernel UD polygons


batchUD <- function(DataGroup, Scale = 50, UDLev = 50, Res=500, polyOut=FALSE)
    {
    require(sp)
    #require(maptools)
    #require(rgdal)
    require(adehabitatHR)   ### adehabitat is deprecated, switched to adehabitatHR on 23 Dec 2016
    require(geosphere)  ## needed for centroid
    require(tidyverse)

    if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
    if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
    if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")

    ##~~~~~~~~~~~~~~~~~~~~~~~~##
    ###### DATA PREPARATION ####
    ##~~~~~~~~~~~~~~~~~~~~~~~~##
    ## kernelUD requires a SpatialPointsDataFrame in a projected equal-area CRS as input
    ## we create a SPDF called 'TripCoords' in a projected CRS either from raw data or an existing SPDF
    ## only column in @dat should be the 'ID' field
  
    
    if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
    {
      ## set the minimum fields that are needed
      DataGroup <- DataGroup %>%
        dplyr::select(ID, Latitude, Longitude,DateTime) %>%
        arrange(ID, DateTime)
      mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
      DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=proj.UTM )
      TripCoords <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
      TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
    }else{  ## if data are already in a SpatialPointsDataFrame then check for projection
      if(is.projected(DataGroup)){
        TripCoords <- DataGroup
        TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
      }else{ ## project data to UTM if not projected
        mid_point<-data.frame(centroid(cbind(DataGroup@data$Longitude, DataGroup@data$Latitude)))
        proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
        TripCoords <- spTransform(DataGroup, CRS=proj.UTM)
        TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
      }
  
    }

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### REMOVING IDs WITH TOO FEW LOCATIONS ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  ###### ENSURE THAT kernelUD does not fail by retaining ONLY TRACKS WITH 5 OR MORE LOCATIONS
  validIDs <- names(which(table(TripCoords$ID)>5))
  TripCoords<-TripCoords[(TripCoords@data$ID %in% validIDs),]
  TripCoords@data$ID<-droplevels(as.factor(TripCoords@data$ID))        ### encountered weird error when unused levels were retained (27 Feb 2017)
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### SETTING PARAMETERS FOR kernelUD : THIS NEEDS MORE WORK!! ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ### deprecated as of 8 March 2019
  
  ## kernelUD requires a grid and extent parameter to ensure that all kernel boundaries are contained in the spatial grid
    # Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
    # if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5} #changed from 3 to 5 on 23 Dec 2016 to avoid 'too small extent' error

    ### CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)
  ### NEED TO DO: replace 500 with a 'resolution' parameter that users can choose
  ### NEED TO DO: link resolution of grid to H-parameter ('Scale')
  
  minX<-min(coordinates(TripCoords)[,1]) - Scale*2000
  maxX<-max(coordinates(TripCoords)[,1]) + Scale*2000
  minY<-min(coordinates(TripCoords)[,2]) - Scale*2000
  maxY<-max(coordinates(TripCoords)[,2]) + Scale*2000
  xrange<-seq(minX,maxX, by = diff(range(coordinates(TripCoords)[,1]))/Res)   ### if Res should be provided in km we need to change this
  yrange<-seq(minY,maxY, by = diff(range(coordinates(TripCoords)[,2]))/Res)   ### if Res should be provided in km we need to change this
  grid.locs<-expand.grid(x=xrange,y=yrange)
  INPUTgrid<-SpatialPixels(SpatialPoints(grid.locs), proj4string=proj4string(TripCoords))
  #  plot(INPUTgrid)
    

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### ESTIMATING KERNEL DISTRIBUTION  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## may need to insert extent=BExt, but hopefully avoided by custom-specified grid
    ## switched from same4all=T to =F because we provide a fixed input grid
    
  KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=INPUTgrid,same4all=F)
    
    
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### OPTIONAL POLYGON OUTPUT ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  if(polyOut==TRUE){
    require(sf)
    KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")
    HR_sf <- st_as_sf(KDE.Sp) %>%
      st_transform(4326) ### convert to longlat CRS
    
    plot(HR_sf["id"])
    return(list(KDE.Surface=KDE.Surface, UDPolygons=HR_sf))
  }else{
    return(KDE.Surface)
    }  ## changed from KDE.Spdf to replace with cleaned version
}
