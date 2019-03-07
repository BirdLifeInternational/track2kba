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
    require(tidyverse)
    require(sf)
    require(lwgeom)

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
      DataGroup=tracks
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
  UIDs <- unique(TripCoords$ID)
  
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### SETTING PARAMETERS FOR kernelUD ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## kernelUD requires a grid and extent parameter to ensure that all kernel boundaries are contained in the spatial grid
    Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
    if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5} #changed from 3 to 5 on 23 Dec 2016 to avoid 'too small extent' error

    ### NEED TO EXPLORE: CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)  
    
    
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### ESTIMATING KERNEL DISTRIBUTION  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    #### ESTIMATING THE UD AND EXTRACTING THE CONTOURS 
    KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), same4all=TRUE)		##grid=1000, extent=BExt, 
    KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")
    
    
    
    
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### CONVERT, PLOT, AND RETURN OUTPUT ###
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    HR_sf <- st_as_sf(KDE.Sp) %>% st_buffer(Scale) %>%
                st_transform(4326) %>% ### convert to longlat CRS
                st_set_precision(100000) %>%
                lwgeom::st_make_valid()
    plot(HR_sf["id"])

    
    iba = st_intersection(HR_sf) # all intersections
    plot(iba["n.overlaps"])
    sf::st_is_valid(HR_sf)
    HR_sf_valid <- HR_sf %>% st_set_precision(100000) %>% lwgeom::st_make_valid()
    #### CONVERT THE OUTPUT INTO POLYGONS
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
