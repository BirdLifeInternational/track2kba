### findScale ###############################################################################
#########################################################################################

#' Find an appropriate smoothing parameter
#'
#' \code{findScale} takes a tracking data set and outputs a series of candidate smoothing parameter values. Additionally, it compares the scale of movement resolved by the sampling resolution of the data set, to a grid of desired resolution.
#'
#' The purpose of this function is to provide guidance regarding the two most sensitive steps in the track2KBA analysis: specification of the (1) smoothing parameter and the (2) grid cell size for kernel density estimation (KDE). Specifically, the goal is to allow for exploration of the effect of these parameters and their inter-relatedness, so that an informed decision may be made regarding their specification in subsequent track2KBA steps.
#'
#' Kernel density estimation has been identified as particularly sensitive to the specification of the smoothing parameter (AKA bandwidth, or 'H' value), that is, the parameter that defines the width of the normal distribution around each location. Many techniques for identifying 'optimal' smoothing parameters have been proposed (see Gitzen, Millspaugh, and Kernohan for a classic review; see Fleming and Calabreses 2017 for a later implementation) and many of these techniques have their merits; however, in the track2KBA implementation of KDE we have opted for simplicity.
#'
#' In the context of the track2KBA analysis, the smoothing parameter ought to represent the relevant scale at which the animal interacts with the environment. Therefore, when selecting a \emph{Scale} value for subsequent analysis, the user must take into account the movement ecology of the study species. For species which use Area-Restricted Search (ARS) behavior when foraging, First Passage Time analysis may be used to identify the scale of interaction with the environment (Fauchald and Tveraa 2003), however not all species use ARS when foraging and therefore different techniques must be used.
#'
#' What minimum spatial scales are detectable by the data also depends on the sampling resolution. Therefore, when applying First Passage Time analysis, \code{findScale} sets the range of scales at which movements are analyzed based on the distribution of forward, between-point displacements in the data.
#'
#' The grid cell size also affects the output of kernel density-based space use analyses. Therefore, by specifying the \emph{Res} parameter you can check whether your desired grid cell size is reasonable, given the scale of movement resolved by your data.
#'
#' @param DataGroup SpatialPointsDataFrame or data.frame of animal relocations as formatted by \code{\link{formatFields}}. Must include 'ID' field. If input is a data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param ARSscale logical scalar (TRUE/FALSE). Do you want to calculate the scale of area-restricted search using First Passage Time analysis? NOTE: does not allow for duplicate date-time stamps.

#' @param Res numeric. The desired grid cell resolution (square kilometers) for subsequent kernel analysis (NOT performed in this function). If this is not specified, the scale of movement is compared to a 500-cell grid, with spatial extent determined by the latitudinal and longitudinal extent of the data.
#' @param Trip_summary data.frame. Output of \code{\link{tripSummary}} function. If not specified, \code{\link{tripSummary}} will be called within the function.
#' @param FPTscales numeric vector. Set of spatial scales at which to calculate First Passage Time. If not specified, the distribution of between-point distances will be used to derive a set. 
#' @param peakWidth numeric. How many scale-steps either side of focal scale used to identify a peak. Default is 1, whereby a peak is defined as any scale at which the variance in log FPT increases from the previous scale, and decreases for the following one.
#' @param findPeak character. Which method should be used to select the focal peak for each ID. Options are "first", "max", and "steep". "steep" is a FPTscale at which the variance in log FPT changes the most compared to the surrounding scale(s).
#'
#' @return This function returns a one-row dataframe with the foraging range in the first column (i.e. 'med_max_distance') calculated by \code{\link{tripSummary}}, and the median step length (i.e. between point distance) for the data set. The subsequent columns contain various candidate smoothing parameter ('h') values calculated in the following ways:
#' \enumerate{
#'   \item 'mag' - log of the foraging range, calculated as the median maximum trip distance
#'   \item 'href' - reference bandwidth a simple, data-driven method which takes into account the number of points, and the variance in X and Y directions.
#'
#'    \eqn{sqrt((X + Y)* (n^(-1/6)))}; where X=Longitude/Easting, Y=Latitude/Northing, and n=number of relocations
#'   \item 'scaleARS' - spatial scale of area-restricted Search behavior as estimated using First Passage Time analysis (see \code{\link[adehabitatLT]{fpt}})
#' }
#' 
#' All values are in kilometers.
#'
#' @examples
#' \dontrun{HVALS <- findScale(DataGroup, ARSscale = T, Trip_summary = trip_distances)}
#'
#' @export
#' @import dplyr
#' @import sp
#'


findScale <- function(DataGroup, ARSscale=TRUE, Res=NULL, Trip_summary=NULL, FPTscales = NULL, peakWidth=1, findPeak="first") {

  ##################################################################
  ### CREATE PROJECTED DATAFRAME ###  ***** NEED TO ADD CLEAN TRACKS BIT
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## set the minimum fields that are needed
    mid_point <- data.frame(geosphere::centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(DataGroup$Longitude) < -170 &  max(DataGroup$Longitude) > 170) {
      longs <- ifelse(DataGroup$Longitude < 0, DataGroup$Longitude + 360, DataGroup$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
    
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=proj.UTM )
    Trips.Projected <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    
  }else{   ## if data are already in a SpatialPointsDataFrame then check for projection
    if(is.projected(DataGroup)){
      Trips.Projected <- DataGroup

    }else{ ## if not projected, project data to custom laea
      mid_point <- data.frame(geosphere::centroid(cbind(DataGroup@data$Longitude, DataGroup@data$Latitude)))
      
      ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
      if (min(DataGroup@data$Longitude) < -170 &  max(DataGroup@data$Longitude) > 170) {
        longs <- ifelse(DataGroup@data$Longitude < 0, DataGroup@data$Longitude + 360,DataGroup@data$Longitude)
        mid_point$lon<-ifelse(median(longs) > 180, median(longs)-360, median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      Trips.Projected <- sp::spTransform(DataGroup, CRS=proj.UTM)
    }
  }
  
  #### prep data frame to fill ####
  HVALS <- data.frame(
    href=0,
    ARSscale=0,
    stringsAsFactors=F
  )

  ##################################################################
  #### Calculate href for each ID, then get average for dataset ####
  ##################################################################

  IDs <- unique(Trips.Projected$ID)
  href_list <- vector(mode="list", length(IDs))
  
  href_list <- lapply(split(Trips.Projected, Trips.Projected$ID), function(x)
    {
    xy <- coordinates(x)
    
    varx <- stats::var(xy[, 1])
    vary <- stats::var(xy[, 2])
    sdxy <- sqrt(0.5 * (varx + vary))
    n <- nrow(xy)
    ex <- (-1/6)
    href <- sdxy * (n^ex)
    return(href)  
   }
  )
  hrefs <- do.call(rbind, href_list)
  href <- mean(hrefs)

  ##################################################################
  ##### calculate mean foraging range ####
  ##################################################################

  ### Use tripSummary
  if(is.null(Trip_summary)){
    message("As no 'Trip_summary' was supplied, the foraging range, and mag, cannot be calculated.")
    ForRangeH <- data.frame(med_max_dist = NA, mag = NA)
    max_dist <- 0
    
  } else {
    ForRangeH <- Trip_summary %>%
      ungroup() %>%
      summarise(med_max_dist = round(median(.data$max_dist), 2),
        mag = round(log(.data$med_max_dist), 2)
        )
    max_dist <- max(Trip_summary$max_dist)
        }

  ## Calculate median step length in data ## 
  poss_dist <- purrr::possibly(geosphere::distm, otherwise = NA)
  
  if("trip_id" %in% names(DataGroup)){
    grouped <- as.data.frame(Trips.Projected@data) %>%
      tidyr::nest(coords=c(.data$Longitude, .data$Latitude)) %>%
      group_by(.data$ID, .data$trip_id)
  } else {
    grouped <- as.data.frame(Trips.Projected@data) %>%
      tidyr::nest(coords=c(.data$Longitude, .data$Latitude)) %>%
      group_by(.data$ID)
  }
  
  ## all summary in one pipe
  med_displace <- grouped %>% 
    mutate(prev_coords = dplyr::lag(.data$coords)) %>%
    mutate(Dist = purrr::map2_dbl(.data$coords, .data$prev_coords, poss_dist)) %>%
    dplyr::summarise(value = round(median(na.omit(.data$Dist)), 2) / 1000) ## convert to km
  
  ##################################################################
  ##### calculate scale of ARS ####
  ##################################################################

  if(ARSscale == T){
    
    Trips.Projected$X <- Trips.Projected@coords[,1]
    Trips.Projected$Y <- Trips.Projected@coords[,2]

    Tripslt <- adehabitatLT::as.ltraj(data.frame(Trips.Projected$X, Trips.Projected$Y), date=Trips.Projected$DateTime, id=Trips.Projected$ID, typeII = TRUE)
    
    ##################################################
    ### Determination of FPTscales ###
    ##################################################
    
    # Relating the scale of movement in data to the user's desired Res value
    minX <- min(coordinates(Trips.Projected)[,1])
    maxX <- max(coordinates(Trips.Projected)[,1])
    minY <- min(coordinates(Trips.Projected)[,2])
    maxY <- max(coordinates(Trips.Projected)[,2])

    if(is.null(Res)){Res <- (max(abs(minX - maxX) / 500, abs(minY - maxY) / 500)) / 1000
    message(sprintf("No 'Res' was specified. Movement scale in the data was compared to a 500-cell grid with cell size of %s km squared.", round(Res, 3)))}

    minScale <- max(0.5, quantile(med_displace$value, 0.25))
    
    if(minScale > 20) {minScale <- 20
    message("The average step length in your data is greater than 20km. Data at this resolution are likely inappropriate for identifying Area-Restricted Search behavior. Consider using another Scale parameter method (e.g. href).")
    } ## set 20km as absolute minimum start point for FPTscales
    if (minScale < Res*0.1228){warning("Your chosen 'Res' is very large compared to the scale of movement in your data. To avoid encompassing space use patterns in very few cells later on, consider reducing 'Res'.")}

    if(!is.null(Trip_summary) & (!is.na(max_dist) & max_dist < 200) ) {maxScale <- max(Trip_summary$max_dist)} else {maxScale <- 200}


    ### FPTscales NEED TO BE SET DEPENDING ON DATASET - THIS CAN FAIL IF MAXDIST <100 so we need to set this vector conditional on maxdist
    ### Setting the end of one seq() call the same number as the start of another, creates two of this value. However, Steffen changed to this in response to an error
    
    if(is.null(FPTscales)){
      if(maxScale<20){FPTscales <- c(seq(minScale, maxScale, by = max(0.5, quantile(med_displace$value, 0.25))))} 
      if(maxScale>=20 & maxScale<50){FPTscales <- c(seq(minScale, 20,
        by = max(0.5, quantile(med_displace$value, 0.25))),
        seq(20, maxScale,
          by = max(1, quantile(med_displace$value, 0.5))))}
      if(maxScale>=50 & maxScale<100){FPTscales <- c(seq(minScale, 20,
        by = max(0.5, quantile(med_displace$value, 0.25))),
        seq(21, 50,
          by = max(1, quantile(med_displace$value, 0.5))),
        seq(50, maxScale,
          by = max(5, quantile(med_displace$value, 0.75))))}
      if(maxScale>100){FPTscales <- c(seq(minScale, 20,
        by = max(0.5, quantile(med_displace$value, 0.25))),
        seq(21, 50,
          by = max(1, quantile(med_displace$value, 0.5))),
        seq(55, 100,
          by = max(5, quantile(med_displace$value, 0.75))),
        seq(100, maxScale,
          by = max(10, quantile(med_displace$value, 0.9))))}
      FPTscales <- unique(FPTscales) ## remove duplicated values
      
    }

    ## FPT analysis
    fpt.out <- adehabitatLT::fpt(Tripslt, radii = FPTscales, units = "seconds")
    out_scales <- adehabitatLT::varlogfpt(fpt.out, graph = FALSE)

    ### working on a replacement for old scaleARS peak identification code ## 
    # turn rows into list of single row dfs
    out_scales_list <- split(out_scales, seq(nrow(out_scales)))
    out_scales_list <- setNames(split(out_scales, seq(nrow(out_scales))), rownames(out_scales))
    # find all peaks for each individual (m is # of pnts either side of peak that has lower or equal value to focal point)
    pks_list <- lapply(out_scales_list, function (x, m = peakWidth){
      # x <- unlist(out_scales_list[[4]]) # single-row df to vector
      x <- unlist(x, use.names = F) # single-row df to vector
      shape <- diff(sign(diff(x, na.pad = FALSE)))
      pks <- unlist( lapply(which(shape < 0), function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(NULL)
      }) )
      
      # steepness around peak
      steepness <- unlist( lapply(pks, function(i){
        if(length(i) > 0) {
          s <- sum(diff(x[(i - m) : i]), abs(diff(x[i : (i + m)]))) }
        else { s <- NULL }
        return(s)
      }) )   
      
      pks       <- unlist(pks)    # all peaks
      firstpeak <- pks[1]         # first peak
      maxpeak   <- pks[pks %in% which( x == suppressWarnings(max(x[pks])) )] # max peak
      steeppeak <- pks[which.max(steepness)] # steepest peak
      
      pk_list <- list(allpeaks=pks, first=firstpeak, max=maxpeak, steep=steeppeak)
      return(pk_list)
    })
    
    nulls <- unlist( lapply( pks_list, function(x) { return(is.null(x$allpeaks)) } ) )
    if(length(nulls) > 0) {
      message("No peak found for ID(s):",  paste(names(pks_list[nulls]), collapse=" ") )
    }
    # select only peak type chosen by fxn argument findPeak 
    ars.scales <- unlist( lapply( pks_list, function(x) {
      return(x[[findPeak]])
    } ) )
    
    AprScale <- round(median(ars.scales), 2) 
    HVALS$ARSscale <- AprScale
    
   } else {HVALS$ARSscale <- NA}

  ######### Compile dataframe
  HVALS$href <- round(href/1000, 2)
  HVALS <- data.frame(med_max_dist = ForRangeH$med_max_dist, step_length = round(median(med_displace$value),2), mag = ForRangeH$mag, href = HVALS$href, ARSscale = HVALS$ARSscale)

  return(HVALS)
}
