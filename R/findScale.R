### findScale ###############################################################################
#########################################################################################

#' Find an appropriate smoothing parameter
#'
#' \code{findScale} takes a tracking data set and outputs a series of candidate smoothing parameter values. Additionally, it compares the scale of movement resolved by the sampling resolution of the data set, to a grid of desired resolution.
#'
#' The purpose of this function is to provide guidance regarding the two most sensitive steps in the track2KBA analysis: specification of the (1) smoothing parameter and the (2) grid cell size for kernel density estimation (KDE). Specifically, the goal is to allow for exploration of the effect of these parameters and their inter-relatedness, so that an informed decision may be made regarding their specification in subsequent track2KBA steps.
#'
#' Kernel density estimation has been identified as particularly sensitive to the specification of the smoothing parameter ((AKA bandwidth, or 'H' value), that is, the parameter that defines the width of the normal distribution around each location. Many techniques for identifying 'optimal' smoothing parameters have been proposed (see Gitzen, Millspaugh, and Kernohan for a classic review; see Fleming and Calabreses 2017 for a later implementation) and many of these techniques have their merits; however, in the track2KBA implementation of KDE we have opted for simplicity.
#'
#' In the context of the track2KBA analysis, the smoothing parameter ought to represent the relevant scale at which the animal interacts with the environment. Therefore, when selecting a \emph{Scale} value for subsequent analysis, the user must take into account the movement ecology of the study species. For species which use Area-Restricted Search (ARS) behavior when foraging, First Passage Time analysis may be used to identify the scale of interaction with the environment (Fauchald and Tveraa 2003), however not all species use ARS when foraging and therefore different techniques must be used.
#'
#' What minimum spatial scales are detectable by the data also depends on the sampling resolution. Therefore, when applying First Passage Time analysis, \code{findScale} sets the range of scales at which movements are analyzed based on the distribution of forward, between-point displacements in the data.
#'
#' The grid cell size also affects the output of kernel density-based space use analyses. Therefore, by specifying the \emph{Res} parameter you can check whether your desired grid cell size is reasonable, given the scale of movement resolved by your data.
#'
#' @param Trips SpatialPointsDataFrame or data.frame of animal relocations as formatted by \code{\link{formatFields}}. Must include 'ID' field. If input is a data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param ARSscale logical scalar (TRUE/FALSE). Do you want to calculate the scale of area-restricted search using First Passage Time analysis? NOTE: does not allow for duplicate date-time stamps.

#' @param Res numeric (in kilometers). The desired grid cell resolution (square kilometers) for subsequent kernel analysis (NOT performed in this function). If this is not specified, the scale of movement is compared to a 500-cell grid, with spatial extent determined by the latitudinal and longitudinal extent of the data.
#' @param Trip_summary data.frame. Output of \code{\link{tripSummary}} function. If not specified, \code{\link{tripSummary}} will be called within the function.
#' @param FPTscales data.frame. Output of \code{\link{tripSummary}} function. If not specified, \code{\link{tripSummary}} will be called within the function.
#' @param Trip_summary data.frame. Output of \code{\link{tripSummary}} function. If not specified, \code{\link{tripSummary}} will be called within the function.
#'
#' @return Returns a one-row dataframe with smoothing parameter ('H') values and the foraging range estimated from the data set.
#'
#' This function returns a one-row dataframe with the foraging range in the first column (i.e. 'med_max_distance') calculated by \code{\link{tripSummary}}. The following columns contain various candidate smoothing parameter ('h') values calculated in the following ways:
#' \enumerate{
#'   \item 'mag' - log of the foraging range, calculated as the median maximum trip distance
#'   \item 'scaled_mag' - \eqn{med_max_distance / mag}
#'   \item 'href' - reference bandwidth a simple, data-driven method which takes into account the number of points, and the variance in X and Y directions.
#'
#'    \eqn{sqrt((X + Y)* (n^(-1/6)))}; where X=Longitude/Easting, Y=Latitude/Northing, and n=number of relocations
#'   \item 'scaleARS' - spatial scale of area-restricted Search behavior as estimated using First Passage Time analysis (see \code{\link[adehabitatLT]{fpt}})
#' }
#'
#' @examples
#' \dontrun{HVALS <- findScale(Trips, ARSscale = T, Trip_summary = trip_distances)}
#'
#' @export
#' @import dplyr
#' @import sp
#'


findScale <- function(Trips, ARSscale=T, Res=100, Trip_summary=NULL, FPTscales = NULL, plotPeaks = FALSE, Peak = "Flexible") {

  ### Packages
  # pkgs <- c('sp', 'tidyverse', 'geosphere', 'adehabitatHR')
  # for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}

  ### Warning
  if(!"Returns" %in% names(Trips)) warning("Your data may not have been split into trips AND filtered of non-trip periods. Many of these points may be stationary periods, and therefore skew the resulting ARSscale parameter. If unsure, run tracks through tripSplit() and specify that rmColLocs=T.")

  ##################################################################
  ### CREATE PROJECTED DATAFRAME ###  ***** NEED TO ADD CLEAN TRACKS BIT
  if(class(Trips)!= "SpatialPointsDataFrame") {
    Trips.wgs <- SpatialPointsDataFrame(SpatialPoints(data.frame(Trips$Longitude, Trips$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = Trips, match.ID=F)
    mid_point<-data.frame(geosphere::centroid(cbind(Trips.wgs$Longitude, Trips.wgs$Latitude)))

    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(Trips.wgs$Longitude) < -170 &  max(Trips.wgs$Longitude) > 170) {
      longs=ifelse(Trips.wgs$Longitude<0,Trips.wgs$Longitude+360,Trips.wgs$Longitude)
      mid_point$lon<-ifelse(median(longs)>180,median(longs)-360,median(longs))
    }

    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    Trips.Projected <- spTransform(Trips.wgs, CRS=proj.UTM)

  } else { Trips.Projected <- Trips}

  ### If dataset are SPDF but user has NOT projected them,  do so!
  if(sp::is.projected(Trips.Projected) != TRUE) {
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
  ##### calculate mean foraging range ####
  ##################################################################

  ### Use tripSummary
  if(is.null(Trip_summary)){
    warning("As no Trip_summary was supplied, the foraging range, mag, and scaled_mag cannot be calculated.")
    ForRangeH <- data.frame(med_max_dist = NA, mag = NA, scaled_mag = NA)
    max_dist <- 0
    
  } else {
    ForRangeH <- Trip_summary %>%
      ungroup() %>%
      summarise(med_max_dist = round(median(.data$max_dist), 2),
        mag = round(log(med_max_dist), 2)) %>%
      #mutate(mag=ifelse(mag<1,1,mag)) %>%
      mutate(scaled_mag = round(.data$med_max_dist / .data$mag, 2)
      )
    max_dist <- max(Trip_summary$max_dist)
    
    }

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

    Tripslt <- adehabitatLT::as.ltraj(data.frame(Trips.Projected$X, Trips.Projected$Y), date=as.POSIXct(Trips.Projected$TrackTime, origin="1970/01/01", tz="GMT"), id=Trips.Projected$ID, typeII = TRUE)

    ##################################################
    ### Determination of FPTscales ###
    ##################################################

    ## helper function to calculate distance unless no previous location
    poss_dist <- purrr::possibly(geosphere::distm, otherwise = NA)

    ## all summary in one pipe
    med_displace <- as.data.frame(Trips@data) %>%
      tidyr::nest(.data$Longitude, .data$Latitude, .key = "coords") %>%
      group_by(.data$trip_id) %>%
      mutate(prev_coords = dplyr::lag(.data$coords)) %>%
      mutate(Dist = purrr::map2_dbl(.data$coords, .data$prev_coords, poss_dist)) %>%
      dplyr::summarise(value = round(median(na.omit(.data$Dist)), 2) / 1000) ## convert to km

    # Relating the scale of movement in data to the user's desired Res value
    minX <- min(coordinates(Trips)[,1])
    maxX <- max(coordinates(Trips)[,1])
    minY <- min(coordinates(Trips)[,2])
    maxY <- max(coordinates(Trips)[,2])

    if(Res > 99){Res <- (max(abs(minX - maxX) / 500,
      abs(minY - maxY) / 500)) / 1000
    warning(sprintf("No grid resolution ('Res') was specified, or the specified resolution was >99 km and therefore ignored. Movement scale in the data was compared to a 500-cell grid with cell size of %s km squared.", round(Res, 3)), immediate. = TRUE)}

    minScale <- max(0.5, quantile(med_displace$value, 0.25))
    if(minScale > 20) {minScale <- 20
    warning("The average between-point displacement in your data is greater than 20km. Data at this resolution are likely inappropriate for identifying Area-Restricted Search behavior. Consider using another Scale parameter method (e.g. Href).")
    } ## set 20km as absolute minimum start point for FPTscales
    if (minScale < Res*0.1228){warning("Your chosen grid cell size (i.e. 'Res') is very large compared to the scale of movement in your data. To avoid encompassing space use patterns in very few cells later on, consider reducing 'Res'.", immediate. = TRUE)}

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
    Temp <- as.double(out_scales[1,])
    #plot(Scales, Temp, type="l", ylim=c(0, max(out_scales, na.rm=T)))
    Tripslt<-NULL
    fpt.out<-NULL

    ars.scales <- NULL
    UIDs <- unique(Trips.Projected$ID)
    for(i in 1:length(UIDs))
    {
      if(length(FPTscales) == length(which(is.na(out_scales[i,])))) {print(paste("Warning: ID", UIDs[i], "is smaller than smallest Scale and will be ignored")); next}
      Temp <- as.double(out_scales[i,])
      # lines(FPTscales,Temp)
      if(plotPeaks == TRUE){
        plot(FPTscales, Temp, type="l", main=paste("ID:", UIDs[i]))
      }

      q <- which(!is.na(Temp))
      p <- 2
      while(!is.na(Temp[q[p]]) & Temp[q[p]] < Temp[q[p-1]] & q[p] != length(Temp)) {p <- p + 1}
      while(!is.na(Temp[q[p]]) & Temp[q[p]] > Temp[q[p-1]]) {p <- p + 1}

      rfpt <- FPTscales[q[p-1]]
      if(suppressWarnings(min(which(is.na(Temp))) == p)) {print(paste("No peak was found for:", "ID", UIDs[i])); next}
      FirstPeak <- FPTscales[q[p-1]]
      MaxPeak <- FPTscales[which(Temp == max(Temp[q[p-1]:length(Temp)], na.rm=T))]
      
      if( (FirstPeak==FPTscales[length(FPTscales)-1]) && (FirstPeak == MaxPeak) ) {{print(paste("No peak was found for:", "ID", UIDs[i])); next}}

      if(Peak == "Flexible")  #MB# what is this part doing?
      {
        if(FirstPeak < MaxPeak[1])
        {
          MaxPeak <- MaxPeak[MaxPeak >= FirstPeak]
          ifelse(MaxPeak[1] < FirstPeak + (max(FPTscales)/3), ars.sc <- MaxPeak[1], ars.sc <- FirstPeak)
        }  else  {ars.sc <- FirstPeak}
      } else if(Peak == "First"){
        ars.sc <- FirstPeak
      } else if(Peak == "User"){
        print("Select Peak on Graph")
        N <- identify(FPTscales, Temp, n=1)
        ars.sc <- FPTscales[N]
      }
      
      if(plotPeaks == TRUE){
        abline(v=ars.sc, col="red", lty=2)
      }
      ars.scales <- c(ars.scales, ars.sc)
      #print(ars.sc)
      #readline("proceed?")
    }

    AprScale <- round(median(ars.scales), 2)            ### changed from mean to median to make output less susceptible to choice of input scales
    plot((FPTscales), Temp, type="l", ylim=c(0, max(out_scales, na.rm=T)), xlab="Scales (km)", ylab="Variance in first passage time")
    for(i in 1:length(UIDs))
    {
      Temp <- as.double(out_scales[i,])
      lines((FPTscales),Temp)
    }
    abline(v=ars.scales, col="red", lty=2)
    abline(v=AprScale, col="darkred", lty=1, lwd=3)
    text(max(FPTscales)/2, 1, paste(AprScale, "km"), col="darkred", cex=3)

    HVALS$ARSscale <- AprScale ## add ARS scale to data frame
  } else {HVALS$ARSscale <- NA}


  ##################################################################
  ######### Compile dataframe
  HVALS$href <- round(href/1000, 2)
  HVALS <- data.frame(med_max_dist = ForRangeH$med_max_dist, mag = ForRangeH$mag, scaled_mag = ForRangeH$scaled_mag, href = HVALS$href, ARSscale = HVALS$ARSscale)
  # HVALS <- cbind.data.frame(HVALS, ForRangeH) %>%
  #   dplyr::select(.data$med_max_dist, .data$mag, .data$scaled_mag, .data$href, .data$ARSscale)

  return(HVALS)

}
