### findScale ##################################################################

#' Find an appropriate smoothing parameter
#'
#' \code{findScale} takes a tracking data set and outputs a series of candidate 
#' smoothing parameter values. Additionally, it compares the scale of movement 
#' resolved by the sampling resolution of the data set, to a grid of desired 
#' resolution.
#'
#' The purpose of this function is to provide guidance regarding the two most 
#' sensitive steps in the track2KBA analysis: specification of the (1) smoothing
#'  parameter and the (2) grid cell size for kernel density estimation (KDE). 
#'  Specifically, the goal is to allow for exploration of the effect of these 
#'  parameters and their inter-relatedness, so that an informed decision may be 
#'  made regarding their specification in subsequent track2KBA steps.
#'
#' Kernel density estimation has been identified as particularly sensitive to 
#' the specification of the smoothing parameter (AKA bandwidth, or 'H' value), 
#' that is, the parameter that defines the width of the normal distribution 
#' around each location. Many techniques for identifying 'optimal' smoothing 
#' parameters have been proposed (see Gitzen, Millspaugh, and Kernohan for a 
#' classic review; see Fleming and Calabreses 2017 for a later implementation) 
#' and many of these techniques have their merits; however, in the track2KBA 
#' implementation of KDE we have opted for simplicity.
#'
#' In the context of the track2KBA analysis, the smoothing parameter ought to 
#' represent the relevant scale at which the animal interacts with the 
#' environment. Therefore, when selecting a \emph{Scale} value for subsequent 
#' analysis, the user must take into account the movement ecology of the study 
#' species. For species which use Area-Restricted Search (ARS) behavior when 
#' foraging, First Passage Time analysis may be used to identify the scale of 
#' interaction with the environment (Fauchald and Tveraa 2003), however not all 
#' species use ARS when foraging and therefore different techniques must be 
#' used.
#'
#' What minimum spatial scales are detectable by the data also depends on the 
#' sampling resolution. Therefore, when applying First Passage Time analysis, 
#' \code{findScale} sets the range of scales at which movements are analyzed 
#' based on the distribution of forward, between-point displacements in the 
#' data.
#'
#' The grid cell size also affects the output of kernel density-based space use 
#' analyses. Therefore, by specifying the \emph{res} parameter you can check 
#' whether your desired grid cell size is reasonable, given the scale of 
#' movement resolved by your data.
#'
#' @param tracks SpatialPointsDataFrame or data.frame of animal relocations as 
#' formatted by \code{\link{formatFields}}. Must include 'ID' field. If input is
#'  a data.frame or unprojected SpatialPointsDF, must also include 'Latitude' 
#'  and 'Longitude' fields.
#' @param scaleARS logical scalar (TRUE/FALSE). Do you want to calculate the 
#' scale of area-restricted search using First Passage Time analysis? NOTE: does
#'  not allow for duplicate date-time stamps.

#' @param res numeric. The desired grid cell resolution (square kilometers) for 
#' subsequent kernel analysis (NOT performed in this function). If this is not 
#' specified, the scale of movement is compared to a 500-cell grid, with spatial
#'  extent determined by the latitudinal and longitudinal extent of the data.
#' @param sumTrips data.frame. Output of \code{\link{tripSummary}} function. If 
#' not specified, \code{\link{tripSummary}} will be called within the function.
#' @param scalesFPT numeric vector. Set of spatial scales at which to calculate 
#' First Passage Time. If not specified, the distribution of between-point 
#' distances will be used to derive a set. 
#' @param peakWidth numeric. How many scale-steps either side of focal scale 
#' used to identify a peak. Default is 1, whereby a peak is defined as any scale
#'  at which the variance in log FPT increases from the previous scale, and 
#'  decreases for the following one.
#' @param peakMethod character. Which method should be used to select the focal 
#' peak for each ID. Options are "first", "max", and "steep". "steep" is a value
#'  of scalesFPT at which the variance in log FPT changes the most compared to 
#'  the surrounding scale(s).
#'
#' @return This function returns a one-row dataframe with the foraging range in 
#' the first column (i.e. 'med_max_distance') calculated by 
#' \code{\link{tripSummary}}, and the median step length 
#' (i.e. between point distance) for the data set. The subsequent columns 
#' contain various candidate smoothing parameter ('h') values calculated in the 
#' following ways:
#' \enumerate{
#'   \item 'mag' - log of the foraging range (i.e. median maximum trip distance)
#'   \item 'href' - reference bandwidth a simple, data-driven method which takes
#'    into account the number of points, and the variance in X and Y directions.
#'
#'    \eqn{sqrt((X + Y)* (n^(-1/6)))}; where X=Longitude/Easting, 
#'    Y=Latitude/Northing, and n=number of relocations
#'   \item 'scaleARS' - spatial scale of area-restricted Search behavior as 
#'   estimated using First Passage Time analysis 
#'   (see \code{\link[adehabitatLT]{fpt}})
#' }
#' 
#' If the scaleARS option is used, a diagnostic plot is shown which illustrates 
#' the change in variance of log-FPT  values calculated at each FPT scale. Grey
#'  vertical lines indicate the peaks identified for each individual using 
#'  peakMethod method chosen, and the red line is the median of these, and the 
#'  resulting scaleARS in the output table.
#' 
#' All values are in kilometers.
#'
#' @examples
#' \dontrun{HVALS <- findScale(
#' tracks, scaleARS = TRUE, sumTrips = trip_distances)
#' }
#'
#' @export
#' @import dplyr
#' @import sp
#' @importFrom stats setNames


findScale <- function(
  tracks, scaleARS=TRUE, res=NULL, sumTrips=NULL, 
  scalesFPT = NULL, peakWidth=1, peakMethod="first") {
  
  ### prep data frame to fill -------------------------------------------------
  HVALS <- data.frame(
    href=0,
    scaleARS=0,
    stringsAsFactors=FALSE
  )

  ### Calculate href for each ID, then get average for dataset ----------------

  IDs <- unique(tracks$ID)
  href_list <- vector(mode="list", length(IDs))
  
  href_list <- lapply(split(tracks, tracks$ID), function(x)
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

  ##### calculate mean foraging range -----------------------------------------

  ### Use tripSummary ---------------------------------------------------------
  if(is.null(sumTrips)){
    message(
      "As no 'sumTrips' was supplied, the foraging range and mag, 
      cannot be calculated."
      )
    ForRangeH <- data.frame(med_max_dist = NA, mag = NA)
    max_dist <- 0
    
  } else {
    ForRangeH <- sumTrips %>%
      ungroup() %>%
      summarise(med_max_dist = round(median(.data$max_dist), 2),
        mag = round(log(.data$med_max_dist), 2)
        )
    max_dist <- max(sumTrips$max_dist)
        }

  ## Calculate median step length in data -------------------------------------
  poss_dist <- purrr::possibly(geosphere::distm, otherwise = NA)
  
  if("trip_id" %in% names(tracks)){
    grouped <- as.data.frame(tracks@data) %>%
      tidyr::nest(coords=c(.data$Longitude, .data$Latitude)) %>%
      group_by(.data$ID, .data$trip_id)
  } else {
    grouped <- as.data.frame(tracks@data) %>%
      tidyr::nest(coords=c(.data$Longitude, .data$Latitude)) %>%
      group_by(.data$ID)
  }
  
  med_displace <- grouped %>% 
    mutate(prev_coords = dplyr::lag(.data$coords)) %>%
    mutate(Dist = purrr::map2_dbl(
      .data$coords, .data$prev_coords, poss_dist)
      ) %>%
    dplyr::summarise(value = round(median(na.omit(.data$Dist)), 2) / 1000) 
  
  ### calculate scale of ARS --------------------------------------------------

  if(scaleARS == TRUE){
    
    tracks$X <- tracks@coords[,1]
    tracks$Y <- tracks@coords[,2]

    Tripslt <- adehabitatLT::as.ltraj(
      data.frame(tracks$X, tracks$Y),
      date=tracks$DateTime, id=tracks$ID, 
      typeII=TRUE
      )
    
    ### Determination of scalesFPT --------------------------------------------

    # Relating the scale of movement to the user's desired res value ----------
    minX <- min(coordinates(tracks)[,1])
    maxX <- max(coordinates(tracks)[,1])
    minY <- min(coordinates(tracks)[,2])
    maxY <- max(coordinates(tracks)[,2])

    if(is.null(res)){res <- (max(
      abs(minX - maxX) / 500, abs(minY - maxY) / 500)
      ) / 1000
    message(
      sprintf("No 'res' was specified. Movement scale in the data was compared 
        to a 500-cell grid with cell size of %s km squared.", round(res, 3)
        )
      )
    }

    minScale <- max(0.5, quantile(med_displace$value, 0.25))
    
    ### set 20km as absolute minimum start point for scalesFPT ----------------
    if(minScale > 20) {minScale <- 20
    message("The average step length in your data is greater than 20km. Data at 
      this resolution are likely inappropriate for identifying 
      Area-Restricted Search behavior. Consider using another scale parameter 
      (e.g. href).")
    } 
    ### check if cell size is small enough to resolve KDE ---------------------
    if (minScale < res*0.1228){
      message(
        "Your chosen 'res' is very large compared to the scale of movement in 
        your data. To avoid encompassing space use patterns in very few cells 
        later on, consider reducing 'res'."
        )
      }

    ### set max FPTscale to 200 km, or max foraging range ---------------------
    if(!is.null(sumTrips) & (!is.na(max_dist) & max_dist < 200) ) {
      maxScale <- max(sumTrips$max_dist)} else {maxScale <- 200}

    ### generate set of FPTscales based on distrib. of point2point displacement
    if(is.null(scalesFPT)){
      if(maxScale<20){scalesFPT <- c(
        seq(
          minScale, maxScale, by = max(0.5, quantile(med_displace$value, 0.25))
          )
        )
      } 
      if(maxScale>=20 & maxScale<50){
        scalesFPT <- c(
          seq(minScale, 20, 
            by = max(0.5, quantile(med_displace$value, 0.25)) ),
          seq(20, maxScale,
            by = max(1, quantile(med_displace$value, 0.5)))
          )
        }
      if(maxScale>=50 & maxScale<100){
        scalesFPT <- c(
          seq(minScale, 20,
            by = max(0.5, quantile(med_displace$value, 0.25))),
          seq(21, 50,
            by = max(1, quantile(med_displace$value, 0.5))),
          seq(50, maxScale,
            by = max(5, quantile(med_displace$value, 0.75)))
        )
      }
      if(maxScale>100){
        scalesFPT <- c(
          seq(minScale, 20,
            by = max(0.5, quantile(med_displace$value, 0.25))),
          seq(21, 50,
            by = max(1, quantile(med_displace$value, 0.5))),
          seq(55, 100,
            by = max(5, quantile(med_displace$value, 0.75))),
          seq(100, maxScale,
            by = max(10, quantile(med_displace$value, 0.9)))
        )
      }
      
      scalesFPT <- unique(scalesFPT) ## remove duplicated values --------------
      
    }

    ### FPT analysis -----------------------------------------------------------
    ## Calc. FPT for each FPTscale --------------------------------------------
    fpt.out <- adehabitatLT::fpt(Tripslt, radii = scalesFPT, units = "seconds")
    # Calc. variance in log FPT for each FPTscale ----------------------------
    out_scales <- adehabitatLT::varlogfpt(fpt.out, graph = FALSE)

    # turn rows into list of single row dfs ---------------------------------
    out_scales_list <- split(out_scales, seq(nrow(out_scales)))
    out_scales_list <- setNames(
      split(out_scales, seq(nrow(out_scales))), rownames(out_scales)
      )
    # find all peaks for each individual  ----------------------------------
    # m is # of pnts either side of pk with lower or equal value to focal pnt
    pks_list <- lapply(out_scales_list, function (x, m = peakWidth){
      # x <- unlist(out_scales_list[[4]]) # single-row df to vector
      x <- unlist(x, use.names = FALSE) # single-row df to vector
      shape <- diff(sign(diff(x, na.pad = FALSE)))
      pks <- unlist( lapply(which(shape < 0), function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) 
        else return(NULL)
      }) )
      
      # steepness around peak --------------------------------------------
      steepness <- unlist( lapply(pks, function(i){
        if(length(i) > 0) {
          s <- sum(diff(x[(i - m) : i]), abs(diff(x[i : (i + m)]))) }
        else { s <- NULL }
        return(s)
      }) )   
      
      pks       <- unlist(pks)               # all peaks 
      firstpeak <- pks[1]                    # first peak
      maxpeak   <- pks[pks %in% which(       # max peak
        x == suppressWarnings(max(x[pks])) 
        )]
      steeppeak <- pks[which.max(steepness)] # steepest peak
      
      pk_list <- list(
        allpeaks=pks, first=firstpeak, max=maxpeak, steep=steeppeak
        )
      return(pk_list)
    })
    
    nulls <- unlist( 
      lapply( pks_list, function(x) { return(is.null(x$allpeaks)) } ) 
      )
    if(length(nulls) > 0) {
      message("No peak found for ID(s):", 
        paste(names(pks_list[nulls]), collapse=" ") 
        )
    }
    # select only peak type chosen by fxn argument peakMethod -----------
    ars.scales <- unlist( lapply( pks_list, function(x) {
      return(x[[peakMethod]])
    } ) )
    AprScale <- round(median(ars.scales), 2) 
    HVALS$scaleARS <- AprScale
    
    # print diagnostic plot --------------------------------------------
    plot(scalesFPT, out_scales_list[[1]],
      ylim=c(0, max(na.omit(out_scales))), type="l",
      ylab="var(log FPT)")
    lapply(out_scales_list, function(x){
      points(scalesFPT, x, type="l")
    })
    abline(v = ars.scales, col="grey")
    abline(v = AprScale, col="red", lwd=2)
    
   } else {HVALS$scaleARS <- NA}

  ######### Compile dataframe -------------------------------------------------
  HVALS$href <- round(href/1000, 2)
  HVALS <- data.frame(
    med_max_dist = ForRangeH$med_max_dist, 
    step_length = round(median(med_displace$value),2), 
    mag = ForRangeH$mag, 
    href = HVALS$href, 
    scaleARS = HVALS$scaleARS
    )

  return(HVALS)
}

