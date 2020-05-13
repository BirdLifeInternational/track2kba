## findKBA  #####################################################################################################################

#' Delineating sites of potential importance to conservation
#'
#'
#' \code{findKBA} uses the utilization distributions of individual animals to 
#' identify areas of aggregation (i.e. where a large proportion of individuals' 
#' core areas overlap).
#'
#' The function first calculates a minimum threshold number of animals necessary
#'  for an area to be important, based on the representativeness of the tracking
#'   data (as quantified by \code{\link{repAssess}}). \code{findKBA} then 
#'   summarises the number of individual core UDs overlapping in an area and 
#'   compares that number against the thresholds.
#' The areas identified are POTENTIAL Key Biodiversity Areas. That is, they are 
#' areas of ecological relevance to the species, but must yet be assessed 
#' against global criteria (conservation status and global population size of 
#' the species) to determine whether they achieve global (or regional) KBA 
#' status. To assess the areas against the criteria, a population estimate is 
#' needed, from which a maximum number of animals in the population that may use
#'  each area is calculated, again adjusted based on the sample 
#'  representativeness. 
#' 
#' The criteria for site assessment are published in the KBA standard, which may
#'  be found here: \url{http://www.keybiodiversityareas.org/what-are-kbas}.
#' 
#' @param KDE estUDm or  SpatialPixels/GridDataFrame. If estUDm, as created by \code{\link{estSpaceUse}} or 
#' \code{adehabitatHR::kernelUD}, if Spatial*, each column should correspond to 
#' the Utilization Distribution of a single individual or track.
#' @param represent Numeric (between 0-1). Output value provided by 
#' \code{\link{repAssess}} which assesses how representative the tracking data 
#' are for characterising the space use of the wider population.
#' @param popSize Numeric, the number of individuals breeding or residing at the
#'  origin location from where animals were tracked, quantifying the population 
#'  that the tracking data represent. This number will be used to calculate how 
#'  many animals use the delineated areas of aggregation. If no value for 
#'  \code{popSize} is provided then output will be as the proportion of the 
#'  population.
#' @param levelUD Numeric (percentage). Specifies the quantile used for 
#' delineating the core use areas of individuals based on the kernel density 
#' distribution. Default set to 50\% based on Lascelles et al. (2016). For 
#' penguins higher values can be accepted, see Dias et al. (2018).
#' @param polyOut Logical. (Default TRUE) Should the output be a polygon dataset
#'  (TRUE) or grid of animal densities (FALSE). See 'Value' below for more 
#'  details.
#'
#' @return if \code{polyOut = TRUE} function returns an object of class 
#' \code{sf} containing polygon data with three data columns:
#'   Column \code{N_IND} indicates the number of tracked individuals whose core 
#'   use area (at \code{levelUD}) overlapped with this polygon.
#'
#'   Column \code{N_animals} estimates the number of animals from the 
#'   represented population that regularly use the polygon area. If no value for
#'    (at \code{popSize}) is provided, this number is the proportion of the 
#'    represented population using the area.
#'
#'   Column \code{potentialKBA} indicates whether the polygon can be considered 
#'   a potential KBA (TRUE) or not (FALSE).
#'
#' if \code{polyOut = F} function returns a gridded density surface of class 
#' SpatialPixelsDataFrame, with the same three aforementioned columns as cell 
#' values. 
#'
#'  If \code{polyOut = TRUE} the user may choose to automatically produce a plot
#'   of the result using \code{plot=TRUE}. The map produced displays the areas 
#'   which hold aggregations above a certain threshold proportion of the 
#'   population. If there are no areas displayed on the map, then either the 
#'   species doesn't aggregate, the Scale is too small to identify aggregations 
#'   in this species, or the tracked sample aren't representative enough to meet
#'   the thresholds.
#'
#' @examples
#' \dontrun{
#' findKBA(KDE, represent=represent$out)
#' }
#' @export
#' @import dplyr
#' @import sf

findKBA <- function(
  KDE, represent, popSize = NULL, levelUD = 50, polyOut = TRUE){

  if(!class(KDE) %in% c(
    "estUDm", "SpatialPixelsDataFrame", "SpatialGridDataFrame")
    ) {
    stop("KDE should be of class 'estUDm' provided by adehabitatHR::kernelUD or 
      track2kba::estSpaceUse, or an sp class-SpatialPixelsDataFrame or 
      SpatialGridDataFrame.")
  }
  
  # deal with class of KDE input ----------------------------------------------
  if(class(KDE) == "estUDm") {
  KDEpix <- adehabitatHR::estUDm2spixdf(KDE)
  }
  if(class(KDE) %in% c("SpatialPixelsDataFrame", "SpatialGridDataFrame")) {
    KDEpix <- KDE
  }
  SampSize <- ncol(KDEpix)

  # ensure proportion value
  represent <- ifelse(represent > 1, represent/100, represent)


  ### CALCULATE THRESHOLD PROP OF *POPULATION* NEEDED ------------------------
  if (represent < 0.5) message("UNREPRESENTATIVE SAMPLE: below 50% 
    representativeness it sites of importance cannot be identified with 
    confidence")

  thresh <- ifelse(represent <= 0.7, 0.5,
                 ifelse(represent < 0.8, 0.2,
                        ifelse(represent < 0.9, 0.125, 0.1)))
  ### 'correcting' estimates of the proportion of the population in each cell
  corr <- represent

  if(SampSize < 10) {
    warning("LOW SAMPLE SIZE: identifying sites based on <10 tracked individuals
      is not recommended. You may use indEffectTest() to test whether 
      individuals are site faithful between foraging trips 
      (if animal is central-place forager), and if NOT consider using 'tripID' 
      as independent samples instead of individual.")
    if(SampSize < 5) {
      ### if sample size tiny, make it impossible to identify potential KBA --- 
      thresh <- SampSize + 1 
    }
  }

  ###### CONVERTING OUTPUT TO PROPORTIONAL UD FOR EACH INDIVIDUAL--------------
  ## create SpatialPixelsDataFrame
  if(sp::is.projected(KDEpix) != TRUE) stop("Please re-calculate your kernel UD
    after projecting the data into an equal-area projection")

  ## calculate area of each pixel
  pixArea <- KDEpix@grid@cellsize[[1]]
  ## output reported by kernelUD is intensity / m2. This intensity is multiplied
  # by pixel area and sums to 1 for each individual (exceptions near borders)
  ## we sort this calc cumulative sum --> this is effectively the "% UD"

  # if the input was from adehabitatHR (estUDm) convert cell values to 0-1 ----
  if(class(KDE) == "estUDm"){

  KDEpix@data <- KDEpix@data %>%
    mutate(rowname = seq_len(nrow(KDEpix@data))) %>%
    tidyr::gather(key = "ID", value = "UD", -.data$rowname) %>%
    mutate(usage = .data$UD * (pixArea^2)) %>%
    arrange(.data$ID, desc(.data$usage)) %>%
    group_by(.data$ID) %>%
    mutate(cumulUD = cumsum(.data$usage)) %>%
    dplyr::select(.data$rowname, .data$ID, .data$cumulUD) %>%
    arrange(.data$rowname) %>%
    tidyr::spread(key = .data$ID, value = .data$cumulUD) %>%
    dplyr::select(-.data$rowname)
  
  }

  ### COUNT THE NUMBER OF OVERLAPPING UD KERNELS >levelUD ---------------------
  ## convert pixels to 1 if they are below levelUD and 0 otherwise  
  # then sum the number of overlapping 1s (i.e individuals)
  . <- NULL # makes R CMD Check happy
  Noverlaps <- KDEpix
  Noverlaps@data <- as.data.frame(
    ifelse(Noverlaps@data < (levelUD / 100), 1, 0)
    ) %>%
    mutate(N_IND = rowSums(.)) %>%
    dplyr::select(.data$N_IND)

  ### Classify each cell as POTENTIAL (KBA) or not based on thres and corr ----
  potentialKBA <- Noverlaps

  ### Introduce population size -----------------------------------------------
  if(is.null(popSize)){
    potentialKBA@data$N_animals <- (corr * (potentialKBA@data$N_IND / SampSize))
    message(
    "No value for colony size provided. Output for N_animals is in % of colony 
      size"
      )
    potentialKBA@data <- potentialKBA@data %>%
      mutate(potentialKBA = ifelse( .data$N_animals >= thresh, TRUE, FALSE) )
    } else {   ## provide the number of ind expected if colony size is given
    potentialKBA@data$N_animals <- corr * popSize * (potentialKBA@data$N_IND / SampSize)
    potentialKBA@data <- potentialKBA@data %>%
      mutate(potentialKBA = ifelse(
        (.data$N_animals/popSize) >= thresh, TRUE, FALSE) 
        )
    } 
  Noverlaps <- NULL

  if(polyOut==TRUE){
      
      #### CONVERT OUTPUT INTO POLYGONS WITH KBA ASSESSMENT INFO --------------
      # slow conversion
      KBApoly <- as(potentialKBA, "SpatialPolygonsDataFrame")

      potentialKBA <- NULL
    
      ### aggregate all pixel-polygons with the same number of animals
      OUTMAP <- raster::aggregate(
        KBApoly, 
        c('N_animals','N_IND','potentialKBA')
        )
      KBApoly <- NULL
    
      ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
      KBA_sf <- sf::st_as_sf(OUTMAP) %>%
        sf::st_union(by_feature = TRUE) %>%
        smoothr::smooth(method = "densify") %>%
        sf::st_transform(4326) %>%
        arrange(.data$N_IND)
      OUTMAP <- NULL

    return(KBA_sf)
      
  } else {
    return(potentialKBA)
    }
  
} 



