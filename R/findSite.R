## findSite  ##################################################################

#' Delineating sites of potential importance to conservation
#'
#'
#' \code{findSite} uses the core areas (based on utilization distributions) of
#' individual animals to identify areas used regularly used by a significant
#' portion of the local source population (i.e. the tracked population).
#'
#' \code{findSite} estimates the proportion of the local source population using
#' an area based on the proportion of overlap among individual core areas and
#' the degree of representativeness as quantified by \code{\link{repAssess}}).
#' This value is then compared to a threshold of importance (i.e. a certain % of
#'  the population) to delineate areas as 'potentialSites'. Thresholds area
#' either set automatically set on the representativenss of the sample
#' (lower rep==higher threshold), or set manually by the user.
#'
#' The areas identified are sites of ecological relevance to the populations,
#' which may be significant for the wider region or entire species, which cane
#'  be assessed using global (or regional) criteria, such as those of the Key
#' Biodiversity Area program.
#'
#' The KBA criteria for site assessment are published in the KBA standard, which
#'  may be found here: \url{http://www.keybiodiversityareas.org/}.
#'
#' If grid used for estimating core areas (i.e. KDE) is very memory-heavy
#' (e.g. >10,000 cells) use \code{polyOut = FALSE} to speed things up.
#'
#' @param KDE estUDm or SpatialPixels/GridDataFrame. If estUDm, as created by
#' \code{\link{estSpaceUse}} or \code{adehabitatHR::kernelUD}, if Spatial*,
#' each column should correspond to the Utilization Distribution of a single
#' individual or track.
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
#' delineating the core use (or home range) areas of individuals based on the
#' kernel density estimation (e.g core area=50, home range=95).
#' @param thresh Numeric (percentage). Threshold percentage of local source
#' population needed to be found using a location for it to be considered part
#' of a 'potentialSite'. Default is set based on degree of representativeness.
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
#'   represented population that predictably use the polygon area during the
#'   tracked season. If no value for (at \code{popSize}) is provided, this
#'   number is the proportion of the represented population using the area.
#'
#'   Column \code{potentialSite} indicates whether the polygon can be considered
#'    a potential Site (TRUE) or not (FALSE).
#'
#' if \code{polyOut = FALSE} function returns a gridded surface of class
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
#' KDE <- track2KBA::KDE_example
#' 
#' ## identify potential sites
#' pot_site <- findSite(KDE, represent = 90, levelUD = 50)
#' 
#' @export
#' @importFrom adehabitatHR estUDm2spixdf
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom sf st_union

findSite <- function(
  KDE, represent, popSize = NULL, levelUD, thresh, polyOut = FALSE) {

  classKDE <- class(KDE)

  if (!classKDE %in% c(
    "estUDm", "SpatialPixelsDataFrame", "SpatialGridDataFrame")
    ) {
    stop("KDE should be of class 'estUDm' provided by adehabitatHR::kernelUD or
      track2KBA::estSpaceUse, or an sp class-SpatialPixelsDataFrame or
      SpatialGridDataFrame.")
  }

  # deal with class of KDE input ----------------------------------------------
  if (classKDE == "estUDm") {
  KDE <- adehabitatHR::estUDm2spixdf(KDE)
  }

  SampSize <- ncol(KDE)

  # ensure proportion value
  represent <- ifelse(represent > 1, represent / 100, represent)

  ### CALCULATE THRESHOLD PROP OF *POPULATION* NEEDED ------------------------
  if (represent < 0.5) message(
  "UNREPRESENTATIVE SAMPLE: sample below 50% representativeness. Sites of
  importance cannot be identified with confidence")

  if (missing("thresh")) {
    thresh <- ifelse(represent <= 0.7, 0.5,
                     ifelse(represent < 0.8, 0.2,
                            ifelse(represent < 0.9, 0.125, 0.1)))
  } else { thresh <- thresh / 100
      if (1 / ncol(KDE) > thresh) {message(
        "NOTE: Selected 'thresh' is lower than the 1/sample size, which means an
        area could be delineated as important although only visited by one
        tracked individual."
      ) }
    }
  ### 'correcting' estimates of the proportion of the population in each cell
  corr <- represent

  if (SampSize < 10) {
    warning(
    "LOW SAMPLE SIZE: identifying sites based on <10 tracked individuals
    is not recommended. You may use indEffectTest() to test whether
    individuals are site faithful between foraging trips (if animal is a
    central-place forager), and if NOT consider using 'tripID' as independent
     samples instead of individual.")
    if (SampSize < 5) {
      # if sample size tiny, make it impossible to identify potential Site
      thresh <- SampSize + 1
    }
  }

  ###### CONVERTING OUTPUT TO PROPORTIONAL UD FOR EACH INDIVIDUAL--------------
  ## create SpatialPixelsDataFrame
  if (sp::is.projected(KDE) != TRUE) stop("Please re-calculate your kernel UD
    after projecting the data into an equal-area projection")

  ## calculate area of each pixel
  pixArea <- KDE@grid@cellsize[[1]]
  ## output reported by kernelUD is intensity / m2. This intensity is multiplied
  # by pixel area and sums to 1 for each individual (exceptions near borders)
  ## we sort this calc cumulative sum --> this is effectively the "% UD"

  # if the input was from adehabitatHR (estUDm) convert cell values to 0-1 ----
  if (classKDE == "estUDm") {
    KDE@data <- KDE@data %>%
      mutate(rowname = seq_len(nrow(KDE@data))) %>%
      tidyr::pivot_longer(!.data$rowname, names_to = "ID", values_to = "UD") %>%
      mutate(usage = .data$UD * (pixArea^2)) %>%
      arrange(.data$ID, desc(.data$usage)) %>%
      group_by(.data$ID) %>%
      mutate(cumulUD = cumsum(.data$usage)) %>%
      dplyr::select(.data$rowname, .data$ID, .data$cumulUD) %>%
      arrange(.data$rowname) %>%
      tidyr::pivot_wider(names_from = .data$ID, values_from = .data$cumulUD) %>%
      dplyr::select(-.data$rowname)
  }

  ### COUNT THE NUMBER OF OVERLAPPING UD KERNELS >levelUD ---------------------
  ## convert pixels to 1 if they are below levelUD and 0 otherwise
  # then sum the number of overlapping 1s (i.e individuals)
  . <- NULL # makes R CMD Check happy
  Noverlaps <- KDE
  KDE <- NULL
  Noverlaps@data <- as.data.frame(
    ifelse(Noverlaps@data < (levelUD / 100), 1, 0)
  ) %>%
    mutate(N_IND = rowSums(.))

  cols <- colnames(Noverlaps@data)[-which(colnames(Noverlaps@data) == "N_IND")]

  Noverlaps@data$ID_IND <- apply(Noverlaps@data, 1, function(x) {
    paste(na.omit(
      cols[which(x == 1)]
    ), collapse = " ")
  })

  Noverlaps@data <- Noverlaps@data %>% dplyr::select(.data$N_IND, .data$ID_IND)

  ### Classify each cell as POTENTIAL (Site) or not based on thres and corr ----
  potentialSite <- Noverlaps
  Noverlaps     <- NULL

  ### Introduce population size -----------------------------------------------
  if (is.null(popSize)) {
    potentialSite@data$N_animals <- corr * (potentialSite@data$N_IND / SampSize)
    message(
    "No value for population size provided. Output for N_animals is in % of pop
      size"
      )
    potentialSite@data <- potentialSite@data %>%
      mutate(potentialSite = ifelse(.data$N_animals >= thresh, TRUE, FALSE))
    } else {   ## provide the number of ind expected if colony size is given
    potentialSite@data$N_animals <- (
      corr * popSize * (potentialSite@data$N_IND / SampSize)
    )
    potentialSite@data <- potentialSite@data %>%
      mutate(potentialSite = ifelse(
        (.data$N_animals / popSize) >= thresh, TRUE, FALSE)
        )
    }

  if (polyOut == TRUE) {

      #### CONVERT OUTPUT INTO POLYGONS WITH Site ASSESSMENT INFO --------------
      # slow conversion
      ### aggregate all pixel-polygons with the same number of animals
      OUTMAP <- raster::aggregate(
        as(potentialSite, "SpatialPolygonsDataFrame"),
        c("N_animals", "N_IND", "potentialSite")
        )
      potentialSite <- NULL

      ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
      Site_sf <- sf::st_as_sf(OUTMAP) %>%
        sf::st_union(by_feature = TRUE) %>%
        sf::st_transform(4326) %>%
        arrange(.data$N_IND)
      OUTMAP <- NULL

    return(Site_sf)

  } else {
    return(potentialSite)
    }
}
