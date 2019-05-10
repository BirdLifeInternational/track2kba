## findKBA  #####################################################################################################################

#' Delineating areas of aggregation of tracked animals to identify potential Key Areas for Biodiversity (KBA).
#'
#' \code{findKBA} uses the utilization distributions of individual animals to identify areas of aggregation (i.e. where a large proportion of individuals' core areas overlap).
#'
#' The function first calculates a minimum threshold number of animals necessary for an area to be important, based on the representativeness of the tracking data (as quantified by \code{\link{repAssess}}). \code{findKBA} then summarises the number of individual core UDs overlapping in an area and compares that number against the thresholds.
#' The areas identified are POTENTIAL Key Biodiversity Areas. That is, they are areas of ecological relevance to the species, but must yet be assessed against global criteria (conservation status and global population size of the species) to determine whether they achieve global (or regional) KBA status. To assess the areas against the criteria, a population estimate is needed, from which a maximum number of animals in the population that may use each area is calculated, again adjusted based on the sample representativeness. 
#' 
#' The criteria for site assessment are published in the KBA standard, which may be found here: \url{http://www.keybiodiversityareas.org/what-are-kbas}.
#' 
#' @param KDE.Surface estUDm or SpatialPixels/GridDataFrame. If estUDm, as created by \code{\link{estSpaceUse}} or \code{adehabitatHR::kernelUD}. If Spatial*, each column should correspond to the Utilization Distribution of a single individual or track. Only accepted if the utilization distribution was calculated in a projected coordinate reference system.
#' @param Represent Numeric (between 0-1). Output value provided by \code{\link{repAssess}} which assesses how representative the tracking data are for characterising the space use of the wider population. If this value is <0.7 then a warning will be issued as the data do not meet the representativeness criteria for a KBA.
#' @param Col.size Numeric, the number of individuals breeding or residing at the origin location from where animals were tracked, quantifying the population that the tracking data represent. This number will be used to calculate how many animals use the delineated areas of aggregation. If no value for \code{Col.size} is provided then output will be as the proportion of the population.
#' @param UDLev Numeric (percentage). Specifies the quantile used for delineating the core use areas of individuals based on the kernel density distribution. Default set to 50\% based on Lascelles et al. (2016). For penguins higher values can be accepted, see Dias et al. (2018).
#' @param polyOut Logical. (Default TRUE) Should the output be a polygon dataset (TRUE) or grid of animal densities (FALSE). See 'Value' below for more details.
#' @param plotit Logical. If TRUE then a map of identified areas will be drawn. NOTE: this only works if \code{polyOut = TRUE}
#'
#' @return if \code{polyOut = TRUE} function returns an object of class \code{sf} containing polygon data with three data columns:
#'   Column \code{N_IND} indicates the number of tracked individuals whose core use area (at \code{UDLev}) overlapped with this polygon.
#'
#'   Column \code{N_animals} estimates the number of animals from the represented population that regularly use the polygon area. If no value for (at \code{Col.size}) is provided, this number is the proportion of the represented population using the area.
#'
#'   Column \code{potentialKBA} indicates whether the polygon can be considered a potential KBA (TRUE) or not (FALSE).
#'
#' if \code{polyOut = F} function returns a gridded density surface of class SpatialPixelsDataFrame, with the same three aforementioned columns as cell values. 
#'
#'  If \code{polyOut = TRUE} the user may choose to automatically produce a plot of the result using \code{plotit=TRUE}. The map produced displays the areas which hold aggregations above a certain threshold proportion of the population. If there are no areas displayed on the map, then either the species doesn't aggregate, the Scale is too small to identify aggregations in this species, or the tracked sample aren't representative enough to meet the thresholds.
#'
#' @examples
#' \dontrun{
#' findKBA(KDE.Surface, Represent=Represent$out)
#' }
#' @export
#' @import dplyr
#' @import ggplot2
#' @import sf

findKBA <- function(KDE.Surface, Represent, Col.size = NULL, UDLev = 50, polyOut = TRUE, plotit = TRUE){

  #### LOAD PACKAGES ####
  # pkgs <- c('sp', 'sf','smoothr','raster','tidyverse', 'geosphere', 'adehabitatHR')
  # for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE))}

  if(class(KDE.Surface) == "list") { KDE.Surface <- KDE.Surface$KDE.Surface } 
  if(!class(KDE.Surface) %in% c("estUDm", "SpatialPixelsDataFrame", "SpatialGridDataFrame")) {
    stop("KDE.Surface should be of class 'estUDm' provided by adehabitatHR::kernelUD or track2kba::estSpaceUse, or an sp class-SpatialPixelsDataFrame or SpatialGridDataFrame.")
  }
  
  # if estUDm, convert to SPixDF
  if(class(KDE.Surface) == "estUDm") {
  KDEpix <- adehabitatHR::estUDm2spixdf(KDE.Surface)
  }
  
  SampSize <- ncol(KDEpix)
  
  if(SampSize < 10) warning("LOW SAMPLE SIZE: identifying a KBA based on <10 tracked individuals is not recommended. You may use IndEffectTest() to test whether individuals are site faithful between foraging trips, and if NOT consider using 'tripID' as independent samples instead of individual.")


  #### CALCULATING THRESHOLD OF PROP OF TRACKED ANIMALS NEEDED FROM LASCELLES ET AL. 2016 ####
  Represent <- ifelse(Represent > 1, Represent/100, Represent)   ## convert to proportion if people enter percent value

  #threshlkup<-data.frame(rep=c(0.9,0.8,0.7),thresh=c(10,12.5,20),corr=c(0.9,0.75,0.5))
  if (Represent < 0.7) warning("UNREPRESENTATIVE SAMPLE: you either did not track a sufficient number of birds to characterise the colony's space use or your species does not lend itself to KBA identification due to its dispersed movement")

  thresh <- ifelse(Represent <= 0.7, SampSize * 0.5, # length(KDE.Surface) is number of individuals in dataset
                 ifelse(Represent < 0.8, SampSize * 0.2,
                        ifelse(Represent < 0.9, SampSize * 0.125, SampSize * 0.1)))
  ## Correction factor: 'correcting' estimates of the proportion of the population in a cell, based on the representativeness of the sample
  corr <- ifelse(Represent <= 0.7, 0.25,
                 ifelse(Represent < 0.8, 0.5,
                        ifelse(Represent < 0.9, 0.75, 0.9)))

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### CONVERTING OUTPUT TO PROPORTIONAL UD FOR EACH INDIVIDUAL  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## create SpatialPixelsDataFrame
  # UDLev=50
  # KDEpix <- adehabitatHR::estUDm2spixdf(KDE.Surface)
  if(sp::is.projected(KDEpix) != TRUE) stop("Please re-calculate your kernel UD after projecting the data into a coordinate reference system where units are identical on x- and y-axis")

  ## calculate area of each pixel
  # pixArea <- KDE.Surface[[1]]@grid@cellsize[1]
  pixArea <- KDEpix@grid@cellsize[[1]]
  ## output reported by kernelUD is intensity / m2
  ## this intensity is multiplied by pixel area (in m2)
  ## this usage sums to 1 for each individual (some individuals bordering the grid may not sum to 1)
  ## we sort this usage and calculate the cumulative sum -> this is effectively the "% UD"
  ## to find the 50% UD for an individual, simply use all grid cells where the output value is <0.5

  KDEpix@data <- KDEpix@data %>%
    mutate(rowname = 1:nrow(KDEpix@data)) %>%
    tidyr::gather(key = "ID", value = "UD", -.data$rowname) %>%
    mutate(usage = .data$UD * (pixArea^2)) %>%
    arrange(.data$ID, desc(.data$usage)) %>%
    group_by(.data$ID) %>%
    mutate(cumulUD = cumsum(.data$usage)) %>%
    dplyr::select(.data$rowname, .data$ID, .data$cumulUD) %>%
    arrange(.data$rowname) %>%
    tidyr::spread(key = .data$ID, value = .data$cumulUD) %>%
    dplyr::select(-.data$rowname)


  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #### COUNT THE NUMBER OF OVERLAPPING UD KERNELS ABOVE THE UDLev==50
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

  ## convert pixels to 1 if they are below UDLev and 0 if they are outside this quantile and sums the number of individuals with a '1' in each cell (i.e. number of overlapping individuals)
  Noverlaps <- KDEpix
  Noverlaps@data <- as.data.frame(ifelse(Noverlaps@data < (UDLev / 100), 1, 0)) %>%
    mutate(N_IND = rowSums(.)) %>%
    dplyr::select(.data$N_IND)

  # KDEpix <- NULL

  ## Classify each cell as POTENTIAL (KBA) or not based on threshold
  potentialKBA <- Noverlaps
  potentialKBA@data <- potentialKBA@data %>%
    mutate(potentialKBA = ifelse(.data$N_IND >= thresh, TRUE, FALSE))
  ## 'Correct' N animal estimates in cells of POTENTIAL KBA status by the representativeness-set correction factor
  if(!is.null(Col.size)){
    potentialKBA@data$N_animals <- corr * Col.size * (potentialKBA@data$N_IND / SampSize)
    }else{   ## provide the number of ind expected if colony size is given
    potentialKBA@data$N_animals <- (corr * 100 * (potentialKBA@data$N_IND / SampSize)) / 100
    warning("No value for colony size provided. Output for N_animals is in % of colony size")}   ## if no colony size is given then provide output in proportion of population
  # KDE.Surface <- NULL
  Noverlaps <- NULL

  if(polyOut==TRUE){
      
      ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
      #### CONVERT OUTPUT INTO POLYGONS WITH KBA ASSESSMENT INFO
      ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ### the first step is very slow
      KBApoly <- as(potentialKBA, "SpatialPolygonsDataFrame")
      # KBApot <- subset(KBApoly, KBA=="potential")
      #if(dim(KBApoly@data)[1]==0) stop("No areas are used by a sufficient proportion of individuals to qualify as potential KBA.")
      potentialKBA <- NULL
    
        ### aggregate all pixel-sized polygons into big polygons with the same number of birds
      OUTMAP <- raster::aggregate(KBApoly, c('N_animals','N_IND','potentialKBA'))
      KBApoly <- NULL
    
        ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
      KBA_sf <- sf::st_as_sf(OUTMAP) %>%
        sf::st_union(by_feature = TRUE) %>%
        smoothr::smooth(method = "densify") %>%
        #drop_crumbs(threshold = units::set_units(100, km^2)) %>%
        #fill_holes(threshold = units::set_units(100, km^2)) %>%
        sf::st_transform(4326) %>% 
        arrange(.data$N_IND)
      OUTMAP <- NULL
    
      ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
      #### RETURN SIMPLE FEATURE WITH KBA INFO AS OUTPUT AND PLOT
      ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
      # ### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE
    
    if(plotit == TRUE) {
      coordsets <- sf::st_bbox(KBA_sf)
    
      KBAPLOT <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>%
        ggplot() +
        geom_sf(mapping = aes(fill=N_animals, colour=N_animals)) +
        coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE) +
        borders("world", fill="dark grey", colour="grey20") +
        # geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
        theme(panel.background=element_blank(),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          axis.text=element_text(size=16, colour="black"),
          axis.title=element_text(size=16),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        guides(colour=FALSE) +
        scale_fill_continuous(name = "N animals") +
        ylab("Longitude") +
        xlab("Latitude")
      if(is.null(Col.size)) { ## make legend title percent
        KBAPLOT <- KBAPLOT + scale_fill_continuous(name = "Prop. of animals")
      }
      print(KBAPLOT)
    } ## end plotit=T loop
    return(KBA_sf)
      
  } else {
    return(potentialKBA)
    }
  
} ### end findKBA function loop



