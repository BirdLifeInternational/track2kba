## mapSite  #####################################################################

#' Make simple maps of aggregation and important sites 
#'
#' \code{mapSite} uses output from \code{findSite} to create maps illustrating 
#' density of animals in space, and borders of potentially important areas for 
#' the population. 
#'
#' If the input is simple features polygons (i.e. \code{polyOut = TRUE} in 
#' \code{findSite}), areas which meet threshold of importance are displayed 
#' (in red) on top of of the estimated density of animals in space. Black 
#' borders are political and coastline borders.If there are no red borders areas
#'  displayed on the map, then either the species doesn't aggregatee enough to 
#'  meet the threshold, or the tracked sample aren't representative enough to 
#'  identify significant aggregations.
#' 
#' If input is SpatialPixelsDataFrame (i.e. \code{polyOut = FALSE} in 
#' \code{findSite}), a simple density surface map is plotted. 
#'
#' @param Site Simple feature MULTIPOLYGON object or SpatialPixelsDataFrame. Must
#'  be output of \code{\link{findSite}} function).
#' @param colony data.frame. Optional. Must contain columns named 'Latitude' and
#'  'Longitude', with coordinate locations to display reference point of, for 
#'  example, a breeding or tagging site.
#' @param show logical. show plot, or just save it. Note, saving plot only works
#'  for Simple Features input. Default is TRUE. 
#' @return Returns a figure of either single map with all core ranges displayed 
#' together, or a series of facetted maps, each of which shows a utilization 
#' distribution corresponding to a level of ID in \emph{KDE}.
#'
#' @seealso \code{\link{estSpaceUse}}
#'
#' @examples
#' \dontrun{kde_maps <- mapKDE(Trips)}
#' @export
#' @importFrom sf st_bbox
#' @importFrom ggplot2 geom_sf coord_sf borders ggplot theme element_rect
#' @importFrom ggplot2 ylab xlab scale_fill_continuous scale_colour_continuous
#' @importFrom ggplot2 geom_point guides aes element_line

mapSite <- function(Site, colony=NULL, show=TRUE) {
  
  if(class(Site)[1] == "sf"){
    if (!requireNamespace("maps", quietly = TRUE)) { stop(
      "Packages \"maps\"is needed for this function to work. Please install it.", 
      call. = FALSE)
    }
    coordsets <- sf::st_bbox(Site)
    
    csf <- ggplot2::coord_sf(
      xlim = c(coordsets$xmin, coordsets$xmax), 
      ylim = c(coordsets$ymin, coordsets$ymax), 
      expand = FALSE
      )
    csf$default <- TRUE
    
    if(any(Site$N_animals > 1)) {
      label <- "N animals"
    } else {
      label <- "Prop. animals"
    }
    
    denseplot <- Site %>% filter(.data$N_animals > 0) %>% ggplot() +
      geom_sf(mapping = aes(fill=.data$N_animals, colour=.data$N_animals)) +
      borders("world", colour="black", fill = NA) +
      csf +
      scale_fill_continuous(high = "#132B43", low = "#56B1F7", name = label) +
      scale_colour_continuous(high = "#132B43", low = "#56B1F7") + 
      theme(panel.background=element_rect(colour = NA, fill="white"),
        panel.grid.major=element_line(colour="transparent"),
        panel.grid.minor=element_line(colour="transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      ylab("Latitude") +  xlab("Longitude") + guides(colour=FALSE)
    # if any areas are potentialSites, add red border
    if(any(Site$potentialSite == TRUE)) {
      potSitearea <- Site %>% group_by(.data$potentialSite) %>% 
        summarise(N_animals = max(.data$N_animals)) %>% 
        filter(.data$potentialSite==TRUE)
      
      denseplot <- denseplot + 
        geom_sf( data=potSitearea, colour="red", fill=NA, size=1.1 ) + 
        csf
    }
    if(!is.null(colony)){ 
      denseplot <- denseplot +
        geom_point(
          data=colony, 
          aes(x=.data$Longitude, y=.data$Latitude), 
          fill='dark orange', color='black', pch=21, size=2.5,
          )
      }
    if(show == TRUE){
      print(denseplot)
    } else { return(denseplot) }
    
  } else if(class(Site) == "SpatialPixelsDataFrame") {
    plot(Site[Site$N_IND>0,])
  }
  
}
