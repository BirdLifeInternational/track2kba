## mapKBA  ##################################################################

#' Make simple maps of aggregation and important sites 
#'
#' \code{mapKBA} uses output from \code{findKBA} to create maps illustrating density of animals in space, and borders of potentially important areas for the population. 
#'
#' If the input is simple features polygons (i.e. \code{polyOut = TRUE} in \code{findKBA}), areas which meet threshold of importance are displayed (in red) on top of of the estimated density of animals in space. Black borders are political and coastline borders.If there are no red borders areas displayed on the map, then either the species doesn't aggregatee enough to meet the threshold, or the tracked sample aren't representative enough to identify significant aggregations.
#' 
#' If input is SpatialPixelsDataFrame (i.e. \code{polyOut = FALSE} in \code{findKBA}), a simple density surface map is plotted. 
#'
#' @param KBA Simple feature MULTIPOLYGON object or SpatialPixelsDataFrame. Must be output of \code{\link{findKBA}} function).
#' @param Colony data.frame. Optional. Must contain columns named 'Latitude' and 'Longitude', with coordinate locations to display reference point of, for example, a breeding or tagging site.
#' @param Show logical. Show plot, or just save it. Note, saving plot only works for Simple Features input. Default is TRUE. 
#' @return Returns a figure of either single map with all core ranges displayed together, or a series of facetted maps, each of which shows a utilization distribution corresponding to a level of ID in \emph{KDE}.
#'
#' @seealso \code{\link{estSpaceUse}}
#'
#' @examples
#' \dontrun{kde_maps <- mapKDE(Trips)}
#' @export
#' @importFrom sf st_bbox
#' @importFrom ggplot2 geom_sf
#' @importFrom ggplot2 coord_sf
#' @importFrom ggplot2 borders
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 scale_colour_continuous
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 aes

mapKBA <- function(KBA, Colony=NULL, Show=TRUE) {
  
  if(class(KBA)[1] == "sf"){
    ###

    coordsets <- sf::st_bbox(KBA)
    
    csf <- ggplot2::coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE)
    csf$default <- TRUE
    
    if(any(KBA$N_animals > 1)) {
      label <- "N animals"
    } else {
      label <- "Prop. animals"
    }
    
    denseplot <- KBA %>% filter(N_animals > 0) %>% ggplot() +
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
    # if any areas are potentialKBAs, add red border
    if(any(KBA$potentialKBA == TRUE)) {
      potKBAarea <- KBA %>% group_by(potentialKBA) %>% summarise(N_animals = max(N_animals)) %>% filter(potentialKBA==TRUE)
      denseplot <- denseplot + geom_sf(data=potKBAarea, colour="red", fill=NA, size=1.1) + 
        csf
    }
    if(!is.null(Colony)){ 
      denseplot <- denseplot +
        geom_point(data=Colony, aes(x=.data$Longitude, y=.data$Latitude), col='dark orange', shape=16, size=2)
      }
    if(Show == TRUE){
      print(denseplot)
    } else { return(denseplot) }
    
  } else if(class(KBA) == "SpatialPixelsDataFrame") {
    plot(KBA)
  }
  
}
