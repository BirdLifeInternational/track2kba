## mapKDE  ##################################################################

#' Make simple maps of Kernel Density Estimates
#'
#' \code{mapKDE} uses output from \code{estSpaceUse} to create maps illustrating
#'  utilization distributions for each ID.
#'
#' If the input is simple features polygons, these will be displayed for all IDs
#'  on same map. If input estUDm utilization distribution surface, each ID level
#'  gets its own facet displaying the full UD.
#'
#'
#' @param KDE Simple feature MULTIPOLYGON or estUDm object. Must be output of
#' \code{\link{estSpaceUse}} function).
#' @param colony data.frame. Optional.'Latitude' and 'Longitude' locations to
#' display reference point of, for example, a breeding or tagging site.
#' @param show logical. show plot, or just save it. Note, saving plot only works
#'  for Simple Features input. Default is TRUE.
#' @return Returns a figure of either single map with all core ranges displayed
#' together, or a series of facetted maps, each of which shows a utilization
#' distribution corresponding to a level of ID in \emph{KDE}.
#'
#' @seealso \code{\link{estSpaceUse}}
#'
#' @examples
#' ## make some play data
#'dataGroup <- data.frame(Longitude = c(1, 1.01, 1.02, 1.04, 1.05, 1.03, 1), 
#'                        Latitude =  c(1, 1.01, 1.02, 1.03, 1.021, 1.01, 1),
#'                        ID = rep("A", 7),
#'                        DateTime = format(
#'                          lubridate::ymd_hms("2021-01-01 00:00:00") + 
#'                          lubridate::hours(0:6)
#'                          )
#' )
#' ## project tracks
#'tracks <- projectTracks(dataGroup, projType = "azim", custom = TRUE)
#' ## estimate utilization distributions for each track
#'KDE <- estSpaceUse(tracks, scale=10, levelUD = 50)
#' ## map it
#'mapKDE(KDE)
#'
#' @export
#' @importFrom maps map
#' @importFrom sf st_bbox
#' @importFrom ggplot2 geom_sf coord_sf borders ggplot theme element_rect
#' @importFrom ggplot2 theme element_rect ylab xlab
#' @importFrom graphics image

mapKDE <- function(KDE, colony=NULL, show=TRUE) {
  ### Polygon data ###
  if (inherits(KDE, "sf")) {

    coordsets <- sf::st_bbox(KDE)
    UDPLOT <- ggplot(KDE) + geom_sf(data = KDE, aes(col = .data$id), fill = NA) +
      coord_sf(
        xlim = c(coordsets$xmin, coordsets$xmax),
        ylim = c(coordsets$ymin, coordsets$ymax), expand = TRUE) +
      borders("world", colour = "black", fill = NA) +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
      ylab("Latitude") +
      xlab("Longitude")
    if (!is.null(colony)) {
      UDPLOT <- UDPLOT +
        geom_point(
          data = colony,
          aes(x = .data$Longitude, y = .data$Latitude),
          fill = "dark orange", color = "black", pch = 21, size = 2.5,
        )
    }
    if (show == TRUE) {
      print(UDPLOT)
    } else { return(UDPLOT) }
  } else if (class(KDE)[1] == "estUDm") {
    print(image(KDE))
  } else { stop(
    "KDE must be either KDE.Surface or UDPolygons layer of estSpaceUse output"
    )
  }

}
