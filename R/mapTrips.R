## mapTrips  ##################################################################

#' Make simple maps of foraging trips 
#'
#' \code{mapTrips} uses output from \code{tripSplit} to create maps illustrating movements for each ID. 
#'
#' This function only works with the output of \code{tripSplit}.
#'
#'
#' @param Trips SpatialPointsDataFrame. Must be output of \code{\link{tripSplit}} function).
#' @param Colony data.frame. Containing 'Latitude' and 'Longitude' fields specifying the central location(s) from which trips begin. If more than one location, each row should correspond to an appropriate location (Lat/Lon) for each ID value in \emph{Trips}.
#' @return Returns a figure of facetted maps, each of which corresponds to a level of ID in \emph{Trips}.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @examples
#' \dontrun{Trips <- mapTrips(Trips, Colony)}
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 vars
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank


mapTrips <- function(Trips, Colony){
  
  if(length(unique(Trips@data$ID))>25){
    selectIDs <- unique(Trips@data$ID)[1:25]
    plotdat <-  Trips@data %>% dplyr::filter(.data$ID %in% selectIDs)
    warning("Too many individuals to plot. Only the first 25 ID's will be shown")
  }else{plotdat <- Trips@data}
  
  TRACKPLOT <- plotdat %>% dplyr::mutate(complete=ifelse(.data$Returns=="No","No","Yes")) %>%
    dplyr::arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
    ggplot(aes(.data$., x=.data$Longitude, y=.data$Latitude, col=.data$complete)) +
    geom_path() +
    geom_point(data=Colony, aes(x=.data$Longitude, y=.data$Latitude), col='red', shape=16, size=2) +
    facet_wrap(ggplot2::vars(.data$ID)) +
    theme(panel.background=element_rect(fill="white", colour="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      panel.border = element_blank())
  
  ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ###
  if (min(tracks$Longitude) < -170 &  max(tracks$Longitude) > 170) {
    plotdat <-  Trips@data %>%
      mutate(Longitude=ifelse(.data$Longitude<0, .data$Longitude+360, .data$Longitude))
    
    Colony$Longitude  <- ifelse(Colony$Longitude<0, Colony$Longitude+360, Colony$Longitude)
    longlimits <- c(min(plotdat$Longitude)-2, max((plotdat$Longitude)+2))
    longbreaks <- round(seq(longlimits[1], longlimits[2], by=10)/10,0)*10
    longlabels <- ifelse(longbreaks>180, longbreaks-360, longbreaks)
    
    TRACKPLOT <- plotdat %>% mutate(complete=ifelse(.data$Returns=="No", "No", "Yes")) %>%
      arrange(.data$ID, .data$TrackTime) %>% # filter(ifelse... include condition to only show 20 Ind
      ggplot(.data$., aes(x=.data$Longitude, y=.data$Latitude, col=.data$complete)) +
      geom_path() +
      geom_point(data=Colony, aes(x=.data$Longitude, y=.data$Latitude), col='red', shape=16, size=2) +
      scale_x_continuous(limits = longlimits,
        breaks = longbreaks,
        labels = longlabels) +
      facet_wrap(ggplot2::vars(.data$ID)) +
      theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_blank())}
  
  base::print(TRACKPLOT)
  
}

