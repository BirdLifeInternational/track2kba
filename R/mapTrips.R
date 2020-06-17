## mapTrips  ##################################################################

#' Make simple maps of foraging trips 
#'
#' \code{mapTrips} uses output from \code{tripSplit} to create maps illustrating
#'  movements for each ID. 
#'
#' This function only works with the output of \code{tripSplit}.
#'
#'
#' @param trips SpatialPointsDataFrame. Must be output of 
#' \code{\link{tripSplit}} function).
#' @param colony data.frame. Containing 'Latitude' and 'Longitude' fields 
#' specifying the central location(s) from which trips begin. If more than one 
#' location, each row should correspond to an appropriate location (Lat/Lon) for
#'  each ID value in \emph{trips}.
#' @return Returns a figure of facetted maps, each of which corresponds to a 
#' level of ID in \emph{trips}.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @examples
#' \dontrun{trips <- mapTrips(trips, colony)}
#' @export
#' @importFrom ggplot2 ggplot aes scale_x_continuous geom_path geom_point
#' @importFrom ggplot2 facet_wrap vars theme element_rect element_blank
#' @importFrom rlang .data

mapTrips <- function(trips, colony){
  
  if(length(unique(trips@data$ID))>25){
    selectIDs <- unique(trips@data$ID)[1:25]
    plotdat <- trips@data %>% dplyr::filter(.data$ID %in% selectIDs)
    message("Too many individuals to plot. Only first 25 ID's will be shown")
  } else {plotdat <- trips@data}
  TRACKPLOT <- plotdat %>% 
    dplyr::mutate(complete=ifelse(.data$Returns=="No","No","Yes")) %>%
    dplyr::arrange(.data$ID, .data$TrackTime) %>%
    ggplot(
      aes(.data$., x=.data$Longitude, y=.data$Latitude, col=.data$complete)
      ) +
    geom_path() +
    geom_point(
      data=colony, aes(x=.data$Longitude, y=.data$Latitude), 
      fill='dark orange', color='black', pch=21, size=2
    ) +
    facet_wrap(~ID) +
    theme(panel.background=element_rect(fill="white", colour="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      panel.border = element_blank())
  
  ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ------------------------
  if (min(trips$Longitude) < -170 &  max(trips$Longitude) > 170) {
    plotdat <-  trips@data %>%
      mutate(
        Longitude=ifelse(
          .data$Longitude<0, .data$Longitude+360, .data$Longitude
          )
        )
    
    colony$Longitude  <- ifelse(
      colony$Longitude<0, colony$Longitude+360, colony$Longitude
      )
    
    longlimits <- c(min(plotdat$Longitude) - 2, max((plotdat$Longitude) + 2))
    longbreaks <- round(seq(longlimits[1], longlimits[2], by=10) / 10, 0) * 10
    longlabels <- ifelse(longbreaks > 180, longbreaks - 360, longbreaks)

    TRACKPLOT <- plotdat %>% 
      mutate(complete = ifelse(.data$Returns=="No", "No", "Yes")) %>%
      arrange(.data$ID, .data$TrackTime) %>% 
      ggplot(.data$., 
        aes(x = .data$Longitude, y = .data$Latitude, col = .data$complete)
        ) +
      geom_path() +
      geom_point(data = colony, 
        aes(x=.data$Longitude, y=.data$Latitude), 
        fill='dark orange', color='black', pch=21, size=2
        ) +
      scale_x_continuous(limits = longlimits,
        breaks = longbreaks,
        labels = longlabels) +
      facet_wrap(~.data$ID)  +
      theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_blank())
    }
  
  base::print(TRACKPLOT)
  
}

