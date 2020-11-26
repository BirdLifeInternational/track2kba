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
#' @param IDs numeric vector. Sequence of numeric indices for the IDs you wish
#' to map. Max of 25. 
#' @param colorBy character string. Either "complete" if trips are to be coloured as complete or incomplete, or "trip" if trips are to be coloured by trip ID. 
#' @return Returns a figure of facetted maps, each of which corresponds to a 
#' level of ID in \emph{trips}.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @examples
#' \dontrun{
#' trips <- mapTrips(trips, colony)
#' trips <- mapTrips(trips, colony, IDs = 2:10) # show IDs  #2 through #10
#' }
#' @export
#' @importFrom ggplot2 ggplot aes scale_x_continuous geom_path geom_point
#' @importFrom ggplot2 facet_wrap vars theme element_rect element_blank
#' @importFrom ggplot2 aes_string labs
#' @importFrom rlang .data

mapTrips <- function(trips, colony, IDs=NULL,  colorBy = c("complete", "trip")){
  
  if(is.null(IDs)){
    IDs <- 1:dplyr::n_distinct(trips@data$ID)
  }
  
  N_ceiling <- ifelse(
    n_distinct(trips@data$ID) > 25, 25, n_distinct(trips@data$ID)
    )
  
  if(length(IDs) > N_ceiling) {
    IDs <- IDs[1:N_ceiling]
    message("Too many individuals to plot. Only 25 IDs will be shown.")
  }
  
  selectIDs <- unique(trips@data$ID)[IDs]
  plotdat <- trips@data %>% dplyr::filter(.data$ID %in% selectIDs)
  if(length(colorBy)==2){
    coldat <- "complete"
    message("Trips colored by completeness. Indicate colBy=='trip' to color by
    trips.")
  } else if(colorBy == "complete"){coldat <- "complete"
  } else if(colorBy == "trip") {coldat <- "colID"}

  ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ------------------------
  if(min(plotdat$Longitude) < -170 & max(plotdat$Longitude) > 170) {
    plotdat <- plotdat@data %>%
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
    
    plotdat %>% 
      mutate(complete = ifelse(.data$Returns=="No", "No", "Yes"),
             colID = as.character(
               x = factor(
                 x = tripID, 
                 labels = seq_len(length.out = n_distinct(x = tripID))
                 ))) %>%
      arrange(.data$ID, .data$DateTime) -> forplot 
      
    TRACKPLOT <- ggplot(
      data = forplot,
      aes_string(x = "Longitude", y = "Latitude", col = coldat)
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
            panel.border = element_blank()) + 
      labs(col = if_else(coldat == "colID",  "Trip", "Complete"))
  } else {
    plotdat %>% 
      dplyr::arrange(.data$ID, .data$DateTime)  %>%
      dplyr::group_by(.data$ID) %>% 
      dplyr::mutate(complete=ifelse(.data$Returns=="No","No","Yes"),
                    colID = as.character(x = factor(
        x = tripID, 
        labels = seq_len(length.out = n_distinct(x = tripID))
      ))) -> forplot
      
    TRACKPLOT <- ggplot(data = forplot, 
                        aes_string(x = "Longitude", y = "Latitude", col = coldat)
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
            panel.border = element_blank()) + 
      labs(col = if_else(coldat == "colID",  "Trip", "Complete"))
  }
  
  base::print(TRACKPLOT)
  
}

