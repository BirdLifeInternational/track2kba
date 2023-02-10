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
#' @param colorBy character string. Either "complete" if trips are to be
#' coloured as complete or incomplete, or "trip" if trips are to be coloured by
#' trip ID.
#' @return Returns a figure of facetted maps, each of which corresponds to a
#' level of ID in \emph{trips}.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @examples
#' ## make some play data
#' dataGroup <- data.frame(Longitude = rep(c(1:10, 10:1), 2), 
#'                         Latitude =  rep(c(1:10, 10:1), 2),
#'                         ID = c(rep("A", 20), rep("B", 20)),
#'                         DateTime = format(
#'                         lubridate::ymd_hms("2021-01-01 00:00:00") +
#'                         lubridate::hours(0:19))
#' )
#' colony <- data.frame(
#'  Longitude = dataGroup$Longitude[1], Latitude = dataGroup$Latitude[1]
#' )
#' Trips <- tripSplit(dataGroup,
#'                    colony=colony,
#'                    innerBuff=2,
#'                    returnBuff=20,
#'                    duration=1,
#'                    nests = FALSE,
#'                    rmNonTrip = TRUE
#' )
#' ## Visualize trips
#' mapTrips(Trips, colony)                   # add colony location to each facet
#' mapTrips(Trips, colony, colorBy = "trip") # color trips by their order
#' 
#' @export
#' @importFrom ggplot2 ggplot aes scale_x_continuous geom_path geom_point
#' @importFrom ggplot2 facet_wrap vars theme element_rect element_blank
#' @importFrom ggplot2 aes_string labs
#' @importFrom rlang .data
#' @importFrom dplyr filter

mapTrips <- function(trips, colony, IDs=NULL, colorBy = c("complete", "trip")) {

  if (is.null(IDs)) {
    IDs <- 1:dplyr::n_distinct(trips@data$ID)
  }

  N_ceiling <- ifelse(
    n_distinct(trips@data$ID) > 25, 25, n_distinct(trips@data$ID)
    )

  if (length(IDs) > N_ceiling) {
    IDs <- IDs[1:N_ceiling]
    message("Too many individuals to plot. Only 25 IDs will be shown.")
  }
  if (nrow(colony) > 1) {
    if (!"ID" %in% colnames(colony)) {
      message("ID column missing from colony object")
      }
    colony <- dplyr::filter(colony, .data$ID %in% IDs)
  }
  selectIDs <- unique(trips@data$ID)[IDs]
  plotdat <- trips@data %>% dplyr::filter(.data$ID %in% selectIDs)
  if (length(colorBy) == 2) {
    coldat <- "complete"
    message("Trips colored by completeness. Indicate colorBy=='trip' to color by
    trips.")
  } else if (colorBy == "complete") {coldat <- "complete"
  } else if (colorBy == "trip") {coldat <- "colID"}

  ##### DIFFERENT PLOT FOR BIRDS CROSSING THE DATELINE ------------------------
  if (min(plotdat$Longitude) < -170 & max(plotdat$Longitude) > 170) {
    plotdat <- plotdat %>%
      mutate(
        Longitude = ifelse(
          .data$Longitude < 0, .data$Longitude + 360, .data$Longitude
        )
      )

    colony$Longitude <- ifelse(
      colony$Longitude < 0, colony$Longitude + 360, colony$Longitude
    )

    longlimits <- c(min(plotdat$Longitude) - 2, max((plotdat$Longitude) + 2))
    longbreaks <- round(seq(longlimits[1], longlimits[2], by = 10) / 10, 0) * 10
    longlabels <- ifelse(longbreaks > 180, longbreaks - 360, longbreaks)

    plotdat %>%
      mutate(complete = ifelse(.data$Returns == "No", "No", "Yes"),
             colID = format(
               x = factor(
                 x = .data$tripID,
                 labels = seq_len(length.out = n_distinct(.data$tripID))
                 ))) %>%
      arrange(.data$ID, .data$DateTime) -> forplot

    TRACKPLOT <- ggplot(
      data = forplot,
      aes_string(x = "Longitude", y = "Latitude", col = coldat)
      ) +
      geom_path() +
      geom_point(data  = colony,
                 aes(x = .data$Longitude, y = .data$Latitude),
                 fill  = "dark orange", color = "black", pch = 21, size = 2
      ) +
      scale_x_continuous(limits = longlimits,
                         breaks = longbreaks,
                         labels = longlabels) +
      facet_wrap(~.data$ID)  +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(colour = "black", fill = "white"),
            panel.border = element_blank()) +
      labs(col = if_else(coldat == "colID",  "Trip", "Complete"))
  } else {
    plotdat %>%
      dplyr::arrange(.data$ID, .data$DateTime) %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::mutate(complete = ifelse(.data$Returns == "No", "No", "Yes"),
                    colID = format(x = factor(
        x = .data$tripID,
        labels = seq_len(length.out = n_distinct(.data$tripID))
      ))) -> forplot
    TRACKPLOT <- ggplot(data = forplot,
                        aes_string(x = "Longitude", y = "Latitude",
                                   col = coldat)
      ) +
      geom_path() +
      geom_point(
        data = colony, aes(x = .data$Longitude, y = .data$Latitude),
        fill = "dark orange", color = "black", pch = 21, size = 2
      ) +
      facet_wrap(~.data$ID) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(colour = "black", fill = "white"),
            panel.border = element_blank()) +
      labs(col = if_else(coldat == "colID",  "Trip", "Complete"))
  }

  base::print(TRACKPLOT)

}
