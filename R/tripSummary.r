## tripSummary   ###############################################################

#' Summary of trip movements
#'
#' \code{tripSummary} provides a simple summary of foraging trip distances,
#' durations, and directions performed by central place foraging animals.
#'
#' \emph{nests}=T may be used if it is desired, for example, to use specific
#' nest locations instead of one central location for all individuals/tracks.
#'
#' @param trips SpatialPointsDataFrame, as produced by \code{\link{tripSplit}}.
#' @param colony data.frame with 'Latitude' and 'Longitude' columns specifying
#' the locations of the central place (e.g. breeding colony). If
#' \code{nests=TRUE}, \code{colony} should have a third column, 'ID'
#' with corresponding character values in the 'ID' field in \emph{trips}.
#' @param nests logical scalar (TRUE/FALSE). Were central place
#' (e.g. deployment) locations used in \code{tripSplit} specific to each unique
#'  'ID'? If so, each place must be matched with an 'ID' value in both
#'  \code{trips} and \code{colony} objects.
#' @param extraDist logical scalar (TRUE/FALSE). If TRUE, the distance between
#' the first and last points of each trip and the \code{colony} will be added to
#' the 'total_dist' (total distance travelled) for each trip.
#'
#' @return Returns a tibble data.frame grouped by ID. Trip characteristics
#' included are trip duration (in hours), maximum distance and cumulative
#' distance travelled (in kilometers), direction (in degrees, measured from
#' origin to furthest point of track), start and end times as well as a unique
#' trip identifier ('tripID') for each trip performed by each individual in the
#'  data set. Distances are calculated on a great circle.
#'
#' If the beginning of a track is starts out on a trip which is followed by only
#'  one point within \emph{InnerBuff}, this is considered an 'incomplete' trip
#'  and will have an NA for duration. If an animal leaves on a trip but does not
#'  return within the \emph{ReturnBuff} this will be also classified an
#'  'incomplete trip'.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @examples 
#' ## make some play data
#' dataGroup <- data.frame(Longitude = rep(c(1:10, 10:1), 2), 
#'                Latitude =  rep(c(1:10, 10:1), 2),
#'                ID = c(rep("A", 20), rep("B", 20)),
#'                DateTime = format(
#'                lubridate::ymd_hms("2021-01-01 00:00:00") + 
#'                lubridate::hours(0:19))
#' )
#' 
#' colony <- data.frame(
#' Longitude = dataGroup$Longitude[1], Latitude = dataGroup$Latitude[1]
#' )
#' ## split tracks into trips
#' trips <- tripSplit(dataGroup, colony=colony, 
#'                 innerBuff = 1, 
#'                 returnBuff = 1, 
#'                 duration = 0.5,
#'                 rmNonTrip = FALSE
#' )
#' ## summarise trip characteristics
#' sumTrips <- tripSummary(trips, colony)
#'
#' @export
#' @importFrom dplyr first last group_by ungroup if_else left_join mutate

tripSummary <- function(trips, colony=NULL, nests=FALSE, extraDist=FALSE)
  {

  if (!"Latitude" %in% names(colony)) stop("colony missing Latitude field")
  if (!"Longitude" %in% names(colony)) stop("colony missing Longitude field")

  ### helper function to calculate distance unless no previous location -------
  poss_dist <- purrr::possibly(geosphere::distm, otherwise = NA)

  if(inherits(trips, "SpatialPointsDataFrame")) {
    trips <- as.data.frame(trips@data)
  } else { trips <- trips }

  ## summaries ----------------------------------------------------------------
  trip_distances <- trips %>%
    tidyr::nest(coords = c(.data$Longitude, .data$Latitude)) %>%
    group_by(.data$tripID) %>%
    mutate(prev_coords = dplyr::lag(.data$coords)) %>%
    ungroup() %>%
    mutate(Dist = purrr::map2_dbl(
      .data$coords, .data$prev_coords, poss_dist)
      ) %>%
    mutate(Dist = if_else(is.na(.data$Dist), .data$ColDist, .data$Dist)) %>%
    mutate(count = 1) %>%
    group_by(.data$ID, .data$tripID) %>%
    dplyr::summarise(n_locs = sum(.data$count),
                     departure = min(.data$DateTime),
                     return = max(.data$DateTime),
                     duration = ifelse("No" %in% unique(.data$Returns), NA,
                                        as.numeric(
                                          difftime(max(.data$DateTime),
                                                   min(.data$DateTime),
                                                   units = "hours")
                                          )
                ),
                total_dist = sum(.data$Dist, na.rm = TRUE) / 1000,
                max_dist = max(.data$ColDist) / 1000) %>%
    mutate(
      direction = 0,
      duration  = ifelse(.data$duration == 0, NA, .data$duration),
      complete  = ifelse(
        is.na(.data$duration), "incomplete trip", "complete trip"
        ),
      complete = ifelse(.data$tripID == "-1", "non-trip", .data$complete)
      )
  if (extraDist == TRUE) { # add dist btwn colony+1st+last pnts of trips to tot
    extra_dist <- trips %>%
      group_by(.data$tripID) %>%
      dplyr::summarise(
        firstlast_dist = (first(.data$ColDist) + last(.data$ColDist)) / 1000
    )

    trip_distances <- trip_distances %>%
      left_join(extra_dist, by = "tripID") %>%
      mutate(
        total_dist = .data$total_dist + .data$firstlast_dist
    ) %>%
      dplyr::select(-.data$firstlast_dist)
  }

  ### LOOP OVER EACH TRIP TO CALCULATE DIRECTION TO FURTHEST POINT FROM COLONY
  for (i in unique(trip_distances$tripID)) {
    if (i == "-1") {
      trip_distances$direction <- rep(NA)
    } else {
      x <- trips[trips$tripID == i, ]

      maxdist <- cbind(
        x$Longitude[x$ColDist == max(x$ColDist)],
        x$Latitude[x$ColDist == max(x$ColDist)]
      )
      if (dim(maxdist)[1] > 1) {maxdist <- maxdist[1, ]}

      if (nests == TRUE) {
        origin <- colony[match(unique(x$ID), colony$ID), ] %>%
          dplyr::select(.data$Longitude, .data$Latitude)
      } else {
        origin <- colony
        b <- geosphere::bearing(c(origin$Longitude, origin$Latitude), maxdist)
      }

      ## great circle (ellipsoidal) bearing of trip ---------------------------
      b <- geosphere::bearing(c(origin$Longitude, origin$Latitude), maxdist)
      ## convert the azimuthal bearing to a compass direction -----------------
      trip_distances$direction[trip_distances$tripID == i] <- (b + 360) %% 360
    }
  }

if ("incomplete trip" %in% trip_distances$complete) warning(
  "Some trips did not return to the specified returnBuffer distance from the
  colony. The return DateTime given for these trips refers to the last location
  of the trip, and NOT the actual return time to the colony.")

return(trip_distances)
}
