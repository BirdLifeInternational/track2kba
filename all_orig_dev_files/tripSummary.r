## tripSummary   #####################################################################################################################

#' Summary of trip movements
#'
#' \code{tripSummary} provides a simple summary of foraging trips made by central place foraging animals.
#'
#' \emph{Nests}=T may be used if it is desired, for example, to use specific nest locations instead of one central location for all individuals/tracks.
#'
#' @param Trips projected SpatialPointsDataFrame. Specifically, as produced by \code{\link{tripSplit}}.
#' @param Colony data.frame with 'Latitude' and 'Longitude' columns specifying the locations of the central place (e.g. breeding colony). If Nests=TRUE, Colony should have a third column, 'ID' with corresponding character values in the 'ID' field in \emph{Trips}.
#' @param Nests logical scalar (TRUE/FALSE). Were central place locations used in \code{tripSplit} specific to each unique 'ID'? If so, each place must be matched with an 'ID' value in both \emph{Trips} and \emph{Colony} objects.
#'
#' @return a data.frame providing trip duration, distances, and directions for each trip performed by each individual in the data set. Distances are great circle as calculated by \code{\link[geosphere]{distGeo}}.
#'
#' @seealso \code{\link{tripSplit}}
#'
#' @export

tripSummary <- function(Trips, Colony=Colony, Nests=FALSE)
{

  pkgs <-c('sp', 'dplyr', 'geosphere', 'lubridate', 'purrr')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}


  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")


  ### SUMMARISE MAX DIST FROM COLONY AND TRIP TRAVELLING TIME FOR EACH TRIP

  ## helper function to calculate distance unless no previous location
  poss_dist <- possibly(geosphere::distm, otherwise = NA)

  ## all summary in one pipe
  trip_distances <- as.data.frame(Trips@data) %>%
    dplyr::filter(trip_id!=-1) %>%   ### this removes the non-trip locations
    tidyr::nest(Longitude, Latitude, .key = "coords") %>%
    group_by(trip_id) %>%
    mutate(prev_coords = dplyr::lag(coords)) %>%
    ungroup() %>%
    mutate(Dist = map2_dbl(coords, prev_coords, poss_dist)) %>%
    mutate(Dist=if_else(is.na(Dist),ColDist,Dist)) %>%
    mutate(count=1) %>%
    group_by(ID, trip_id) %>%
    summarise(n_locs=sum(count),
      departure=min(DateTime),
      return=max(DateTime),
      duration=ifelse("N" %in% unique(Returns),NA,((max(TrackTime)-min(TrackTime))/3600)),
      total_dist=sum(Dist, na.rm=T)/1000,
      max_dist=max(ColDist)/1000) %>%
    mutate(direction=0) %>%
    mutate(complete=ifelse(is.na(duration),"incomplete trip","complete trip"))


  ### LOOP OVER EACH INDIVIDUAL TRIP TO CALCULATE DIRECTION TO FURTHEST POINT FROM COLONY
  for (i in unique(trip_distances$trip_id)){			### removed as.numeric as this only works with numeric ID
    x <- Trips@data[Trips@data$trip_id==i,]
    maxdist <- cbind(x$Longitude[x$ColDist==max(x$ColDist)], x$Latitude[x$ColDist==max(x$ColDist)])
    if(dim(maxdist)[1]>1){maxdist <- maxdist[1, ]}

    if(Nests == TRUE) {origin <- Colony[match(unique(x$ID), Colony$ID),] %>% dplyr::select(Longitude,Latitude)}else{origin <- Colony}
    b <- bearing(origin, maxdist)			## great circle (ellipsoidal) bearing of trip
    trip_distances$direction[trip_distances$trip_id==i] <- (b + 360) %% 360  ## convert the azimuthal bearing to a compass direction
    #trip_distances$bearingRhumb[trip_distances$trip_id==i]<-bearingRhumb(origin,maxdist) 	## constant compass bearing of trip
  }
  if("incomplete trip" %in% trip_distances$complete) warning("Some trips did not return to the specified return buffer around the colony. The return date given for these trips refers to the last location of the trip, and NOT the actual return time to the colony.")
  return(trip_distances)
}
