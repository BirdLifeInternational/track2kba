## projectTracks ###############################################################

#' Project tracking data
#'
#' \code{projectTracks} is a convenience function to project tracking data to a
#' an equal-area projection for use in kernel density analysis.
#'
#' Data are transformed to either a World Cylindrical Equal Area, or a Lambert
#' equal-area projection. Cylindrical projections generally appear better for
#' data that are distributed more along one axis, while azimuthal appear better
#' for data that is distributed evenly within a radius. The most important thing
#'  is that the data are in an equal-area projection for Kernel Density Analysis
#'  (e.g. \code{\link{estSpaceUse}}).
#'
#' If \code{custom=TRUE}, the projection will be centered on the data. This
#' is particularly preferable for data that cross the international dateline,
#' or near the poles. However, it is important to recognize that this projection
#' is specific to inpute dataset (i.e. \code{dataGroup}) so if
#' \code{projectTracks} is run again with even slightly different data, the
#' projections will differ, which may cause issues down the line if merging
#' spatial datasets again.
#'
#' NOTE that these projections may not be the most appropriate for your data and
#' it is almost certainly better to manually identify a projection appropriate
#' for your study region. Custom projections are centered on the centroid of
#' the tracking locations, which is biased for locations close to the poles. In
#' this case it would be better identify an appropriate polar projection for
#' your study are instead of relying on \code{projectTracks}. So it is not
#' strictly necessary for \code{projectTracks} to be used in track2KBA analysis,
#'  what is important is that an equal-area projection of some kind is used when
#'  constructing utilization distributions.
#'
#' @param dataGroup data.frame or SpatialPointsDataFrame. Tracking data, with
#' fields as named by \code{\link{formatFields}}. Must contain 'Latitude' and
#' 'Longitude' columns.
#' @param projType character. Select type of equal-area projection to use. Two
#' options are are available: 'cylin' projects to a World Cylindrical Equal Area
#' projection, and 'azim' projects to a Lambert Azimuthal EA.
#' @param custom logical (TRUE/FALSE). Choose whether projection will use
#' default centering parameters or whether to set projection center on centroid
#' of latitude and longitude in dataGroup.
#'
#' Input data can be tracks split into trips (i.e. output of
#' \code{\link{tripSplit}})
#'
#'
#' @return Returns a SpatialPointsDataFrame, which can be used for the following
#'  functions: \code{\link{findScale}}, \code{\link{estSpaceUse}},
#'  \code{\link{indEffectTest}}, \code{\link{repAssess}}
#'
#' @seealso \code{\link{tripSummary}}
#'
#' @examples
#'
#' tracks_raw <- track2KBA::boobies
#'
#' ## format data
#' tracks_formatted <- formatFields(
#'   dataGroup = tracks_raw,
#'   fieldID   = "track_id",
#'   fieldLat  ="latitude",
#'   fieldLon  ="longitude",
#'   fieldDate ="date_gmt",
#'   fieldTime ="time_gmt"
#'   )#'
#'   
#' tracks_prj <- projectTracks(
#' tracks_formatted, 
#' projType = "azim", 
#' custom = "TRUE"
#' )
#'
#'
#' @export
#' @importFrom geosphere centroid
#' @importFrom sp CRS SpatialPoints SpatialPointsDataFrame spTransform

projectTracks <- function(dataGroup, projType, custom) {

  if (!"Latitude" %in% names(dataGroup)) {
    stop("dataGroup missing Latitude field")}
  if (!"Longitude" %in% names(dataGroup)) {
    stop("dataGroup missing Latitude field")}

  mid_point <- data.frame(
    geosphere::centroid(cbind(dataGroup$Longitude, dataGroup$Latitude))
  )

  if (projType == "azim") {
    if (custom == TRUE) {
      proj <- sp::CRS(
        paste(
          "+proj=laea +lon_0=", mid_point$lon,
          " +lat_0=", mid_point$lat, sep = ""
        ))
      message("NOTE: projection is data specific")
    } else {
      proj <- sp::CRS("+proj=laea +lat_0=0 +lon_0=-0")
      message("NOTE: projection center default used, which is at 0 Lat 0 Lon. If
        your data fall far from this location, the shape of your data will be
        highly distorted. Either set 'custom=T' or use a region-specific EA
        projection of your choosing (best option).")
    }
  } else if (projType == "cylin") {
    if (custom == TRUE) {
      proj <- sp::CRS(
        paste(
          "+proj=cea +lon_0=", mid_point$lon,
          " +lat_ts=", mid_point$lat,  # latitude of true scale
          " +x_0=0 +y_0=0", sep = ""
        ))
      message("NOTE: projection is data specific")
    } else {
      proj <- sp::CRS("+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0")
      message("NOTE: projection center default used, which is at 0 Lat 0 Lon. If
        your data fall far north or south of this location (e.g. near the poles)
        the shape of your data will be highly distorted, either set 'custom=T'
        or use a region-specific EA projection of your choosing (best option).")
    }
  }

  if (class(dataGroup) != "SpatialPointsDataFrame") {

    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE ----------------
    if (min(dataGroup$Longitude) < -170 & max(dataGroup$Longitude) > 170) {
      longs <- ifelse(
        dataGroup$Longitude < 0, dataGroup$Longitude + 360, dataGroup$Longitude
        )
      mid_point$lon <- ifelse(
        median(longs) > 180, median(longs) - 360, median(longs)
        )
      }

    dataGroup.Wgs <- sp::SpatialPoints(
      data.frame(dataGroup$Longitude, dataGroup$Latitude),
      proj4string = sp::CRS(SRS_string = "EPSG:4326")
      )
    Tracks_prj <- sp::spTransform(dataGroup.Wgs, CRSobj = proj)
    Tracks_prj <- sp::SpatialPointsDataFrame(Tracks_prj, data = dataGroup)

  } else {
    ## if SPDF and not projected, project -------------------------------------
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE ---------------
    if (
      min(dataGroup@data$Longitude) < -170 & max(dataGroup@data$Longitude) > 170
      ) {
      longs <- ifelse(
        dataGroup@data$Longitude < 0,
        dataGroup@data$Longitude + 360, dataGroup@data$Longitude
        )
      mid_point$lon <- ifelse(
        median(longs) > 180, median(longs) - 360, median(longs)
        )
      }

    Tracks_prj <- sp::spTransform(dataGroup, CRSobj = proj)
  }
  return(Tracks_prj)
}
