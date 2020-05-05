## projectTracks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Project tracking data
#'
#' \code{projectTracks} Projects tracking data to a custom equal-area projection for use in track2KBA analysis.
#'
#' @param Tracks data.frame. Tracking data, with fields named as named by \code{\link{tripSummary}}.
#' @param Reproject logical (TRUE/FALSE). If your Tracks dataframe is already projected, would you like to reproject these to a custom equal-area projection?
#' Data are transformed to a lambert equal-area projection with it's center determined by the data. Data must contain 'Latitude' and 'Longitude' columns. It is not strictly necessary for this projection to be used in track2KBA analysis, what is important is that an equal-area projection, of any kind, is used when constructing kernel density estimates.
#' 
#' @return Returns a SpatialPointsDataFrame, which can be used for the following functions: \code{\link{findScale}}, \code{\link{estSpaceUse}},  \code{\link{IndEffectTest}},  \code{\link{repAssess}} .
#'
#'@seealso \code{\link{tripSummary}}
#'
#' @examples
#' \dontrun{
#' 
#' data(boobies)
#' 
#' tracks <- formatFields(boobies, BLformat=TRUE)
#' 
#' tracks_prj <- project(tracks)
#' 
#' }
#'
#' @export

projectTracks <- function(Tracks, Reproject=FALSE){
  
  mid_point <- data.frame(geosphere::centroid(cbind(Tracks$Longitude, Tracks$Latitude)))
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
  
  if(class(Tracks)!= "SpatialPointsDataFrame") {

    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(Tracks$Longitude) < -170 &  max(Tracks$Longitude) > 170) {
      longs <- ifelse(Tracks$Longitude < 0, Tracks$Longitude + 360, Tracks$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
    
    Tracks.Wgs <- SpatialPoints(data.frame(Tracks$Longitude, Tracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    Tracks_prj <- spTransform(Tracks.Wgs, CRSobj=proj.UTM )
    Tracks_prj <- SpatialPointsDataFrame(Tracks_prj, data = Tracks)
    
  } else if(is.projected(Tracks) & Reproject == FALSE){
    Tracks_prj <- Tracks
    message("if you wish to reproject these data to custom equal-area projection, set Reproject=TRUE")
  } else { ## if SPDF and not projected, project 
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(Tracks@data$Longitude) < -170 &  max(Tracks@data$Longitude) > 170) {
      longs <- ifelse(Tracks@data$Longitude < 0, Tracks@data$Longitude + 360,Tracks@data$Longitude)
      mid_point$lon<-ifelse(median(longs) > 180, median(longs)-360, median(longs))}
    
    Tracks_prj <- sp::spTransform(Tracks, CRSobj=proj.UTM)
  }
  return(Tracks_prj)
}
