## cropData   ################################################################
#' crop tracking data to a specified geographic area
#'
#' \code{cropData} employs \code{splitSingleID} to reduce tracking data to locations
#' contained within a user-defined geographic area (e.g. a country, an EEZ, or any other geography).
#'
#' This function removes all locations outside the specified geography
#' and re-labels IDs to split tracks that may have intermittently left the geography.
#'
#'
#' @param dataGroup data.frame. Must contain 'Latitude', 'Longitude', 'ID' and
#' 'DateTime' columns (correct format may be assured using
#' \code{\link{formatFields}} function).
#' @param targetArea character or object of class \code{\link{sf}}. If character, must be a valid country name from Containing 'Latitude' and 'Longitude' fields
#' \code{rnaturalearth::ne_countries()$name}. Please note that this is a convenience function and we are not able to 
#' cater for all possible spellings of certain countries, so if this does not work, please provide a valid 
#' \code{sf} object.
#' @param timeGap numeric (in hours). After how many hours outside the \code{targetArea} will an animal's journey
#' be split into separate trips. When animals leave the geography and return, some locations will be outside
#' the geography and not considered for important site identification. This can lead to problems during data 
#' interpolation to a regular time interval. the \code{timeGap} should be longer than the regular time gap between
#' locations of an animal's trajectory. Defaults to 3 days (72 hours).
#'
#'
#' @examples
#' ## make some play data
#' dataGroup <- data.frame(Longitude = rep(c(1:40, 40:1), 2), 
#'                         Latitude =  rep(c(40:1, 1:40), 2),
#'                         ID = c(rep("A", 80), rep("B", 80)),
#'                         DateTime = as.character(
#'                         lubridate::ymd_hms("2021-01-01 00:00:00") +
#'                         lubridate::days(0:79))
#' )
#' 
#' ## crop data to country of 
#' Trips <- tripSplit(dataGroup,
#'                    targetArea="Libya",
#'                    timeGap=72
#' )
#'                    
#' @export
#' @importFrom rlang .data
#' @importFrom rnaturalearth ne_countries
#' @importFrom rnaturalearthdata
#' @importFrom sf st_join
#' @importFrom sf st_as_sf
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr lag
#' @importFrom dplyr arrange
#'
cropData <- function(
  dataGroup, targetArea, timeGap=72) {
  
  if ("sf" %in% class(targetArea)) {
    message(
      "Cropping data to user-specified geography.")
  } else {
  if (is.character(targetArea)) {
    library(rnaturalearth)
    library(rnaturalearthdata)
    world <- ne_countries(scale = "medium", returnclass = "sf")
    if (targetArea %in% world$name) {
      targetCountry<-targetArea
      targetArea<-world %>% filter(name==targetArea)
      message(
        "Cropping data to user-specified country.")
    }else{
      stop("Country name does not exist. Check spelling of country names in rnaturalearth::ne_countries()$name")
    }
  }else{
    stop("targetArea input should be either an sf object or a country name (as character).")
  }
}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # EXTRACT SUBSET OF DATA FOR SPATIAL UNIT
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # convert into sf object
  DATA_sf <- dataGroup %>%
    dplyr::filter(!is.na(Latitude)) %>%
    dplyr::filter(!is.na(Longitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"),crs = 4326)        # EPSG code for WGS84
  
  # extract locations from a geographic area
  TD<-st_join(DATA_sf, targetArea, join = st_within) %>%
    #dplyr::filter(name==targetCountry)
    dplyr::filter(!is.na(name))
  

  # create function to re-define 'trips' (because Bird_ID becomes discontinuous when birds leave and re-enter the target)
  TDdf<-as.data.frame(TD) %>% dplyr::arrange(ID, DateTime) %>%
    dplyr::mutate(Latitude= unlist(map(TD$geometry,2)),
           Longitude = unlist(map(TD$geometry,1))) %>%
    dplyr::select(ID,Latitude,Longitude, DateTime) %>%
    group_by(ID) %>%
    dplyr::mutate(tdiff=as.numeric(difftime(DateTime,dplyr::lag(DateTime), unit="hours")))
  
  # for each time gap of >72 hrs assign a new tripID
  anim<-unique(TDdf$ID)
  TDdf<-TDdf %>% dplyr::mutate(TripID=1)
  for (a in anim){							
    x<-TDdf %>% dplyr::filter(ID==a) %>%
      dplyr::arrange(DateTime)
    splits<-x$DateTime[x$tdiff>timeGap]
    if(length(splits)>1){
      for(l in 2:length(splits)){
        TDdf$TripID[TDdf$ID==a & TDdf$DateTime>=splits[l]]<-l
      }}
  }
  
  cropped_tracks<-TDdf %>% dplyr::mutate(tripID=paste(ID,TripID)) %>%
    dplyr::select(ID,tripID,Latitude,Longitude, DateTime) %>%
    dplyr::arrange(ID, tripID, DateTime) %>%
    dplyr::ungroup()
  
  return(cropped_tracks)
}