## move2kba   ##################################################################

#' Import Movebank data sets for track2KBA analysis
#'
#' \code{move2KBA} imports data from Movebank repository and re-formats them to
#' fit track2KBA functions.
#'
#' This is a wrapper function for functions in \code{move} package to import and
#'  format tracking data from Movebank. It also attains study site location data
#'   (lat/lons).
#'
#' @param movebankID character or numeric. Character: full name of the study, as
#'  stored on Movebank. Numeric: Movebank ID of the study. Both can be obtained 
#'  on the Study Details page on Movebank (\url{www.movebank.org}) or with 
#'  \code{\link[move]{getMovebankID}}.
#' @param user Username associated with your Movebank account.
#' @param password password associated with your Movebank username.
#' @param filename character. File path to .csv downloaded 
#' from \url{www.movebank.org}.
#'
#' @return Returns a list object of length two, containing tracking data 
#' (accessed using: \code{dataset$data}) and study site location information 
#' (accessed using: \code{dataset$site}) .
#'
#' @seealso \code{\link[move]{getMovebankData}} for data download, 
#' \code{\link[move]{getMovebank}} for study metadata, 
#' \code{\link[move]{getMovebankID}} for getting study ID number
#'
#' @examples
#' \dontrun{
#'
#' dataset <- move2KBA(movebankID=xxx, user="myusername", password="mypassword")
#'
#' tracks <- dataset$data  ## access tracking data
#' site <- dataset$site    ## access study site coordinates
#'}
#'
#' @export
#' @importFrom magrittr %>%

move2KBA <- function(movebankID=NULL, user=NULL, password=NULL, filename=NULL)
{
  if (!requireNamespace("move", quietly = TRUE)) {
    stop("Package \"move\" needed for this function to work. Please install.",
      call. = FALSE)
  }
  
  ### IMPORT FROM MOVEBANK IF CREDENTIALS SUPPLIED ----------------------------
  if (any(!is.null(movebankID), !is.null(user), !is.null(password))) {
    loginStored <- move::movebankLogin(username=user, password=password)

    if(is.character(movebankID)) # convert study name to numeric movebankID
    {movebankID <- move::getMovebankID(movebankID, loginStored)}

    ### IMPORT FROM MOVEBANK
    input <- move::getMovebankData(study=movebankID, login=loginStored)
    try(
      study <- move::getMovebank(entity_type="study", 
      login=loginStored, study_id=movebankID), 
      silent=FALSE
      )

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME ------------------------------
    tracks <- input@data %>% 
      dplyr::select(
        .data$deployment_id, .data$timestamp, 
        .data$location_lat, .data$location_long) %>%
      rename(ID=.data$deployment_id, DateTime=.data$timestamp,
        Latitude=.data$location_lat, Longitude=.data$location_long) %>%
      arrange(.data$ID, .data$DateTime)

    ### IMPORT FROM FILE IF NO LOGIN IS PROVIDED ------------------------------
  } else {
    if(is.null(filename)) stop(
    "No filename provided, credentials for Movebank login are missing. Please 
    provide a numeric Movebank ID, user and password")
    input <- utils::read.csv(filename, stringsAsFactors = FALSE)

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input %>% 
      dplyr::select(.data$`individual-local-identifier`, .data$timestamp,
        .data$`location-lat`, .data$`location-long`) %>%
      rename(ID=.data$`individual-local-identifier`, DateTime=.data$timestamp, 
        Latitude=.data$`location-lat`, Longitude=.data$`location-long`) %>%
      arrange(.data$ID, .data$DateTime)

  }

  ### EXTRACT COLONY INFORMATION FROM EITHER TRACKS OR FIRST LINE IN TRACKS ---
  if('study' %in% ls()){
    Colony <- study %>% 
      dplyr::select(.data$main_location_long, .data$main_location_lat) %>% 
      rename(
        Longitude=.data$main_location_long, Latitude=.data$main_location_lat 
        )
  }else{
    Colony <- tracks[1,] %>% dplyr::select(.data$Longitude, .data$Latitude)
    message("Study had no information on colony location. 
      Used first locations as colony location. 
      If a different location desired, please manually specify.")}

  tracks_list <- list(data=tracks, site=Colony)

  return(tracks_list)
  #return(list(tracks=tracks, Colony=Colony))
}
