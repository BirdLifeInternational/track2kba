## move2kba   #####################################################################################################################

#' Import Movebank data sets for track2KBA analysis
#'
#' \code{move2KBA} imports data from Movebank repository and re-formats them to fit track2KBA functions.
#'
#' This is a wrapper function for functions in \code{move} package to import and format tracking data from Movebank. It also attains study site location data (lat/lons).
#'
#' @param MovebankID character or numeric. Character: full name of the study, as stored on Movebank. Numeric: Movebank ID of the study. Both can be obtained on the Study Details page on Movebank (\url{www.movebank.org}) or with \code{\link[move]{getMovebankID}}.
#' @param User Username associated with your Movebank account.
#' @param Password Password associated with your Movebank username.
#' @param filename character. File path to .csv downloaded from \url{www.movebank.org}.
#'
#' @return A list object of length two, containing tracking data (accessed using: \code{dataset$data}) and study site location information (accessed using: \code{dataset$site}) .
#'
#' @seealso \code{\link[move]{getMovebankData}} for data download, \code{\link[move]{getMovebank}} for study metadata, \code{\link[move]{getMovebankID}} for getting study ID number
#'
#' @examples
#' \dontrun{
#'
#' dataset <- move2KBA(MovebankID=xxxxxxxxx, User="myusername", Password="mypassword")
#'
#' tracks <- dataset$data  ## access tracking data
#' site <- dataset$site    ## access study site coordinates
#'}
#'
#' @export

move2KBA <- function(MovebankID=NULL, User=NULL, Password=NULL, filename=NULL)
  {

  pkgs <- c('move', 'tidyverse', 'lubridate','data.table')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE))}

  ### IMPORT FROM MOVEBANK IF CREDENTIALS SUPPLIED
  if (any(!is.null(MovebankID), !is.null(User), !is.null(Password))) {
    loginStored <- move::movebankLogin(username=User, password=Password)

    if(is.character(MovebankID)) # convert character study name to numeric MovebankID
    {MovebankID <- move::getMovebankID(MovebankID, loginStored)}

    ### IMPORT FROM MOVEBANK
    input <- move::getMovebankData(study=MovebankID, login=loginStored)
    try(study <- move::getMovebank(entity_type="study", login=loginStored, study_id=MovebankID), silent=F)

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input@data %>% dplyr::select(deployment_id, timestamp, location_lat, location_long) %>%
      rename(ID=deployment_id, DateTime=timestamp, Latitude=location_lat, Longitude=location_long) %>%
      arrange(ID, DateTime)
    head(tracks)

    ### IMPORT FROM FILE IF NO LOGIN IS PROVIDED ###
  } else {
    if(is.null(filename)) stop("No filename provided, and one of the credentials for Movebank login is missing. Please provide a numeric Movebank ID, User and Password")
    input <- fread(filename)

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input %>% dplyr::select(`individual-local-identifier`, timestamp, `location-lat`, `location-long`) %>%
      rename(ID=`individual-local-identifier`, DateTime=timestamp, Latitude=`location-lat`, Longitude=`location-long`) %>%
      arrange(ID, DateTime)
    head(tracks)

 }

  ### EXTRACT COLONY INFORMATION FROM EITHER TRACKS OR FIRST LINE IN TRACKS ###
  if('study' %in% ls()){
    Colony <- study %>% dplyr::select(main_location_long, main_location_lat) %>% rename(Longitude=main_location_long, Latitude=main_location_lat )
  }else{
    Colony <- tracks[1,] %>% dplyr::select(Longitude, Latitude)
    warning("Study had no information on colony location. Used first location of tracks as colony location. If a different location is desired, please manually specify.")}

  tracks_list <- list(data=tracks, site=Colony)

  return(tracks_list)
#return(list(tracks=tracks, Colony=Colony))
}
