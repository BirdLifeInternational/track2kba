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
#' @return Returns a list object of length two, containing tracking data (accessed using: \code{dataset$data}) and study site location information (accessed using: \code{dataset$site}) .
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
  if (!requireNamespace("move", quietly = TRUE)) {
    stop("Package \"move\" needed for this function to work. Please install it.",
      call. = FALSE)
  }
  
  ### IMPORT FROM MOVEBANK IF CREDENTIALS SUPPLIED
  if (any(!is.null(MovebankID), !is.null(User), !is.null(Password))) {
    loginStored <- move::movebankLogin(username=User, password=Password)

    if(is.character(MovebankID)) # convert character study name to numeric MovebankID
    {MovebankID <- move::getMovebankID(MovebankID, loginStored)}

    ### IMPORT FROM MOVEBANK
    input <- move::getMovebankData(study=MovebankID, login=loginStored)
    try(study <- move::getMovebank(entity_type="study", login=loginStored, study_id=MovebankID), silent=FALSE)

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input@data %>% dplyr::select(.data$deployment_id, .data$timestamp, .data$location_lat, .data$location_long) %>%
      rename(ID=.data$deployment_id, DateTime=.data$timestamp, Latitude=.data$location_lat, Longitude=.data$location_long) %>%
      arrange(.data$ID, .data$DateTime)

    ### IMPORT FROM FILE IF NO LOGIN IS PROVIDED ###
  } else {
    if(is.null(filename)) stop("No filename provided, and one of the credentials for Movebank login is missing. Please provide a numeric Movebank ID, User and Password")
    input <- utils::read.csv(filename, stringsAsFactors = F)

    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input %>% dplyr::select(.data$`individual-local-identifier`, .data$timestamp, .data$`location-lat`, .data$`location-long`) %>%
      rename(ID=.data$`individual-local-identifier`, DateTime=.data$timestamp, Latitude=.data$`location-lat`, Longitude=.data$`location-long`) %>%
      arrange(.data$ID, .data$DateTime)

  }

  ### EXTRACT COLONY INFORMATION FROM EITHER TRACKS OR FIRST LINE IN TRACKS ###
  if('study' %in% ls()){
    Colony <- study %>% dplyr::select(.data$main_location_long, .data$main_location_lat) %>% rename(Longitude=.data$main_location_long, Latitude=.data$main_location_lat )
  }else{
    Colony <- tracks[1,] %>% dplyr::select(.data$Longitude, .data$Latitude)
    warning("Study had no information on colony location. Used first location of tracks as colony location. If a different location is desired, please manually specify.")}

  tracks_list <- list(data=tracks, site=Colony)

  return(tracks_list)
  #return(list(tracks=tracks, Colony=Colony))
}
