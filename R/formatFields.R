## formatFields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Format tracking data for track2BKA analysis
#'
#' \code{formatFields} formats the column names of a data frame so that they are accepted by track2KBA functions.
#'
#' By matching up the names of your existing columns with those recognized by track2KBA functions, \code{formatFields} re-formats the data frame, and converts the date/date-time fields into a singe date-time field of class POSIXct.
#'
#' If date-time is combined in a single column, please use \emph{field_DateTime} instead of \emph{field_Date} and \emph{field_Time}.
#'
#' @param tracks data.frame or data.table.
#' @param field_ID character. Unique identifier; e.g. for individuals or tracks.
#' @param field_Lat numeric. Name of column corresponding to the LATITUDINAL positions.
#' @param field_Lon numeric. Name of column corresponding to the LONGITUDINAL positions.
#' @param field_DateTime character. If existing, this is the name of the column corresponding to the combined DATE & TIME.
#' @param field_Date character. Name of column corresponding to the DATE only.
#' @param field_Time character. Name of column corresponding to the TIME only.
#' @param format_DT character. What is the format of the data in your DateTime, Date, and Time columns (e.g. "ymd_HMS")? Specify the format following the standard in \code{\link[lubridate]{parse_date_time}}.
#'
#' @return Returns a data.frame, with 'ID', 'Latitude', 'Longitude', and 'DateTime' (class POSIXct) columns.
#'
#' @examples
#' \dontrun{
#' ## using data as formatted on \url{www.seabirdtracking.org},
#' i.e. with separate Date and Time fields
#'  tracks_formatted <- formatFields(tracks, 
#'  field_ID = "track_id", 
#'  field_Lat="latitude", 
#'  field_Lon="longitude", 
#'  field_Date="date_gmt", 
#'  field_Time="time_gmt"
#'  )
#'
#' ## using data with only single Date field
#' tracks_formatted <- formatFields(tracks, 
#' field_Lat="lat", 
#' field_Lon="lon", 
#' field_Date="Date", 
#' format_DT = "dmy")
#' }
#'
#' @export
#'
formatFields <- function(tracks, field_ID, field_Lat, field_Lon,  field_DateTime=NULL, field_Date=NULL, field_Time=NULL, format_DT=NULL) {

  #### INPUT CHECKS

  ## check that df is either a data.frame or a data.table
  if (! "data.frame" %in% class(tracks)) {
    if( ! "data.table" %in% class(tracks)){
    stop("Object is not a data.table or data.frame. Please try again.")
    }
  }

  # If they don't specify field_ID and there is an ID field in dataframe, use that.
  if(missing(field_ID)){ # if field_ID missing
    if("ID" %in% colnames(tracks)) { # AND there isn't already an ID field present in the dataframe
      warning("No field_ID was specified, so the pre-existing column named 'ID' was used. If another field is desired as identifier, please specify it in the field_id argument.")
      field_ID <- "ID"
    } else { stop("No field_ID was specified. Please specify the desired IDentifier in the field_ID argument.") } # if no field_ID AND no pre-existing 'ID'
  }
  
  # IF  they specify field_ID and another field with name "ID" then rename old field ID
  if(field_ID!="ID" & "ID" %in% names(tracks)){
    tracks<-tracks %>% dplyr::rename(oldID=ID)}  
  
  # ## check that they have supplied all 3 of field_ID, field_Lat, field_Lon
  # if (missing(field_ID) | missing(field_Lat) | missing(field_Lon)){
  #   stop("Please supply field_ID, field_Lat and field_Lon to specify which columns of the data sheet to use for ID, Latitude and Longitude.")
  # } else if (class(field_ID) != "character" | class(field_Lat) != "character" | class(field_Lon) != "character"){
  #   stop("Please supply the field_ID, field_Lat and field_Lon as indivIDual character strings.")
  # }  ## check inputs are right (character) format

  ## ==== OR USE THE FOLLOWING SYNTAX, one for each field =====
  ## print( ifelse( is.null(a), 'a not specified', paste('a =',a) ) ) ## only useful for OPTIONAL arguments

  #### convert df to data.frame instead of data.table
  tracks <- as.data.frame(tracks)

  ## check that field_ID, field_Lat, field_Lon supplied are actually names of fields in the dataframe.
  if (! field_ID %in% colnames(tracks)){
    stop("The field_ID supplied does not exist in the data frame.")
  } else if (! field_Lat %in% colnames(tracks)){
    stop("The field_Lat supplied does not exist in the data frame.")
  } else if (! field_Lon %in% colnames(tracks)){
    stop("The field_Lon supplied does not exist in the data frame.")
  }

  #### Date and Time, or DateTime?
  ## check: have the supplied both Date and Time AND DateTime?
  if ( ( ( (! is.null(field_Date)) & (! is.null(field_Time)) ) & (! is.null(field_DateTime)) ) ) { ## IF there IS BOTH DateTime AND Date and Time
    warning("You have supplied a DateTime field as well as a Date field and a Time field. Using the DateTime field, discarding Date and Time fields. If you wish to use the Date and Time fields instead, please remove the DateTime fields from the function input argument field_DateTime.")
    field_Date <- NULL
    field_Time <- NULL
  }


  ## ============= OPTION 1 - DateTime supplied ===============
  if (! is.null(field_DateTime)) {
    if(! lubridate::is.POSIXct(tracks[, field_DateTime])) {
      warning("Column supplied to the 'field_DateTime' is not of class POSIXct, the function will attempt to convert it.")
      if(is.null(format_DT)){
        warning("No format was supplied for the input DateTime field, a default format ('ymd_HMS') was attempted. If an error is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
        format_DT <- "ymd_HMS"
        tracks$DateTime <- lubridate::parse_date_time(tracks[, field_DateTime], format_DT, tz = "UTC")
        }
      tracks$DateTime <- lubridate::parse_date_time(tracks[, field_DateTime], format_DT, tz = "UTC")
    } else {
    tracks <- tracks %>% dplyr::rename(DateTime = field_DateTime)
    }
  }


  ## =========== OPTION 2 - Date and Time, or only Date supplied ============
  ## MB ## added conditions to handle user input of format of Date, or Date and Time columns (if not set, a default is tried)
  if (is.null(field_DateTime) & ! is.null(field_Date) & ! is.null(field_Time)) { # if both Date and Time supplied

    if(is.null(format_DT)){     # if format of DateTime/(Date + Time) field(s) not supplied, set to default "ymd_HMS"
      warning("No format was supplied for the input Date and Time fields, a default format ('ymd_HMS') was attempted when combining the fields. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")
      format_DT <- "ymd_HMS"
      tracks$DateTime <- lubridate::parse_date_time(paste(tracks[, field_Date], tracks[, field_Time]), format_DT, tz = "UTC")
    } else {
    tracks$DateTime <- lubridate::parse_date_time(paste(tracks[, field_Date], tracks[, field_Time]), format_DT, tz = "UTC")
    }
  } else {                                                                       # if only Date supplied (and Date column not missing)
    if(! is.null(field_Date)){
      warning("Only a Date column (field_Date) was supplied, this will be used to create the DateTime column. If you have a Time column, please indicate which it is in the 'field_Time' argument. ")
      if(is.null(format_DT)){   # if format of Date field not supplied, set to default "ymd"
        warning("No format was supplied for the input Date field, a default format ('ymd') was attempted. If an warning that 'no formats are found' is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
        format_DT <- "ymd"
        tracks$DateTime <- lubridate::parse_date_time(tracks[, field_Date], format_DT, tz = "UTC")
      } else {
     tracks$DateTime <- lubridate::parse_date_time(tracks[, field_Date], format_DT, tz = "UTC")
      }
    }
  }

  ## ============= checks complete ================

  #### FORMAT COLUMNS

  tracks <- tracks %>% dplyr::rename(ID=field_ID, Latitude=field_Lat, Longitude=field_Lon) %>%
    dplyr::mutate(ID = as.character(.data$ID))

  return(tracks)
}
