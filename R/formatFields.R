## formatFields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Format tracking data for track2BKA analysis
#'
#' \code{formatFields} formats the column names of a data frame so that they are accepted by track2KBA functions.
#'
#' If data are already in format of BirdLife Seabird tracking database (\url{www.seabirdtracking.org}), use \code{BL_format == TRUE} and formatting conversion will occur automatically. 
#'
#' By matching up the names of your existing columns with those recognized by track2KBA functions, \code{formatFields} re-formats the data frame, and converts the date/date-time fields into a singe date-time field of class POSIXct.
#' 
#'
#' If date-time is combined in a single column, please use \emph{field_DateTime} instead of \emph{field_Date} and \emph{field_Time}.
#'
#' @param tracks data.frame or data.table.
#' @param BL_format logical. Is data set already in format of BirdLife Seabird tracking database? If so, indicate \emph{field_Time} and ignore following arguments. 
#' @param field_ID character. Unique identifier; e.g. for individuals or tracks.
#' @param field_Lat numeric. Name of column corresponding to the LATITUDINAL positions.
#' @param field_Lon numeric. Name of column corresponding to the LONGITUDINAL positions.
#' @param field_DateTime character. If existing, this is the name of the column corresponding to the combined DATE & TIME.
#' @param field_Date character. Name of column corresponding to the DATE only.
#' @param field_Time character. Name of column corresponding to the TIME only.
#' @param format_DT character. What is the format of the data in your DateTime, Date, and Time columns (e.g. "ymd_HMS")? Specify the format following the standard in \code{\link[lubridate]{parse_date_time}}.
#' @param cleanDF logical scalar (T/F). Should columns which are non-essential for track2KBA analysis be removed from dataframe, or not? Removal will speed analysis up a bit. 
#'
#' @return Returns a data.frame, with 'ID', 'Latitude', 'Longitude', and 'DateTime' (class POSIXct) columns.
#'
#' @examples
#' \dontrun{
#' ## using data as formatted on \url{www.seabirdtracking.org},
#' tracks_formatted <- formatFields(tracks=tracks_raw, BL_format=TRUE)
#' 
#' ## using data with user-custom format
#' i.e. with separate Date and Time fields
#'  tracks_formatted <- formatFields(
#'  tracks=tracks_raw, 
#'  field_ID = "ID", 
#'  field_Lat="lat", 
#'  field_Lon="long", 
#'  field_Date="dateGMT", 
#'  field_Time="timeGMT"
#'  )
#'
#' ## using data with only single Date field
#' tracks_formatted <- formatFields(
#' tracks=tracks_raw, 
#' field_Lat="lat", 
#' field_Lon="lon", 
#' field_Date="Date", 
#' format_DT = "dmy")
#' }
#'
#' @export
#'
formatFields <- function(tracks, BL_format=FALSE, field_ID, field_Lat, field_Lon,  field_DateTime=NULL, field_Date=NULL, field_Time=NULL, format_DT=NULL, cleanDF=FALSE) {

  #### INPUT CHECKS
  
  ## check that df is either a data.frame or a data.table
  if (! "data.frame" %in% class(tracks)) {
    if( ! "data.table" %in% class(tracks)){
    stop("Object is not a data.table or data.frame. Please try again.")
    }
  }
  #### convert df to data.frame instead of data.table
  tracks <- as.data.frame(tracks)
  
  #### if data are already in Seabird Database format use that 
  
  if (BL_format == TRUE){
    field_ID   <- "track_id"
    field_Lat  <- "latitude"
    field_Lon  <- "longitude"
    field_Date <- "date_gmt"
    field_Time <- "time_gmt"
  }
  
  # If user doesn't specify field_ID and there is an ID field in dataframe, use that.
  if(missing(field_ID)){ # if field_ID missing
    if("ID" %in% colnames(tracks)) { # AND there isn't already an ID field present in the dataframe
      warning("No field_ID was specified, so the pre-existing column named 'ID' was used. If another field is desired as identifier, please specify it in the field_id argument.")
      field_ID <- "ID"
    } else { stop("No field_ID was specified. Please specify the desired IDentifier in the field_ID argument.") } # if no field_ID AND no pre-existing 'ID'
  }

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
  # IF there IS BOTH DateTime AND Date and Time
  if ( ( ( (! is.null(field_Date)) & (! is.null(field_Time)) ) & (! is.null(field_DateTime)) ) ) {
    message("Both a DateTime, and separate Date and Time fields supplied. If you wish to use the separate Date and Time fields, please remove the DateTime field from the function input argument 'field_DateTime'.")
    field_Date <- NULL
    field_Time <- NULL
  }

  ## ============= OPTION 1 - DateTime supplied ===============
  if (! is.null(field_DateTime)) {
    if(! lubridate::is.POSIXct(tracks[, field_DateTime])) {
      message("Column supplied to the 'field_DateTime' is not of class POSIXct, the function will attempt to convert it.")
      if(is.null(format_DT)){
        message("No format was supplied for the input DateTime field, a default format ('ymd_HMS') was attempted. If an error is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
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

    if( is.null(format_DT) ){     # if format of DateTime/(Date + Time) field(s) not supplied, set to default "ymd_HMS"
      if(BL_format == FALSE) {message("No format was supplied for the Date and Time fields, a default format ('ymd_HMS') was attempted when combining the fields. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")}
      format_DT <- "ymd_HMS"
      tracks$DateTime <- lubridate::parse_date_time(paste(tracks[, field_Date], tracks[, field_Time]), format_DT, tz = "UTC")
    } else {
    tracks$DateTime <- lubridate::parse_date_time(paste(tracks[, field_Date], tracks[, field_Time]), format_DT, tz = "UTC")
    }
  } else {                                                                       # if only Date supplied (and Date column not missing)
    if(! is.null(field_Date)){
      message("Only a Date column (field_Date) was supplied, this will be used to create the DateTime column. If you have a Time column, please indicate which it is in the 'field_Time' argument. ")
      if(is.null(format_DT)){   # if format of Date field not supplied, set to default "ymd"
        message("No format was supplied for the input Date field, a default format ('ymd') was attempted. If an warning that 'no formats are found' is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
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

  if(cleanDF==TRUE){
    tracks <- tracks %>%
      dplyr::select(.data$ID, .data$Latitude, .data$Longitude, .data$DateTime) %>%
      arrange(.data$ID, .data$DateTime)
  } 
  return(tracks)
}
