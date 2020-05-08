## formatFields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Format tracking data
#'
#' \code{formatFields} formats the column names of a data frame so that they are accepted by track2KBA functions.
#'
#' If data are already in format of BirdLife Seabird tracking database (\url{www.seabirdtracking.org}), use \code{formatBL = TRUE} and formatting conversion will occur automatically. 
#'
#' By matching up the names of your existing columns with those recognized by track2KBA functions, \code{formatFields} re-formats the data frame, and converts the date/date-time fields into a singe date-time field of class POSIXct.
#' 
#'
#' If date-time is combined in a single column, please use \emph{fieldDateTime} instead of \emph{fieldDate} and \emph{fieldTime}.
#'
#' @param dataGroup data.frame or data.table.
#' @param formatBL logical. Is data set already in format of BirdLife Seabird tracking database? If so, indicate \emph{fieldTime} and ignore following arguments. 
#' @param fieldID character. Unique identifier; e.g. for individuals or dataGroup.
#' @param fieldLat numeric. Name of column corresponding to the LATITUDINAL positions.
#' @param fieldLon numeric. Name of column corresponding to the LONGITUDINAL positions.
#' @param fieldDateTime character. If existing, this is the name of the column corresponding to the combined DATE & TIME.
#' @param fieldDate character. Name of column corresponding to the DATE only.
#' @param fieldTime character. Name of column corresponding to the TIME only.
#' @param formatDT character. What is the format of the data in your DateTime, Date, and Time columns (e.g. "ymd_HMS")? Specify the format following the standard in \code{\link[lubridate]{parse_date_time}}.
#' @param cleanDF logical scalar (T/F). Should columns which are non-essential for track2KBA analysis be removed from dataframe, or not? Removal will speed analysis up a bit. 
#'
#' @return Returns a data.frame, with 'ID', 'Latitude', 'Longitude', and 'DateTime' (class POSIXct) columns.
#'
#' @examples
#' \dontrun{
#' ## using data as formatted on \url{www.seabirdtracking.org},
#' tracks_formatted <- formatFields(dataGroup=tracks_raw, formatBL=TRUE)
#' 
#' ## using data with user-custom format
#' i.e. with separate Date and Time fields
#'  tracks_formatted <- formatFields(
#'  dataGroup=tracks_raw, 
#'  fieldID = "ID", 
#'  fieldLat="lat", 
#'  fieldLon="long", 
#'  fieldDate="dateGMT", 
#'  fieldTime="timeGMT"
#'  )
#'
#' ## using data with only single Date field
#' tracks_formatted <- formatFields(
#' dataGroup=tracks_raw, 
#' fieldLat="lat", 
#' fieldLon="lon", 
#' fieldDate="Date", 
#' formatDT = "dmy")
#' }
#'
#' @export
#'
formatFields <- function(dataGroup, formatBL=FALSE, fieldID, fieldLat, fieldLon,  
  fieldDateTime=NULL, fieldDate=NULL, fieldTime=NULL, formatDT=NULL, cleanDF=FALSE) {

  #### INPUT CHECKS
  
  ## check that df is either a data.frame or a data.table
  if (! "data.frame" %in% class(dataGroup)) {
    if( ! "data.table" %in% class(dataGroup)){
    stop("Object is not a data.table or data.frame. Please try again.")
    }
  }
  #### convert df to data.frame instead of data.table
  dataGroup <- as.data.frame(dataGroup)
  
  #### if data are already in Seabird Database format use that 
  
  if (formatBL == TRUE){
    fieldID   <- "track_id"
    fieldLat  <- "latitude"
    fieldLon  <- "longitude"
    fieldDate <- "date_gmt"
    fieldTime <- "time_gmt"
  }
  
  # If user doesn't specify fieldID and there is an ID field in dataframe, use that.
  if(missing(fieldID)){ # if fieldID missing
    if("ID" %in% colnames(dataGroup)) { # AND there isn't already an ID field present in the dataframe
      warning("No fieldID was specified, so the pre-existing column named 'ID' was used. If another field is desired as identifier, please specify it in the fieldID argument.")
      fieldID <- "ID"
    } else { stop("No fieldID was specified. Please specify the desired IDentifier in the fieldID argument.") } # if no fieldID AND no pre-existing 'ID'
  }

  ## check that fieldID, fieldLat, fieldLon supplied are actually names of fields in the dataframe.
  if (! fieldID %in% colnames(dataGroup)){
    stop("The fieldID supplied does not exist in the data frame.")
  } else if (! fieldLat %in% colnames(dataGroup)){
    stop("The fieldLat supplied does not exist in the data frame.")
  } else if (! fieldLon %in% colnames(dataGroup)){
    stop("The fieldLon supplied does not exist in the data frame.")
  }

  #### Date and Time, or DateTime?
  ## check: have the supplied both Date and Time AND DateTime?
  # IF there IS BOTH DateTime AND Date and Time
  if ( ( ( (! is.null(fieldDate)) & (! is.null(fieldTime)) ) & (! is.null(fieldDateTime)) ) ) {
    message("Both a DateTime, and separate Date and Time fields supplied. If you wish to use the separate Date and Time fields, please remove the DateTime field from the function input argument 'fieldDateTime'.")
    fieldDate <- NULL
    fieldTime <- NULL
  }

  ## ============= OPTION 1 - DateTime supplied ===============
  if (! is.null(fieldDateTime)) {
    if(! lubridate::is.POSIXct(dataGroup[, fieldDateTime])) {
      message("Column supplied to the 'fieldDateTime' is not of class POSIXct, the function will attempt to convert it.")
      if(is.null(formatDT)){
        message("No format was supplied for the input DateTime field, a default format ('ymd_HMS') was attempted. If an error is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
        formatDT <- "ymd_HMS"
        dataGroup$DateTime <- lubridate::parse_date_time(dataGroup[, fieldDateTime], formatDT, tz = "UTC")
        }
      dataGroup$DateTime <- lubridate::parse_date_time(dataGroup[, fieldDateTime], formatDT, tz = "UTC")
    } else {
    dataGroup <- dataGroup %>% dplyr::rename(DateTime = fieldDateTime)
    }
  }

  ## =========== OPTION 2 - Date and Time, or only Date supplied ============
  ## MB ## added conditions to handle user input of format of Date, or Date and Time columns (if not set, a default is tried)
  if (is.null(fieldDateTime) & ! is.null(fieldDate) & ! is.null(fieldTime)) { # if both Date and Time supplied

    if( is.null(formatDT) ){     # if format of DateTime/(Date + Time) field(s) not supplied, set to default "ymd_HMS"
      if(formatBL == FALSE) {message("No format was supplied for the Date and Time fields, a default format ('ymd_HMS') was attempted when combining the fields. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")}
      formatDT <- "ymd_HMS"
      dataGroup$DateTime <- lubridate::parse_date_time(
        paste(dataGroup[, fieldDate], dataGroup[, fieldTime]),
        formatDT, tz = "UTC"
        )
    } else {
    dataGroup$DateTime <- lubridate::parse_date_time(
      paste(dataGroup[, fieldDate], dataGroup[, fieldTime]),
      formatDT, tz = "UTC"
      )
    }
  } else {                                                                       # if only Date supplied (and Date column not missing)
    if(! is.null(fieldDate)){
      message("Only a Date column (fieldDate) was supplied, this will be used to create the DateTime column. If you have a Time column, please indicate which it is in the 'fieldTime' argument. ")
      if(is.null(formatDT)){   # if format of Date field not supplied, set to default "ymd"
        message("No format was supplied for the input Date field, a default format ('ymd') was attempted. If an warning that 'no formats are found' is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
        formatDT <- "ymd"
        dataGroup$DateTime <- lubridate::parse_date_time(dataGroup[, fieldDate], formatDT, tz = "UTC")
      } else {
     dataGroup$DateTime <- lubridate::parse_date_time(dataGroup[, fieldDate], formatDT, tz = "UTC")
      }
    }
  }

  ## ============= checks complete ================
  #### FORMAT COLUMNS
  dataGroup <- dataGroup %>% dplyr::rename(ID=fieldID, Latitude=fieldLat, Longitude=fieldLon) %>%
    dplyr::mutate(ID = as.character(.data$ID))

  if(cleanDF==TRUE){
    dataGroup <- dataGroup %>%
      dplyr::select(.data$ID, .data$Latitude, .data$Longitude, .data$DateTime) %>%
      arrange(.data$ID, .data$DateTime)
  } 
  return(dataGroup)
}
