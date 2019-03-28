## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## mIBA R package, formatCols function VERSION 3: USER SPECIFIES INPUT COLUMNS
## function to find and format the fields of a data frame, ready to input into tripSplit
##
## tripSplit inputs needed:
## Colony data frame, fields: Latitude, Longitude
## tracks data frame, fields: ID, DateTime, TrackTime (made from DateTime), Latitude, Longitude
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## testing
# df <- tracks
# class(df)
# head(df)
# 
# field_id <- "track_id"
# field_lat <- "latitude"
# field_lon <- "longitude"
# field_date <- "date_gmt"
# field_time <- "time_gmt"
# field_datetime <- NULL

formatFields <- function(df, field_id, field_lat, field_lon, field_date=NULL, field_time=NULL, field_datetime=NULL, format_DT=NULL) {
  
  #### INPUT CHECKS

  ### NEW VERSION: Require people to supply a datetime column in posixCT format.

  ## check that df is either a data.frame or a data.table
  if (! "data.frame" %in% class(df)) {
    if( ! "data.table" %in% class(df)){
    stop("Object is not a data.table or data.frame. Please try again.")
    }
  }
  ## check that they have supplied all 3 of field_id, field_lat, field_lon
  if (missing(field_id) | missing(field_lat) | missing(field_lon)){
    stop("Please supply field_id, field_lat and field_lon to specify which columns of the data sheet to use for ID, Latitude and Longitude.")
  } else if (class(field_id) != "character" | class(field_lat) != "character" | class(field_lon) != "character"){
    stop("Please supply the field_id, field_lat and field_lon as individual character strings.")
  } 

  ## ==== OR USE THE FOLLOWING SYNTAX, one for each field =====
  ## print( ifelse( is.null(a), 'a not specified', paste('a =',a) ) ) ## only useful for OPTIONAL arguments
  
  #### convert df to data.frame instead of data.table
  df <- as.data.frame(df)
  
  ## check that field_id, field_lat, field_lon supplied are actually names of fields in the dataframe.
  if (! field_id %in% colnames(df)){
    stop("The field_id supplied does not exist in the data frame.")
  } else if (! field_lat %in% colnames(df)){
    stop("The field_lat supplied does not exist in the data frame.")
  } else if (! field_lon %in% colnames(df)){
    stop("The field_lon supplied does not exist in the data frame.")
  }

  #### Date and Time, or DateTime?
  ## check: have the supplied both Date and Time AND DateTime?
  if ( ( ( (! is.null(field_date)) & (! is.null(field_time)) ) & (! is.null(field_datetime)) ) ) { ## IF there IS BOTH datetime AND date and time
    warning("You have supplied a DateTime field as well as a Date field and a Time field. Using the DateTime field, discarding Date and Time fields. If you wish to use the Date and Time fields instead, please remove the DateTime fields from the function input arguments.")
    field_date <- NULL
    field_time <- NULL
  }
  
  
  ## ============= OPTION 1 - DateTime supplied ===============
  if (! is.null(field_datetime)) {
    if(! is.POSIXct(df[, field_datetime])) {
      warning("Column supplied to the 'field_datetime' is not of class POSIXct, the function will attempt to convert it.")
      if(is.null(format_DT)){
        warning("No format was supplied for the input DateTime field, a default format ('ymd_HMS') was attempted. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")
        format_DT <- "ymd_HMS"
        df$DateTime <- lubridate::parse_date_time(df[, field_datetime], format_DT)
        }
      df$DateTime <- lubridate::parse_date_time(df[, field_datetime], format_DT)
    } else {
    ## TODO: Attempt to format DateTime correctly!
    df <- df %>% dplyr::rename(DateTime=field_datetime)
    }
  }
  
  
  ## =========== OPTION 2 - Date and Time, or only Date supplied ============
  ## MB ## added conditions to handle user input of format of Date, or Date and Time columns (if not set, a default is tried)
  if (is.null(field_datetime) & ! is.null(field_date) & ! is.null(field_time)) { # if both date and time supplied
    
    if(is.null(format_DT)){     # if format of datetime/(date + time) field(s) not supplied, set to default "ymd_HMS"
      warning("No format was supplied for the input date and time fields, a default format ('ymd_HMS') was attempted when combining the fields. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")
      format_DT <- "ymd_HMS"
      df$DateTime <- lubridate::parse_date_time(paste(df[, field_date], df[, field_time]), format_DT)
    } else {
    df$DateTime <- lubridate::parse_date_time(paste(df[, field_date], df[, field_time]), format_DT)
    }
  } else {                                                                       # if only date supplied (and date column not missing)
    if(! is.null(field_date)){
      warning("Only a date column (field_date) was supplied, this will be used to create the DateTime column. If you have a time column, please indicate which it is in the 'field_time' argument. ")
      if(is.null(format_DT)){   # if format of date field not supplied, set to default "ymd"
        warning("No format was supplied for the input date field, a default format ('ymd') was attempted. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")
        format_DT <- "ymd"
        df$DateTime <- lubridate::parse_date_time(df[, field_date], format_DT)
      } else {
     df$DateTime <- lubridate::parse_date_time(df[, field_date], format_DT)
      }
    }
  }
  
  ## ============= checks complete ================
  
  #### FORMAT COLUMNS

  df <- df %>% dplyr::rename(ID=field_id, Latitude=field_lat, Longitude=field_lon)
  
  
  ## if (DateTime column does not exist){
    ## make DateTime column
    ## remove Date and Time columns
    ## }

  ## Make TrackTime column.
  

  ## ============== outputs to user ==================
  ## output a confirmation message saying that 'ID field is track_id' etc.
  
  return(df)
}













