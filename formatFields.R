## formatFields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## mIBA R package, formatFields function VERSION 4: USER SPECIFIES INPUT COLUMNS
## function to find and format the fields of a data frame, ready to input into tripSplit
##
## tripSplit inputs needed:
## Colony data frame, fields: Latitude, Longitude
## tracks data frame, fields: ID, DateTime, TrackTime (made from DateTime), Latitude, Longitude
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



formatFields <- function(df, field_ID, field_Lat, field_Lon, field_Date=NULL, field_Time=NULL, field_DateTime=NULL, format_DT=NULL) {
  
  #### INPUT CHECKS

  ### NEW VERSION: Require people to supply a DateTime column in posixCT format.

  ## check that df is either a data.frame or a data.table
  if (! "data.frame" %in% class(df)) {
    if( ! "data.table" %in% class(df)){
    stop("Object is not a data.table or data.frame. Please try again.")
    }
  }
  
  # If they don't specify field_ID and there is an ID field in dataframe, use that. 
  if(missing(field_ID)){ # if field_ID missing
    if("ID" %in% colnames(df)) { # AND there isn't already an ID field present in the dataframe
      warning("No field_ID was specified, so the pre-existing column named 'ID' was used. If another field is desired as identifier, please specify it in the field_id argument.")
      field_ID <- "ID"
    } else { stop("No field_ID was specified. Please specify the desired IDentifier in the field_ID argument.") } # if no field_ID AND no pre-existing 'ID'
  }

  # ## check that they have supplied all 3 of field_ID, field_Lat, field_Lon
  # if (missing(field_ID) | missing(field_Lat) | missing(field_Lon)){
  #   stop("Please supply field_ID, field_Lat and field_Lon to specify which columns of the data sheet to use for ID, Latitude and Longitude.")
  # } else if (class(field_ID) != "character" | class(field_Lat) != "character" | class(field_Lon) != "character"){
  #   stop("Please supply the field_ID, field_Lat and field_Lon as indivIDual character strings.")
  # }  ## check inputs are right (character) format

  ## ==== OR USE THE FOLLOWING SYNTAX, one for each field =====
  ## print( ifelse( is.null(a), 'a not specified', paste('a =',a) ) ) ## only useful for OPTIONAL arguments
  
  #### convert df to data.frame instead of data.table
  df <- as.data.frame(df)
  
  ## check that field_ID, field_Lat, field_Lon supplied are actually names of fields in the dataframe.
  if (! field_ID %in% colnames(df)){
    stop("The field_ID supplied does not exist in the data frame.")
  } else if (! field_Lat %in% colnames(df)){
    stop("The field_Lat supplied does not exist in the data frame.")
  } else if (! field_Lon %in% colnames(df)){
    stop("The field_Lon supplied does not exist in the data frame.")
  }

  #### Date and Time, or DateTime?
  ## check: have the supplied both Date and Time AND DateTime?
  if ( ( ( (! is.null(field_Date)) & (! is.null(field_Time)) ) & (! is.null(field_DateTime)) ) ) { ## IF there IS BOTH DateTime AND Date and Time
    warning("You have supplied a DateTime field as well as a Date field and a Time field. Using the DateTime field, discarding Date and Time fields. If you wish to use the Date and Time fields instead, please remove the DateTime fields from the function input arguments.")
    field_Date <- NULL
    field_Time <- NULL
  }
  
  
  ## ============= OPTION 1 - DateTime supplied ===============
  if (! is.null(field_DateTime)) {
    if(! is.POSIXct(df[, field_DateTime])) {
      warning("Column supplied to the 'field_DateTime' is not of class POSIXct, the function will attempt to convert it.")
      if(is.null(format_DT)){
        warning("No format was supplied for the input DateTime field, a default format ('ymd_HMS') was attempted. If an error is produced, see help page ('?lubridate::parse_date_Time') for information on Date formats.")
        format_DT <- "ymd_HMS"
        df$DateTime <- lubridate::parse_date_Time(df[, field_DateTime], format_DT)
        }
      df$DateTime <- lubridate::parse_date_Time(df[, field_DateTime], format_DT)
    } else {
    ## TODO: Attempt to format DateTime correctly!
    df <- df %>% dplyr::rename(DateTime=field_DateTime)
    }
  }
  
  
  ## =========== OPTION 2 - Date and Time, or only Date supplied ============
  ## MB ## added conditions to handle user input of format of Date, or Date and Time columns (if not set, a default is tried)
  if (is.null(field_DateTime) & ! is.null(field_Date) & ! is.null(field_Time)) { # if both Date and Time supplied
    
    if(is.null(format_DT)){     # if format of DateTime/(Date + Time) field(s) not supplied, set to default "ymd_HMS"
      warning("No format was supplied for the input Date and Time fields, a default format ('ymd_HMS') was attempted when combining the fields. If an error is produced, see help page ('?lubridate::parse_date_time') for information on date formats.")
      format_DT <- "ymd_HMS"
      df$DateTime <- lubridate::parse_date_time(paste(df[, field_Date], df[, field_Time]), format_DT)
    } else {
    df$DateTime <- lubridate::parse_date_time(paste(df[, field_Date], df[, field_Time]), format_DT)
    }
  } else {                                                                       # if only Date supplied (and Date column not missing)
    if(! is.null(field_Date)){
      warning("Only a Date column (field_Date) was supplied, this will be used to create the DateTime column. If you have a Time column, please indicate which it is in the 'field_Time' argument. ")
      if(is.null(format_DT)){   # if format of Date field not supplied, set to default "ymd"
        warning("No format was supplied for the input Date field, a default format ('ymd') was attempted. If an warning that 'no formats are found' is produced, see help page ('?lubridate::parse_date_time') for information on Date formats.")
        format_DT <- "ymd"
        df$DateTime <- lubridate::parse_date_time(df[, field_Date], format_DT)
      } else {
     df$DateTime <- lubridate::parse_date_time(df[, field_Date], format_DT)
      }
    }
  }
  
  ## ============= checks complete ================
  
  #### FORMAT COLUMNS

  df <- df %>% dplyr::rename(ID=field_ID, Latitude=field_Lat, Longitude=field_Lon) %>% 
    mutate(ID = as.character(ID))
  
  
  ## if (DateTime column does not exist){
    ## make DateTime column
    ## remove Date and Time columns
    ## }

  ## Make TrackTime column.
  

  ## ============== outputs to user ==================
  ## output a confirmation message saying that 'ID field is track_ID' etc.
  
  return(df)
}













