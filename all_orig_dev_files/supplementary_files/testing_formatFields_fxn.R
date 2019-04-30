#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### testing Lizzie's formatCols fxn ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("formatFields.R")

############## Trying different data structures ##################

### for a dataset with an already defined DateTime field
formatCols(tracks, field_id = "ID", field_lat="latitude", field_lon="longitude", field_datetime="date_time")

formatCols(tracks, field_id = "ID", field_lat="latitude", field_lon="longitude", field_datetime="date_time", field_date = "date_time", field_time="date_time")

## does not work this way! Function requires character strings (i.e. colnames)
# formatCols(tracks, field_id = tracks[1], field_lat=tracks[4], field_lon=tracks[3], field_datetime=tracks[2])

##### both date and time field supplied #####
head(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time="time_gmt"))
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time="time_gmt"))

## only a date column supplied
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time=NULL))

## supplying a date column to the field_datetime argument DOES NOT WORK...
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_datetime="date_gmt"))
## unless the format_DT is set correctly!
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_datetime="date_gmt", format_DT = "ymd"))


##### Attempting to specificy the format of the date and time fields date and time field supplied #####
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time="time_gmt", format_DT="ymd_HMS"))
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", format_DT="ymd"))
str(formatCols(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt"))