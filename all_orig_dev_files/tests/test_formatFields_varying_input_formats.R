
###### Chize-formatted GLS data (light-mantled albatross)
tracks <- fread("C:/Users/Martim Bill/Documents/political_connectivity/data/extra_data/GLS_Chizé_Birdlife_190211/lmsa/lmsa_locs.csv") 


source("formatFields.R")

### for a dataset with an already defined DateTime field
# tracks <- formatFields(tracks, field_id = "ID", field_lat="latitude", field_lon="longitude", field_datetime="date_time")
str(tracks)
### for a dataset with both date and time fields #####
x <- formatFields(tracks, field_Lat="lat", field_Lon="lon", field_Date="Date", format_DT = "dmy")
