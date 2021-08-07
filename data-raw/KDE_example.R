## code to prepare `KDE_example` dataset goes here

tracks_raw <- track2KBA::boobies

## format data
tracks_formatted <- formatFields(
  dataGroup = tracks_raw,
  fieldID   = "track_id",
  fieldLat  ="latitude",
  fieldLon  ="longitude",
  fieldDate ="date_gmt",
  fieldTime ="time_gmt"
)


library(dplyr)
tracks_formatted <- dplyr::filter(
  tracks_formatted, ID %in% c("69324", "69302", "69343", "69304")) %>% 
  dplyr::filter(row_number() %% 40 == 1)


## project dataset
tracks_prj <- projectTracks(
  tracks_formatted,
  projType = "azim",
  custom = "TRUE"
)

## get utilization distributions
KDE_example <- estSpaceUse(tracks_prj, scale = 30, res = 10, levelUD = 50)

usethis::use_data(KDE_example, overwrite = TRUE)
