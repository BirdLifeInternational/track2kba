## move2kba   #####################################################################################################################

## STEFFEN OPPEL, 2019

## this script allows to pull data from Movebank and format them for the track2kba functions

## User must provide a Movebank study name and a username and password if the data are restricted
## Output is a data frame in a specified format with ID, DateTime, Latitude and Longitude fields

## NEED TO DO: extract colony/deployment location
                                       
move2kba <- function(MovebankID, User, Password, filename)
  {
  
  pkgs <-c('move', 'tidyverse', 'lubridate','data.table')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}
  
  
  ### IMPORT FROM MOVEBANK IF LOGIN IS PROVIDED ###
  if (c('MovebankID','Username','Password') %in% ls()) {
    loginStored <- movebankLogin(username=User, password=Password)
    input <- getMovebankData(study=MovebankID, login=loginStored)
    try(deploy <- getMovebank(entity_type="deployments", login=loginStored,study_id=MovebankID), silent=T)


    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input@data %>% dplyr::select(deployment_id,timestamp,location_lat,location_long) %>%
      rename(ID=deployment_id,DateTime=timestamp,Latitude=location_lat,Longitude=location_long) %>%
      arrange(ID, DateTime)
    head(tracks)
    
    ### EXTRACT COLONY INFORMATION FROM EITHER TRACKS OR FIRST LINE IN TRACKS ###
    # if('deploy' in ls()){
    #   Colony <- deploy %>% dplyr::select(lon,lat)
    # }else{
    # Colony<- tracks[1,] %>% dplyr::select(Longitude, Latitude)
    # warning("Study had no information on colony location. Used first location of tracks as colony location.")}
    # 
  
  ### IMPORT FROM FILE IF NO LOGIN IS PROVIDED ###
  } else {
    input <- fread(filename)
    
    ### EXTRACT THE IMPORTANT COLUMNS AND RENAME
    tracks <- input@data %>% dplyr::select(individual-local-identifier,timestamp,location_lat,location_long) %>%
      rename(ID=individual-local-identifier,DateTime=timestamp,Latitude=location_lat,Longitude=location_long) %>%
      arrange(ID, DateTime)
    head(tracks)
    
    # ### EXTRACT COLONY INFORMATION FROM EITHER TRACKS OR FIRST LINE IN TRACKS ###
    # if('deploy' in ls()){
    #   Colony <- deploy %>% dplyr::select(lon,lat)
    # }else{
    #   Colony<- tracks[1,] %>% dplyr::select(Longitude, Latitude)
    #   warning("Study had no information on colony location. Used first location of tracks as colony location.")}
    # 
  }

  return(tracks)  
#return(list(tracks=tracks, Colony=Colony))
}