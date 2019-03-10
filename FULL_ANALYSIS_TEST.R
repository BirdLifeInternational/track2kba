############################################################################################################
####### IDENTIFYING MARINE KBAs FROM SEABIRD TRACKING DATA #################################################
############################################################################################################

### full analysis for package development
### started by steffen oppel on 1 March 2019


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#library(maptools)
#require(geosphere)
#require(sp)
#library(rgdal)
library(tidyverse)
library(data.table)
#library(maps)
#library(rgeos)
#library(adehabitatHR)
library(lubridate)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OUR track2iba FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\track2iba")
#setwd("C:/Users/Martim Bill/Documents/track2iba")
# source("tripSplit.r")
# source("tripSummary.r")
# source("scaleARS.r")
# source("bootstrap.r")
# source("batchUD.r")
# source("varianceTest.r")
#source("findIBA.r")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORT DATA FROM MOVEBANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("move2kba.r")
### either from Movebank account
### working only if password is provided
#tracks<-move2kba(MovebankID=654043458,User="Steffen",Password="Fr1gateb1rd30")  ## Grey Seal data

## Galapagos tortoises: 174165487

### or from csv file downloaded from Movebank
tracks<-move2kba(filename="example_data/GreySeal_McConnell_UK.csv")

### CREATE COLONY DATA FRAME
Colony<- tracks[1,] %>% dplyr::select(Longitude,Latitude)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND PREPARE SAMPLE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # tracks <- fread("example_data/Dataset_1004_2019-03-01.csv")     ## MUPE
tracks <- fread("example_data/Dataset_1012_2019-03-01.csv")     ## MABO St Helena
# tracks <- fread("example_data/Dataset_1151_2019-03-01.csv")     ## SHAG
# tracks <- fread("example_data/Dataset_1245_2019-03-01.csv")       ## RAZO
# tracks <- fread("example_data/R56Data.csv")       ## Luke Halpin dateline crossing data set


### CREATE COLONY DATA FRAME

Colony<- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  rename(Longitude=lon_colony,Latitude=lat_colony)


### Convert Dates and Times

tracks <- tracks %>%
  mutate(DateTime = ymd_hms(paste(date_gmt,time_gmt, sep = " "))) %>%
  #mutate(TrackTime = as.double(DateTime)) %>%
  dplyr::select(track_id, latitude, longitude,DateTime) %>%
  rename(ID=track_id,Latitude=latitude,Longitude=longitude)
head(tracks)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSplit FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(tracks)
source("tripSplit.r")
Trips <- tripSplit(tracks, Colony=Colony, InnerBuff=15, ReturnBuff=20, Duration=2, plotit=T, nests = F, rmColLocs = T)
dim(Trips)

# Trips <- Trips[!Trips$trip_id %in% names(which(table(Trips$trip_id) < 5)), ] # remove trips with less than 5 points


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSummary FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("tripSummary.r")
trip_distances <- tripSummary(Trips, Colony = Colony, nests = F)
trip_distances
dim(trip_distances)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN findScale FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trips <- spTransform(Trips, CRS=CRS("+proj=longlat + datum=wgs84")) ## test that it handles non-projected SPDF data

source("findScale_wip.r")
dev.new()
HVALS <- findScale(Trips, 
  ARSscale = T,
  Colony = Colony)
HVALS
HVALS1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN batchUD FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("batchUD_clean.r")
KDE.Surface <- batchUD(DataGroup=Trips[Trips$trip_id != "-1",], Scale = HVALS$mag, UDLev = 50, polyOut=F, Res=2)

plot(KDE.Surface[[1]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE bootstrap FUNCTION AND ASSIGN THRESHOLD FOR IBA IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("bootstrap_NEW.r")

before <- Sys.time()
test_NEW <- bootstrap(Trips, Scale=HVALS$mag, Iteration=100, Res=2, BootTable = F)
Sys.time() - before

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE findIBA FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("findIBA_clean.r")
IBAs <- findIBA(KDE.Surface, representativity=test_NEW$out, Col.size = 2000) ## error here if smoothr not installed!



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE IndEffectTest function to assess whether all trips of an individual should be included
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~













