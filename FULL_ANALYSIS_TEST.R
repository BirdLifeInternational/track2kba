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
setwd("C:/Users/Martim Bill/Documents/track2iba")
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
tracks<-move2kba(MovebankID=114336340,User="Steffen",Password="xxx")

### or from csv file downloaded from Movebank
tracks<-move2kba(filename="example_data/MovebankExampleData.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND PREPARE SAMPLE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # tracks <- fread("example_data/Dataset_1004_2019-03-01.csv")     ## MUPE
tracks <- fread("example_data/Dataset_1012_2019-03-01.csv")     ## MABO St Helena
# tracks <- fread("example_data/Dataset_1151_2019-03-01.csv")     ## SHAG
# tracks <- fread("example_data/Dataset_1245_2019-03-01.csv")       ## RAZO
# tracks <- fread("example_data/R56Data.csv")       ## Luke Halpin dateline crossing data set


### CREATE COLONY DATA FRAME

Colony <- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  rename(Longitude=lon_colony,Latitude=lat_colony)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### formatFields - get data formatted for the other functions ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("formatFields.R")

### for a dataset with an already defined DateTime field
# tracks <- formatFields(tracks, field_id = "ID", field_lat="latitude", field_lon="longitude", field_datetime="date_time")
str(tracks)
### for a dataset with both date and time fields #####
tracks <- formatFields(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time="time_gmt")
## for a dataset with only a date column
# tracks <- formatFields(tracks, field_id = "track_id", field_lat="latitude", field_lon="longitude", field_date="date_gmt", field_time=NULL)

str(tracks)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSplit FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(tracks)
source("tripSplit.r")
Trips <- tripSplit(tracks, Colony=Colony, InnerBuff=2, ReturnBuff=20, Duration=1, plotit=T, nests = F, rmColLocs = T)
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

source("findScale.r")
##### Feeding the function the tripSummary output via Trips_summary argument speeds up computation #####
before <- Sys.time()
HVALS <- findScale(Trips, 
  ARSscale = T,
  Colony = Colony,
  Trips_summary = trip_distances)
HVALS
Sys.time() - before

# HVALS_RAZO <- HVALS
# HVALS_MUPE <- HVALS
# HVALS_MABO <- HVALS
# HVALS_EUSH <- HVALS

# print(c(HVALS_EUSH$ARSscale, HVALS_MABO$ARSscale, HVALS_MUPE$ARSscale, HVALS_RAZO$ARSscale))

head(tracks[order(tracks$DateTime), ])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN (MB edited) IndEffectTest function to assess whether individuals are site faithful 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("IndEffectTest_wip.R")  ## MBs edited version (27MAR19)

## second change: tell function which variable is the inGroupVar= (vs. GroupVar) (i.e. trip_id[inGroupVar] w/in indID[GroupVar]) 
indEffect <- IndEffectTest(Trips, GroupVar="ID", inGroupVar="trip_id", method="BA", Scale=HVALS$ARSscale, nboots=500)
indEffect$`Kolmogorov-Smirnov`


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN batchUD FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("estSpaceUse.r")
KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = HVALS$ARSscale, UDLev = 50, polyOut=F)

plot(KDE.Surface[[1]])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE bootstrap FUNCTION AND ASSIGN THRESHOLD FOR IBA IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("repAssess.r")

before <- Sys.time()
repr <- repAssess(Trips, Scale=HVALS$ARSscale, Iteration=50, BootTable = F, Res=10)
Sys.time() - before
repr


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE findKBA FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("findKBA.r")
KBAs <- findKBA(KDE.Surface, represent=repr$out, Col.size = 2000) ## error here if smoothr not installed!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Original IndEffectTest (and re-formatting necessary to make it work)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("IndEffectTest.R")

IETtrips <- Trips

IETtrips$indID <- IETtrips$ID
IETtrips$ID <- IETtrips$trip_id

indEffect <- IndEffectTest(IETtrips@data, Grouping_var="indID", method="BA", Scale=HVALS$ARSscale, nboots=500)



