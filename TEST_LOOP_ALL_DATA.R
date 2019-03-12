############################################################################################################
####### IDENTIFYING MARINE KBAs FROM ANIMAL TRACKING DATA #################################################
############################################################################################################

### test loop across all datasets for package development
### started by steffen oppel on 11 March 2019


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(data.table)
library(lubridate)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OUR track2iba FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\track2iba")
#setwd("C:/Users/Martim Bill/Documents/track2iba")
source("tripSplit.r")
source("tripSummary.r")
source("findScale.r")
source("repAssess.r")
source("estSpaceUse.r")
source("findKBA.r")

## troubleshoot crashes by assessing objects
obj_summary<-data.frame()
for (obj in ls()) {obj_summary<-data.frame(obj=obj,size=as.numeric(object.size(get(obj)))) %>% bind_rows(obj_summary) }
obj_summary %>% arrange(desc(size))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND PREPARE SAMPLE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\track2iba\\example_data")
alldata<-list.files(pattern="Dataset")

OUTPUT<-data.frame()

for (f in 1:length(alldata)) {
  tracks <- fread(alldata[f])
  

### CREATE COLONY DATA FRAME
#if('lon_colony' %in% names(tracks){
Colony<- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  rename(Longitude=lon_colony,Latitude=lat_colony)
#}
  
### Convert Dates and Times
tracks <- tracks %>%
  mutate(DateTime = ymd_hms(paste(date_gmt,time_gmt, sep = " "))) %>%
  dplyr::select(track_id, latitude, longitude,DateTime) %>%
  rename(ID=track_id,Latitude=latitude,Longitude=longitude)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSplit FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:\\STEFFEN\\track2iba\\tripSplit.r")
Trips <- tripSplit(tracks, Colony=Colony, InnerBuff=5, ReturnBuff=25, Duration=5, plotit=T, nests = F, rmColLocs = T)
tracks<-NULL
tripSplit<-NULL
splitSingleID<-NULL
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSummary FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:\\STEFFEN\\track2iba\\tripSummary.r")
trip_distances <- tripSummary(Trips, Colony = Colony, nests = F)
OUTPUT<-trip_distances %>% mutate(Dataset=alldata[f]) %>% bind_rows(OUTPUT)
fwrite(OUTPUT,"test_run_trip_sums.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN findScale FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:\\STEFFEN\\track2iba\\findScale.r")
HVALS <- findScale(Trips, ARSscale = T, Colony = Colony,Trips_summary=trip_distances)
trip_distances<-NULL
tripSummary<-NULL
findScale<-NULL
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE repAssess FUNCTION BEFORE estSpaceUse, because the loop will then be lighter (without the massive KDE.Surface object)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:\\STEFFEN\\track2iba\\repAssess.r")
represent <- repAssess(Trips, Scale=HVALS$ARSscale, Iteration=10, BootTable = F)
repAssess<-NULL
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN estSpaceUse FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:\\STEFFEN\\track2iba\\estSpaceUse.r")
KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = HVALS$ARSscale, polyOut=F)
Trips<-NULL
estSpaceUse<-NULL
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE findIBA FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IBAs <- findKBA(KDE.Surface, represent=represent$out)
rm(list=setdiff(ls(), c("OUTPUT","alldata","findKBA")))
gc()

} ### end loop over all data


