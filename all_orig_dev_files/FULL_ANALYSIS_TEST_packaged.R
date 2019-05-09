###### track2KBA steps ######

library(dplyr)
library(track2KBA)

## 1a.####
### move2KBA (Download and format Movebank data) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dataset <- move2KBA(MovebankID=621703893, User="bealhammar", Password="xxx")

### Movebank data (from move2KBA)
tracks <- dataset[["data"]]
Colony <- dataset[["site"]]
# 
# head(tracks)
# head(Colony)


## 1b. ####
### formatFields (upload data in own or STDB format, and re-format) ~~~~~~~~~~~~~~~~~~~

tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1012_2019-03-01.csv")
# tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1151_2019-03-01.csv")
# tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1245_2019-03-01.csv")

## MABO St Helena

Colony <- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  rename(Longitude=lon_colony,Latitude=lat_colony)

tracks <- formatFields(tracks, field_ID = "track_id", field_Lat="latitude", field_Lon="longitude", field_Date="date_gmt", field_Time="time_gmt")


## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~

Trips <- tripSplit(tracks, Colony=Colony, InnerBuff=2, ReturnBuff=20, Duration=1, plotit=T, Nests = F, rmColLocs = T)


## 2b. ####
### tripSummary (summary of trip movements, by individual) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trip_distances <- tripSummary(Trips, Colony = Colony, Nests = F)
trip_distances


## 3. ####
### findScale (get average foraging range, a list of H-value options, and test whether desired grid cell for kernel estimation makes sense given movement scale/tracking resolution) ~~~~~~~~~~~~~~~

HVALS <- findScale(Trips,
  ARSscale = T,
  Colony = Colony,
  Trips_summary = trip_distances
  )
HVALS


## 4. ####
### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

# indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag, nboots=500)
# indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~

KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = HVALS$half_mag, UDLev = 50, polyOut=T)

# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
repr <- repAssess(Trips, listKDE=KDE.Surface$KDE.Surface, Iteration=20, BootTable = F)
Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = F) ## error here if smoothr not installed!
KBAs
# plot the area meeting a certain percentage threshold (e.g. areas used by >75% of population)
plot(KBAs[KBAs$N_animals > 60, 1] )

# or, if there is a population estimate, the absolute number of individuals using the area
KBAs <- findKBA(KDE.Surface, Represent=repr$out, Col.size = 1000) ## error here if smoothr not installed!
plot(KBAs[KBAs$N_animals > 100, 1] )

