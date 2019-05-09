###### track2KBA steps ######

library(dplyr)
library(track2KBA)

tracks <- data.table::fread("C:/Users/Martim Bill/Documents/mIBA_package/test_data/green_turtles_cmnbd.csv")

## 1
### formatFields
## this works!

## Green turtles, Bijagos Archipelago

tracks <- formatFields(tracks, field_ID = "Tag_ID", field_Lat="Latitude", field_Lon="Longitude", field_Date="UTC_Date", field_Time="UTC_Time", format_DT = "dmy_HMS")


## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~
## Inappropriate for sea turtles, since they aren't moving to and from a central place
## 2b. ####
### tripSummary 
## Inappropriate for sea turtles


## 3. ####
### findScale 
# This step utterly fails for non-central place moving animals (i.e. those for which TripSplit is not appropriate to run and which do not have central "Colony" place)

HVALS <- findScale(tracks,
  Colony = Colony
)
HVALS


## 4. ####
### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

# indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag, nboots=500)
# indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~
## works!

KDE.Surface <- estSpaceUse(DataGroup=tracks, Scale = 5, UDLev = 50, polyOut=T)

# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
repr <- repAssess(tracks, listKDE=KDE.Surface$KDE.Surface, Iteration=5, BootTable = F)
Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = T, plotit = T) ## error here if smoothr not installed!
KBAs
# plot the area meeting a certain percentage threshold (e.g. areas used by >75% of population)
plot(KBAs[KBAs$N_animals > 60, 1] )

# or, if there is a population estimate, the absolute number of individuals using the area
KBAs <- findKBA(KDE.Surface, Represent=repr$out, Col.size = 1000) ## error here if smoothr not installed!
plot(KBAs[KBAs$N_animals > 100, 1] )

