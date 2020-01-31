###### track2KBA steps ######

library(dplyr)
library(track2KBA)
library(stringr)

tracksALL <- data.table::fread("C:/Users/Martim Bill/Documents/mIBA_package/test_data/green_turtles_cmnbd.csv")
tracks <- filter(tracksALL, tracksALL$Tag_ID !="6524:60868")
# tracks <- tracksALL[tracksALL$Tag_ID %in% c("6524:60898", "6524:60892", "6524:60891", "6524:60889", "6524:60887", "6524:60886", "6524:60865"), ]

## 1
### formatFields
## this works!

## Green turtles, Bijagos Archipelago

tracks <- formatFields(tracks, field_ID = "Tag_ID", field_Lat="Latitude", field_Lon="Longitude", field_Date="UTC_Date", field_Time="UTC_Time", format_DT = "dmy_HMS")

tracks$ID <- str_replace(tracks$ID, pattern = ":", replacement = ".")

## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~

nests <- tracks %>% group_by(ID) %>% 
  summarise(
    Longitude = first(Longitude),
    Latitude  = first(Latitude)
  )

Trips <- tripSplit(tracks, Colony=nests, InnerBuff = 10, ReturnBuff = 1, Nests = T, rmColLocs = F)

## 2b. ####
### tripSummary 
tripSum <- tripSummary(Trips, nests, Nests = T)
tripSum


## 3. ####
### findScale 
# This step utterly fails for non-central place moving animals (i.e. those for which TripSplit is not appropriate to run and which do not have central "Colony" place)

HVALS <- findScale(Trips,
  Colony = nests,
  Trips_summary = tripSum,
  ARSscale = F
)
HVALS


## 4. ####
### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

# indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag, nboots=500)
# indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~
## works!

KDE.Surface <- estSpaceUse(DataGroup=tracks, Scale = 3, UDLev = 50, polyOut=F)

# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
repr <- repAssess(tracks, KDE=KDE.Surface, Iteration=50, BootTable = F)
Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = T, plotit = T) ## error here if smoothr not installed!
KBAs

# KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = F) ## error here if smoothr not installed!


#### Plot density map/KBA plot
library(ggmap)

xmin <- min(tracks$Longitude)
xmax <- max(tracks$Longitude)
ymin <- min(tracks$Latitude) 
ymax <- max(tracks$Latitude) 

gmap <- ggmap::get_map(location=c(xmin, ymin, xmax, ymax), zoom=8, maptype = "satellite")

# ggmap(gmap)

# expBB <- c(0.5, 0.5, 0.5, 0.5)
expBB <- c(0, 0, 0, 0)

# plot(st_transform(KBAs[KBAs$potentialKBA==T, ], crs = 3857)[3], bgMap = gmap, col=scales::alpha("red", 0.3),  border=scales::alpha("red", 0), key.length=1, expandBB=expBB)
# raster::scalebar(10)
# plot(st_transform(KBAs[KBAs$N_animals > 0.099, ], crs = 3857)[1], bgMap = gmap, key.length=1, border=scales::alpha("red", 0))
plot(st_transform(KBAs[KBAs$N_animals > 0, ], crs = 3857)[1], bgMap = gmap, key.length=1, border=scales::alpha("red", 0), expandBB=expBB)
plot(st_transform(KBAs[KBAs$N_animals > 0, ], crs = 3857)[2], bgMap = gmap, key.length=1, border=scales::alpha("red", 0), expandBB=expBB)
raster::scalebar(40)
