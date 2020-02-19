###### track2KBA steps ######

library(dplyr)
library(track2KBA)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OUR track2iba FUNCTIONS - when they have recently changed after package installation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\track2iba")
setwd("C:/Users/Martim Bill/Documents/track2iba")
source("R/tripSplit.r")
source("R/tripSummary.r")
source("R/findScale.r")
source("R/repAssess.r")
source("R/estSpaceUse.r")
source("R/IndEffectTest.r")
source("R/findKBA.r")




## 1a.####
### move2KBA (Download and format Movebank data) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dataset <- move2KBA(MovebankID=621703893, User="bealhammar", Password="xxx")

### Movebank data (from move2KBA)
tracks <- dataset[["data"]]
colony <- dataset[["site"]]
# 
# head(tracks)
# head(colony)


## 1b. ####
### formatFields (upload data in own or STDB format, and re-format) ~~~~~~~~~~~~~~~~~~~

### Masked Booby
tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1004_2019-03-01.csv")   # Masked Booby
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1008_2019-03-01.csv")   # Masked Booby - NOT IN PACKAGE DATA
# tracks3 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1009_2019-03-01.csv")   # Masked Booby - NOT IN PACKAGE DATA
# tracks4 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1010_2019-03-01.csv")   # Masked Booby
# tracks5 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1011_2019-03-01.csv")   # Masked Booby
# tracks6 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1012_2019-03-01.csv")   # Masked Booby
# tracks7 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1013_2019-03-01.csv")   # Masked Booby
# tracks <- rbind.data.frame(tracks1, tracks2, tracks3, tracks4, tracks5, tracks6, tracks7) 
tracks <- rbind.data.frame(tracks1, tracks2, tracks3) # for README and vignette


# tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1151_2019-03-01.csv") # Black-legged kittiwake

# tracks1 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1245_2019-03-01.csv") # Razorbill
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1246_2019-03-01.csv") # Razorbill

### Shag
# tracks1 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1219_2019-03-01.csv") # Eur. Shag
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1218_2019-03-01.csv") # Eur. Shag
# 
# tracks <- rbind.data.frame(tracks1, tracks2) # combine two EUSH datasets

## MABO St Helena

data("boobies")

tracks <- boobies

colony <- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  rename(Longitude=lon_colony,Latitude=lat_colony)

tracks <- formatFields(tracks, field_ID = "track_id", field_Lat="latitude", field_Lon="longitude", field_Date="date_gmt", field_Time="time_gmt")


## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~

Trips <- tripSplit(tracks, Colony=colony, InnerBuff=5, ReturnBuff=10, Duration=1, plotit=T, Nests = F, rmColLocs = T, cleanDF=T)


## 2b. ####
### tripSummary (summary of trip movements, by individual) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trip_distances <- tripSummary(Trips, Colony = colony, Nests = F)
trip_distances


## 3. ####
### findScale (get average foraging range, a list of H-value options, and test whether desired grid cell for kernel estimation makes sense given movement scale/tracking resolution) ~~~~~~~~~~~~~~~

HVALS <- findScale(Trips,
  ARSscale = T,
  Trip_summary = trip_distances
  )
HVALS

# HVALS <- findScale(Trips,
#   ARSscale = T,
#   Trip_summary = NULL,
#   FPTscales = c(seq(1, 25, 1), seq(30, 50, 5), 75, 100),
#   plotPeaks = T,
#   Peak = "Flexible"
# )
# HVALS
## 4. ####
### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

# indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag, nboots=500)
# indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~

KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = HVALS$mag, UDLev = 50, polyOut=T, plotIt = F)
# KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = 0.5, Res=0.1, UDLev = 50, polyOut=F)

# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
# repr <- repAssess(Trips, KDE=KDE.Surface, Iteration=50, BootTable = F, avgMethod="weighted", Ncores = 11)
repr <- repAssess(Trips, KDE=KDE.Surface, Iteration=2, UDLev=50, avgMethod="mean", Ncores = 5)

Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = F) ## error here if smoothr not installed!
KBAs <- findKBA(KDE.Surface, Represent=repr$out, polyOut = T, plotIt=T) ## error here if smoothr not installed!
KBAs

KBA_sp <- KBAs
KBA_sf <- KBAs

### 8. ###
### plot KBA ### ~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(ggplot2)
#### if output sf (polygons) ####
KBA_sf <- KBAs

coordsets <- sf::st_bbox(KBA_sf)

KBAPLOT <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>%
  ggplot() +
  geom_sf(mapping = aes(fill=N_animals, colour=N_animals)) +
  borders("world", fill="dark grey", colour="grey20") +
  coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE) +
  theme(panel.background=element_blank(),
    panel.grid.major=element_line(colour="transparent"),
    panel.grid.minor=element_line(colour="transparent"),
    axis.text=element_text(size=16, colour="black"),
    axis.title=element_text(size=16),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour=FALSE) +
  scale_fill_continuous(name = "N animals") +
  ylab("Latitude") +
  xlab("Longitude")

KBAPLOT

KBAPLOT <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE)
KBAPLOT

plot(KBAs)
# plot the area meeting a certain percentage threshold (e.g. areas used by >75% of population)
plot(KBAs[KBAs$N_animals > 0, 3] )
plot(KBAs[KBAs$N_IND > 0, 1] )

# or, if there is a population estimate, the absolute number of individuals using the area
KBAs <- findKBA(KDE.Surface, Represent=repr$out, Col.size = 1000) ## error here if smoothr not installed!
plot(KBAs[KBAs$N_animals > 100, 1] )


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add in some background maps for context 

library(ggmap)

xmin <- min(Trips$Longitude)
xmax <- max(Trips$Longitude)
ymin <- min(Trips$Latitude) 
ymax <- max(Trips$Latitude) 

gmap <- ggmap::get_map(location=c(xmin, ymin, xmax, ymax), zoom=8, maptype = "satellite")

# ggmap(gmap)

# expBB <- c(0.5, 0.5, 0.5, 0.5)
expBB <- c(0, 0, 0, 0)

plot(st_transform(KBAs[KBAs$potentialKBA==T, ], crs = 3857)[3], bgMap = gmap, col=scales::alpha("red", 0.3),  border=scales::alpha("red", 0), key.length=1, expandBB=expBB)
raster::scalebar(10)
plot(st_transform(KBAs[KBAs$N_animals > 0.099, ], crs = 3857)[1], bgMap = gmap, key.length=1, border=scales::alpha("red", 0))

colony_sf <- st_transform(st_as_sf(colony, coords = c("Longitude", "Latitude"), 
  crs = "+proj=laea +lon_0=-6.442550477651 +lat_0=56.0611517263499 +ellps=WGS84"), 3857)


#######
potKBA <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>% 
  summarise(
    max_animals = max(na.omit(N_animals)),
    min_animals = min(na.omit(N_animals))
  )

potKBA
