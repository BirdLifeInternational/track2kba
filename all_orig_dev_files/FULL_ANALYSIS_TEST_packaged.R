###### track2KBA steps ######

library(dplyr)
library(track2KBA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OUR track2iba FUNCTIONS - when they have recently changed after package installation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setwd("C:\\STEFFEN\\track2iba")
# setwd("C:/Users/Martim Bill/Documents/track2iba")


## 1a.####
### move2KBA (Download and format Movebank data) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dataset <- move2KBA(MovebankID=621703893, User="bealhammar", Password="xxx")

### Movebank data (from move2KBA)
# tracks <- dataset[["data"]]
# colony <- dataset[["site"]]
# 
# head(tracks)
# head(colony)


## 1b. ####
### formatFields (upload data in own or STDB format, and re-format) ~~~~~~~~~~~~~~~~~~~

### Masked Booby
# tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1004_2019-03-01.csv")   # Masked Booby
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1008_2019-03-01.csv")   # Masked Booby - NOT IN PACKAGE DATA
# tracks3 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1009_2019-03-01.csv")   # Masked Booby - NOT IN PACKAGE DATA
# tracks4 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1010_2019-03-01.csv")   # Masked Booby
# tracks5 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1011_2019-03-01.csv")   # Masked Booby
# tracks6 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1012_2019-03-01.csv")   # Masked Booby
# tracks7 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1013_2019-03-01.csv")   # Masked Booby
# tracks <- rbind.data.frame(tracks1, tracks2, tracks3, tracks4, tracks5, tracks6, tracks7) 
# tracks <- rbind.data.frame(tracks1, tracks2, tracks3) # for README and vignette


# tracks <- data.table::fread("all_orig_dev_files/example_data/Dataset_1151_2019-03-01.csv") # Black-legged kittiwake

# tracks1 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1245_2019-03-01.csv") # Razorbill
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1246_2019-03-01.csv") # Razorbill

### Shag
# tracks1 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1219_2019-03-01.csv") # Eur. Shag
# tracks2 <- data.table::fread("all_orig_dev_files/example_data/Dataset_1218_2019-03-01.csv") # Eur. Shag
# 
# tracks <- rbind.data.frame(tracks1, tracks2) # combine two EUSH datasets

## MABO St Helena (using package example data)

data("boobies")

tracks <- boobies

colony <- tracks[1,] %>% dplyr::select(lon_colony,lat_colony) %>%
  dplyr::rename(Longitude=lon_colony,Latitude=lat_colony)

tracks <- formatFields(tracks, field_ID = "track_id", field_Lat="latitude", field_Lon="longitude", field_Date="date_gmt", field_Time="time_gmt")


## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~

Trips <- tripSplit(tracks, Colony=colony, InnerBuff=2, ReturnBuff=10, Duration=1, plot=T, Nests = F, rmNonTrip = T)

## 2b. ####
### tripSummary (summary of trip movements, by individual) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- subset(Trips, Trips$Returns == "Yes" )

TripSum <- tripSummary(Trips, Colony = colony, Nests = F)
TripSum

frange <- median(TripSum$max_dist)
frange
c(min(TripSum$max_dist), max(TripSum$max_dist))
fduration <- median(TripSum$duration)
fduration
c(min(TripSum$duration), max(TripSum$duration))

## 3. ####
### findScale (get average foraging range, a list of H-value options, and test whether desired grid cell for kernel estimation makes sense given movement scale/tracking resolution) ~~~~~~~~~~~~~~~

HVALS <- findScale(Trips,
  ARSscale = F,
  Trip_summary = TripSum
  )
HVALS

# HVALS <- findScale(Trips,
#   ARSscale = T,
#   Trip_summary = TripSum,
#   FPTscales = seq(1, frange),
#   plotPeaks = T,
#   findPeak = "Flexible"
# )
# HVALS

## 4. ####
Trips <- Trips[Trips$ColDist > 2, ] # remove trip start and end points near colony

### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag, nboots=100)
indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~
h <- HVALS$mag
KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = h, UDLev = 50, polyOut=T, plot = T)
# KDE.Surface <- estSpaceUse(DataGroup=Trips, Scale = 0.5, Res=0.1, UDLev = 50, polyOut=F)
n <- length(KDE.Surface$KDE.Surface)

ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/masked_boobys/indcores_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6)


# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()
repr <- repAssess(Trips, KDE=KDE.Surface$KDE.Surface, Iteration=10, UDLev=50, avgMethod="mean", Ncores = 2)
Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
popsize <- 500 # pairs or individuals? Oppel et al. 2015

KBAs <- findKBA(KDE.Surface, Represent=repr$out, popSize = popsize, polyOut = T, plot=T) 
KBAs

KBA_sp <- findKBA(KDE.Surface, Represent=repr$out, polyOut = F, plot=T) 

KBA_sf <- KBAs


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add in some background maps for context 

# aggs <- KBAs[KBAs$N_animals > 0, ]
aggs <- KBAs[KBAs$N_animals >= (0.1 * popsize), ]

aggs_3857 <- st_transform(aggs, crs = 3857)
# aggs_3857 <- st_transform(aggs, crs = 3857)

### 8. ###
### plot KBA ### ~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(ggplot2)
#### if output sf (polygons) ####

coordsets <- sf::st_bbox(aggs)

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
island <- raster::shapefile("C:/Users/Martim Bill/Documents/geodata/saint_helena_border/polbnda_shn.shp")
island2 <- st_as_sf(aggregate(island, dissolve=T))

colony_sf <- st_as_sf(colony, coords = c("Longitude", "Latitude"), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

KBAPLOT <- ggplot() +
  geom_sf(data=aggs, mapping = aes(fill=N_animals, colour=N_animals)) +
  geom_sf(data=island2, fill="dark grey", colour="grey20") +
  geom_sf(data=colony_sf, color="red", size=4) +
  # borders(rnaturalearth::ne_countries("countries"), fill="dark grey", colour="grey20") +
  coord_sf(xlim = c(coordsets$xmin-.1, coordsets$xmax+.1), ylim = c(coordsets$ymin-.1, coordsets$ymax+.1), expand = FALSE) +
  theme(panel.background=element_blank(),
    panel.grid.major=element_line(colour="transparent"),
    panel.grid.minor=element_line(colour="transparent"),
    axis.text=element_text(size=11, colour="black"),
    axis.title=element_text(size=16),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  guides(colour=FALSE) +
  scale_fill_continuous(name = "N animals") +
  ylab("Latitude") +
  xlab("Longitude")
KBAPLOT

ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/masked_boobys/findKBAplot_", "h", round(h), "_", "n",n, ".png"), width = 10, height=8)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add in some background maps for context 

library(ggmap)
xmin <- st_bbox(aggs)[[1]] - 0.1
xmax <- st_bbox(aggs)[[3]] + 0.1
ymin <- st_bbox(aggs)[[2]] - 0.1
ymax <- st_bbox(aggs)[[4]] + 0.1

gmap <- ggmap::get_map(location=c(xmin, ymin, xmax, ymax), zoom=10, maptype = "satellite")

colony_sf <- st_transform(colony_sf, 3857)
 
# dev.off()


# Use the function:
source("C:/Users/Martim Bill/Documents/R/source_scripts/ggmap_bbox.r")

gmap <- ggmap_bbox(gmap)

ggmap(gmap) + 
  coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  geom_sf(data = aggs_3857, inherit.aes = FALSE, aes(fill=N_animals), color=NA) +
  scale_fill_gradientn(colours=alpha(sf.colors(n=3), 0.9)) +
  geom_sf(data=colony_sf, inherit.aes = FALSE, color="red", size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"),
    axis.text=element_text(size=9, color="black"),
    axis.title=element_text(size=14),
    legend.direction = "horizontal",
    legend.position=c(0.1, -0.08),
    legend.title=element_text(size=14),
    legend.text = element_text(size = 12) ) +
  ylab("Latitude") +
  xlab("Longitude")

ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/masked_boobys/potKBA_", "h", round(h), "_", "n",n, ".png"), width = 10, height=8 )



#######
potKBA <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>% 
  summarise(
    max_animals = max(na.omit(N_animals)),
    min_animals = min(na.omit(N_animals))
  )

potKBA
