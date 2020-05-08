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

tracks <- formatFields(tracks, fieldID = "track_id", fieldLat="latitude", fieldLon="longitude", fieldDate="date_gmt", fieldTime="time_gmt")

## 2a. ####
### tripSplit (split tracks in to discrete trips [and optionally filter]) ~~~~~~~~~~~~~

trips <- tripSplit(tracks, colony=colony, innerBuff=2, returnBuff=10, duration=1, nests = F, rmNonTrip = T)

tripmaps <- mapTrips(trips, colony=colony)

## 2b. ####
### tripSummary (summary of trip movements, by individual) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trips <- subset(trips, trips$Returns == "Yes" )

tripSum <- tripSummary(trips, colony = colony, nests = F)
tripSum

frange <- median(tripSum$max_dist)
frange
c(min(tripSum$max_dist), max(tripSum$max_dist))
fduration <- median(tripSum$duration)
fduration
c(min(tripSum$duration), max(tripSum$duration))

## 3. ####
### project trips to equal-area projection ### 

trips_prj <- projectTracks(trips)

## 4. ####
### findScale (get average foraging range, a list of H-value options, and test whether desired grid cell for kernel estimation makes sense given movement scale/tracking resolution) ~~~~~~~~~~~~~~~

HVALS <- findScale(trips_prj,
  scaleARS = F,
  sumTrips = tripSum
  )
HVALS

# HVALS <- findScale(trips,
#   scaleARS = T,
#   sumTrips = tripSum,
#   scalesFPT = seq(1, frange),
#   plotPeaks = T,
#   peakMethod = "first"
# )
# HVALS

## 4. ####
trips_prj <- trips_prj[trips_prj$ColDist > 2, ] # remove trip start and end points near colony

### IndEffectTest (test whether individuals are site-faithful across trips) ~~~~~~~~~~~

indEffect <- indEffectTest(trips_prj, groupVar="ID", tripID="trip_id", method="BA", scale=HVALS$mag, iterations=10)
indEffect$`Kolmogorov-Smirnov`


## 5. ####
### estSpaceUse (Produce utilization distributions for each individual) ~~~~~~~~~~~~~~~
h <- HVALS$mag
KDE.Surface <- estSpaceUse(tracks=trips_prj, scale = h, levelUD = 50, polyOut=T)
kde_map <- mapKDE(KDE.Surface$UDPolygons, show=F)
ud_map <- mapKDE(KDE.Surface$KDE.Surface, show=F)

n <- length(KDE.Surface$KDE.Surface)

# ggsave( paste0("C:/Users/Martim Bill/Documents/mIBA_package/figures/masked_boobys/indcores_", "h", round(h), "_", "n",n, ".png"), width = 8, height=6)


# plot(KDE.Surface$KDE.Surface[[4]]) # if polyOut=T
# plot(KDE.Surface[[1]])             # if polyOut=F


## 6. ####
### repAssess (Assess representativeness of tracked sample ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

before <- Sys.time()s
repr <- repAssess(trips_prj, KDE=KDE.Surface$KDE.Surface, iteration=1, levelUD=50, avgMethod="mean", nCores = 2)
Sys.time() - before


## 7. ####
### findKBA (Identify areas of significant aggregation) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
popsize <- 500 # pairs or individuals? Oppel et al. 2015

KBA_sf <- findKBA(KDE.Surface$KDE.Surface, represent=repr$out, popSize = popsize, polyOut = T) 
KBA_sp <- findKBA(KDE.Surface$KDE.Surface, represent=repr$out, polyOut = F) 

mapKBA(KBA_sf, colony = colony)
mapKBA(KBA_sp)

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
