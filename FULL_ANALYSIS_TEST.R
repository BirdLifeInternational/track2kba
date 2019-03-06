############################################################################################################
####### IDENTIFYING MARINE KBAs FROM SEABIRD TRACKING DATA #################################################
############################################################################################################

### full analysis for package development
### started by steffen oppel on 1 March 2019


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(maptools)
require(geosphere)
require(sp)
library(rgdal)
library(tidyverse)
library(data.table)
library(maps)
library(rgeos)
library(adehabitatHR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OUR track2iba FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 setwd("C:\\STEFFEN\\track2iba")
#setwd("C:/Users/Martim Bill/Documents/track2iba")
source("tripSplit.r")
source("tripSummary.r")
source("scaleARS.r")
source("bootstrap.r")
source("batchUD.r")
source("varianceTest.r")
source("findIBA.r")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND PREPARE SAMPLE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tracks <- fread("example_data/Dataset_1004_2019-03-01.csv")
# tracks <- fread("example_data/Dataset_1012_2019-03-01.csv")
# tracks <- fread("example_data/Dataset_1151_2019-03-01.csv")
# tracks <- fread("example_data/Dataset_1245_2019-03-01.csv")


### Convert Dates and Times

tracks <- tracks %>%
  mutate(DateTime = ymd_hms(paste(date_gmt,time_gmt, sep = " "))) %>%
  mutate(TrackTime = as.double(DateTime)) %>%
  dplyr::select(track_id, latitude, longitude,DateTime, TrackTime) %>%
  rename(ID=track_id,Latitude=latitude,Longitude=longitude)
head(tracks)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE PROJECTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
proj.UTM <- CRS(paste("+proj=laea +lon_0=", mean(tracks$Longitude), " +lat_0=", mean(tracks$Latitude), sep=""))
DataGroup <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")), data = tracks, match.ID=F)
DataGroup.Projected <- spTransform(DataGroup, CRS=proj.UTM)
plot(DataGroup)
plot(DataGroup[1,],pch=16, col='red', add=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSplit FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(tracks)
source("tripSplit.r")
Trips<-tripSplit(tracks, Colony=tracks[1,3:2], InnerBuff=15, ReturnBuff=50, Duration=5, plotit=T, nests = F)
dim(Trips)


# Trips <- Trips[Trips$trip_id != "-1",]
# Trips <- Trips[!Trips$trip_id %in% names(which(table(Trips$trip_id) < 5)), ] # remove trips with less than 5 points
# Trips$ID <- Trips$trip_id #reset the ID field to individual trips rather than individual birds (optional!)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN tripSummary FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("tripSummary.r")
trip_distances <- tripSummary(Trips, Colony = tracks[1,3:2], nests = F)
trip_distances
dim(trip_distances)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN scaleARS FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ScaleOut <- scaleARS(Trips, Scales = c(seq(0, 250, 0.5)), Peak="Flexible")
# FAME_summary$ARS[FAME_summary$DG==DG] <- ScaleOut


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN batchUD FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UD <- 50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUD(Trips[Trips$trip_id != "-1",], Scale = 50, UDLev = UD)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE polyCount FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gridres <- ifelse(FAME_summary$species[DG]=="EUSH",0.007,0.0286)	### Resolution should be 0.5 km for EUSH and 2km for others, converted to degree based on 70 km width of 1 degree grid cell	
IBA<-polyCount(Output, Res=gridres)		
spoldf <- rasterToPolygons(IBA, n=4)
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\FAME\\OUTPUT_EWAN")
outname<-paste("IBA_PolyCount",FAME_summary$species[DG], FAME_summary$site[DG],sep="_")
unlink(outname, recursive = T, force = T)								### deletes existing folder of the same name to avoid function crashing on re-run
writeOGR(spoldf, dsn=outname, layer=outname,  driver="ESRI Shapefile")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE bootstrap FUNCTION AND ASSIGN THRESHOLD FOR IBA IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(length(unique(Trips@data$trip_id))>2){						## FAILS WITH <2! less than 15 trips does not qualify for IBA

test <- bootstrap(Trips, Scale=ScaleOut, Iteration=1)
dev.off()
sumname<-paste("BootstrapOutput",FAME_summary$species[DG], FAME_summary$site[DG],"csv",sep=".")
write.table(test, sumname, row.names=F, sep=",")
representativity<-test$RepresentativeValue[1]
FAME_summary$representativ[DG]<-representativity
thresh<-ifelse(representativity>0.9,10,ifelse(representativity>0.8,12.5,ifelse(representativity>0.7,20,100)))	## set threshold depending on representativity value




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE thresholdRaster FUNCTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(representativity>0.69){
final<-thresholdRaster(IBA, Threshold=thresh)
dev.off()
final@data$Species<-FAME_summary$species[DG]
final@data$Site<-FAME_summary$site[DG]
final@data$Col_size<-FAME_summary$Col_size[DG]
mult<-ifelse(representativity>0.9,0.9,ifelse(representativity>0.8,0.75,ifelse(representativity>0.7,0.5,0)))	### set threshold depending on representativity value
final@data$IBA<-final@data$Col_size*(final@data$MaxPerc/100)*mult		### THIS IS THE NUMBER OF BIRDS USING THE POLYGON



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT OUTPUT TO A SHAPEFILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(FAME_summary, "FAME_IBA_analysis_output.csv", row.names=F, sep=",")
outname<-paste("IBA_Polygon",FAME_summary$species[DG], FAME_summary$site[DG],sep="_")
unlink(outname, recursive = T, force = T)								### deletes existing folder of the same name to avoid function crashing on re-run
writeOGR(final, dsn=outname, layer=outname,  driver="ESRI Shapefile")
}		## close if loop for low representativity
}		## close if loop for no IBA if n trips too small
rm(summary, input, tracks, DataGroup.Wgs, DataGroup.Projected, final, test, thresh, Trips, DataGroup, output, outname, mult)
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN THE variance Test to assess whether all trips of an individual should be included
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output <- batchUD(DataGroup, Scale = ScaleOut/2, UDLev = UD)


## variance test
bird_string<-as.character(Output@data$ID)
vt<-varianceTest(Output, bird_string, Iteration=10)
vt
### to choose randomly just one trip per individual
if (vt < 0.25)
{
bird_idtrip=datagroupsUDd@data
birds=unique(bird_idtrip$originalID)
trips=numeric()
set.seed(1)
for (x in 1:length(birds)) trips=c(trips, as.character(sample(bird_idtrip[bird_idtrip$originalID==(birds[x]),]$ID,1)))
DataGroupTrips2=DataGroupTrips[DataGroupTrips$ID%in%trips,]
datagroupsUDd2=batchUD(DataGroupTrips2, Scale=fpt.scales/2, UDLev=50)
datagroupsUDd=datagroupsUDd2
DataGroupTrips=DataGroupTrips2
}












