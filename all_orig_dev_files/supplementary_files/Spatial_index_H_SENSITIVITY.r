############################################################################################################
####### SPATIAL AGGREGATION INDEX ANALYSIS FOR SEABIRDS   ##################################################
############################################################################################################
## SENSITIVITY OF BA OVERLAP ANALYSIS TO H VALUE in kernelUD
## raised by Ewan Wakefield
## script by Steffen Oppel
## 27 Feb 2018

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc)
library(data.table)
require(maps)
require(mapdata)
require(adehabitatHR)
require(foreign)
require(lme4)
require(geosphere)
require(sp)
library(rgdal)
require(rgeos)
library(raster)
library(trip)
library(vegan)
library(adehabitatLT)
require(foreach)
require(doParallel)
require(parallel)
library(lubridate)
library(tidyverse)
library(dplyr)
library(ggplot2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD CUSTOM FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index6.r")
source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index6.r")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v16.csv", header=T, sep=',')

overview<- overview %>%
  mutate(ReturnBuff=as.numeric(ReturnBuff)) %>%
  mutate(Duration=as.numeric(Duration)) %>%
  filter(n_individuals>4) %>%
  arrange(n_individuals) %>%
  dplyr::select(DataGroup, Family, scientific_name,site_name,colony_name,LATITUDE,LONGITUDE,breed_stage,n_individuals,n_tracks,device,InnerBuff,ReturnBuff,Duration)

head(overview)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD SAVED RESULTS FROM DATA GROUP SPECIFIC ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this analysis is based on Spatial_index_FINAL_ANALYSIS.r
## output compiled by SeabirdPrioritisation_data_aggregation_for_paper.R

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AggIndex_compiled_output_v16.RData")
load("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AggIndex_compiled_output_v16.RData")


### use breed stage as brood-guard or other
ORIG$chick<-ifelse(ORIG$breed_stage=="brood-guard",1,0)
TRIPS$chick<-ifelse(TRIPS$breed_stage=="brood-guard",1,0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA FROM OUTSIDE THE ATLANTIC OCEAN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

overview<-overview %>% filter(!(scientific_name %in% c("Fregata minor","Pterodroma ultima")))
ORIG<-ORIG %>% filter(DataGroup %in% overview$DataGroup)
TRIPS<-TRIPS %>% filter(DataGroup %in% overview$DataGroup)
						



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCALE H PARAMETER TO 0.5 - 50 km
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(TRIPS)
hist(log(TRIPS$max_dist))

HVALS<-TRIPS %>% group_by(DataGroup) %>%
  summarise(med_dist=median(max_dist), mag=log(max(max_dist))) %>%
  #mutate(mag=ifelse(mag<1,1,mag)) %>%
  mutate(H=med_dist/mag) %>%
  mutate(H1=ifelse(H>15,scales::rescale(H,to = c(15, 50)),H))%>%
  mutate(H2=scales::rescale(mag,to = c(0.5, 50)))
summary(HVALS)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START OF SPECIES-SPECIFIC BA CALCULATION FOR EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
#setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
ORIG<-fread("BA_sensitivity.csv")
dgs<-ORIG$DataGroup[is.na(ORIG$BAalt1)]
dgs<-dgs[!(dgs %in% c(77,226,181,172,247,76,242,167))]
dgs<-rev(dgs)
dgs

for (dg in dgs){        ### START OVERALL LOOP OVER EACH DATA GROUP (deliberately kept serial to facilitate troubleshooting)

#### SELECT DATA FOR ANALYSIS ##########################
tracksname<-sprintf("CLEAN_SeabirdTracks_DataGroup%s.csv",dg)
#setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")
tracks<-read.table(tracksname, header=T, sep=",")
tracks$DateTime <- ymd_hms(tracks$DateTime, tz="GMT")
try(tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID),silent=T)


#### REPLACE NON-NUMERIC IDs ##########################
names(tracks)[match("BirdId",names(tracks))]<-"bird_id"
if(is.numeric(tracks$bird_id)==F){tracks$bird_id<-as.numeric(as.factor(tracks$bird_id))}
if(is.numeric(tracks$ID)==F){tracks$ID<-as.numeric(as.factor(tracks$ID))}


#### CREATE COLONY LOCATION AND APPROPRIATE COORDINATE REFERENCE SYSTEM FOR PROJECTION ##########################
## FOR GPS DATA - take the first location for each Individual

if(overview$device[overview$DataGroup==dg]=="GPS"){
loc <- tracks %>%
	group_by(ID) %>%
	summarise(Latitude=mean(Latitude[1:3], na.rm=T), Longitude=mean(Longitude[1:3], na.rm=T))
loc<-as.data.frame(loc)

}else{

loc <- overview %>%
	filter(DataGroup==dg) %>%
	mutate(Latitude=LATITUDE, Longitude=LONGITUDE)%>%
	dplyr::select(Latitude,Longitude)
loc<-as.data.frame(loc)
}


### PROJECT COORDINATES FOR SPATIAL ANALYSES - this now needs WGS84 in CAPITAL LETTERS (not wgs84 anymore!!)
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tracks, match.ID=F)
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude[1], " +lat_0=", loc$Latitude[1], sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)
input <- DataGroup.Projected



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- input[input@data$Latitude>90,]
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	  Trip <- tripSplit(Track=Temp, Colony=loc,
					InnerBuff=overview$InnerBuff[overview$DataGroup==dg],
					ReturnBuff=as.numeric(overview$ReturnBuff[overview$DataGroup==dg]),
					Duration = overview$Duration[overview$DataGroup==dg],
					plotit=F, nests=ifelse(overview$device[overview$DataGroup==dg]=="GPS",T,F))
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE BA WITH ALTERNATIVE H
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls()[!ls() %in% c("ORIG", "batchUDOL", "Trips","dg","dgs","HVALS","KDE.Sp","overview","tripSplit")])
gc()
Output <- batchUDOL(DataGroup=Trips, Scale = HVALS$H1[HVALS$DataGroup==dg], UDLev = 50)
ORIG$BAalt1[ORIG$DataGroup==dg]<-Output$OverlapIndex
rm(list = ls()[!ls() %in% c("ORIG", "batchUDOL", "Trips","dg","dgs","HVALS","KDE.Sp","overview","tripSplit")])
gc()
Output <- batchUDOL(DataGroup=Trips, Scale = HVALS$H2[HVALS$DataGroup==dg], UDLev = 50)
ORIG$BAalt2[ORIG$DataGroup==dg]<-Output$OverlapIndex

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
#setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
fwrite(ORIG,"BA_sensitivity.csv")

rm(list = ls()[!ls() %in% c("ORIG", "batchUDOL", "dg","dgs","HVALS","KDE.Sp","overview","tripSplit")])
gc()
} ### end of loop across all data groups

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF SPECIES-SPECIFIC INDEX CALCULATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
#setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
ORIG<-fread("BA_sensitivity.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN AVIAN BODY MASS DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

mass<-fread("AvianBodyMass.csv")
overview<-merge(overview, mass, by="scientific_name", all.x=T)
head(overview)

ORIG<-ORIG %>% mutate(mass=as.numeric(overview$Mean_mass[match(DataGroup,overview$DataGroup)]))
TRIPS<-TRIPS %>% mutate(mass=as.numeric(overview$Mean_mass[match(DataGroup,overview$DataGroup)]))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN SamplingRate DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

SR<-fread("Seabird_DataGroup_SamplingRates.csv")
overview<-merge(overview, SR[,c(1,10)], by="DataGroup", all.x=T)
head(overview)

ORIG<-ORIG %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))
TRIPS<-TRIPS %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN COLONY SIZE DATA AND MODIFY TO COMMON CURRENCY (pairs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

colonies<-fread("Seabird_priority_colony_sizes.csv")
colonies<-colonies %>% select(DataGroup,scientific_name,PopYearStart,PopYearEnd,PopMin,PopMax,PopUnits,MostRecent)
head(colonies)
unique(colonies$PopUnits)

colonies<-colonies %>%
  dplyr::arrange(desc(MostRecent)) %>%
  mutate(COL_SIZE=ifelse(PopUnits %in% c('breeding pairs','chicks'), PopMin, PopMin*0.5)) %>%
  group_by(DataGroup,scientific_name) %>%
  summarise(COL_SIZE=dplyr::first(COL_SIZE))

overview<-merge(overview, colonies, by=c("DataGroup","scientific_name"), all.x=T)
head(overview)


### MISSING COLONY SIZE DATA ###

misscol<-overview %>% filter(is.na(COL_SIZE))
misscol<-misscol %>% group_by(DataGroup,Family,scientific_name,site_name,colony_name,LATITUDE,LONGITUDE) %>%
  summarise(N_tracks=sum(n_tracks))
#fwrite(misscol,"ColonySizes_NEEDED.csv")




ORIG<-ORIG %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))
TRIPS<-TRIPS %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL ANALYSIS OF FAMILY EFFECT FOR BA OVERLAP INDEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(ORIG %>% arrange(BA))
head(ORIG %>% arrange(desc(BA)))

### FOR ORIGINAL DATA ###
m0BA<-lmer(BA~chick+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
m1BA<-lmer(BA~Family+chick+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
anova(m0BA, m1BA)
summary(m1BA)

### FOR ALTERNATIVE 1 ###
m0BA<-lmer(BAalt1~chick+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
m1BA<-lmer(BAalt1~Family+chick+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
anova(m0BA, m1BA)
summary(m1BA)

### FOR ALTERNATIVE 2 ###
m0BA<-lmer(BAalt2~chick+mass+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
m1BA<-lmer(BAalt2~Family+chick+mass+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG)
mBAaov<-anova(m0BA, m1BA)
summary(m1BA)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE FIGURE 2 BA BOXPLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(ORIG$BA)

ORIG %>% arrange(BA) %>% select(scientific_name,breed_stage,BA,Scale)

### CREATE SORT ORDER OF FAMILIES
sortfam<-ORIG %>%
  filter(breed_stage %in% c('brood-guard')) %>%
  group_by(Family) %>%
  summarise(range=median(BA))%>%
  arrange(desc(range)) %>%
  select(Family)


#pdf("Fig2.pdf", height=8, width=8)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Fig1.jpg", height=780, width=680, quality=100)

ORIG %>%
  filter(breed_stage %in% c('incubation','brood-guard')) %>%
  mutate(breed_stage=as.factor(breed_stage)) %>%
  mutate(breed_stage=factor(breed_stage, levels=levels(breed_stage)[c(2,1)])) %>%
  mutate(Family=as.factor(Family)) %>%
  mutate(Family=fct_relevel(Family,sortfam$Family)) %>%
  
  ggplot(aes(x=Family, y=BA, width=1), size=1)+geom_boxplot()+
  facet_wrap("breed_stage", ncol=1, scales = "fixed")+
  xlab("Seabird family") +
  ylab("Bhattacharyya's affinity index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=90, vjust=0.5, hjust=0.99), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()



