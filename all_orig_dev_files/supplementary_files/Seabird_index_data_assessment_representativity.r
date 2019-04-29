############################################################################################################
####### SEABIRD PRIORITISATION BASED ON SPATIAL AGGREGATION INDEX ANALYSIS #################################
############################################################################################################
## this analysis is based on Spatial_index_COMPLETE_REANALYSIS.r
## output compiled by Seabird_index_OUTPUT_COMPILATION.r
## developed by steffen.oppel@rspb.org.uk in December 2016

## v4 modified 23 August 2017
## based on complete new reanalysis of cleaned data from 23 Aug 2017

## v5modified on 7 Sept 2017 after inclusion of storm petrel data
## removed the less interesting figures, exploratory plots comparing lots of indices (see v1-4 for code)
## removed latitudinal variation plot (most families have very little latitudinal variation)

## MEETING ON 8 SEPT 2017: Maria requested to eliminate all data with <10 individuals, and check representativity for all data sets


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
require(foreign)
library(readr)
library(data.table)
 require(foreach)
  require(doParallel)
require(parallel)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD SAVED RESULTS FROM DATA GROUP SPECIFIC ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\AggIndex_compiled_output_v11.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSESS REPRESENTATIVITY OF ALL SIMULATED DATASETS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(SIMUL)





cl<-makeCluster(detectCores())
registerDoParallel(cl)

Representativity <- foreach(dg=unique(SIMUL$DataGroup), .combine = rbind, .packages=c("dplyr","tidyverse")) %dopar% {
rm(RepresentativeValue)
test<-SIMUL %>%
	filter(DataGroup==dg) %>%
	select(DataGroup,SampleSize,InclusionMean,n_tracks)

try(M1 <- nls((test$InclusionMean ~ (a*test$SampleSize)/(1+b*test$SampleSize)), data=test, start=list(a=1,b=0.1)), silent = TRUE)
 if ('M1' %in% ls()){       ### run this only if nls was successful
  Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
  RepresentativeValue <- test %>%
	mutate(pred=predict(M1)) %>%
	group_by(DataGroup,SampleSize) %>%
	summarise(out=max(pred)/Asymptote*100) %>%
	filter(out==max(out)) %>%
	select(DataGroup,out) %>%
	mutate(type='asymptote')%>%
	mutate(asym=Asymptote)
}else{
  RepresentativeValue <- test %>%
	filter(SampleSize==max(SampleSize)) %>%
	group_by(DataGroup,SampleSize) %>%
	summarise(out=mean(InclusionMean)) %>%
	select(DataGroup,out)%>%
	mutate(type='inclusion')%>%
	mutate(asym=out)
}
rm(M1,test)
return(RepresentativeValue)

}			## close loop over all data groups

Representativity


ORIG<-merge(ORIG,Representativity, by="DataGroup")
fwrite(ORIG,"SeabirdPrioritisation_SpatIndex_OrigData_all.csv")


dim(ORIG[ORIG$out>70,])





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT ORIGINAL DATA SET SO THAT GRADIENT IS VISIBLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(ORIG)
library(RColorBrewer)

colourpalette<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#000120')


pdf("Seabird_prioritisation_BA_IBA_scatter_representative.pdf", height=7, width=9)

ORIG %>%
  mutate(BreedingStage=ifelse(breed_stage=="brood-guard","brood-guard","other")) %>%
  filter(Family!="") %>%
  filter(type=="asymptote") %>%
  filter(out>70) %>%
  
  ggplot(aes(x=log(IBA20), y=BA,colour=Family, pch=BreedingStage))+
  geom_point(size=2.5) +
  #scale_colour_brewer(type = "qual" , palette = "Set1") +
  scale_colour_manual(values = colourpalette)+
  xlab("log(size) of marine IBA") +
  ylab("Bhattacharrya's Affinity index") +
  geom_segment(aes(x=log(MPAsizes[5]),y=0.5,xend=log(MPAsizes[5]), yend=1),colour="red", size=0.5, linetype=2) +
  geom_segment(aes(x=0,y=0.5,xend=log(MPAsizes[5]), yend=0.5),colour="red", size=0.5, linetype=2) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT CLUMPEDNESS RATIO ON ORIGINAL DATA SET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(ORIG)


pdf("IBA_MMA_ratio.pdf", height=7, width=9)

ORIG %>% 
  mutate(Clumpedness=log(IBA20)/log(MMA)) %>%
  #filter(Clumpedness<(-20)) %>%
  filter(Family!="") %>%
  filter(type=="asymptote") %>%
  filter(out>70) %>%
  
  
  ggplot(aes(x=Family, y=Clumpedness, width=1), size=1)+geom_boxplot()+
  xlab("Seabird family") +
  ylab("IBA/MMA ratio") +
  geom_hline(aes(yintercept=1),colour="red", size=1) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()









