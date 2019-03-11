## findKBA  #####################################################################################################################

## written by Steffen Oppel, Martin Beal, Lizzie Pearmain and Jonathan Handley, 2019
## partly based on now deprecated functions polyCount and thresholdRaster written by Phil Taylor and Mark Miller in 2011

## findKBA uses the Utilisation Distribution of tracked individuals to identify areas where core UDs overlap.
## The function first calculates thresholds based on representativeness of the tracking data (as quantified by 'repAssess')
## findKBA summarises the number of individual core UDs overlapping in an area and compares that number against thresholds.
## Output is a map and a simple feature with polygons of potential Key Areas for Biodiversity.
## NOTE: the potential KBAs must be assessed against global criteria (conservation status and global population size of the species) to determine whether they meet KBA criteria: http://www.keybiodiversityareas.org/what-are-kbas


## Polys must be an 'estUDm' object created by 'estSpaceUse' or 'adehabitatHR::kernelUD'
## representativity is the output value provided by 'repAssess' (number between 0-1).
## Col.size (optional) is the number of individuals breeding or residing at the origin location from where animals were tracked
## NOTE: if colony size is not provided the output will be as the proportion (0-100%) of colony size

findIBA <- function(KDE.Surface, representativity, Col.size = NA, UDLev=50, plotit=TRUE){
  
  #### LOAD PACKAGES ####
  # require(adehabitatHR)
  # require(raster)
  # require(sp)
  # require(sf)
  # require(smoothr)
  # require(rnaturalearth)
  pkgs <-c('sp', 'sf','smoothr','raster','rnaturalearth','raster','tidyverse', 'geosphere', 'adehabitatHR')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}
  
  
  #### ERROR CHECKING ####
  if(class(KDE.Surface) %in% c("list")) {KDE.Surface<-KDE.Surface$KDE.Surface}  ### for users who specified polyOut=T in the batchUD function
  if(!class(KDE.Surface) %in% c("estUDm")) stop("KDE.Surface should be of class 'estUDm' provided by adehabitatHR::kernelUD or track2iba::batchUD")
  if(length(KDE.Surface)<10) warning("LOW SAMPLE SIZE: identifying an IBA based on <10 tracked individuals is not recommended")
  

  #### CALCULATING THRESHOLD OF PROP OF TRACKED ANIMALS NEEDED FROM LASCELLES ET AL. 2016 #### 
  #threshlkup<-data.frame(rep=c(0.9,0.8,0.7),thresh=c(10,12.5,20),corr=c(0.9,0.75,0.5))
  if (representativity<0.7) warning("UNREPRESENTATIVE SAMPLE: you either did not track a sufficient number of birds to characterise the colony's space use or your species does not lend itself to IBA identification due to its dispersed movement")
  thresh<-ifelse(representativity<=0.7,length(KDE.Surface)*0.5,
                 ifelse(representativity<0.8,length(KDE.Surface)*0.2,
                        ifelse(representativity<0.9,length(KDE.Surface)*0.125,length(KDE.Surface)*0.1))) 
  corr<-ifelse(representativity<=0.7,0.25,
                 ifelse(representativity<0.8,0.5,
                        ifelse(representativity<0.9,0.75,0.9))) 
  
  

  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### CONVERTING OUTPUT TO PROPORTIONAL UD FOR EACH INDIVIDUAL  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## create SpatialPixelsDataFrame
  #UDLev=50
  KDEpix <- estUDm2spixdf(KDE.Surface)
  if(is.projected(KDEpix)!=TRUE) stop("Please re-calculate your kernel UD after projecting the data into a coordinate reference system where units are identical on x- and y-axis")
  
  ## calculate area of each pixel
  pixArea<-KDE.Surface[[1]]@grid@cellsize[1]
  
  
  
  ## output reported by kernelUD is intensity / m2
  ## this intensity multiplied by pixel area (in m2) 
  ## this usage sums to 1 for each individual (some individuals bordering the grid may not sum to 1)
  ## we sort this usage and calculate the cumulative sum -> this is effectively the "% UD"
  ## to find the 50% UD for an individual, simply use all grid cells where the output value is <0.5
  
  KDEpix@data <- KDEpix@data %>% 
    mutate(rowname=1:nrow(KDEpix@data)) %>%
    gather(key="ID",value="UD",-rowname) %>%
    mutate(usage=UD*(pixArea^2)) %>%
    arrange(ID, desc(usage)) %>%
    group_by(ID) %>%
    mutate(cumulUD = cumsum(usage)) %>%
    dplyr::select(rowname,ID,cumulUD) %>%
    arrange(rowname) %>%
    spread(key=ID, value=cumulUD) %>%
    dplyr::select(-rowname)
  
  
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #### COUNT THE NUMBER OF OVERLAPPING UD KERNELS ABOVE THE UDLev==50
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  
  
  ## convert to a 0/1 pixel depending on whether UDLev=50 was exceeded
  Noverlaps<-KDEpix
  Noverlaps@data<-as.data.frame(ifelse(Noverlaps@data<(UDLev/100),1,0)) %>%
    mutate(N_IND = rowSums(.)) %>%
    dplyr::select(N_IND)
  #plot(Noverlaps)
  
  ## CONVERT TO POTENTIAL IBA BASED ON THRESHOLDS
  potentialIBA<-Noverlaps
  potentialIBA@data<-potentialIBA@data %>% 
    mutate(IBA = ifelse(N_IND>= thresh,"potential","no"))
  if(!is.na(Col.size)){potentialIBA@data$N_animals<- corr*Col.size*(potentialIBA@data$N_IND/length(KDE.Surface))}else{   ## provide the number of ind expected if colony size is given
    potentialIBA@data$N_animals<- corr*100*(potentialIBA@data$N_IND/length(KDE.Surface))
    warning("No value for colony size provided. Output for N_animals is in % of colony size")}   ## if no colony size is given then provide output in percent

  
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#### CONVERT OUTPUT INTO POLYGONS WITH IBA ASSESSMENT INFO
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  ### the first step is very slow
IBApoly <- as(potentialIBA, "SpatialPolygonsDataFrame") ## this gives a warning that needs to be sorted
IBApoly <- subset(IBApoly, IBA=="potential")

if(dim(IBApoly@data)[1]==0) stop("No areas are used by a sufficient proportion of individuals to qualify as potential KBA.")

  ### aggregate all pixel-sized polygons into big polygons with the same number of birds 
OUTMAP <- aggregate(IBApoly, c('N_animals','N_IND','IBA'))
dim(OUTMAP@data)

  ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
IBA_sf <- st_as_sf(OUTMAP) %>%
  st_union(by_feature=T) %>%
  smoothr::smooth(method = "densify") %>%
  #drop_crumbs(threshold = units::set_units(100, km^2)) %>%
  #fill_holes(threshold = units::set_units(100, km^2)) %>%
  st_transform(4326)



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#### RETURN SIMPLE FEATURE WITH IBA INFO AS OUTPUT AND PLOT
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# ### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE


if(plotit == TRUE) {
  coordsets<-st_bbox(IBA_sf)

  IBAPLOT<-ggplot() +  
    geom_sf(data = IBA_sf, mapping = aes(fill=N_animals),colour="transparent") +
    coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE) +
    borders("world",fill="black",colour="black") +
    geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          axis.text=element_text(size=16, color="black"), 
          axis.title=element_text(size=16),
          panel.border = element_blank()) +
    ylab("Longitude") +
    xlab("Latitude")
  print(IBAPLOT)
  
} ## end plotit=T loop

return(IBA_sf)

} ### end findIBA function loop



