## findKBA  #####################################################################################################################

## written by Steffen Oppel, Martin Beal, Lizzie Pearmain and Jonathan Handley, 2019
## partly based on now deprecated functions polyCount and thresholdRaster written by Phil Taylor and Mark Miller in 2011

## represent is the output value provided by 'repAssess' (number between 0-1).
## Col.size (optional) is the number of individuals breeding or residing at the origin location from where animals were tracked
## NOTE: if colony size 

#####
#' Delineating areas of aggregation of tracked animals to identify potential Key Areas for Biodiversity (KBA).
#'
#' \code{findKBA} uses the utilisation distribution of tracked individuals to identify areas where core UDs overlap.
#' The function first calculates thresholds based on representativeness of the tracking data (as quantified by \code{\link{repAssess}}) 
#' \code{findKBA} then summarises the number of individual core UDs overlapping in an area and compares that number against thresholds.
#' Output is a map and a simple polygon feature of potential Key Areas for Biodiversity.
#' The potential KBAs must be assessed against global criteria (conservation status and global population size of the species) to determine whether they meet KBA criteria. See \url{http://www.keybiodiversityareas.org/what-are-kbas} for more details. 
#'
#' @param Polys Must be an \code{estUDm} object created by \code{\link{estSpaceUse}} or \code{adehabitatHR::kernelUD}. Only accepted if the utilisation distribution was calculated in a projected coordinate reference system.
#' @param represent Numeric, between 0-1. Output value provided by \code{\link{repAssess}} which assesses how representative the tracking data are for characterising the space use of the wider population. If this value is <0.7 then a warning will be issued as the data do not meet the representativeness criteria for a KBA.
#' @param Col.size Numeric, the number of individuals breeding or residing at the origin location from where animals were tracked, quantifying the population that the tracking data represent. This number will be used to calculate how many animals use the delineated areas of aggregation. If no value for \code{Col.size} is provided then output will be as the proportion (0-100%) of colony size. 
#' @param UDLev Numeric, in % (0-100). Specifies the quantile used for delineating the core use areas of individuals based on the kernel density distribution. Default set to 50% based on Lascelles et al. (2016). For penguins higher values can be accepted, see Dias et al. (2018).
#' @param plotit Logical. If TRUE then a map of identified areas will be drawn.
#' 
#' 
#' @return An object of class \code{sf} containing polygon data with three data columns:
#'   Column \code{N_IND} indicates the number of tracked individuals whose core use area (at \code{UDLev}) overlapped with this polygon.
#'   Column \code{N_animals} estimates the number of animals from the represented population that regularly use the polygon area. If no value for (at \code{Col.size}) was provided, this number is in % of the size of the represented population.
#'   Column \code{KBA} indicates whether the polygon can be considered a potential KBA.
#'   
#' @examples
#' data(seabird)
#' Trips <- tripSplit(seabird, Colony=seabird[1,3:4], InnerBuff=5, ReturnBuff=15, Duration=2, rmColLocs = T)
#' HVALS <- findScale(Trips, ARSscale = T, Colony = seabird[1,3:4])
#' KDE.Surface <- batchUD(DataGroup=Trips, Scale = HVALS$ARSscale, polyOut=F)
#' represent <- repAssess(Trips, Scale=HVALS$ARSscale, Iteration=100)
#' findKBA(KDE.Surface, represent=represent$out)

#' 
#' \dontrun{
#' batchUD(DataGroup=Trips, Scale = HVALS$ARSscale, polyOut=F)
#' }



findKBA <- function(KDE.Surface, represent, Col.size = NA, UDLev=50, plotit=TRUE){
  
  #### LOAD PACKAGES ####
  pkgs <-c('sp', 'sf','smoothr','rnaturalearth','raster','tidyverse', 'geosphere', 'adehabitatHR')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}
  
  
  #### ERROR CHECKING ####
  if(class(KDE.Surface) %in% c("list")) {KDE.Surface<-KDE.Surface$KDE.Surface}  ### for users who specified polyOut=T in the batchUD function
  if(!class(KDE.Surface) %in% c("estUDm")) stop("KDE.Surface should be of class 'estUDm' provided by adehabitatHR::kernelUD or track2kba::estSpaceUse")
  if(length(KDE.Surface)<10) warning("LOW SAMPLE SIZE: identifying a KBA based on <10 tracked individuals is not recommended")
  

  #### CALCULATING THRESHOLD OF PROP OF TRACKED ANIMALS NEEDED FROM LASCELLES ET AL. 2016 ####
  represent<-ifelse(represent>0,represent/100,represent)   ## convert to proportion if people enter percent value
  #threshlkup<-data.frame(rep=c(0.9,0.8,0.7),thresh=c(10,12.5,20),corr=c(0.9,0.75,0.5))
  if (represent<0.7) warning("UNREPRESENTATIVE SAMPLE: you either did not track a sufficient number of birds to characterise the colony's space use or your species does not lend itself to KBA identification due to its dispersed movement")
  thresh<-ifelse(represent<=0.7,length(KDE.Surface)*0.5,
                 ifelse(represent<0.8,length(KDE.Surface)*0.2,
                        ifelse(represent<0.9,length(KDE.Surface)*0.125,length(KDE.Surface)*0.1))) 
  corr<-ifelse(represent<=0.7,0.25,
                 ifelse(represent<0.8,0.5,
                        ifelse(represent<0.9,0.75,0.9))) 
  
  

  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### CONVERTING OUTPUT TO PROPORTIONAL UD FOR EACH INDIVIDUAL  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## create SpatialPixelsDataFrame
  #UDLev=50
  KDEpix <- estUDm2spixdf(KDE.Surface)
  if(is.projected(KDEpix)!=TRUE) stop("Please re-calculate your kernel UD after projecting the data into a coordinate reference system where units are identical on x- and y-axis")
  
  ## calculate area of each pixel
  pixArea<-KDE.Surface[[1]]@grid@cellsize[1]
  KDE.Surface<-NULL
  gc()
  ## output reported by kernelUD is intensity / m2
  ## this intensity is multiplied by pixel area (in m2) 
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
  
  KDEpix<-NULL
  gc()
  
  ## CONVERT TO POTENTIAL KBA BASED ON THRESHOLDS
  potentialKBA<-Noverlaps
  potentialKBA@data<-potentialKBA@data %>% 
    mutate(KBA = ifelse(N_IND>= thresh,"potential","no"))
  if(!is.na(Col.size)){potentialKBA@data$N_animals<- corr*Col.size*(potentialKBA@data$N_IND/length(KDE.Surface))}else{   ## provide the number of ind expected if colony size is given
    potentialKBA@data$N_animals<- corr*100*(potentialKBA@data$N_IND/length(KDE.Surface))
    warning("No value for colony size provided. Output for N_animals is in % of colony size")}   ## if no colony size is given then provide output in percent
  Noverlaps<-NULL
  gc()
  
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#### CONVERT OUTPUT INTO POLYGONS WITH KBA ASSESSMENT INFO
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  ### the first step is very slow
KBApoly <- as(potentialKBA, "SpatialPolygonsDataFrame")
#KBApoly <- subset(KBApoly, KBA=="potential")
#if(dim(KBApoly@data)[1]==0) stop("No areas are used by a sufficient proportion of individuals to qualify as potential KBA.")
potentialKBA<-NULL
gc()

  ### aggregate all pixel-sized polygons into big polygons with the same number of birds 
OUTMAP <- aggregate(KBApoly, c('N_animals','N_IND','KBA'))
KBApoly<-NULL
gc()

  ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
KBA_sf <- st_as_sf(OUTMAP) %>%
  st_union(by_feature=T) %>%
  smoothr::smooth(method = "densify") %>%
  #drop_crumbs(threshold = units::set_units(100, km^2)) %>%
  #fill_holes(threshold = units::set_units(100, km^2)) %>%
  st_transform(4326)
OUTMAP<-NULL
gc()


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#### RETURN SIMPLE FEATURE WITH KBA INFO AS OUTPUT AND PLOT
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# ### CREATE MULTIPANEL PLOT OF FORAGING TRIPS WITH INCOMPLETE TRIPS SHOWN AS DASHED LINE

if(plotit == TRUE) {
  coordsets<-st_bbox(KBA_sf)

  KBAPLOT<- KBA_sf %>% filter(KBA=="potential") %>%
    ggplot() +  
    geom_sf(mapping = aes(fill=N_animals),colour="transparent") +
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
  print(KBAPLOT)
  
} ## end plotit=T loop

return(KBA_sf)

} ### end findKBA function loop



