
#### NEW FUNCTION TO OVERLAP ALL INDIVIDUAL HOME RANGE CENTRES AND IDENTIFY AREAS WHERE >THRESHOLD OCCUR

### NEW FUNCTION
## combines previous polyCount and rasterThresh functions into one
## requires col size and representativity as input
## first calculates threshold based on representativity
## overlays all individual UDs and finds areas of intersection where required number of individual UDs overlap

findIBA <- function(KDE.Surface, representativity=0.8, Col.size = NA, UDLev=50, plotit=TRUE){
  
  #### LOAD PACKAGES ####
  require(adehabitatHR)
  require(raster)
  require(sp)
  require(sf)
  require(smoothr)
  require(rnaturalearth)
  
  #### ERROR CHECKING ####
  if(class(KDE.Surface) %in% c("list")) {KDE.Surface<=KDE.Surface$KDE.Surface}  ### for users who specified polyOut=T in the batchUD function
  if(!class(KDE.Surface) %in% c("estUDm")) stop("KDE.Surface should be of class 'estUDm' provided by adehabitatHR::kernelUD or track2iba::batchUD")
  if(length(KDE.Surface)<10) warning("LOW SAMPLE SIZE: identifying an IBA based on <10 tracked individuals is not recommended")
  

  #### CALCULATING THRESHOLD OF PROP OF TRACKED ANIMALS NEEDED FROM LASCELLES ET AL. 2016 #### 
  #threshlkup<-data.frame(rep=c(0.9,0.8,0.7),thresh=c(10,12.5,20),corr=c(0.9,0.75,0.5))
  if (representativity<0.7) warning("UNREPRESENTATIVE SAMPLE: you either did not track a sufficient number of birds to characterise the colony's space use or your species does not lend itself to IBA identification due to its dispersed movement")
  thresh<-ifelse(representativity<=0.7,length(KDE.Surface)*0.5,
                 ifelse(representativity<0.8,length(KDE.Surface)*0.2,
                        ifelse(representativity<0.7,length(KDE.Surface)*0.125,length(KDE.Surface)*0.1)))
  corr<-ifelse(representativity<=0.7,0.25,
                 ifelse(representativity<0.8,0.5,
                        ifelse(representativity<0.7,0.75,0.9)))
  
  
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  #### COUNT THE NUMBER OF OVERLAPPING UD KERNELS ABOVE THE UDLev==50
  ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  
  ## create SpatialPixelsDataFrame
  #UDLev=50
  KDEpix<-estUDm2spixdf(KDE.Surface)
  if(is.projected(KDEpix)!=TRUE) stop("Please re-calculate your kernel UD after projecting the data into a coordinate reference system where units are identical on x- and y-axis")

  ## find the sum of area usage for each individual to calculate what 50% is
  thresholdUD<-KDEpix@data %>% gather(key="ID",value="UD") %>%
    group_by(ID) %>%
    summarise(thresh=(sum(UD)*(UDLev/100))/nrow(KDEpix@data))

  ## convert to a 0/1 pixel depending on whether UDLev=50 was exceeded
  Noverlaps<-KDEpix
  Noverlaps@data<-Noverlaps@data %>% 
    mutate(rowname=1:nrow(KDEpix@data)) %>%
    gather(key="ID",value="UD",-rowname) %>%
    left_join(thresholdUD, by="ID") %>%
    mutate(value=ifelse(UD<thresh,0,1)) %>% 
    group_by(rowname) %>%
    summarise(N_IND=sum(value)) %>%   ### change that to max for bootstrap function
    dplyr::select(N_IND)
  
  ## CONVERT TO POTENTIAL IBA BASED ON THRESHOLDS
  potentialIBA<-Noverlaps
  potentialIBA@data<-potentialIBA@data %>% 
    mutate(IBA = ifelse(N_IND>= thresh,"potential","no"))
  if(!is.na(Col.size)){potentialIBA@data$N_birds<- corr*Col.size*(potentialIBA@data$N_IND/length(KDE.Surface))}else{   ## provide the number of ind expected if colony size is given
    potentialIBA@data$N_birds<- corr*100*(potentialIBA@data$N_IND/length(KDE.Surface))
    warning("No value for colony size provided. Output for N_birds is in % of colony size")}   ## if no colony size is given then provide output in percent

  
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#### CONVERT OUTPUT INTO POLYGONS WITH IBA ASSESSMENT INFO
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  ### the first step is very slow
IBApoly <- as(potentialIBA, "SpatialPolygonsDataFrame")
IBApoly <- subset(IBApoly, IBA=="potential")

  ### aggregate all pixel-sized polygons into big polygons with the same number of birds 
OUTMAP <- aggregate(IBApoly, c('N_birds','N_IND','IBA'))
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
  #world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))
  # todo: fix how the background world map can be plotted to IBA_sf extent

  IBAPLOT<-ggplot() +  
    geom_sf(data = IBA_sf, mapping = aes(fill=N_birds),colour="transparent") +
    #geom_sf(data = world1, mapping = aes(fill = ID), lwd = 0) +
    geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          panel.border = element_blank())
  print(IBAPLOT)
  
} ## end plotit=T loop

return(IBA_sf)

} ### end findIBA function loop



