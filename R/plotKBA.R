

plotKBA <- function(KBAobj, Colony=NULL){
  
  if(any(class(KBAobj)=="sf")){
    
    KBA_sf <- KBAobj
    
    coordsets <- sf::st_bbox(KBA_sf)
    
    KBAPLOT <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>%
      ggplot() +
      geom_sf(mapping = aes(fill=N_animals), colour=NA, alpha=1) +
      borders("world", fill="dark grey", colour="grey20") +
      coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE) +
      theme(panel.background=element_blank(),
        panel.grid.major=element_line(colour="transparent"),
        panel.grid.minor=element_line(colour="transparent"),
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      guides(colour=FALSE) +
      scale_fill_continuous(name = "N animals") +
      ylab("Latitude") +
      xlab("Longitude")
    
    if(max(na.omit(KBA_sf$N_animals) < 1.1)) { ## make legend title proportion
      KBAPLOT <- KBA_sf %>% dplyr::filter(.data$potentialKBA==TRUE) %>%
        ggplot() +
        borders("world", fill="dark grey", colour="grey20") +
        geom_sf(mapping = aes(fill=N_animals, colour="transparent")) +
        coord_sf(xlim = c(coordsets$xmin, coordsets$xmax), ylim = c(coordsets$ymin, coordsets$ymax), expand = FALSE) +
        theme(panel.background=element_blank(),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          axis.text=element_text(size=14, colour="black"),
          axis.title=element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        guides(colour=FALSE) +
        scale_fill_continuous(name = "Prop. of animals") +
        ylab("Latitude") +
        xlab("Longitude")
    }
    
    if(!is.null(Colony)){
      KBAPLOT <- KBAPLOT + 
        geom_point(data=Colony, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2)
    }
    
    print(KBAPLOT)
    return(KBAPLOT)
    
  } else if(class(KBAobj)==c("SpatialPixelsDataFrame", "SpatialGridDataFrame")){
    
    KBA_sp <- KBAobj
    
    plot(KBA_sp[KBA_sp$N_animals>0, ])
    
  }
  
}

