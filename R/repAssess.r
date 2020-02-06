#' Assess sample representativeness
#'
#' \code{repAssess_} estimates the degree to which the space use of a tracked sample of animals represents that of the larger population. 
#'
#' Representativeness is assessed by calculating the proportion of out-sample points included in in-sample space use areas.
#'
#' First, if no list of utilization distributions is supplied, estSpaceUse is called to estimate UDs for each ID in the tracking dataset. Then, the set of IDs is iteratively sub-sampled, and in each iteration the individual UDs are pooled and the points of the un-selected (outsample) IDs are overlayed on the 50\% contour area of the KDE. The proportion of these outsample points which overlap the pooled KDE area (i.e. the inclusion rate) are taken as the estimate of representativeness at that sample size. This process is iterated as many times as desired at each sample size. A representative dataset would approach an inclusion rate of 0.5.
#' 
#' By fitting a trend line to the relationship between sample size inclusion rate we can identify the sample size at which the curve approaches an asymptote, signifying that any new data would simply add to existing knowledge.
#'
#' @param DataGroup SpatialPointsDataFrame or data.frame of animal relocations. Must include 'ID' field. If input is data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param KDE Several input options: an estUDm, a SpatialPixels/GridDataFrame, or a list object. If estUDm, as created by \code{\link{estSpaceUse}} or \code{adehabitatHR::kernelUD}, if Spatial*, each column should correspond to the Utilization Distribution of a single individual or track, and if a list it should be output from \code{\link{estSpaceUse}} when the argument \code{polyOut = TRUE}. If \code{KDE} is not supplied, then UDs will be produced by applying DataGroup to \code{estSpaceUse} using the input \code{Scale} value.
#' @param Iteration numeric. Number of times to repeat sub-sampling procedure. The higher the iterations, the more robust the result. 
#' @param Scale numeric. This value sets the smoothing (h) parameter for Kernel Density Estimation. Only needs to be set if nothing is supplied to \code{KDE}.
#' @param Res numeric. Grid cell resolution (in square kilometers) for kernel density estimation. Default is a grid of 500 cells, with spatial extent determined by the latitudinal and longitudinal extent of the data. Only needs to be set if nothing is supplied to \code{KDE}.
#' @param BootTable logical (TRUE/FALSE). Do you want to save the full results table to the working directory?
#'  
#' @return A single-row data.frame, with columns '\emph{SampleSize}' signifying the maximum sample size in the data set, '\emph{out}' signifying the percent representativeness of the sample, '\emph{type}' is the type  of asymptote value used to calculate the '\emph{out}' value, and '\emph{asym}' is the asymptote value used.
#'
#' There are three potential values for '\emph{type}': 'asymptote' is the ideal, where the asymptote value is calculated from the parameter estimates of the successful nls model fit. Sometimes nls fit is successful but the curve either flips or does not approach 0.5. In these cases the 'type' is 'asymptote_adj' and the inclusion rate is compared to an artificially adjusted asymptote value of 0.5. Lastly, when nls fails to converge at all, then the mean inclusion rate is taken for the largest sample size; 'type'=='inclusion.'
#'
#' @examples
#' \dontrun{repr <- repAssess(Trips, Scale=10, Iteration=1, BootTable = F, n.cores = 1)}
#'
#' @export
#' @importFrom foreach %do%


repAssess <- function(DataGroup, KDE=NULL, Iteration=50, Scale=NULL, Res=NULL, BootTable=FALSE){
  
  pkgs <- c('sp', 'geosphere', 'adehabitatHR','foreach','dplyr','data.table', 'raster')
  for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}
  
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
    if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
    ## set the minimum fields that are needed
    CleanDataGroup <- DataGroup %>%
      dplyr::select(.data$ID, .data$Latitude, .data$Longitude, .data$DateTime) %>%
      arrange(.data$ID, .data$DateTime)
    mid_point <- data.frame(geosphere::centroid(cbind(CleanDataGroup$Longitude, CleanDataGroup$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanDataGroup$Longitude) < -170 &  max(CleanDataGroup$Longitude) > 170) {
      longs = ifelse(CleanDataGroup$Longitude < 0, CleanDataGroup$Longitude + 360, CleanDataGroup$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
    
    DataGroup.Wgs <- SpatialPoints(data.frame(CleanDataGroup$Longitude, CleanDataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=proj.UTM )
    TripCoords <- SpatialPointsDataFrame(DataGroup.Projected, data = CleanDataGroup)
    TripCoords@data <- TripCoords@data %>% dplyr::select(.data$ID)
    
  }else{  ## if data are already in a SpatialPointsDataFrame then check for projection
    if(is.projected(DataGroup)){
      TripCoords <- DataGroup
      TripCoords@data <- TripCoords@data %>% dplyr::select(.data$ID)
    }else{ ## project data to UTM if not projected
      if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
      if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
      
      mid_point <- data.frame(geosphere::centroid(cbind(DataGroup@data$Longitude, DataGroup@data$Latitude)))
      
      ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
      if (min(DataGroup@data$Longitude) < -170 &  max(DataGroup@data$Longitude) > 170) {
        longs = ifelse(DataGroup@data$Longitude < 0, DataGroup@data$Longitude + 360,DataGroup@data$Longitude)
        mid_point$lon<-ifelse(median(longs) > 180, median(longs)-360, median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      TripCoords <- sp::spTransform(DataGroup, CRS=proj.UTM)
      TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
    }
  }
  
  UIDs <- unique(TripCoords$ID)
  NIDs <- length(UIDs)

  if(NIDs<50){Nloop <- seq(1, (NIDs - 1), 1)}
  if(NIDs>=50 & NIDs<100){Nloop <- c(seq(1, 19, 1), seq(20, (NIDs - 1), 3))}
  if(NIDs>=100){Nloop <- c(seq(1, 20, 1), seq(21, 49, 3), seq(50, (NIDs - 1), 6))}
  
  DoubleLoop <- data.frame(SampleSize = rep(Nloop, each=Iteration), Iteration=rep(seq(1:Iteration), length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  
  # Determine class of KDE, and convert to Raster
  if(is.null(KDE)){
    if(is.null(Res)) { Res <- 100 }
    KDE.Surface <- estSpaceUse(DataGroup=TripCoords, Scale = Scale, Res = Res, UDLev = 50, polyOut=F)
    KDEraster <- stack(lapply(KDE.Surface, function(x) raster::raster(x, values=T)))
    
  } else if(class(KDE) == "list") { 
    KDE.Surface <- KDE$KDE.Surface 
    KDEraster <- stack(lapply(KDE.Surface, function(x) raster::raster(x, values=T)))
    
  } else { 
      KDE.Surface <- KDE 
      KDEraster <- stack(lapply(KDE, function(x) raster(as(x,"SpatialPixelsDataFrame"), values=T)))
      # KDEraster <- raster::stack(KDE.Surface)
  }
  
  # convert estSpaceUse output (list of estUDs) to RasterLayer list
  # KDEraster <- lapply(KDE.Surface, function(x) raster::raster(x, values=T))
  
  ###
  
  Result <- data.frame()
  
  Result <- foreach::foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp", "dplyr", "raster")) %do% {
    
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)
    
    RanNum <- sample(UIDs, N, replace=F)
    NotSelected <- TripCoords[!TripCoords$ID %in% RanNum,]
    SelectedTracks <- TripCoords[TripCoords$ID %in% RanNum,]
    if(all(stringr::str_detect(names(KDEraster), pattern = "^X"))){
      Selected <- KDEraster[[paste("X", RanNum, sep = "")]]
    } else {
      Selected <- KDEraster[[ RanNum ]]
    }

    KDEstack <- raster::stack(Selected)  # list of RasterLayers to RasterStack
    # KDEstack <- Selected
    weights <- as.vector(unname(table(SelectedTracks$ID))) # see how uneven sample sizes are per individual

    KDEcmbnd <- raster::weighted.mean(KDEstack, w=weights)
    # KDEcmbnd <- raster::calc(KDEstack, mean)  # average together individual UDs
    
    ### Calculating inclusion value, using Kernel surface ######
    KDElev <- KDEcmbnd
    pixArea <- res(KDElev)[1]
    
    ### original ## 
    df <- data.frame(UD = getValues(KDElev)) %>%
      mutate(rowname = 1:length(getValues(KDElev))) %>%
      mutate(usage = .data$UD * (pixArea^2)) %>%
      arrange(desc(.data$usage)) %>%
      mutate(cumulUD = cumsum(.data$usage)) %>%
      mutate(INSIDE = ifelse(.data$cumulUD < 0.5, 1, NA)) %>%
      arrange(.data$rowname) %>%
      dplyr::select(.data$INSIDE)
    
    ### volume UD change ### 
    
    # maxVol <- max(getValues(KDElev))
    # if(maxVol > 1){ values(KDElev) <- values(KDElev) / 100 } # if pixel values are from adehabitat, change them to 0-1
    
    # df <- data.frame(UD = getValues(KDElev)) %>%
    #   mutate(rowname = 1:length(getValues(KDElev))) %>% 
    #   mutate(INSIDE = ifelse(.data$UD < 0.5, 1, NA)) %>%
    #   arrange(.data$rowname) %>% 
    #   dplyr::select(.data$INSIDE)
    # 
    KDElev[] <- df$INSIDE
    
    # plot(KDElev)
    Overlain_Raster <- raster::extract(KDElev, NotSelected)
    
    Output$InclusionMean <- length(which(!is.na(Overlain_Raster)))/nrow(NotSelected)
    
    return(Output)
  }

  if(BootTable==T){
    write.csv(Result,"bootout_temp.csv", row.names=F)
  }
  
  try(M1 <- stats::nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
  if ('M1' %in% ls()){       ### run this only if nls was successful
    Asymptote <- (base::summary(M1)$coefficients[1] / summary(M1)$coefficients[2])
    Result$pred <- stats::predict(M1)
    
    ## Calculate RepresentativeValue 
    RepresentativeValue <- Result %>%
      group_by(SampleSize) %>%
      summarise(out = max(pred) / ifelse(Asymptote < 0.45, 0.5, Asymptote)*100) %>%
      dplyr::filter(out == max(.data$out)) %>%
      mutate(type = ifelse(Asymptote < 0.45, 'asymptote_adj', 'asymptote')) %>%
      mutate(asym = Asymptote) 
    
    if(Asymptote < 0.45) {
      RepresentativeValue$asym_adj <- 0.5 }
    
    ## Plot
    P2 <- Result %>%
      group_by(.data$SampleSize) %>%
      dplyr::summarise(
        meanPred = mean(na.omit(.data$pred)),
        sdInclude = sd(.data$InclusionMean))
    
    yTemp <- c(P2$meanPred + 0.5 * P2$sdInclude, rev(P2$meanPred - 0.5 * P2$sdInclude))
    xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
    
    # pdf("track2kba_repAssess_output.pdf", width=6, height=5)  ## avoids the plotting margins error
    print(plot(InclusionMean ~ SampleSize,
      data = Result, pch = 16, cex = 0.2, col="darkgray", ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))), ylab = "Inclusion", xlab = "Sample Size"))
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
    lines(P2, lty=1,lwd=2)
    text(x=0, y=0.99, paste(round(RepresentativeValue$out, 2), "%", sep=""), cex=2, col="gray45", adj=0)
    # dev.off()
    
  }else{ ### if nls is unsuccessful then use mean output for largest sample size
    RepresentativeValue <- Result %>%
      dplyr::filter(SampleSize == max(.data$SampleSize)) %>%
      group_by(.data$SampleSize) %>%
      summarise(out = mean(.data$InclusionMean)) %>%
      mutate(type = 'inclusion') %>%
      mutate(asym = .data$out)
  }
  
  print(ifelse(exists("M1"),"nls (non linear regression) successful, asymptote estimated for bootstrap sample.",
    "WARNING: nls (non linear regression) unsuccessful, likely due to 'singular gradient', which means there is no asymptote. Data may not be representative, output derived from mean inclusion value at highest sample size. Check bootstrap output csv file"))
  
  return(as.data.frame(RepresentativeValue))
  
}


