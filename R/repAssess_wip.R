#####
################### EDIT to SPEED THINGS UP!#########

### Instead producing a pooled Kernel Estimation for each iteration, estSpaceUse is run once producing individual KDEs, which are then sub-sampled (boostrapped) and averaged together to produce the group-level kernel in each iteration.

#####

#' Assess sample representativeness
#'
#' \code{repAssess} estimates the degree to which the space use of a tracked sample of animals represents that of the larger population.
#'
#' Representativeness is assessed by calculating the proportion of out-sample points included in in-sample space use areas.
#'
#' First, if no list of utilization distributions is supplied, estSpaceUse is called to estimate UDs for each ID in the tracking dataset. Then, the set of IDs is iteratively sub-sampled, and in each iteration the individual UDs are pooled and the points of the un-selected (outsample) IDs are overlayed on the 50\% contour area of the KDE. The proportion of these outsample points which overlap the pooled KDE area (i.e. the inclusion rate) are taken as the estimate of representativeness at that sample size. A representative dataset would approach an inclusion rate of 0.5.
#' By fitting a trend line to the relationship between sample size inclusion rate we can identify the sample size at which the curve approaches an asymptote, signifying that any new data would simply add to existing knowledge.
#'
#' @param DataGroup SpatialPointsDataFrame or data.frame of animal relocations. Must include 'ID' field. If input is data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param Scale numeric. This value sets the smoothing (h) parameter for Kernel Density Estimation.
#' @param Iteration numeric. Number of times to repeat sub-sampling procedure.
#' @param Res numeric. Grid cell resolution (in square kilometers) for kernel density estimation. Default is a grid of 500 cells, with spatial extent determined by the latitudinal and longitudinal extent of the data.
#' @param BootTable logical (TRUE/FALSE). Do you want to save the full results table to the working directory?
#' @param Ncores numeric. The number of processing cores to use in parallel processing. Check how many are available with parallel::detectCores().
#' @return A single-row data.frame, with columns '\emph{SampleSize}' signifying the maximum sample size in the data set, '\emph{out}' signifying the percent representativeness of the sample, '\emph{type}' is the type of asymptote value used to calculate the '\emph{out}' value, and '\emph{asym}' is the asymptote value used.
#'
#' There are three potential valules for '\emph{type}': 'asymptote' is the ideal, where the asymptote value is calculated from the parameter estimates of the successful nls model fit. Sometimes nls fit is successful but the curve either flips or does not approach 0.5. In these cases the 'type' is 'asymptote_adj' and the inclusion rate is compared to an artificially adjusted asymptote value of 0.5. Finally, when nls fails, then the mean inclusion rate is taken for the largest sample size; 'type'=='inclusion.'
#'
#' @examples
#' \dontrun{repr <- repAssess(Trips, Scale=10, Iteration=1, BootTable = F, n.cores = 1)}
#'
#' @export
#' @importFrom foreach %dopar%


repAssess_wip <- function(DataGroup, listKDE=NA, Iteration=10, Scale=NA, Res=NA, BootTable=FALSE, Ncores=1){
  
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
    if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
    ## set the minimum fields that are needed
    CleanDataGroup <- DataGroup %>%
      dplyr::select(ID, Latitude, Longitude,DateTime) %>%
      arrange(ID, DateTime)
    mid_point <- data.frame(centroid(cbind(CleanDataGroup$Longitude, CleanDataGroup$Latitude)))
    
    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanDataGroup$Longitude) < -170 &  max(CleanDataGroup$Longitude) > 170) {
      longs = ifelse(CleanDataGroup$Longitude < 0, CleanDataGroup$Longitude + 360, CleanDataGroup$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}
    
    DataGroup.Wgs <- SpatialPoints(data.frame(CleanDataGroup$Longitude, CleanDataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=proj.UTM )
    TripCoords <- SpatialPointsDataFrame(DataGroup.Projected, data = CleanDataGroup)
    TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
    
  }else{  ## if data are already in a SpatialPointsDataFrame then check for projection
    if(is.projected(DataGroup)){
      TripCoords <- DataGroup
      TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
    }else{ ## project data to UTM if not projected
      if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
      if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
      
      mid_point <- data.frame(centroid(cbind(DataGroup@data$Longitude, DataGroup@data$Latitude)))
      
      ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
      if (min(DataGroup@data$Longitude) < -170 &  max(DataGroup@data$Longitude) > 170) {
        longs = ifelse(DataGroup@data$Longitude < 0, DataGroup@data$Longitude + 360,DataGroup@data$Longitude)
        mid_point$lon<-ifelse(median(longs) > 180, median(longs)-360, median(longs))}
      
      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      TripCoords <- spTransform(DataGroup, CRS=proj.UTM)
      TripCoords@data <- TripCoords@data %>% dplyr::select(ID)
    }
  }
  
  UIDs <- unique(TripCoords$ID)
  NIDs <- length(UIDs)
  
  Nloop <- seq(1, (NIDs - 1), 1)
  
  DoubleLoop <- data.frame(SampleSize = rep(Nloop, each=Iteration), Iteration=rep(seq(1:Iteration), length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  
  # first case scenario: no listKDE is supplied
  if(is.na(listKDE)){
    if(is.na(Res)) { Res <- 100 }
    KDE.Surface <- estSpaceUse(DataGroup=TripCoords, Scale = Scale, Res = Res, UDLev = 50, polyOut=F)
  } else { 
    KDE.Surface <- listKDE 
  }
  
  # convert estSpaceUse output (list of estUDs) to RasterLayer list
  KDEraster <- lapply(KDE.Surface, function(x) raster(x, values=T))
  
  ###
  
  Ncores <- ifelse(Ncores==1, detectCores()/2, Ncores) ## use user-specified value if provided to avoid computer crashes by using only half the available cores
  cl <- makeCluster(n.cores)  
  registerDoParallel(cl)
  Result <- data.frame()
  
  Result <- foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp", "dplyr", "raster")) %dopar% {
    
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)
    
    RanNum <- sample(UIDs, N, replace=F)
    NotSelected <- TripCoords[!TripCoords$ID %in% RanNum,]
    Selected <- KDEraster[names(KDEraster) %in% RanNum]
    
    KDEstack <- stack(Selected)  # list of RasterLayers to RasterStack
    KDEcmbnd <- calc(KDEstack, mean)  # average together individual UDs
    
    ### Calculating inclusion value, using Kernel surface ######
    KDElev <- KDEcmbnd
    pixArea <- res(KDElev)[1]
    
    df <- data.frame(UD = getValues(KDElev)) %>%
      mutate(rowname = 1:length(getValues(KDElev))) %>%
      mutate(usage = UD*(pixArea^2)) %>%
      arrange(desc(usage)) %>%
      mutate(cumulUD = cumsum(usage)) %>%
      mutate(INSIDE = ifelse(cumulUD < 0.5, 1, NA)) %>%
      arrange(rowname) %>%
      dplyr::select(INSIDE)
    
    KDElev[] <- df$INSIDE
    
    # plot(KDElev)
    Overlain_Raster <- raster::extract(KDElev, NotSelected)
    
    Output$InclusionMean <- length(which(!is.na(Overlain_Raster)))/nrow(NotSelected)
    
    return(Output)
  }
  ## stop the cluster
  on.exit(stopCluster(cl))
  
  if(BootTable==T){
    data.table::fwrite(Result,"bootout_temp.csv", row.names=F, sep=",")
  }
  
  try(M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
  if ('M1' %in% ls()){       ### run this only if nls was successful
    Asymptote <- (base::summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
    Result$pred <- stats::predict(M1)
    
    ## Calculate RepresentativeValue 
    RepresentativeValue <- Result %>%
      group_by(SampleSize) %>%
      summarise(out = max(pred) / ifelse(Asymptote < 0.45, 0.5, Asymptote)*100) %>%
      dplyr::filter(out == max(out)) %>%
      mutate(type = ifelse(Asymptote < 0.45, 'asymptote_adj', 'asymptote')) %>%
      mutate(asym = Asymptote) 
    
    if(Asymptote < 0.45 | Asymptote > 60) {
      RepresentativeValue$asym_adj <- 0.5 }
    
    ## Plot
    P2 <- Result %>% 
      group_by(SampleSize) %>% 
      dplyr::summarise(
        meanPred = mean(na.omit(pred)),
        sdInclude = sd(InclusionMean))
    yTemp <- c(P2$meanPred + 0.5 * P2$sdInclude, rev(P2$meanPred - 0.5 * P2$sdInclude))
    xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
    pdf("track2kba_repAssess_output.pdf",width=6, height=5)  ## avoids the plotting margins error
    plot(InclusionMean ~ SampleSize, 
      data = Result, pch = 16, cex = 0.2, col="darkgray", ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))), ylab = "Inclusion", xlab = "SampleSize")
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
    lines(P2, lty=1,lwd=2)
    text(x=0, y=1,paste(round(RepresentativeValue$out, 2), "%", sep=""), cex=2, col="gray45", adj=0)  
    dev.off()
  }else{ ### if nls is unsuccessful then use mean output for largest sample size
    RepresentativeValue <- Result %>%
      filter(SampleSize==max(SampleSize)) %>%
      group_by(SampleSize) %>%
      summarise(out=mean(InclusionMean)) %>%
      mutate(type='inclusion')%>%
      mutate(asym=out)
  }
  
  print(ifelse(exists("M1"),"nls (non linear regression) successful, asymptote estimated for bootstrap sample.",
    "WARNING: nls (non linear regression) unsuccessful, likely due to 'singular gradient', which means there is no asymptote. Data may not be representative, output derived from mean inclusion value at highest sample size. Check bootstrap output csv file"))
  
  return(RepresentativeValue)
  
}


