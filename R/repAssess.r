#' Assess sample representativeness
#'
#' \code{repAssess} estimates the degree to which the space use of a tracked sample of animals represents that of the larger population. 
#'
#' Representativeness is assessed by calculating the proportion of out-sample points included in in-sample space use areas.
#'
#' First, the set of IDs is iteratively sub-sampled, and in each iteration a set of individual UDs are pooled and the points of the un-selected (outsample) IDs are overlaid on the % contour area ('UDLev') of the KDE. The proportion of these outsample points which overlap the pooled KDE area is known as the inclusion rate, and represents an estimate of representativeness at each sample size. Then, a non-linear function is fit to the relationship between the inclusion rate and sample size (i.e. number of tracks/animals) in order to estimate the point at which the relationship reaches an information asymptote (i.e. no more information added per track). By estimating this asymptote, \code{repAssess} provides a estimate of representativeness by dividing the inclusion rate estiamted at the maximum sample size, by this asymptote. Finally, using this relationship, a minimum representative sample sizes are also calculated.
#' 
#' \code{\link{repAssess}} accepts UDs calculated outside of \code{track2KBA}, if they have been converted to class \code{RasterStack} or \code{SpatialPixelsDataFrame}. However, one must make sure that the cell values in these cases are probability densities (i.e. <=0) not volumne contour values (i.e. 0-1, or 0-100).
#' 
#' When setting \code{avgMethod} care must be taken. If the input are classic KDEs (e.g. from \code{\link{estSpaceUse}}) then the weighted mean is likely the optimal way to pool individual UDs. However, if any other method (for example AKDE, auto-correlated KDE) was used to estimate UDs, then the arithmetic mean is the safer option. 
#'
#' @param DataGroup SpatialPointsDataFrame or data.frame of animal relocations. Must include 'ID' field. If input is data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param KDE Several input options: an estUDm, a SpatialPixels/GridDataFrame, a list object, or a RasterStack. If estUDm, as created by \code{\link{estSpaceUse}} or \code{adehabitatHR::kernelUD}, if Spatial*, each column should correspond to the Utilization Distribution of a single individual or track, and if a list it should be output from \code{\link{estSpaceUse}} when the argument \code{polyOut = TRUE}. If a RasterStack, each layer must be an individual UD. If \code{KDE} is not supplied, then UDs will be produced by applying DataGroup to \code{estSpaceUse} using the input \code{Scale} value.
#' @param Iteration numeric. Number of times to repeat sub-sampling procedure. The higher the iterations, the more robust the result. 
#' @param Scale numeric. This value sets the smoothing (h) parameter for Kernel Density Estimation. Only needs to be set if nothing is supplied to \code{KDE}.
#' @param Res numeric. Grid cell resolution (in square kilometers) for kernel density estimation. Default is a grid of 500 cells, with spatial extent determined by the latitudinal and longitudinal extent of the data. Only needs to be set if nothing is supplied to \code{KDE}.
#' @param UDLev numeric. Specify which contour of the Utilization Distribution the home range estimate (\code{KDE}) represented (e.g. 50, 95). 
#' @param avgMethod character. Choose whether to use the arithmetic or weighted mean when combining individual IDs. Options are :'mean' arithmetic mean, or 'weighted', which weights each UD by the numner of points per level of ID.
#' @param BootTable logical (TRUE/FALSE). Do you want to save the full results table to the working directory?
#' @param Ncores numeric. The number of processing cores to use. For heavy operations, the higher the faster. NOTE: CRAN sets a maximum at 2 cores. If using the git-hub version of the package, this is set a maximum of one fewer than the maximum cores in your computer.
#'  
#' @return A single-row data.frame, with columns '\emph{SampleSize}' signifying the maximum sample size in the data set, '\emph{out}' signifying the percent representativeness of the sample, '\emph{type}' is the type  of asymptote value used to calculate the '\emph{out}' value, and '\emph{asym}' is the asymptote value used.
#'
#' There are two potential values for '\emph{type}': 'type' == 'asymptote' is the ideal, where the asymptote value is calculated from the parameter estimates of the successful nls model fit. Secondly, when nls fails to converge at all, then the mean inclusion rate is taken for the largest sample size; 'type'=='inclusion.'minRep' signifies the sample size at which ~70% of the space use information is encompassed, and 'fullrep' signifies the sample size approaching the asymptote, i.e. representing 99% of the space use.
#'
#' @examples
#' \dontrun{repr <- repAssess(Trips, Scale=10, Iteration=1, BootTable = F, n.cores = 1)}
#'
#' @export
#' @importFrom foreach %dopar%
#' @importFrom graphics abline identify lines points polygon text

repAssess <- function(DataGroup, KDE=NULL, Iteration=50, Scale=NULL, Res=NULL, UDLev=50, avgMethod = "mean", Ncores=1, BootTable=FALSE){

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
      TripCoords@data <- TripCoords@data %>% dplyr::select(.data$ID)
    }
  }
  
  UIDs <- unique(TripCoords$ID)
  
  # Determine class of KDE, and convert to Raster
  if(is.null(KDE)){
    stop("No Utilization Distributions supplied to KDE argument. See track2KBA::estSpaceUse()")
  } else if(class(KDE) == "estUDm"){ 
    KDEraster <- raster::stack(lapply(KDE, function(x) raster::raster(as(x,"SpatialPixelsDataFrame"), values=T)))
  } else if(class(KDE) == "SpatialPixelsDataFrame") {
    KDEraster <- raster::stack(KDE)
  } else if(class(KDE) == "RasterStack") {
    KDEraster <- KDE
  }
  
  # assure only IDs with UDs are in tracking data
  UDnames <- names(KDEraster) 
  if( length(UDnames) != length(UIDs) ) {
    if( all(stringr::str_detect(names(KDEraster), pattern = "^X")) ){
      UDnames <- substring(UDnames, 2) # remove "X" from raster names
      TripCoords <- TripCoords[TripCoords$ID %in% UDnames, ]
    } else {
      TripCoords <- TripCoords[TripCoords$ID %in% UDnames, ]
    }
  }
  
  UIDs <- unique(TripCoords$ID) # get IDs, after any filtered out
  NIDs <- length(UIDs)

  if(NIDs<50){Nloop <- seq(1, (NIDs - 1), 1)}
  if(NIDs>=50 & NIDs<100){Nloop <- c(seq(1, 19, 1), seq(20, (NIDs - 1), 3))}
  if(NIDs>=100){Nloop <- c(seq(1, 20, 1), seq(21, 49, 3), seq(50, (NIDs - 1), 6))}
  
  DoubleLoop <- data.frame(SampleSize = rep(Nloop, each=Iteration), Iteration=rep(seq(1:Iteration), length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  
  ###
  if(Ncores > 1){
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package \"parallel\" needed for this function to work. Please install it.",
        call. = FALSE)  }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package \"doParallel\" needed for this function to work. Please install it.",
        call. = FALSE)  }
    maxCores <- parallel::detectCores()
    Ncores <- ifelse(Ncores == maxCores, Ncores - 1, Ncores) # ensure that at least one core is un-used
    cl <- parallel::makeCluster(Ncores)
    doParallel::registerDoParallel(cl)
  }

  Result <- data.frame()

  Result <- foreach::foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp", "dplyr", "raster")) %dopar% {
    
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
    
    if(avgMethod=="weighted"){
      weights <- as.vector(unname(table(SelectedTracks$ID)))   # get number of points per ID
      KDEcmbnd <- raster::weighted.mean(KDEstack, w=weights)   # weighted mean
    } else if( avgMethod == "mean"){
      KDEcmbnd <- raster::calc(KDEstack, mean)                 # arithmetic mean
    }
    
    ### Calculating inclusion value, using Kernel surface ######
    KDElev <- KDEcmbnd
    pixArea <- raster::res(KDElev)[1]
    
    ### original ## 
    df <- data.frame(UD = raster::getValues(KDElev)) %>%
      mutate(rowname = 1:length(raster::getValues(KDElev))) %>%
      mutate(usage = .data$UD * (pixArea^2)) %>%
      arrange(desc(.data$usage)) %>%
      mutate(cumulUD = cumsum(.data$usage)) %>%
      mutate(INSIDE = ifelse(.data$cumulUD < (UDLev/100), 1, NA)) %>%
      arrange(.data$rowname) %>%
      dplyr::select(.data$INSIDE)
    
    KDElev[] <- df$INSIDE
    
    # plot(KDElev)
    Overlain_Raster <- raster::extract(KDElev, NotSelected)
    
    Output$InclusionMean <- length(which(!is.na(Overlain_Raster)))/nrow(NotSelected)
    
    return(Output)
  }
  ## stop the cluster
  if(Ncores > 1) {on.exit(parallel::stopCluster(cl))}

  try(M1 <- stats::nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
  if ('M1' %in% ls()){       ### run this only if nls was successful
    a <- base::summary(M1)$coefficients[1]
    b <- base::summary(M1)$coefficients[2]
    
    tAsymp <- (UDLev/100)
    Asymptote <- (a / b)
    Result$pred <- stats::predict(M1)
    
    ## Calculate RepresentativeValue 
    RepresentativeValue <- Result %>%
      group_by(.data$SampleSize) %>%
      summarise(
        out      = max(pred) / Asymptote*100
        ) %>%
      dplyr::filter(out == max(.data$out)) %>%
      mutate(
        est_asym = Asymptote,
        tar_asym = (UDLev/100)
        ) 
    
    # Calculate minimum and fully representative sample sizes, and convert inclusions into rep. estimates
    minRep  <- 0.7*Asymptote
    fullRep <- 0.95*Asymptote
    RepresentativeValue$minRep  <- ceiling(minRep / (a - (minRep*b)))
    RepresentativeValue$fullRep <- ceiling(fullRep / (a - (fullRep*b)))
    
    Result <- Result %>% mutate(
      rep_est = .data$pred / Asymptote*100, 
      is_rep  = ifelse(.data$rep_est >= 70, TRUE, FALSE)
    )

    ## Plot
    P2 <- Result %>%
      group_by(.data$SampleSize) %>%
      dplyr::summarise(
        meanPred = mean(na.omit(.data$pred)),
        sdInclude = sd(.data$InclusionMean))
    
    yTemp <- c(P2$meanPred + Asymptote * P2$sdInclude, rev(P2$meanPred - Asymptote * P2$sdInclude))
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
      dplyr::filter(.data$SampleSize == max(.data$SampleSize)) %>%
      group_by(.data$SampleSize) %>%
      summarise(
        out = mean(.data$InclusionMean)
        ) %>%
      mutate(type = 'inclusion') %>%
      mutate(asym = .data$out)
  }
  
  if(BootTable==T){
    utils::write.csv(Result,"bootout_temp.csv", row.names=F)
  }
  
  print(ifelse(exists("M1"),"nls (non linear regression) successful, asymptote estimated for bootstrap sample.",
    "WARNING: nls (non linear regression) unsuccessful, likely due to 'singular gradient' (e.g. small sample), which means there is no asymptote. Data may not be representative, output derived from mean inclusion value at highest sample size. Check bootstrap output csv file"))
  if(Asymptote < (tAsymp - 0.1) ) { warning("Estimated asymptote differs from target; be aware that representativeness value is based on estimated asymptote (i.e. est_asym).") }
  
  return(as.data.frame(RepresentativeValue))
}


