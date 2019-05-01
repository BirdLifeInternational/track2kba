## repAssess  #####################################################################################################################

#' Assess sample representativeness
#'
#' \code{repAssess} estimates the degree to which the space use of a tracked sample of animals represents that of the larger population.
#'
#' Representativeness is assessed by calculating the proportion of out-sample points included in in-sample space use areas.
#'
#' First, the data set is iteratively sub-sampling, producing a pooled utilization distribution (via Kernel Density Estimation), and isolating the 50\% contour. Then, the un-selected points are overlayed on the 50\% contour area and the proportion which overlap is calculated. A representative dataset would approach an inclusion rate of 0.5.
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
#'
repAssess <- function(DataGroup, Scale=10, Iteration=50, Res=100, BootTable=T, Ncores=1)
{
  ## do we need geosphere? PLEASE CHECK!
  # pkgs <- c('sp', 'geosphere', 'adehabitatHR','foreach','doParallel','dplyr','data.table', 'parallel')
  # for(p in pkgs) {suppressPackageStartupMessages(require(p, quietly=TRUE, character.only=TRUE,warn.conflicts=FALSE))}

  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")

  if(class(DataGroup) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
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
      mid_point <- data.frame(centroid(cbind(DataGroup@data$Longitude, DataGroup@data$Latitude)))

      ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
      if (min(DataGroup@data$Longitude) < -170 &  max(DataGroup@data$Longitude) > 170) {
        longs = ifelse(DataGroup@data$Longitude < 0, DataGroup@data$Longitude + 360,DataGroup@data$Longitude)
        mid_point$lon<-ifelse(median(longs) > 180, median(longs)-360, median(longs))}

      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
      TripCoords <- spTransform(DataGroup, CRS=proj.UTM)
      TripCoords@data <- TripCoords@data %>% dplyr::select(.data$ID)
    }

  }

  proj.UTM <- CRS(proj4string(TripCoords))
  UIDs <- unique(TripCoords$ID)
  NIDs <- length(UIDs)

  ### N OF SAMPLE SIZE STEPS NEED TO BE SET DEPENDING ON DATASET - THIS CAN FAIL IF NID falls into a non-existent sequence
  if(NIDs<22){Nloop <- seq(1, (NIDs - 1), 1)}
  if(NIDs>=22 & NIDs<52){Nloop <- c(seq(1, 19, 1), seq(20, (NIDs - 1), 3))}
  if(NIDs>=52 & NIDs<102){Nloop <- c(seq(1, 20, 1), seq(21, 49, 3), seq(50, (NIDs - 1), 6))}
  if(NIDs>=102){Nloop <- c(seq(1, 20, 1), seq(21, 50, 3), seq(51, 99, 6), seq(100, (NIDs - 1), 12))}

  DoubleLoop <- data.frame(SampleSize = rep(Nloop, each=Iteration), Iteration=rep(seq(1:Iteration), length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])
  UDLev <- 50

  ### CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)
  minX <- min(coordinates(TripCoords)[,1]) - Scale*2000
  maxX <- max(coordinates(TripCoords)[,1]) + Scale*2000
  minY <- min(coordinates(TripCoords)[,2]) - Scale*2000
  maxY <- max(coordinates(TripCoords)[,2]) + Scale*2000

  ### if users do not provide a resolution, then split data into ~500 cells
  if(Res > 99){Res <- (max(abs(minX - maxX)/500,
                         abs(minY - maxY)/500))/1000
  warning(sprintf("No grid resolution ('Res') was specified, or the specified resolution was >99 km and therefore ignored. Space use was calculated on a 500-cell grid, with cells of  %s square km.", round(Res,3)),immediate. = TRUE)}


  ### specify sequence of grid cells and combine to SpatialPixels
  xrange <- seq(minX, maxX, by = Res*1000)
  yrange <- seq(minY, maxY, by = Res*1000)
  grid.locs <- expand.grid(x = xrange, y = yrange)
  INPUTgrid <- SpatialPixels(SpatialPoints(grid.locs), proj4string = proj4string(TripCoords))

  #### ERROR CATCH IF PEOPLE SPECIFIED TOO FINE RESOLUTION ####
  if (max(length(xrange), length(yrange)) > 600){warning("Your grid has a pretty large number of cells - this will slow down computation. Reduce 'Res' to speed up the computation.")}
  if (max(length(xrange), length(yrange)) > 1200){stop("Are you sure you want to run this function at this high spatial resolution ('Res')? Your grid is >1 million pixels, computation will take many hours (or days)!")}



  #########~~~~~~~~~~~~~~~~~~~~~~~~~#########
  ### PARALLEL LOOP OVER AL ITERATIONS ######
  #########~~~~~~~~~~~~~~~~~~~~~~~~~#########
  #before <- Sys.time()
  Ncores<-ifelse(Ncores==1, parallel::detectCores()/2, Ncores) ## use user-specified value if provided to avoid computer crashes by using only half the available cores
  cl <- parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl)
  Result <- data.frame()

  Result <- foreach::foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp","adehabitatHR","dplyr")) %dopar% {

    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]

    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)

    RanNum <- sample(UIDs, N, replace=F)
    NotSelected <- TripCoords[!TripCoords$ID %in% RanNum,]
    Selected <- TripCoords[TripCoords$ID %in% RanNum,]
    Selected <- as(Selected, 'SpatialPoints')

    ##### Calculate Kernel

    KDE.Surface <- adehabitatHR::kernelUD(Selected, h=(Scale * 1000), grid=INPUTgrid)

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ### Calculating inclusion value, using Kernel surface ######
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    KDEpix <- as(KDE.Surface, "SpatialPixelsDataFrame")
    pixArea <- KDE.Surface@grid@cellsize[1]

    KDEpix@data <- KDEpix@data %>%
      rename(UD = .data$ud) %>%
      mutate(rowname = 1:nrow(KDEpix@data)) %>%
      mutate(usage = .data$UD * (pixArea^2)) %>%
      arrange(desc(.data$usage)) %>%
      mutate(cumulUD = cumsum(.data$usage)) %>%
      mutate(INSIDE = ifelse(.data$cumulUD < 0.5, 1, NA)) %>%
      arrange(.data$rowname) %>%
      dplyr::select(.data$INSIDE)

    ########

    Overlain <- over(NotSelected, KDEpix)
    Output$InclusionMean <- length(which(!is.na(Overlain$INSIDE)))/nrow(NotSelected)

    return(Output)
    }
  ## stop the cluster
  on.exit(parallel::stopCluster(cl))


  if(BootTable==T){
    data.table::fwrite(Result,"bootout_temp.csv", row.names=F, sep=",")
  }

  try(M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
  if ('M1' %in% ls()){       ### run this only if nls was successful
    Asymptote <- (base::summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
    Result$pred <- stats::predict(M1)

    ## Calculate RepresentativeValue
    RepresentativeValue <- Result %>%
      group_by(.data$SampleSize) %>%
      summarise(out = max(.data$pred) / ifelse(Asymptote < 0.45, 0.5, Asymptote)*100) %>%
      dplyr::filter(out == max(.data$out)) %>%
      mutate(type = ifelse(Asymptote < 0.45 | Asymptote > 0.6, 'asymptote_adj', 'asymptote')) %>%
      mutate(asym = Asymptote)

  if(Asymptote < 0.45 | Asymptote > 0.6) {
    RepresentativeValue$asym_adj <- 0.5 }

    ## Plot
    P2 <- Result %>%
      group_by(.data$SampleSize) %>%
      dplyr::summarise(
        meanPred = mean(na.omit(.data$pred)),
        sdInclude = sd(.data$InclusionMean))

    yTemp <- c(P2$meanPred + 0.5 * P2$sdInclude, rev(P2$meanPred - 0.5 * P2$sdInclude))
    xTemp <- c(P2$SampleSize, rev(P2$SampleSize))

    pdf("track2kba_repAssess_output.pdf", width=6, height=5)  ## avoids the plotting margins error
    plot(InclusionMean ~ SampleSize,
      data = Result, pch = 16, cex = 0.2, col="darkgray", ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))), ylab = "Inclusion", xlab = "SampleSize")
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
    lines(P2, lty=1,lwd=2)
    text(x=0, y=1, paste(round(RepresentativeValue$out, 2), "%", sep=""), cex=2, col="gray45", adj=0)
    dev.off()

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

