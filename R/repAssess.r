#' Assess sample representativeness
#'
#' \code{repAssess} estimates the degree to which the space use of a tracked 
#' sample of animals represents that of the larger population. 
#'
#' Representativeness is assessed by calculating the proportion of out-sample 
#' points included in in-sample space use areas.
#'
#' First, the set of IDs is iteratively sub-sampled, and in each iteration a set
#'  of individual UDs are pooled and the points of the un-selected (outsample) 
#'  IDs are overlaid on the % contour area ('levelUD') of the KDE. The 
#'  proportion of these outsample points which overlap the pooled KDE area is 
#'  known as the inclusion rate, and represents an estimate of 
#'  representativeness at each sample size. Then, a non-linear function is fit 
#'  to the relationship between the inclusion rate and sample size (i.e. number 
#'  of tracks/animals) in order to estimate the point at which the relationship 
#'  reaches an information asymptote (i.e. no more information added per track).
#'   By estimating this asymptote, \code{repAssess} provides an estimate of 
#'   representativeness by dividing the inclusion rate estimated at the maximum 
#'   sample size, by this asymptote. Finally, using this relationship, a minimum
#'    representative sample sizes are also calculated.
#' 
#' \code{\link{repAssess}} accepts UDs calculated outside of \code{track2KBA}, 
#' if they have been converted to class \code{RasterStack} or 
#' \code{SpatialPixelsDataFrame}. However, one must make sure that the cell 
#' values in these cases are probability densities (i.e. <=0) not volumne 
#' contour values (i.e. 0-1, or 0-100).
#' 
#' When setting \code{avgMethod} care must be taken. If the input are classic 
#' KDEs (e.g. from \code{\link{estSpaceUse}}) then the weighted mean is likely 
#' the optimal way to pool individual UDs. However, if any other method (for 
#' example AKDE, auto-correlated KDE) was used to estimate UDs, then the 
#' arithmetic mean is the safer option. 
#'
#' @param tracks SpatialPointsDataFrame or data.frame of animal relocations. 
#' Must include 'ID' field. If input is data.frame or unprojected 
#' SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields.
#' @param KDE Kernel Density Estimates for individual animals. Several input 
#' options: an estUDm, a SpatialPixels/GridDataFrame, or a RasterStack. 
#' If estUDm, must be as created by \code{\link{estSpaceUse}} or 
#' \code{adehabitatHR::kernelUD}, if Spatial* each column should correspond to 
#' the Utilization Distribution of a single individual or track. If a 
#' RasterStack, each layer must be an individual UD. 
#' @param iteration numeric. Number of times to repeat sub-sampling procedure. 
#' The higher the iterations, the more robust the result. 
#' @param res numeric. Grid cell resolution (in square kilometers) for kernel 
#' density estimation. Default is a grid of 500 cells, with spatial extent 
#' determined by the latitudinal and longitudinal extent of the data. Only needs
#'  to be set if nothing is supplied to \code{KDE}.
#' @param levelUD numeric. Specify which contour of the Utilization Distribution
#'  the home range estimate (\code{KDE}) represented (e.g. 50, 95). 
#' @param avgMethod character. Choose whether to use the arithmetic or weighted 
#' mean when combining individual IDs. Options are :'mean' arithmetic mean, or 
#' 'weighted', which weights each UD by the numner of points per level of ID.
#' @param bootTable logical (TRUE/FALSE). If TRUE, output is a list, containing 
#' in the first slot the representativeness results summarized in a table, and 
#' in the second the full results of the iterated inclusion calculations.
#' @param nCores numeric. The number of processing cores to use. For heavy 
#' operations, the higher the faster. NOTE: CRAN sets a maximum at 2 cores. If 
#' using the git-hub version of the package, this can be set to a maximum of one
#'  fewer than the maximum cores in your computer.
#'  
#' @return if \code{bootTable=FALSE} (the default) A single-row data.frame is 
#' returned, with columns '\emph{SampleSize}' signifying the maximum sample size
#'  in the data set, '\emph{out}' signifying the percent 
#' representativeness of the sample, '\emph{type}' is the type  of asymptote 
#' value used to calculate the '\emph{out}' value, and '\emph{asym}' is the 
#' asymptote value used. if \code{bootTable=TRUE}, a list returned with above
#' dataframe in first slot and full iteration results in second slot.
#'
#' There are two potential values for '\emph{type}': 'type' == 'asymptote' is 
#' the ideal, where the asymptote value is calculated from the parameter 
#' estimates of the successful nls model fit. Secondly, when nls fails to 
#' converge at all, then the mean inclusion rate is taken for the largest sample
#'  size; 'type'=='inclusion.'Rep70' signifies the sample size at which ~70% of 
#'  the space use information is encompassed, and 'Rep95' signifies the sample 
#'  size approaching the asymptote, i.e. representing 99% of the space use.
#'
#' @examples
#' \dontrun{repr <- repAssess(
#' Trips, KDE=KDE, iteration=1, bootTable = FALSE, nCores = 1)
#' }
#'
#' @export
#' @importFrom foreach %dopar%
#' @importFrom graphics abline identify lines points polygon text

repAssess <- function(
  tracks, KDE=NULL, iteration=50, res=NULL, levelUD=50, 
  avgMethod = "mean", nCores=1, bootTable=FALSE){
  
  tracks@data <- tracks@data %>% dplyr::select(.data$ID)
  
  UIDs <- unique(tracks$ID)
  
  # Determine class of KDE, and convert to Raster -----------------------------
  if(is.null(KDE)){
    stop("No Utilization Distributions supplied to KDE argument. 
      See track2KBA::estSpaceUse()")
  } else if(class(KDE) == "estUDm"){ 
    KDEraster <- raster::stack(lapply(KDE, function(x) {
      raster::raster(as(x, "SpatialPixelsDataFrame"), values=TRUE)
      } ))
  } else if(class(KDE) == "SpatialPixelsDataFrame") {
    KDEraster <- raster::stack(KDE)
  } else if(class(KDE) == "RasterStack") {
    KDEraster <- KDE
  }
  
  # assure only IDs with UDs are in tracking data -----------------------------
  UDnames <- names(KDEraster) 
  if( length(UDnames) != length(UIDs) ) {
    if( all(stringr::str_detect(names(KDEraster), pattern = "^X")) ){
      UDnames <- substring(UDnames, 2) # remove "X" from raster names
      tracks <- tracks[tracks$ID %in% UDnames, ]
    } else {
      tracks <- tracks[tracks$ID %in% UDnames, ]
    }
  }
  
  UIDs <- unique(tracks$ID) # get IDs, after any filtered out 
  NIDs <- length(UIDs)

  if(NIDs<50){Nloop <- seq(1, (NIDs - 1), 1)}
  if(NIDs>=50 & NIDs<100){Nloop <- c(seq(1, 19, 1), seq(20, (NIDs - 1), 3))}
  if(NIDs>=100){Nloop <- c(seq(1, 20, 1), seq(21, 49, 3), seq(50, (NIDs - 1), 6))}

  
  DoubleLoop <- data.frame(
    SampleSize = rep(Nloop, each=iteration), 
    iteration=rep(seq_len(iteration), length(Nloop))
    )
  LoopNr <- seq_len(dim(DoubleLoop)[1])	
  
  ###
  if(nCores > 1){
    if (!requireNamespace("parallel", quietly = TRUE) | 
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages \"parallel\" and \"doParallel\" needed for this function 
        to work. Please install.",
        call. = FALSE)  }
    maxCores <- parallel::detectCores()
    # ensure that at least one core is un-used 
    nCores <- ifelse(nCores == maxCores, nCores - 1, nCores)
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
  } else { foreach::registerDoSEQ() }

  Result <- data.frame()

  LoopN <- NULL # make R CMD Check happy
  Result <- foreach::foreach(
    LoopN = LoopNr, .combine = rbind, .packages = c("sp", "dplyr", "raster")
    ) %dopar% {
    
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$iteration[LoopN]
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,iteration=i)
    
    RanNum <- sample(UIDs, N, replace=FALSE)
    NotSelected <- tracks[!tracks$ID %in% RanNum,]
    SelectedTracks <- tracks[tracks$ID %in% RanNum,]
    
    # if ID lvls start with number, add X for indexing ------------------------
    if( all(stringr::str_detect(unique(tracks$ID), pattern = "^[0-9]")) ){ 
      Selected <- KDEraster[[paste("X", RanNum, sep = "")]]
    } else {
      Selected <- KDEraster[[ RanNum ]]
    }

    KDEstack <- raster::stack(Selected) # list of RasterLayers to RasterStack

    if(avgMethod=="weighted"){
      # weighted mean - number of points per ID -------------------------------
      weights <- as.vector(unname(table(SelectedTracks$ID)))   
      KDEcmbnd <- raster::weighted.mean(KDEstack, w=weights)   
    } else if( avgMethod == "mean"){
      KDEcmbnd <- raster::calc(KDEstack, mean)                 # arithmetic mean
    }
    
    ### Calculating inclusion value, using Kernel surface ---------------------
    KDElev <- KDEcmbnd
    pixArea <- raster::res(KDElev)[1]
    
    df <- data.frame(UD = raster::getValues(KDElev)) %>%
      mutate(rowname = seq_len(length(raster::getValues(KDElev)))) %>%
      mutate(usage = .data$UD * (pixArea^2)) %>%
      arrange(desc(.data$usage)) %>%
      mutate(cumulUD = cumsum(.data$usage)) %>%
      mutate(INSIDE = ifelse(.data$cumulUD < (levelUD/100), 1, NA)) %>%
      arrange(.data$rowname) %>%
      dplyr::select(.data$INSIDE)
    
    KDElev[] <- df$INSIDE
    
    Overlain_Raster <- raster::extract(KDElev, NotSelected)
    
    Output$InclusionMean <- length(
      which(!is.na(Overlain_Raster)))/nrow(NotSelected
        )
    
    return(Output)
  }
  
  if(nCores > 1) {on.exit(parallel::stopCluster(cl))} # stop the cluster

  try(M1 <- stats::nls(
    Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize), 
    data = Result, start = list(a=1,b=0.1)
    ), silent = TRUE)
  if ('M1' %in% ls()){       # run only if nls was successful
    a <- base::summary(M1)$coefficients[1]
    b <- base::summary(M1)$coefficients[2]
    
    tAsymp <- (levelUD/100)
    Asymptote <- (a / b)
    Result$pred <- stats::predict(M1)
    
    ## Calculate Representativeness value -------------------------------------
    RepresentativeValue <- Result %>%
      group_by(.data$SampleSize) %>%
      summarise(
        out      = max(.data$pred) / Asymptote*100
        ) %>%
      dplyr::filter(.data$out == max(.data$out)) %>%
      mutate(
        est_asym = Asymptote,
        tar_asym = (levelUD/100)
        ) 
    
    # calculate 70% and 95% representative sample sizes -----------------------
    Rep70  <- 0.7*Asymptote
    Rep95 <- 0.95*Asymptote
    # convert inclusions into rep. estimates ----------------------------------
    RepresentativeValue$Rep70  <- ceiling(Rep70 / (a - (Rep70*b)))
    RepresentativeValue$Rep95 <- ceiling(Rep95 / (a - (Rep95*b)))
    
    Result <- Result %>% mutate(
      rep_est = .data$pred / Asymptote*100, 
      is_rep  = ifelse(.data$rep_est >= 70, TRUE, FALSE)
    )

    ### Plot ------------------------------------------------------------------
    P2 <- Result %>%
      group_by(.data$SampleSize) %>%
      dplyr::summarise(
        meanPred = mean(na.omit(.data$pred)),
        sdInclude = sd(.data$InclusionMean))
    
    yTemp <- c(
      P2$meanPred + Asymptote * P2$sdInclude, 
      rev(P2$meanPred - Asymptote * P2$sdInclude)
      )
    xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
    
    plot(InclusionMean ~ SampleSize,
      data = Result, pch = 16, 
      cex = 0.2, col="darkgray", 
      ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))),
      ylab = "Inclusion", xlab = "Sample Size"
      )
    polygon(x = xTemp, y = yTemp, col = "gray93", border = FALSE)
    points(InclusionMean ~ SampleSize, 
      data=Result, pch=16, cex=0.2, col="darkgray"
      )
    lines(P2, lty=1,lwd=2)
    text(
      x=0, y=0.99, 
      paste(round(RepresentativeValue$out, 2), "%", sep=""), 
      cex=2, col="gray45", adj=0
      )

  }else{
    ### if nls is unsuccessful use mean output for largest sample size --------
    RepresentativeValue <- Result %>%
      dplyr::filter(.data$SampleSize == max(.data$SampleSize)) %>%
      group_by(.data$SampleSize) %>%
      summarise(
        out = mean(.data$InclusionMean)
        ) %>%
      mutate(type = 'inclusion') %>%
      mutate(asym = .data$out)
  }
  
  
  if(exists("M1")) {
  message("nls (non linear regression) successful, asymptote estimated for 
    bootstrap sample.")
  } else {  
      warning("nls (non linear regression) unsuccessful, likely due to small 
        sample or iterations. Data may not be representative, output derived 
        from mean inclusion value at highest sample size.") }

  if(Asymptote < (tAsymp - 0.1) ) { warning("Estimated asymptote differs from 
    target; be aware that representativeness value is based on 
    estimated asymptote (i.e. est_asym).") }
  
  if(bootTable==TRUE){
    listedResults <- list(as.data.frame(RepresentativeValue), Result)
    return(listedResults)
  } else {
    return(as.data.frame(RepresentativeValue))
  }
  
}


