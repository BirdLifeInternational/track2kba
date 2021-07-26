#' Assess sample representativeness
#'
#' \code{repAssess} estimates the degree to which the space use of a tracked 
#' sample of animals represents that of the larger population. 
#'
#' Representativeness is assessed by fitting statistical model to the 
#' relationship between sample size and inclusion rate. Incusion rate is the 
#' proportion of out-sample points included in in-sample space use areas.
#'
#' First, the set of IDs is iteratively sub-sampled, and in each iteration a set
#'  of individual Utilization Distributions (UD, 'KDE' argument) are pooled and 
#'  the points of the un-selected (out-sample) IDs are overlaid on the % contour
#'   area ('levelUD') of the UD. The proportion of these outsample points which 
#'  overlap the pooled UD area is known as the inclusion rate, and represents 
#'  an estimate of representativeness at each sample size. Then, a non-linear 
#'  function is fit to the relationship between the inclusion rate and sample 
#'  size (i.e. number of tracks/animals) in order to estimate the point at which
#'   the relationship reaches an asymptote (i.e. no more information added per 
#'   new track). \code{repAssess} then estimates the representativeness of the 
#'   sample by dividing the inclusion rate estimated at the maximum sample size 
#'   minus 3 (for samples where n < 20), 2 (for samples < 50) or 1 (for sample 
#'   >100) by this asymptote. The maximum sample size appearing in the plot
#'   will be different than the true 'n' of the dataset in order to account for 
#'   the possible number of combinations of individuals, thereby ensuring a 
#'   robust result. The maximum sample size reflects the number of KDEs, so if 
#'   any ID has fewer than 5 points, this ID is omitted from the analysis. 
#'   Finally, using this relationship, minimum representative sample sizes 
#'   (70% and 95%) are also calculated.
#' 
#' \code{\link{repAssess}} accepts UDs calculated outside of \code{track2KBA}, 
#' if they have been converted to class \code{RasterStack} or 
#' \code{SpatialPixelsDataFrame}. However, one must make sure that the cell 
#' values represent continuous probability densities (i.e. values >=0 which 
#' integrate to 1 over the raster) and not not discrete probability masses 
#' (i.e. values >=0 which sum to 1), nor home range quantiles (i.e. 0-1, or 
#' 0-100 representing % probability of occurrence).
#' 
#' When setting \code{avgMethod} care must be taken. If the number of points
#' differ greatly among individuals and the UDs are calculated as classic 
#' KDEs (e.g. from \code{\link{estSpaceUse}}) then the weighted mean is likely 
#' the optimal way to pool individual UDs. However, if any other method (for 
#' example AKDE, auto-correlated KDE) was used to estimate UDs, then the 
#' arithmetic mean is the safer option. 
#' 
#' NOTE: this function does not work with fewer than 4 IDs (tracks or 
#' individual animals).
#'
#' @param tracks SpatialPointsDataFrame of spatially projected animal 
#' relocations. Must include 'ID' field.
#' @param KDE Kernel Density Estimates for individual animals. Several input 
#' options: an estUDm, a SpatialPixels/GridDataFrame, or a RasterStack. 
#' If estUDm, must be as created by \code{\link{estSpaceUse}} or 
#' \code{adehabitatHR::kernelUD}, if Spatial* each column should correspond to 
#' the Utilization Distribution of a single individual or track. If a 
#' RasterStack, each layer must be an individual UD. 
#' @param iteration numeric. Number of times to repeat sub-sampling procedure. 
#' The higher the iterations, the more robust the result. 
#' @param levelUD numeric. Specify which contour of the utilization distribution
#'  (\code{KDE}) you wish to filter to (e.g. core area=50, home range=95). 
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
#' returned, with columns '\emph{SampleSize}' signifying the sample size (i.e.,
#' number of KDEs)'\emph{out}' signifying the percent representativeness of the 
#' sample,'\emph{type}' is the type  of asymptote value used to calculate the 
#'  '\emph{out}' value, and '\emph{asym}' is the asymptote value used. 
#'  If \code{bootTable=TRUE}, a list returned with above dataframe in first slot
#'   and full iteration results in second slot.
#'
#' There are two potential values for '\emph{type}':'asymptote' is 
#' the ideal, where the asymptote value is calculated from the parameter 
#' estimates of the successful nls model fit. 'inclusion' is used if the nls 
#' fails to converge, or if the fit model is flipped and the asymptote value 
#' is negative.  In these casess, the mean inclusion rate is taken for the 
#' largest sample size.'Rep70' signifies the sample size which is ~70% 
#' representative, and 'Rep95' signifies the sample 
#'  size which approahces the asymptote.
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
  tracks, KDE=NULL, iteration=1, levelUD, avgMethod = "mean", nCores=1, 
  bootTable=FALSE){
  
  if(!(levelUD >= 1 & levelUD <= 100)) {stop("levelUD must be between 1-100%")}
  
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
    # convert IDs to 'valid' raster names 
    newlvls <- raster::validNames(unique(tracks$ID))
    IDs2keep <- unique(tracks$ID)[newlvls %in% UDnames]
    tracks <- tracks[tracks$ID %in% IDs2keep, ]
    
  }
  
  UIDs <- unique(tracks$ID) # get IDs, after any filtered out 
  NIDs <- length(UIDs)
  
  Nloop <- seq(1, (NIDs - 3), 1)
  if(NIDs >= 50 & NIDs < 100){ # for big samples, skip some sample sizes 
    Nloop <- c(seq(1, 25, 1), seq(26, (NIDs - 2), 3))
    }
  if(NIDs >= 100){
    Nloop <- c(seq(1, 25, 1), seq(26, 49, 3), seq(50, (NIDs - 1), 6))
    }

  DoubleLoop <- data.frame(
    SampleSize = rep(Nloop, each=iteration), 
    iteration=rep(seq_len(iteration), length(Nloop))
    )
  LoopNr <- seq_len(dim(DoubleLoop)[1])	
  
  ### Parallel or sequential processing?----------------------------------------
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
    
    Output <- data.frame(SampleSize = N, InclusionRate = 0,iteration=i)
    
    RanNum <- sample(UIDs, N, replace=FALSE)
    NotSelected <- tracks[!tracks$ID %in% RanNum,]
    SelectedTracks <- tracks[tracks$ID %in% RanNum,]
    
      Selected <- KDEraster[[raster::validNames(RanNum)]]

    KDEstack <- raster::stack(Selected) # list of RasterLayers to RasterStack

    if(avgMethod=="weighted"){
      # weighted mean - number of points per ID -------------------------------
      weights <- as.vector(unname(table(SelectedTracks$ID)))   
      KDEcmbnd <- raster::weighted.mean(KDEstack, w=weights)   
    } else if( avgMethod == "mean"){
      KDEcmbnd <- raster::calc(KDEstack, mean) # arithmetic mean
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
    
    Output$InclusionRate <- length(
      which(!is.na(Overlain_Raster)))/nrow(NotSelected
        )
    
    return(Output)
  }
  
  if(nCores > 1) {on.exit(parallel::stopCluster(cl))} # stop the cluster

  ### if nls is unsuccessful use mean output for largest sample size --------
  RepOutput <- Result %>%
    dplyr::filter(.data$SampleSize == max(.data$SampleSize)) %>%
    group_by(.data$SampleSize) %>%
    summarise(
      out = mean(.data$InclusionRate)
    ) %>%
    mutate(type = 'inclusion') %>%
    mutate(asym = .data$out,
           SampleSize = NIDs)
  
  startparsset <- seq(0.1, 1, 0.1) # set of starting parameters to try
  fit <- FALSE
  i <- 1
  tries <- 100                     # # of times to try re-fitting the model
  while( (fit==FALSE) & (i < tries) ){
    tryCatch({
      startpars <- sample(startparsset, 2)
      M1 <- stats::nls(
        Result$InclusionRate ~ (
          (a * Result$SampleSize)/(1 + b * Result$SampleSize)), 
        data = Result, start = list(a = startpars[1], b = startpars[2])
      )
      fit <- TRUE
      # return(M1)
    }, error=function(e){})
    i <- i + 1
  }
  
  if( exists("M1") ) { # run only if nls was successful

    a <- base::summary(M1)$coefficients[1]
    b <- base::summary(M1)$coefficients[2]
    
    tAsymp <- (levelUD/100)
    Asymptote <- (a / b)
    
    if( Asymptote > 0 ){ # derive rep. if asymptote is positive ----------------
      message("nls (non linear regression) successful, asymptote estimated for 
    bootstrap sample.")
      Result$pred <- stats::predict(M1)
      
      ## Calculate Representativeness value ------------------------------------
      RepOutput <- Result %>%
        group_by(.data$SampleSize) %>%
        summarise(
          out      = max(.data$pred) / Asymptote*100
        ) %>%
        dplyr::filter(.data$out == max(.data$out)) %>%
        mutate(
          SampleSize = NIDs,     
          est_asym = Asymptote,
          tar_asym = (levelUD/100)
        ) 
      
      # calculate 70% and 95% representative sample sizes ----------------------
      Rep70  <- 0.7*Asymptote
      Rep95 <- 0.95*Asymptote
      # convert inclusions into rep. estimates ---------------------------------
      RepOutput$Rep70  <- ceiling(Rep70 / (a - (Rep70*b)))
      RepOutput$Rep95 <- ceiling(Rep95 / (a - (Rep95*b)))
      
      Result <- Result %>% mutate(
        rep_est = .data$pred / Asymptote*100, 
        is_rep  = ifelse(.data$rep_est >= 70, TRUE, FALSE)
      )
      
      ### Plot -----------------------------------------------------------------
      P2 <- Result %>%
        group_by(.data$SampleSize) %>%
        dplyr::summarise(
          meanPred = mean(na.omit(.data$pred)),
          sdInclude = sd(.data$InclusionRate))
      
      yTemp <- c(
        P2$meanPred + Asymptote * P2$sdInclude, 
        rev(P2$meanPred - Asymptote * P2$sdInclude)
      )
      xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
      
      plot(InclusionRate ~ SampleSize,
           data = Result, pch = 16, 
           cex = 0.2, col="darkgray", 
           ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))),
           ylab = "Inclusion", xlab = "Sample Size"
      )
      polygon(x = xTemp, y = yTemp, col = "gray93", border = FALSE)
      points(InclusionRate ~ SampleSize, 
             data=Result, pch=16, cex=0.2, col="darkgray"
      )
      lines(P2, lty=1,lwd=2)
      text(
        x=0, y=0.99, 
        paste(round(RepOutput$out, 1), "%", sep=""), 
        cex=2, col="gray45", adj=0
      )
    } else { # if asymptote is negative 
      message("Model fit was poor, resulting in negative asymptote. Likely due 
      to small sample or few iterations. 'out' derived from mean inclusion value
       at highest sample size.")
      RepOutput <- RepOutput %>% 
        mutate(
          est_asym = Asymptote, tar_asym = tAsymp
        )
    }
    
    if( abs(Asymptote - tAsymp) > 0.1 & Asymptote > 0 ) { 
      message("Estimated asymptote differs from target; be aware that 
  representativeness value is based on estimated asymptote (i.e. est_asym).") }
    
  } else { # if nls fails 
    message("nls (non linear regression) unsuccessful, likely due to small 
        sample or few iterations. Data may not be representative, 'out' derived 
        from mean inclusion value at highest sample size.")
  }
  
  if(bootTable==TRUE){
    listedResults <- list(as.data.frame(RepOutput), Result)
    return(listedResults)
  } else {
    return(as.data.frame(RepOutput))
  }
  
}
