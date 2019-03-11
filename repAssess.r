## bootstrap     ######################################################################################################

## MAIN UPDATE: tidyverse, simple features
## REVISION: use pre-determined H value for species based on Oppel et al. 2018
## explore whether sequence and number of iterations can be reduced to increase speed
## 1-20 at increments of 1, 20-50 at increments of 3, 50-100 at increments of 5, 100-150 at increments of 10, 150-200 at increments of 25, >200 at increments of 50
## max n iterations to n of possible combinations in data - no need to do 100 iterations if only 20 combinations possible
## explore alternative approach of increasing area of 50%UD (may not be much faster though)

## (Based on original by Phil Taylor & Mark Miller, 2012)

#### DESCRIPTION: ##
## This script iteratively sub-samples a dataset of tracking data, investigating the effect of sample size. 
## This is done by estimating the degree to which the space use of the tracked sample of animals is representative of the population's  
## space use. 
## At each iteration the data is split, one half is used as the 'training' data and the 50%UD is calculated from this. The second half is 
## used as 'testing' data and the proportion of points captured within the 50%UD is calculated.
## A perfect dataset would tend towards 0.5. By fitting a trend line to this relationship we can identify the sample size at which the curve
## approaches an asymptote, signifying that any new data would simply add to existing knowledge. This script produces a 
##representativeness value, indicating how close to this point the sample is. 

#### ARGUMENTS: ##
## DataGroup must be a dataframe or SpatialPointsDataFrame with Latitude, Longitude and ID as fields.
## Scale determines the smoothing factor ('h' parameter) used in the kernel analysis.
## Iteration determines the number of times each sample size is iterated.
## Res sets the resolution of grid cells used in kernel analysis (sq. km)

## REVISED BY Steffen Oppel in 2015 to facilitate parallel processing
## updated to adehabitatHR by Steffen Oppel on 27 Dec 2016
## changed to same4all=TRUE on 4 Feb 2017

## REVISED in 2017 to avoid error in nls function of singular gradient
## added mean output for inclusion value even if nls fails

bootstrap <- function(DataGroup, Scale=100, Iteration=50, Res=100, BootTable=T)
{
  
  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)   #### NEED TO FIX
  require(foreach)
  require(doParallel)
  require(parallel)
  
  if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    ## set the minimum fields that are needed
    CleanDataGroup <- DataGroup %>%
      dplyr::select(ID, Latitude, Longitude,DateTime) %>%
      arrange(ID, DateTime)
    mid_point<-data.frame(centroid(cbind(CleanDataGroup$Longitude, CleanDataGroup$Latitude)))
    
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
  
  proj.UTM <- CRS(proj4string(TripCoords))
  
  TripCoords$X <- TripCoords@coords[, 1]
  TripCoords$Y <- TripCoords@coords[, 2]
  BoundBox <- bbox(TripCoords)
  UIDs <- unique(TripCoords$ID)
  NIDs <- length(UIDs)
  Nloop <- seq(1, (NIDs - 1), ifelse(NIDs > 100, 10, 1)) ## change sequence here? (i.e. 1-20  by 1, 20-50 by 3 etc.)
  DoubleLoop <- data.frame(SampleSize = rep(Nloop, each=Iteration), Iteration=rep(seq(1:Iteration), length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- 50
  
  #setup parallel backend to use 4 processors
  before <- Sys.time()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  Result <- data.frame()
  
  Result <- foreach(LoopN = LoopNr, .combine = rbind, .packages = c("sp","adehabitatHR","geosphere","rgdal", "dplyr")) %dopar% {
    
    N <- DoubleLoop$SampleSize[LoopN]
    i <- DoubleLoop$Iteration[LoopN]
    Coverage <- NULL
    Inclusion <- NULL
    History <- NULL
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)
    
    RanNum <- sample(UIDs, N, replace=F)
    SelectedCoords <- coordinates(TripCoords[TripCoords$ID %in% RanNum,])
    NotSelected <- TripCoords[!TripCoords$ID %in% RanNum,]
    Temp <- data.frame(SelectedCoords[,1], SelectedCoords[,2])
    Temp <- SpatialPoints(Temp, proj4string=proj.UTM)      ### added because adehabitatHR requires SpatialPoints object
    
    
    ### CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)
    ### NEED TO DO: link resolution of grid to H-parameter ('Scale')
    
    minX<-min(coordinates(TripCoords)[,1]) - Scale*2000
    maxX<-max(coordinates(TripCoords)[,1]) + Scale*2000
    minY<-min(coordinates(TripCoords)[,2]) - Scale*2000
    maxY<-max(coordinates(TripCoords)[,2]) + Scale*2000
    
    ### if users do not provide a resolution, then split data into ~500 cells
    if(Res>99){Res <- (max(abs(minX-maxX)/500,
      abs(minY-maxY)/500))/1000
    warning(sprintf("No grid resolution ('Res') was specified, or the specified resolution was >99 km and therefore ignored.
      Space use was calculated in square grid cells of %s km", round(Res,3)))}
    
    ### specify sequence of grid cells and combine to SpatialPixels
    xrange<-seq(minX,maxX, by = Res*1000) #diff(range(coordinates(TripCoords)[,1]))/Res)   ### if Res should be provided in km we need to change this
    yrange<-seq(minY,maxY, by = Res*1000) #diff(range(coordinates(TripCoords)[,2]))/Res)   ### if Res should be provided in km we need to change this
    grid.locs<-expand.grid(x=xrange,y=yrange)
    INPUTgrid<-SpatialPixels(SpatialPoints(grid.locs), proj4string=proj4string(TripCoords))
    #  plot(INPUTgrid)
    
    #### ERROR CATCH IF PEOPLE SPECIFIED TOO FINE RESOLUTION ####
    if (max(length(xrange),length(yrange))>600){warning("Your grid has a pretty large number of cells - this will slow down computation. Reduce 'Res' to speed up the computation.")}
    if (max(length(xrange),length(yrange))>1200){stop("Are you sure you want to run this function at this high spatial resolution ('Res')? Your grid is >1 million pixels, computation will take many hours (or days)!")}
    
    
    ##### Calculate Kernel
    
    KDE.Surface <- adehabitatHR::kernelUD(Temp, h=(Scale * 1000), grid=INPUTgrid, same4all=F)		## newer version needs SpatialPoints object and id no longer required in adehabitatHR
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ### Calculating inclusion value, using Kernel surface ######
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    
    KDEpix <- as(KDE.Surface, "SpatialPixelsDataFrame")
    if(is.projected(KDEpix)!=TRUE) stop("Please re-calculate your kernel UD after projecting the data into a coordinate reference system where units are identical on x- and y-axis")
    
    pixArea <- KDE.Surface@grid@cellsize[1]
    
    UDLevCells <- KDEpix
    UDLevCells@data <- KDEpix@data %>% 
      rename(UD = ud) %>% 
      mutate(rowname=1:nrow(KDEpix@data)) %>%
      mutate(usage=UD*(pixArea^2)) %>%
      arrange(desc(usage)) %>%
      mutate(cumulUD = cumsum(usage)) %>%
      mutate(INSIDE = ifelse(cumulUD < 0.5, 1, NA)) %>%
      arrange(rowname) %>%
      dplyr::select(INSIDE) 
    
    
    ########
    
    Overlain <- over(NotSelected, UDLevCells)
    
    Output$InclusionMean <- length(which(!is.na(Overlain$INSIDE)))/nrow(NotSelected)
    
    return(Output)
    }
  ## stop the cluster
  on.exit(stopCluster(cl))
  Sys.time() - before
  
  par(mfrow=c(1,1), mai=c(1,1,1,1))
  #Result <- Output[1:nrow(Output)-1,]
  
  if(BootTable==T){
    write.table(Result,"bootout_temp.csv", row.names=F, sep=",")
  }
  
  try(M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
  if ('M1' %in% ls()){       ### run this only if nls was successful
    Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
    Result$pred <- predict(M1)
    
    ## Calculate RepresentativeValue 
    RepresentativeValue <- Result %>%
      group_by(SampleSize) %>%
      summarise(out = max(pred) / ifelse(Asymptote < 0.45, 0.5, Asymptote)*100) %>%
      filter(out == max(out)) %>%
      mutate(type = ifelse(Asymptote < 0.45, 'asymptote_adj', 'asymptote')) %>%
      mutate(asym = Asymptote)
    
    ## Plot
    P2 <- Result %>% 
      group_by(SampleSize) %>% 
      dplyr::summarise(
        meanPred = mean(na.omit(pred)),
        sdInclude = sd(InclusionMean))
    yTemp <- c(P2$meanPred + 0.5 * P2$sdInclude, rev(P2$meanPred - 0.5 * P2$sdInclude))
    xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
    
    plot(InclusionMean ~ SampleSize, 
      data = Result, pch = 16, cex = 0.2, col="darkgray", ylim = c(0,1), xlim = c(0,max(unique(Result$SampleSize))), ylab = "Inclusion", xlab = "SampleSize")
    polygon(x = xTemp, y = yTemp, col = "gray93", border = F)
    points(InclusionMean ~ SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
    lines(P2, lty=1,lwd=2)
    text(x=0, y=1,paste(round(RepresentativeValue$out, 2), "%", sep=""), cex=2, col="gray45", adj=0)  
    
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

