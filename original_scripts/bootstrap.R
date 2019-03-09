## bootstrap     ######################################################################################################

## MAIN UPDATE: tidyverse, simple features
## REVISION: use pre-determined H value for species based on Oppel et al. 2018
## explore whether sequence and number of iterations can be reduced to increase speed
## 1-20 at increments of 1, 20-50 at increments of 3, 50-100 at increments of 5, 100-150 at increments of 10, 150-200 at increments of 25, >200 at increments of 50
## max n iterations to n of possible combinations in data - no need to do 100 iterations if only 20 combinations possible
## explore alternative approach of increasing area of 50%UD (may not be much faster though)

## Phil Taylor & Mark Miller, 2012

## this script Iteratively subsamples a dataset of tracking data, investigating the
## affect of sample size. At each iteration the data is split, one half is used as the
## 'training' data and the 50%UD is calculated from this. The second half is used as
## 'testing' and the proportion of it, captured within the 50%UD is calculated.
## A perfect dataset would tend towards 0.5.
## By fitting a trend line to this relationship we can establish the sample size at which
## any new data would simply add to the existing knowledge. This script indicates how close to
## this value the inputted data are.

## DataGroup must be a dataframe or SpatialPointsDataFrame with Latitude, Longitude and ID as fields.
## Scale determines the smoothing factor used in the kernel analysis.
## Iteration determines the number of times each sample size is iterated.
    
## REVISED BY Steffen Oppel in 2015 to facilitate parallel processing
## updated to adehabitatHR by Steffen Oppel on 27 Dec 2016
## changed to same4all=TRUE on 4 Feb 2017
    
## REVISED in 2017 to avoid error in nls function of singular gradient
## added mean output for inclusion value even if nls fails

bootstrap <- function(DataGroup, Scale=100, Iteration=50)
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
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    }else{DgProj<-DataGroup@proj4string}

   DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  BoundBox <- bbox(DataGroup)
  UIDs <- unique(DataGroup$ID)
  Ntrips <- length(UIDs)
  Nloop<- seq(1,(Ntrips-1),ifelse(Ntrips>100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each=Iteration), Iteration=rep(seq(1:Iteration),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- 50

#setup parallel backend to use 4 processors
cl<-makeCluster(detectCores())
registerDoParallel(cl)
Result<-data.frame()

Result <- foreach(LoopN=LoopNr, .combine = rbind, .packages=c("sp","adehabitatHR","geosphere","rgdal")) %dopar% {

    N<-DoubleLoop$SampleSize[LoopN]
    i<-DoubleLoop$Iteration[LoopN]
    Coverage <- NULL
    Inclusion <- NULL
    History <- NULL

    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)

     RanNum <- sample(UIDs, N, replace=F)
      SelectedCoords <- coordinates(DataGroup[DataGroup$ID %in% RanNum,])
      NotSelected <- DataGroup[!DataGroup$ID %in% RanNum,]
      Temp <- data.frame(SelectedCoords[,1], SelectedCoords[,2])
      Ext <- (min(Temp[,1]) + 3 * diff(range(Temp[,1])))
      if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(Temp[,1]))))} else {BExt <- 5}
      Temp <- SpatialPoints(Temp,proj4string=DgProj)      ### added because adehabitatHR requires SpatialPoints object

      KDE.Surface <- adehabitatHR::kernelUD(Temp, h=(Scale * 1000), grid=1000, extent=BExt, same4all=TRUE)		## newer version needs SpatialPoints object and id no longer required in adehabitatHR
    try(      ## inserted to avoid function failing if one iteration in bootstrap has crazy extent
      KDE.UD <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent = TRUE)			## syntax differs from older version 
      #KDE.Spl <- kver2spol(KDE.UD)     ## deprecated, newer function would be khr2estUDm, but not required
      #KDE.Spl@proj4string <- DgProj    ## no longer necessary after update to adehabitatHR
      if("KDE.UD" %in% ls()){Overlain <- over(NotSelected, KDE.UD)   ## changed from KDE.Spl
      Output$InclusionMean <- length(which(!is.na(Overlain$area)))/nrow(NotSelected)  ## updated because overlay will not yield predictable number
      }                     ## inserted if statement to avoid function failing when one KDE.UD cannot be calculated
    return(Output)
    }

## stop the cluster
stopCluster(cl)

  par(mfrow=c(1,1), mai=c(1,1,1,1))
  #Result <- Output[1:nrow(Output)-1,]
  write.table(Result,"bootout_temp.csv", row.names=F, sep=",")
  try(M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
 if ('M1' %in% ls()){       ### run this only if nls was successful
  PredData <- data.frame(SampleSize = unique(Result$SampleSize))
  Result$pred<-predict(M1)
  P2 <- aggregate(pred~SampleSize, Result, FUN=mean)
  P2$sd <- aggregate(InclusionMean~SampleSize, Result, FUN=sd)[,2]
  plot(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray", ylim=c(0,1), xlim=c(0,max(PredData)), ylab="Inclusion", xlab="SampleSize")
  yTemp <- c((P2[,2] + 0.5*P2[,3]), rev(P2[,2] - 0.5*P2[,3]))
  xTemp <- c(P2$SampleSize, rev(P2$SampleSize))
  polygon(x=xTemp, y=yTemp, col="gray93", border=F)
  points(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
  lines(P2, lty=1,lwd=2)
  Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
  RepresentativeValue <- max(P2$pred)/Asymptote*100
  print(RepresentativeValue)
  Result$RepresentativeValue <- RepresentativeValue
  text(x=0, y=1,paste(round(RepresentativeValue,2), "%", sep=""), cex=2, col="gray45", adj=0)
}else{RepresentativeValue <- mean(Result$InclusionMean[Result$SampleSize==max(Result$SampleSize)])   ### if nls is unsuccessful then use mean output for largest sample size
  Result$RepresentativeValue <- (RepresentativeValue/(UDLev/100))*100}    ## added by Jono Handley to convert to same scale as nls output
  
  print(ifelse(exists("M1"),"nls (non linear regression) successful, asymptote estimated for bootstrap sample.",
       "WARNING: nls (non linear regression) unsuccessful, likely due to 'singular gradient', which means there is no asymptote. Data may not be representative, output derived from mean inclusion value at highest sample size. Check bootstrap output csv file"))

  return(Result)
  
  }

