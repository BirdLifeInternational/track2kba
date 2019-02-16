## varianceTest    ###########################################################################################

## NEEDS COMPLETE OVERHAUL OR REWRITE
## VIRGINIA MORERA MAY HAVE PROVIDED SOLUTION


## Phil Taylor & Mark Miller, 2012

## this script tests the spatial variance between polygons, investigating
## whether the variance between polygons in a group is significantly different
## to the variance between polygons from different groups.
## groups are determined using the DID field.

## Polys must be a series of polygons of class SpatialPolygonsDataFrame.
## ID must be a vector of equal length the Polys, identifying the group each Polys belongs to.
## Iteration must be an integer indicating the number of times to iterate the random sample.


varianceTest <- function(Polys, DID, Iteration=100)
  {



  require(rgeos)
  require(sp)
  #set.seed(1)
  if(length(Polys) != length(DID)) stop("DID must be the same length as Polys")
  if(class(Iteration) != "numeric") stop("Iteration must be an integer")
  TripTable <- data.frame(Index = 1:length(Polys), DID = DID)
  TripTable$Multiple <- TripTable[,2] %in% names(which(table(TripTable[,2]) > 1))
  if(length(Polys) ==  length(unique(TripTable[,2]))) {stop("Each Poly has a unique DID, Variance is Zero, Polys must be categorised")}
  if(class(Polys) != "SpatialPolygonsDataFrame") stop("Polys must be of class SpatialPolygonsDataFrame")
  Output <- NULL
  for(i in unique(TripTable[TripTable$Multiple == TRUE,2]))
    {
    PValue <- numeric()   #### changed
        for(k in 1:Iteration)
    {
    Selected <- TripTable[TripTable$DID == i,]
    Sample <- Polys[which(TripTable$DID == i),]

    if(nrow(Selected) > length(unique(TripTable[TripTable$DID != i,2])))
      {
      TRand <- sample(nrow(Selected), length(unique(TripTable[TripTable[,2] != i,2])))
      Selected <- Selected[TRand,]
      Sample <- Sample[TRand,]
      }
    Sample <- gSimplify(Sample, tol=100)
    Reps <- 1:nrow(Selected)
    TDist <- NULL
    for(l in Reps)
      {
      Pairs <- Reps[-(which(Reps == l))]
      if(length(which(Pairs > l)) < 1) {break}
      Pairs <- Pairs[which(Pairs > l)]
      for(m in Pairs)
        {
        Temp <- as.double(gDistance(Sample[l,], Sample[m,], hausdorff = T))
        TDist <- rbind(TDist, Temp)
        }
      }
    RanBirds <- sample(unique(TripTable[TripTable[,2] != i,2]), nrow(Selected))
    ITrips <- NULL
    for(l in RanBirds)
      {
      RTrips <- TripTable[TripTable[,2] == l,]
      STrips <- sample(as.character(RTrips[,1]),1)
      ITrips <- c(ITrips, STrips)
      }
    ISample <- Polys[as.double(ITrips),]
    ISample <- gSimplify(ISample, tol=100)
    Reps <- 1:nrow(Selected)
    IDist <- NULL              ### changed location
    for(l in Reps)
      {
      Pairs <- Reps[-(which(Reps == l))]
      if(length(which(Pairs > l)) < 1) {break}
      Pairs <- Pairs[which(Pairs > l)]
      for(m in Pairs)
        {
        Temp <- as.double(gDistance(ISample[l,], ISample[m,], hausdorff = T))
        IDist <- rbind(IDist, Temp)
        }
      }
    MW1 <- suppressWarnings(wilcox.test(TDist, IDist, alternative = "two.sided", paired =F))
    Result <- suppressWarnings(data.frame(DependentValue = TDist[,1], IndependentValue = IDist[,1]))
    boxplot(Result, main= paste("p.value =", round(MW1$p.value, 6)))
    PValue <- c(PValue, round(MW1$p.value, 6))  ### changed

    }

  Temp <- data.frame(DID = i, PVal = mean(PValue))   ### changed
  Output <- rbind(Output, Temp)
  }
  ## The Null Hypothesis is that the within individual DIDs are part of the same
  ## distribution as the between individual.
  ## Rejection level is set at 0.25
  PValueF <- round(mean(Output$PVal, na.rm=T),6)
  par(mfrow=c(1,1), mai=c(0.5,0.5,0.5,0.5))
  hist(Output$PVal, main= paste("p.value =", PValueF), border="grey", breaks=20)
  if(PValueF < 0.25)
    {
    legend("center", PValueF, "Do not reject the Null Hypothesis,
    site fidelity is significant", cex=2, bty="n" )
    } else {legend("center", PValueF, "Reject the Null Hypothesis,
    site fidelity is not significant", cex=2, bty="n" )
    }
  TripTable <- merge(TripTable, Output, all.x = T)
  return(PValueF)
  }