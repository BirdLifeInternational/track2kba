#### IndEffectTest ####

#' Test site fidelity
#'
#' \code{IndEffectTest} tests whether the variance in overlap between space use areas within a group (e.g within individuals) is significant compared to between groups (e.g. between individuals).
#'
#' This function works by producing kernel density areas at a desired contour level (i.e. \emph{UDLEv}) for each level of \emph{tripID} and estimating the degree of overlap between all pairwise comparisons using the desired overlap \emph{method}. Then, comparisons are split into 'within' and 'between' groups, determined by the grouping variable (i.e \emph{GroupVar}) argument.
#'
#' If \emph{conditional=TRUE} then the overlap estimates will range from 0 to \emph{UDLev} (unless \emph{method="HR"}).
#'
#' Then, the empirical distribution of each group is compared in a bootstrapped Kolmogorov-Smirnov test, to check whether differences in the distributions are significant. If so, it indicates that individuals within the \emph{GroupVar} reuse sites more than expected by chance.
#'
#' NOTE: Because \code{IndEffectTes} relies on \code{\link[adehabitatHR]{kerneloverlap}} to estimate overlap, it was not possible to implement a \emph{Res} argument as is done in other track2KBA functions. Therefore, it is advised to either leave the default of 500 cells, or ascertain the number of cells in the grid of chosen \emph{Res} from the output of \link{estSpaceUse}.
#'
#' @param Trips SpatialPointsDataFrame or data.frame. If input is data.frame or unprojected SpatialPointsDF, must include 'Latitude' and 'Longitude' fields.
#' @param tripID character. Column in \emph{Trips} corresponding to the within group ID (e.g. trip-individual combination)
#' @param GroupVar character. Column in \emph{Trips} corresponding to the between group ID (e.g. individual or track)
#' @param plot logical scalar (TRUE/FALSE). Do you want to output a boxplot of the result?
#' @param method character. Which method of overlap estimation to use? See \code{\link[adehabitatHR]{kerneloverlap}} for descriptions of each method.
#' @param conditional logical scalar (T/F). If TRUE, the function sets to 0 the pixels of the grid over which the UD is estimated, outside the home range of the animal estimated at a level of probability equal to percent. Note that this argument has no effect when meth="HR" (from \code{\link[adehabitatHR]{kerneloverlap}}).
#' @param UDLev numeric. The desired contour level of the utilization distribution to be used in overlap estimation. NOTE: this is irrelevant if \emph{conditional=FALSE}.
#' @param Scale numeric (in kilometers). Smoothing ('H') parameter for kernel density estimation.
#' @param grid numeric or SpatialPixels. If numeric, specify the desired number of grid cells over which the utilization distributions will be esimated. A default grid of 500 cells is used.
#' @param nboots numeric. Indicate the desired number of Kolmogorov-Smirnov iterations to run. 500 is an advisable minimum for statistical rigor.
#'
#' @return \code{IndEffectTest} returns a list containing three objects. In the first slot 'Overlap Matrix', the full matrix of overlap comparisons. In the 'Overlap' slot, a dataframe with a column identifying whether each overlap estimate corresponds to a within-group, or a between-group comparison. In the third slot 'Kolmogorov-Smirnov' is the test output of the Kolmogorov-Smirnov test, indicating the D parameter and significance estimates.
#'
#' @examples
#' \dontrun{
#' indEffect <- IndEffectTest(Trips, GroupVar="ID", tripID="trip_id", method="BA", Scale=HVALS$mag)
#'
#' indEffect$`Kolmogorov-Smirnov`}
#'
#' @export
#' @import ggplot2
#' @import sp

IndEffectTest <- function(Trips, tripID, GroupVar, plot=T, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), conditional = TRUE, UDLev=50, Scale, grid = 500, nboots = 1000)
{
  if (!requireNamespace("Matching", quietly = TRUE)) {
    stop("Package \"Matching\" needed for this function to work. Please install it.",
      call. = FALSE)  }
  
  # initial chceks
  if (!"Latitude" %in% names(Trips)) stop("Latitude field does not exist")
  if (!"Longitude" %in% names(Trips)) stop("Longitude field does not exist")
  if (!(tripID) %in% names(Trips)) stop("Within-group field does not exist")
  if (!GroupVar %in% names(Trips)) stop("Group field does not exist")


  # MB # Added this section which converts Trips to spatial dataframe and projects it (and if already is SPDF it accepts this)
  if (class(Trips) != "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    if (!"Latitude" %in% names(Trips)) stop("Latitude field does not exist")
    if (!"Longitude" %in% names(Trips)) stop("Longitude field does not exist")
    ## filter DF to the minimum fields that are needed
    CleanTrips <- Trips %>%
      dplyr::select(.data$GroupVar, .data$tripID, .data$Latitude, .data$Longitude)
    mid_point <- data.frame(geosphere::centroid(cbind(CleanTrips$Longitude, CleanTrips$Latitude)))

    ### PREVENT PROJECTION PROBLEMS FOR DATA SPANNING DATELINE
    if (min(CleanTrips$Longitude) < -170 &  max(CleanTrips$Longitude) > 170) {
      longs <- ifelse(CleanTrips$Longitude < 0, CleanTrips$Longitude + 360, CleanTrips$Longitude)
      mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}

    Trips.Wgs <- SpatialPoints(data.frame(CleanTrips$Longitude, CleanTrips$Latitude), proj4string = CRS("+proj=longlat + datum=wgs84"))
    proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
    Trips.Projected <- spTransform(Trips.Wgs, CRS = proj.UTM )
    TripsSpatial <- SpatialPointsDataFrame(Trips.Projected, data = CleanTrips)
    TripsSpatial@data <- TripsSpatial@data %>% dplyr::select(.data$GroupVar, .data$tripID, .data$Latitude, .data$Longitude)
    Trips.Wgs <- NULL
    Trips.Projected <- NULL

  }else {## if data are already in a SpatialPointsDataFrame then check for projection
    if (is.projected(Trips)) {
      if ("trip_id" %in% names(Trips@data)) {
        TripsSpatial <- Trips }
      TripsSpatial@data <- TripsSpatial@data %>% dplyr::select(.data$GroupVar, .data$tripID, .data$Latitude, .data$Longitude)
    }else {## project data to UTM if not projected
      if (!"Latitude" %in% names(Trips)) stop("Latitude field does not exist")
      if (!"Longitude" %in% names(Trips)) stop("Longitude field does not exist")
      mid_point <- data.frame(geosphere::centroid(cbind(Trips@data$Longitude, Trips@data$Latitude)))

      ### MB  This part prevents projection problems around the DATELINE
      if (min(Trips@data$Longitude) < -170 &  max(Trips@data$Longitude) > 170) {
        longs <- ifelse(Trips@data$Longitude < 0, Trips@data$Longitude + 360, Trips@data$Longitude)
        mid_point$lon <- ifelse(median(longs) > 180, median(longs) - 360, median(longs))}

      proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep = ""))
      TripsSpatial <- spTransform(Trips, CRS = proj.UTM)
      TripsSpatial@data <- TripsSpatial@data %>% dplyr::select(.data$GroupVar, .data$tripID, .data$Latitude, .data$Longitude)
    }
  }

  # remove tripID Trips with < 6 points as they can't be used to calculate kernel
  # MB edit # Changed this step to happen after SPDF set-up. Also added tripID column.
  UIDs <- names(which(table(TripsSpatial@data[, tripID]) > 5))
  TripsSpatial <- TripsSpatial[TripsSpatial@data[, tripID] %in% UIDs, ]
  TripsSpatial@data[ ,tripID] <- droplevels(as.factor(TripsSpatial@data[ ,tripID]))

  # create vector with value of GroupVar for each trip
  gid <- TripsSpatial@data[!duplicated(TripsSpatial@data[, tripID]), ][[GroupVar]]

  # calculate overlap between Trips
  X <- adehabitatHR::kerneloverlap(xy = TripsSpatial[, tripID], method = method, percent = UDLev, conditional = conditional, h = Scale*1000, grid = grid)
  X[lower.tri(X, diag = TRUE)] <- NA

  # assign value of GroupVar to rows and columns
  rownames(X) <- colnames(X) <- gid

  # separate within (WI) and between (BW) group overlaps
  WI <- NULL
  BW <- NULL
  for (i in seq_along(rownames(X))) {
    # i = 1
    x1 <- X[i,]
    x2 <- x1[which(names(x1) == rownames(X)[i])]
    x3 <- x1[which(names(x1) != rownames(X)[i])]
    WI <- c(WI, x2)
    BW <- c(BW, x3)
  }
  BW <- BW[!is.na(BW)]
  WI <- WI[!is.na(WI)]

  ## VMP commented this out since the ks.boot function is robust to ties. Was leftover from when using stats::ks.test function
  # BW <- BW[BW != 0]
  # WI <- WI[WI != 0]

  # organize values in a dataframe for plotting
  Overlaps <- data.frame(Overlap = c(WI, BW), Type = c(rep("Within", length(WI)), rep("Between", length(BW))))

  if(plot==TRUE){
    print(ggplot(data = Overlaps, aes(x = Type, y = Overlap, fill = Type)) + geom_boxplot() + theme_bw())
  }

  # ks <- ks.test(x = WI, y = BW)
  ks <- Matching::ks.boot(WI, BW, alternative = "two.sided", nboots = nboots) # more indicated when data don't come from continuous distr (ours have many 0s)

  # Organise output
  Result <- list()
  Result[1] <- list(X) # overlaps matrix
  Result[2] <- list(Overlaps) # df with overlap values (long format)
  Result[3] <- list(ks) # output from the ks.boot function
  names(Result) <- c("Overlap Matrix", "Overlaps", "Kolmogorov-Smirnov")
  return(Result)
}
