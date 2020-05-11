#### indEffectTest ####

#' Test site fidelity
#'
#' \code{indEffectTest} tests whether the variance in overlap between space use 
#' areas within a group (e.g within individuals) is significant compared to 
#' between groups (e.g. between individuals).
#'
#' This function works by producing kernel density areas at a desired contour 
#' level (i.e. \emph{UDLEv}) for each level of \emph{tripID} and estimating the 
#' degree of overlap between all pairwise comparisons using the desired overlap 
#' \emph{method}. Then, comparisons are split into 'within' and 'between' 
#' groups, determined by the grouping variable (i.e \emph{groupVar}) argument.
#'
#' If \emph{conditional=TRUE} then the overlap estimates will range from 0 to 
#' \emph{levelUD} (unless \emph{method="HR"}).
#'
#' Then, the empirical distribution of each group is compared in a bootstrapped 
#' Kolmogorov-Smirnov test, to check whether differences in the distributions 
#' are significant. If so, it indicates that individuals within the 
#' \emph{groupVar} reuse sites more than expected by chance.
#'
#' NOTE: Because \code{indEffectTest} relies on 
#' \code{\link[adehabitatHR]{kerneloverlap}} to estimate overlap, it was not 
#' possible to implement a \emph{res} argument as is done in other track2KBA 
#' functions. Therefore, it is advised to either leave the default of 500 cells,
#'  or ascertain the number of cells in the grid of chosen \emph{res} from the 
#'  output of \link{estSpaceUse}.
#'
#' @param tracks SpatialPointsDataFrame. Must be in an equal-area projection. 
#' See \code{\link{projectTracks}}
#' @param tripID character. Column in \emph{tracks} corresponding to the within 
#' group ID (e.g. trip-individual combination)
#' @param groupVar character. Column in \emph{tracks} corresponding to the 
#' between group ID (e.g. individual or track)
#' @param plot logical scalar (TRUE/FALSE). Do you want to output a boxplot of 
#' the result?
#' @param method character. Which method of overlap estimation to use? See 
#' \code{\link[adehabitatHR]{kerneloverlap}} for descriptions of each method.
#' @param conditional logical scalar (T/F). If TRUE, the function sets to 0 the
#'  pixels of the grid over which the UD is estimated, outside the home range of
#'   the animal estimated at a level of probability equal to percent. Note that
#'    this argument has no effect when meth="HR" 
#'    (from \code{\link[adehabitatHR]{kerneloverlap}}).
#' @param levelUD numeric. The desired contour level of the utilization 
#' distribution to be used in overlap estimation. NOTE: this is irrelevant if 
#' \emph{conditional=FALSE}.
#' @param scale numeric (in kilometers). Smoothing ('H') parameter for 
#' kernel density estimation.
#' @param grid numeric or SpatialPixels. If numeric, specify the desired number 
#' of grid cells over which the utilization distributions will be esimated. 
#' A default grid of 500 cells is used.
#' @param iterations numeric. Indicate the desired number of Kolmogorov-Smirnov 
#' iterations to run. 500 is an advisable minimum for statistical rigor.
#'
#' @return \code{indEffectTest} returns a list containing three objects. In the 
#' first slot 'Overlap Matrix', the full matrix of overlap comparisons. In the '
#' Overlap' slot, a dataframe with a column identifying whether each overlap 
#' estimate corresponds to a within-group, or a between-group comparison. 
#' In the third slot 'Kolmogorov-Smirnov' is the test output of the 
#' Kolmogorov-Smirnov test, indicating the D parameter and significance 
#' estimates.
#'
#' @examples
#' \dontrun{ indEffect <- indEffectTest(
#' tracks, groupVar="ID", tripID="trip_id", method="BA", scale=HVALS$mag) 
#' }
#'
#' @export
#' @importFrom ggplot2 aes geom_boxplot
#' @import sp

indEffectTest <- function(
  tracks, tripID, groupVar, plot=TRUE, 
  method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"), 
  conditional = TRUE, levelUD=50, scale, grid = 500, iterations = 1000) {
  
  if (!requireNamespace("Matching", quietly = TRUE)) {
    stop("Package \"Matching\" needed for  function to work. Please install.",
      call. = FALSE)  }
  
  if (!(tripID) %in% names(tracks)) stop("Within-group field does not exist")
  if (!groupVar %in% names(tracks)) stop("Group field does not exist")


  tracks@data <- tracks@data %>% 
    dplyr::select({{groupVar}}, {{tripID}}, .data$Latitude, .data$Longitude)

  # remove tripID tracks with < 6 points as can't be used in KDE --------------
  UIDs <- names(which(table(tracks@data[, tripID]) > 5))
  tracks <- tracks[tracks@data[, tripID] %in% UIDs, ]
  tracks@data[ ,tripID] <- droplevels(as.factor(tracks@data[ ,tripID]))

  # create vector with value of groupVar for each trip ------------------------
  gid <- tracks@data[!duplicated(tracks@data[, tripID]), ][[groupVar]]

  # calculate overlap between tracks ------------------------------------------
  X <- adehabitatHR::kerneloverlap(
    xy = tracks[, tripID], 
    method = method, 
    percent = levelUD, 
    conditional = conditional, 
    h = scale*1000, 
    grid = grid
    )
  X[lower.tri(X, diag = TRUE)] <- NA

  # assign value of groupVar to rows and columns -----------------------------
  rownames(X) <- colnames(X) <- gid

  # separate within (WI) and between (BW) group overlaps ---------------------
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

  # organize values in a dataframe for plotting ------------------------------
  Overlaps <- data.frame(
    Overlap = c(WI, BW), 
    Type = c(rep("Within", 
      length(WI)), 
      rep("Between", 
        length(BW)))
    )

  if(plot==TRUE){
    print(ggplot(
      data = Overlaps, 
      aes(x = .data$Type, 
        y = .data$Overlap, 
        fill = .data$Type)
      ) + 
        geom_boxplot())
  }

  # ks <- ks.test(x = WI, y = BW)
  ks <- Matching::ks.boot(
    WI, BW, alternative = "two.sided", nboots = iterations
    )

  # Organise output
  Result <- list()
  Result[1] <- list(X)        # overlaps matrix
  Result[2] <- list(Overlaps) # df with overlap values (long format)
  Result[3] <- list(ks)       # output from the ks.boot function
  names(Result) <- c("Overlap Matrix", "Overlaps", "Kolmogorov-Smirnov")
  
  return(Result)
}
