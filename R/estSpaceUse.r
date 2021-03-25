### estSpaceUse  ###############################################################

#' Estimate the space use of tracked animals using kernel utilization 
#' distribution
#'
#' \code{estSpaceUse} is a wrapper for \code{\link{kernelUD}} which estimates 
#' the utilization distribution (UD) of multiple individuals or tracks in a 
#' tracking dataset.
#'
#' A utilization distribution will be calculated for each unique 'ID'. The data 
#' should be regularly sampled or interpolated (see adehabitatLT package for 
#' functions to this end).
#' 
#' If desired \code{res} results in memory-heavy grid (e.g. >10,000 cells) 
#' use \code{polyOut = FALSE} to speed things up.
#'
#' @param tracks SpatialPointsDataFrame. Must be projected into an equal-area 
#' coordinate system. If not, first run \code{\link{projectTracks}}.
#' @param scale numeric (in kilometers). The smoothing parameter ('H') used in 
#' the kernel density estimation, which defines the width of the normal 
#' distribution around each location. The \code{\link{findScale}} function can 
#' assist in finding sensible scales.
#' @param levelUD numeric (percent). Specify which utilization distribution 
#' contour at which to subset the polygon output. NOTE: This will only affect the output
#'  if \code{polyOut=TRUE}.
#' @param res numeric (in square kilometers). Grid cell resolution 
#' for kernel density estimation. Default is a grid of 
#' 500 cells, with spatial extent determined by the latitudinal and longitudinal
#'  extent of the data.
#' @param polyOut logical scalar (TRUE/FALSE). If TRUE then output will include 
#' a plot of individual UD polygons and a simple feature with kernel UD polygons
#'  for the level of \code{levelUD}. NOTE: creating polygons of UD is both 
#'  computationally slow and prone to errors if the usage included in 
#'  \code{levelUD} extends beyond the specified grid. In this case 
#'  \code{estSpaceUse} will return only the \code{estUDm}object and issue a 
#'  warning.
#' @return Returns an object of class \code{estUDm} which is essentially a list,
#'  with each item representing the utilization distribution of a level of 'ID'.
#'   Values in the output signify the usage probability per unit area for that 
#'   individual in each grid cell. This can be converted into a 
#'   SpatialPixelsDataFrame via the \code{adehabitatHR::estUDm2spixdf} function.
#'
#' If \code{polyOut=TRUE} the output will be a list with two components: 
#' \emph{'KDE.Surface'} is the \code{estUDm} object and \code{UDPolygons} is 
#' polygon object of class \code{sf} (Simple Features) with the UD contour for 
#' each individual at the specified \code{levelUD}.
#'
#' If \code{polyOut=TRUE} but the polygon delineation in 
#' \code{adehabitatHR::getverticeshr} fails, output is an object of class 
#' \code{estUDm} and a warning will be issued.
#'
#' @seealso \code{\link{formatFields}}, \code{\link{tripSplit}}, 
#' \code{\link{findScale}}
#'
#' @examples
#' \dontrun{UDs <- estSpaceUse(tracks=Trips, scale = HVALS$mag, levelUD = 50, 
#' polyOut=TRUE)}
#'
#' @export
#' @import sf
#' @import sp
#' @import adehabitatHR
#' @importFrom stats na.omit quantile sd var
#' @importFrom methods as

estSpaceUse <- function(
  tracks, scale = 50, levelUD, res=NULL, polyOut=FALSE) {
  
  # check for duplicated data
  dup_check <- tracks@data %>% group_by(.data$ID) %>% 
    mutate(duplicates = duplicated(.data$DateTime)) %>% ungroup() %>% 
    summarise(duplicates = sum(.data$duplicates))
  if(dup_check$duplicates > 0){message(
    "WARNING:dataset may contain duplicated data, this will affect KDE"
  )}
  
  tracks@data <- tracks@data %>% dplyr::select(.data$ID)
  
  ### REMOVING IDs WITH TOO FEW LOCATIONS -------------------------------------

  validIDs <- names(which(table(tracks$ID) > 5))
  if(length(validIDs) > n_distinct(tracks$ID) ){
    message(
      paste0("Following ID(s) have too few points for KDE: ", 
             unique(tracks$ID)[!unique(tracks$ID) %in% validIDs])
      )
  }
  tracks <- tracks[(tracks@data$ID %in% validIDs), ]
  tracks@data$ID <- droplevels(as.factor(tracks@data$ID))

  ### SETTING PARAMETERS FOR kernelUD -----------------------------------------

  ### CREATE CUSTOM GRID for kernelUD (instead of same4all=TRUE) --------------
  extendX <- max(diff(range(coordinates(tracks)[,1]))*0.05, scale*2000)
  extendY <- max(diff(range(coordinates(tracks)[,2]))*0.05, scale*2000)
    
  minX <- min(coordinates(tracks)[,1]) - extendX
  maxX <- max(coordinates(tracks)[,1]) + extendX
  minY <- min(coordinates(tracks)[,2]) - extendY
  maxY <- max(coordinates(tracks)[,2]) + extendY
  
  ### if users do not provide a res, split data into ~500 cells ---------------
  if( is.null(res) ){ 
    res <- ( max(abs(minX-maxX)/500, abs( minY-maxY )/500) ) / 1000
  message(
    sprintf(
    "No grid resolution ('res') was specified, or the specified resolution was 
    >99 km and therefore ignored. Space use was calculated on a 500-cell grid, 
    with cells of %s square km", 
    round(res, 3)
      )
    )
  }

  ### specify sequence of grid cells and combine to SpatialPixels -------------
  xrange <- seq(minX,maxX, by = res*1000) 
  yrange <- seq(minY,maxY, by = res*1000) 
  grid.locs <- expand.grid(x=xrange, y=yrange)
  INPUTgrid <- SpatialPixels(
    SpatialPoints(grid.locs), proj4string=proj4string(tracks)
    )

  ### ERROR CATCH IF PEOPLE SPECIFIED TOO FINE RESOLUTION ---------------------
  if (scale < res*0.1228){
  message(
  "Your scale parameter is very small compared to the grid resolution - 
  99.99% of the kernel density for a given location may be within a single grid
  cell, which will limit the amount of overlap of different individual's core 
  use areas. Increase 'scale' or reduce 'res' to avoid this problem.")}
  if ((length(xrange) * length(yrange)) > 100000){
  message(
  "Your grid has a pretty large number of cells - this may slow down 
  computation. Increase 'res' to speed things up.")}
  if ((length(xrange) * length(yrange)) > 1000000){warning("Your grid is 
    >1 million pixels, computation may be VERY slow and may max out R's memory capacity.")}

  ### ESTIMATING KERNEL DISTRIBUTION  -----------------------------------------
  KDE.Surface <- adehabitatHR::kernelUD(
    tracks, h=(scale * 1000), grid=INPUTgrid, same4all=FALSE
    )

  ###### OPTIONAL POLYGON OUTPUT ----------------------------------------------
  if(polyOut==TRUE){
    if(!levelUD >= 1 & levelUD <= 100) {stop("levelUD must be between 1-100%")}
    tryCatch({
          KDE_sp <- adehabitatHR::getverticeshr(
            KDE.Surface, percent = levelUD, unin = "m", unout = "km2"
            )
        }, error=function(e){
        sprintf("Providing individual home range polygons at a UD level of %s 
          percent failed with the following error message: %s. This means that 
          there was estimated space use that extended beyond the grid used for 
          estimating the kernel density. To resolve this, use a lower UD level,
          or a smaller scale parameter.", levelUD,conditionMessage(e))
          }
      )

      if(('KDE_sp' %in% ls())){ ## PROCEED ONLY IF KDE successful

        KDE_sf <- st_as_sf(KDE_sp) %>%
                  st_transform(4326) ### convert to longlat CRS
        
            return(list(KDE.Surface=KDE.Surface, UDPolygons=KDE_sf))

            }else{
            warning(
            sprintf("Providing individual home range polygons at a UD 
              level of %s percent failed. This often means that there was 
              estimated space use that extended beyond the grid used for 
              estimating the kernel density. To resolve this, use a lower UD 
              level, a smaller scale parameter, or a smaller size of grid cells 
              (smaller 'res'). However, the same error may occur if your scale 
              parameter is very small (compared to 'res'), because the grid is 
              extended beyond the bounding box of locations by a distance of 
              2*scale - and with a very small scale parameter that may not 
              actually encompass a full grid cell. ", levelUD), 
              immediate. = TRUE)
            return(KDE.Surface)
              }

  }else{
    return(KDE.Surface)
    }
}

