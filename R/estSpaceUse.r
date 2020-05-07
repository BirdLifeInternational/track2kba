## estSpaceUse  ####

#' Estimate the space use of tracked animals using kernel utilization distribution
#'
#' \code{estSpaceUse} is a wrapper for \code{\link{kernelUD}} which estimates the utilization distribution (UD) of multiple individuals or tracks in a tracking dataset.
#'
#' A utilization distribution will be calculated for each unique 'ID'. The data should be regularly sampled or interpolated (see adehabitatLT package for functions to this end).
#'
#' @param DataGroup data.frame or SpatialPointsDataFrame. Must have an 'ID' field; if input is data.frame or unprojected SpatialPointsDF, must also include 'Latitude' and 'Longitude' fields. For formatting issues see \code{\link{formatFields}}.
#' @param Scale numeric (in kilometers). The smoothing parameter ('H') used in the kernel density estimation, which defines the width of the normal distribution around each location. The \code{\link{findScale}} function can assist in finding sensible scales.
#' @param UDLev numeric (percent). Specify which utilization distribution contour to show in the plotted output. NOTE: this will only affect the output if \code{polyOut=TRUE}.
#' @param Res numeric (in kilometers). Grid cell resolution (in square kilometers) for kernel density estimation. Default is a grid of 500 cells, with spatial extent determined by the latitudinal and longitudinal extent of the data.
#' @param polyOut logical scalar (TRUE/FALSE). If TRUE then output will include a plot of individual UD polygons and a simple feature with kernel UD polygons for the level of \code{UDLev}. NOTE: creating polygons of UD is both computationally slow and prone to errors if the usage included in \code{UDLev} extends beyond the specified grid. In this case \code{estSpaceUse} will return only the \code{estUDm}object and issue a warning.
#' @param plot logical scalar (TRUE/FALSE). If TRUE, map will be produced showing core areas, each level of ID with a different color. NOTE:\code{polyout} must be TRUE for this to work.  
#' @return Returns an object of class \code{estUDm} which is essentially a list, with each item representing the utilization distribution of a level of 'ID'. Values in the output signify the usage probability per unit area for that individual in each grid cell. This can be converted into a SpatialPixelsDataFrame via the \link[adehabitatHR]{estUDm2spixdf} function.
#'
#' If \code{polyOut=TRUE} the output will be a list with two components: \emph{'KDE.Surface'} is the \code{estUDm} object and \code{UDPolygons} is polygon object of class \code{sf} (Simple Features) with the UD contour for each individual at the specified \code{UDLev}.
#'
#' If \code{polyOut=TRUE} but the polygon delineation in \code{adehabitatHR::getverticeshr} fails, output is an object of class \code{estUDm} and a warning will be issued.
#'
#' @seealso \code{\link{formatFields}}, \code{\link{tripSplit}}, \code{\link{findScale}}
#'
#' @examples
#' \dontrun{UDs <- estSpaceUse(DataGroup=Trips, Scale = HVALS$mag, UDLev = 50, polyOut=T)}
#'
#' @export
#' @import sf
#' @import sp
#' @import adehabitatHR
#' @import ggplot2
#' @importFrom stats na.omit quantile sd var
#' @importFrom methods as

estSpaceUse <- function(DataGroup, Scale = 50, UDLev = 50, Res=NULL, polyOut=FALSE, plot=FALSE)
    {

  DataGroup@data <- DataGroup@data %>% dplyr::select(.data$ID)
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### REMOVING IDs WITH TOO FEW LOCATIONS ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### ENSURE THAT kernelUD does not fail by retaining ONLY TRACKS WITH 5 OR MORE LOCATIONS
  validIDs <- names(which(table(DataGroup$ID)>5))
  DataGroup<-DataGroup[(DataGroup@data$ID %in% validIDs),]
  DataGroup@data$ID<-droplevels(as.factor(DataGroup@data$ID))        ### encountered weird error when unused levels were retained (27 Feb 2017)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### SETTING PARAMETERS FOR kernelUD : THIS NEEDS MORE WORK!! ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

  ### CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)
  extendX <- max(diff(range(coordinates(DataGroup)[,1]))*0.05, Scale*2000)
  extendY <- max(diff(range(coordinates(DataGroup)[,2]))*0.05, Scale*2000)
    
  minX <- min(coordinates(DataGroup)[,1]) - extendX
  maxX <- max(coordinates(DataGroup)[,1]) + extendX
  minY <- min(coordinates(DataGroup)[,2]) - extendY
  maxY <- max(coordinates(DataGroup)[,2]) + extendY
  
  ### if users do not provide a resolution, then split data into ~500 cells
  if( is.null(Res) ){ Res <- ( max(abs(minX-maxX)/500, abs( minY-maxY )/500) )/1000
  warning(sprintf("No grid resolution ('Res') was specified, or the specified resolution was >99 km and therefore ignored.
                  Space use was calculated on a 500-cell grid, with cells of %s square km", round(Res,3)),immediate. = TRUE)}

  ### specify sequence of grid cells and combine to SpatialPixels
  xrange <- seq(minX,maxX, by = Res*1000) #diff(range(coordinates(DataGroup)[,1]))/Res)   ### if Res should be provided in km we need to change this
  yrange <- seq(minY,maxY, by = Res*1000) #diff(range(coordinates(DataGroup)[,2]))/Res)   ### if Res should be provided in km we need to change this
  grid.locs <- expand.grid(x=xrange,y=yrange)
  INPUTgrid <- SpatialPixels(SpatialPoints(grid.locs), proj4string=proj4string(DataGroup))
  #  plot(INPUTgrid)

  #### ERROR CATCH IF PEOPLE SPECIFIED TOO FINE RESOLUTION ####
  if (Scale < Res*0.1228){warning("Your scale parameter is very small compared to the grid resolution - 99.99% of the kernel density for a given location will be within a single grid cell, which will limit the amount of overlap of different individual's core use areas. Increase 'Scale' or reduce 'Res' to avoid this problem.",immediate. = TRUE)}
  if (max(length(xrange),length(yrange))>600){warning("Your grid has a pretty large number of cells - this will slow down computation. Increase 'Res' (= make the grid cells larger) to speed up the computation.",immediate. = TRUE)}
  if (max(length(xrange),length(yrange))>1200){stop("Are you sure you want to run this function at this high spatial resolution (= very small grid cell size specified in 'Res')? Your grid is >1 million pixels, computation will take many hours (or days)!")}

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### ESTIMATING KERNEL DISTRIBUTION  ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## may need to insert extent=BExt, but hopefully avoided by custom-specified grid
  ## switched from same4all=T to =F because we provide a fixed input grid
  KDE.Surface <- adehabitatHR::kernelUD(DataGroup, h=(Scale * 1000), grid=INPUTgrid, same4all=F)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ###### OPTIONAL POLYGON OUTPUT ####
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  if(polyOut==TRUE){
    tryCatch({
          KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")
        }, error=function(e){
        sprintf("Providing individual home range polygons at a UD level of %s percent failed with the following error message: %s. This means that there was estimated space use that extended beyond the grid used for estimating the kernel density. To resolve this, use a lower UD level, or a smaller Scale parameter.", UDLev,conditionMessage(e))})

      if(('KDE.Sp' %in% ls())){     ## PROCEED ONLY IF KDE.Sp was successfully created

        HR_sf <- st_as_sf(KDE.Sp) %>%
                  st_transform(4326) ### convert to longlat CRS
        
            return(list(KDE.Surface=KDE.Surface, UDPolygons=HR_sf))

            }else{
            warning(sprintf("Providing individual home range polygons at a UD level of %s percent failed. This often means that there was estimated space use that extended beyond the grid used for estimating the kernel density. To resolve this, use a lower UD level, a smaller Scale parameter, or a smaller size of grid cells (smaller 'Res'). However, the same error may occur if your Scale parameter is very small (compared to 'Res'), because the grid is extended beyond the bounding box of locations by a distance of 2*Scale - and with a very small Scale parameter that may not actually encompass a full grid cell. ", UDLev),immediate. = TRUE)
            return(KDE.Surface)}

  }else{
    return(KDE.Surface)
    }
}

