
#### NEW FUNCTION TO OVERLAP ALL INDIVIDUAL HOME RANGE CENTRES AND IDENTIFY AREAS WHERE >THRESHOLD OCCUR

### NEW FUNCTION
## combines previous polyCount and rasterThresh functions into one
## requires col size and species as input
## first calculates threshold based on species and col size
## overlays all individual UDs and finds areas of intersection where required number of individual UDs overlap

## major problem are invalid geometries: https://www.r-spatial.org/r/2017/03/19/invalid.html
#install.packages("sf", type="source")
## version 1.2    05-04-2012

Polys=Output

polyCount <- function(Polys, spec, col_size)
{
  require(sf)
  require(raster)
  require(maps)
  require(lwgeom)
  
  if(!class(Polys) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")) stop("Polys must be a SpatialPolygonsDataFrame")
  if(is.na(projection(Polys))) stop("Polys must be projected")

  HR_sf <- st_as_sf(Polys) %>% st_transform(4326) ### convert to longlat CRS
plot(HR_sf["ID"])
sf::st_is_valid(HR_sf)
HR_sf_valid <- HR_sf %>% st_set_precision(100000) %>% lwgeom::st_make_valid()

### this simple function throws lots of TopologyException errors
## workaround found here does not work: https://github.com/r-spatial/sf/issues/603

# unique(HR_sf$ID) [1] 69208 69209 69210 69212 69213 69214 69216 69218

## problem with 69216 only exists in combination with certain other polygons
## any 3 of these work together, but all 4 in combination yield TopologyException: side location conflict
## st_buffer SOLVES this problem here but not when in combination with all other polygons!
problem<-HR_sf %>% filter(ID %in% c(69210,69209,69213,69216)) %>% st_buffer(dist=1) %>% st_make_valid()
iba = st_intersection(problem) # all intersections
plot(iba["n.overlaps"])

## works without 69216
test1<-HR_sf %>% filter(ID %in% c(69210,69208,69209,69212,69213,69214,69218,69216)) %>%
  st_buffer(dist=1) %>% st_make_valid()
iba = st_intersection(test1) # all intersections
plot(iba["n.overlaps"])

st_is_valid(HR_sf, reason = TRUE)


########## TRY the gOverlaps function
## worthless because does not return spatially explicit area where overlap occurs

#### described here: http://r-sig-geo.2731867.n2.nabble.com/Counting-overlapping-polygons-in-a-given-area-td7586232.html

 
over(Polys,Polys)
library(rgeos) 
gO <- gOverlaps(Polys, byid=c(TRUE)) 
dim(gO) 

makes a logical matrix with TRUE where scot_BNG[i,] overlaps with 
polys[j]. 

count_by_city <- apply(gO, 2, sum) 








### THIS DOES NOT WORK BECAUSE I CANNOT INSTALL spacebucket
# library(devtools)
# devtools::install_github("mdsumner/spacebucket")
# library(spacebucket)
# geom <- sf::read_sf("C:\\temp\\trouble\\trouble_geom.shp")
# bucket <- spacebucket(geom)
# par(mfrow = c(3, 2))
# for (i in seq(nrow(geom), 1)) {
#   i_overlap <- n_intersections(bucket, i)
#   plot(st_geometry(geom), reset = FALSE, main = sprintf("%i overlaps\n %f", i, sum(st_area(i_overlap))))
#   
#   plot(i_overlap[1], add = TRUE, reset = FALSE)
#   
# }
# 
# 
# https://r-spatial.github.io/sf/reference/geos_binary_ops.html


## polyCount  ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## polyCount calculates the number of overlapping Polygons. The calculation is
## done by overlapping the polygons within a grid and counting the number
## falling within each gridcell. The Res parameter sets the size of the grid
## (in decimal degrees) to be used. Smaller grids will allow for more detailed
## counts but will slow computing time. The function returns a raster object with
## values showing the proportion of Polys (i.e. nPolys/total nPolys) overlapping
## each cell. The raster extent is set at the bounding limits of the Polys or, when
## data crosses the dateline, set to the northern and southern most limits of the
## Polys but longitudinally crossing the circumference of the world.

## Polys must be a SpatialPolygonsDataFrame of the polygons to be counted.
## Res must be a numeric object indicating the resolution in decimal degrees.

## Steffen Oppel revision on 14 Dec 2016 - removed + (Res * 100) from NCol in L. 664
## Steffen Oppel revision on 16 Dec 2016: fixed the point overlay problem by using a proper grid instead
## Steffen Oppel revision on 27 Dec 206: ensured that output was projected in WGS84

## version 1.2    05-04-2012

polyCount <- function(Polys, Res = 0.1)
  {

  require(raster)
  require(maps)

  if(!class(Polys) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")) stop("Polys must be a SpatialPolygonsDataFrame")
  if(is.na(projection(Polys))) stop("Polys must be projected")

  Poly.Spdf <- spTransform(Polys, CRS=CRS("+proj=longlat +ellps=WGS84"))
  DgProj <- Polys@proj4string

  DateLine <- Poly.Spdf@bbox[1,1] < -178 & Poly.Spdf@bbox[1,2] > 178
  if(DateLine == TRUE) {print("Data crosses DateLine")}

  UDbbox <- bbox(Poly.Spdf)
  if(DateLine == TRUE)  {UDbbox[1,] <- c(-180,180)}
  BL <- floor(UDbbox[,1])                 # + (Res/2) - removed on 16 Dec 2016 because it results in some polygons outside the grid
  TR <- ceiling(UDbbox[,2])
  NRow <- ceiling(sqrt((BL[1] - TR[1])^2)/Res)
  NCol <- ceiling(sqrt((BL[2] - TR[2])^2)/Res) #+ (Res * 100)				### THIS LINE CAUSES PROBLEMS BECAUSE IT GENERATES LATITUDES >90 which will cause spTransform to fail
  Grid <- GridTopology(BL, c(Res,Res), c(NRow, NCol))
    newgrid<-SpatialGrid(Grid, proj4string = CRS("+proj=longlat + datum=wgs84"))
    spol <- as(newgrid, "SpatialPolygons")								### this seems to create an orphaned hole
    SpGridProj <- spTransform(spol, CRS=DgProj)
    GridIntersects <- over(SpGridProj, Polys)
    SpGridProj<- SpatialPolygonsDataFrame(SpGridProj, data = data.frame(ID=GridIntersects$ID, row.names=sapply(SpGridProj@polygons,function(x) x@ID)))
    SpGridProj <- subset(SpGridProj, !is.na(SpGridProj@data$ID))
  #SpGrid <- SpatialPoints(Grid, proj4string = CRS("+proj=longlat + datum=wgs84"))
  #SpdfGrid <- SpatialPointsDataFrame(SpGrid, data.frame(Longitude=SpGrid@coords[,1], Latitude=SpGrid@coords[,2]))
  #SpGridProj <- spTransform(SpdfGrid, CRS=DgProj)
  #GridIntersects <- over(SpGridProj, Polys)
  #SpGridProj@data$Intersects$ID <- GridIntersects$ID
  #SpGridProj <- subset(SpGridProj, !is.na(SpGridProj@data$Intersects$ID))   ### SpGridProj[!is.na(SpGridProj@data$Intersects$ID),] 			###
  plot(SpGridProj)

  Count <- 0
  for(i in 1:length(Polys))
    {
    TempB <- Polys[i,]
    Temp <- over(SpGridProj, TempB)[,1]     ### inserted based on Matthew Carroll's advice; MAY NEED TO SWAP arguments in 'over'?
    Temp[is.na(Temp)] <- 0
    Temp[Temp > 0] <- 1
    Count <- Count + Temp
    #Prop <- Count/i                        ### removed to improve efficiency
    }
  Prop <- Count/length(Polys)     ### removed from loop over polys as it only needs to be calculated once
  #GridIntersects$inside<-as.numeric(as.character(GridIntersects$ID))    ### this only works for numeric trip_id!!
  #GridIntersects$Prop <- 0
  #GridIntersects$Prop[!is.na(GridIntersects$inside)] <- Prop    #[,1] removed based on Matthew Carroll's advice, because fixed in L. 678
  SpGridProj@data$Prop <- Prop
  SpGridOUT <- spTransform(SpGridProj, CRS=CRS("+proj=longlat +ellps=WGS84"))   ### show output in WGS84
  SGExtent <- extent(SpGridOUT)
  RT <- raster(SGExtent, ncols=as.double(NCol), nrows=as.double(NRow))
  WgsRas <- (rasterize(x=SpGridOUT,y=RT, field = "Prop"))

  plot(WgsRas, asp=1)
  maps::map("world", add=T, fill=T, col="darkolivegreen3")      ## to avoid conflict with purrr
  projection(WgsRas) <- CRS("+proj=longlat + datum=wgs84")
  return(WgsRas)
  }


## thresholdRaster     ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## thresholdRaster applies a threshold to a raster, and isolates any areas above
## that threshold value. Converting the raster values to polygons is difficult and
## so this part takes some time. The function returns a SpatialPolygonsDataFrame
## containing the polygons that are above threshold, and with an attributes table
## holding each sites Maximum raster value.

## CountRas must be a raster object with values 0 - 1.
## Threshold must be a number indicating the percentage value to be used
## as the threshold.


thresholdRaster <- function(CountRas, Threshold = 10)
    {

    require(raster)
    require(maps)
    require(geosphere)

    plot(CountRas, asp=1)
    maps::map("world", add=T, fill=T, col="darkolivegreen3")
    Threshold <- Threshold/100
    RasSites <- CountRas >= Threshold
    plot(RasSites, asp=1, col=rev(heat.colors(25)))
    maps::map("world", add=T, fill=T, col="darkolivegreen3")      ### to avoid conflict with purrr

    if(length(which(getValues(CountRas) > Threshold)) < 1)
      {
      Mid <- c(bbox(CountRas)[1,1]+(as.numeric(bbox(CountRas)[1,2])- as.numeric(bbox(CountRas)[1,1]))/2,  midPoint(bbox(CountRas)[,2], bbox(CountRas)[,1])[2])
      text(Mid, "No Site Identified", cex=1.25)
      stop("No cells were above the threshold value")
      }

    #Cells <- rasterToPolygons(CountRas, fun=function(x) {x>Threshold})
    Cells <- rasterToPolygons(clump(CountRas>Threshold), dissolve=TRUE)   ## suggested by Ian Cleasby, as more efficient
    DateLine <- Cells@bbox[1,1] < -178 & Cells@bbox[1,2] > 178
    if(DateLine == TRUE)    {Cells <- spTransform(Cells, CRS=DgProj)}

    Sites <- Cells                                                        ## removed based on suggestion from Ian Cleasby: dissolve(Cells)
    ifelse(DateLine == TRUE, projection(Sites) <- DgProj, projection(Sites) <- "+proj=longlat + datum=wgs84")
    Sites <- spTransform(Sites, CRS=CRS("+proj=longlat + datum=wgs84"))
    SiteTable <- data.frame(SiteID = names(Sites), MaxPerc = round(raster::extract(CountRas, Sites, fun=max)*100,2))    ## to avoid conflict with tidyverse
    Sites <- SpatialPolygonsDataFrame(Sites, data=SiteTable)
    print(SiteTable)
    return(Sites)
}







n_intersections <- function(x, n = 2, ...) {
  triangles <- x$index %>%
    dplyr::group_by(.data$triangle_idx) %>%
    dplyr::mutate(nn = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$nn >= n) %>%
    dplyr::transmute(path = .data$path_, .data$triangle_idx)
  gmap <- x$geometry_map %>%
    dplyr::select(.data$object, .data$layer, .data$path)
  ## every unique triangle keeps a record of which path, object, layer
  ## (a bit of redundancy until we get a single path/object index or ...)
  idx <- purrr::map_df(split(triangles, triangles$triangle_idx),
                       function(piece) {
                         ## path joins us to layer + object
                         piece %>% dplyr::inner_join(gmap, "path")
                       }) %>% dplyr::group_by(.data$triangle_idx) %>% tidyr::nest()
  
  ## now build each triangle
  P <- x$primitives$P
  TR <- x$primitives$T
  sf::st_sf(idx = idx, geometry = sf::st_sfc(purrr::map(idx$triangle_idx, ~sf::st_polygon(list(P[TR[.x, ][c(1, 2, 3, 1)], ])))))
}





################## TRY A non-sf SOLUTION ######


library(sp)
library(rgeos)
library(raster)
library(rworldmap)

box <- readWKT("POLYGON((-180 90, 180 90, 180 -90, -180 -90, -180 90))")
proj4string(box) <- CRS("+proj=cea +datum=WGS84")
set.seed(1)
pts <- spsample(box, n=2000, type="random")
pols <- gBuffer(pts, byid=TRUE, width=50) # create circle polys around each point
merge = sample(1:40, 100, replace = T) # create vector of 100 rand #s between 0-40 to merge pols on

Sp.df1 <- gUnionCascaded(pols, id = merge) # combine polygons with the same 'merge' value
# create SPDF using polygons and randomly assigning 1 or 2 to each in the @data df
Sp.df <- SpatialPolygonsDataFrame(Sp.df1,
                                  data.frame(z = factor(sample(1:2, length(Sp.df1), replace = TRUE)),
                                             row.names= unique(merge)))
Sp.df <- crop(Sp.df, box)
colors <- c(rgb(r=0, g=0, blue=220, alpha=50, max=255), rgb(r=220, g=0, b=0, alpha=50, max=255))

land <- getMap()

overlay.map <- spplot(Sp.df, zcol = "z", col.regions = colors,
                      col = NA, alpha = 0.5, breaks=c(0,1),
                      sp.layout = list("sp.polygons", land, fill = "transparent",
                                       col = "grey50"))
overlay.map

# find the count of polygons below each grid cell
GT <- GridTopology(c(-179.5, -89.5), c(0.01, 0.01), c(360, 180))
SG <- SpatialGrid(GT)
SG.proj<-spTransform(SG,proj4string(Output))


o <- over(SG, Polys, returnList=TRUE)
ct <- sapply(o, length)
summary(ct)
SGDF <- SpatialGridDataFrame(SG, data=data.frame(ct=ct))
spplot(SGDF, "ct", col.regions=bpy.colors(20))




