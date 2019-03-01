
#### BULLSHIT TEST

### NEW FUNCTION
## combines previous polyCount and rasterThresh functions into one
## requires col size and species as input
## first calculates threshold based on species and col size
## overlays all individual UDs and finds areas of intersection where required number of individual UDs overlap

set.seed(131)



library(sf)

m = rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))

p = st_polygon(list(m))

n = 100

l = vector("list", n)



for (i in 1:n)
  
  l[[i]] = p + 10 * runif(2)

s = st_sfc(l)

plot(s, col = sf.colors(categorical = TRUE, alpha = .5))

title("overlapping squares")



sf = st_sf(s)

i = st_intersection(sf) # all intersections

plot(i["n.overlaps"])



https://r-spatial.github.io/sf/reference/geos_binary_ops.html


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
