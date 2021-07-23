library(track2KBA)
library(lubridate)

dat <- data.frame(Longitude = c(1, 2, 3, 1), 
                  Latitude =  c(1, 1, 2, 1),
                  ID = rep("A", 4),
                  DateTime = as.character(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:3))
)

expect_silent(projectTracks(dataGroup = dat, projType = "azim", custom = TRUE))
expect_silent(projectTracks(dataGroup = dat, 
                            projType = "cylin", custom = TRUE))
expect_message(projectTracks(dataGroup = dat, projType = "azim", custom = FALSE))
expect_message(projectTracks(dataGroup = dat, 
                            projType = "cylin", custom = FALSE))

xy <- dat[,c(1,2)]
dat_spdf <- sp::SpatialPointsDataFrame(
  data = dat, coords = xy,
  proj4string =
    sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)

expect_silent(projectTracks(dataGroup = dat_spdf,
                            projType = "azim", custom = TRUE))

## can transform projected SPDF
dat_prj <- projectTracks(dataGroup = dat, 
              projType = "cylin", custom = FALSE)

expect_silent(projectTracks(dataGroup = dat_prj,
                           projType = "azim", custom = TRUE))

## needs Lat Long columns to do so
dat_prj <- dat_prj[, -(1:2)]

expect_error(projectTracks(dataGroup = dat_prj,
                           projType = "azim", custom = TRUE), 
             pattern = "missing")

## testing functionality along International Date Line
dat_idl <- data.frame(Longitude = rep(c(-170:-179, 179:170)), 
                      Latitude =  rep(c(1:10, 10:1)),
                      ID = c(rep("A", 20)),
                      DateTime = as.character(
                        ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

expect_silent(
  projectTracks(dataGroup = dat_idl, projType = "azim", custom = TRUE))
expect_silent(
  projectTracks(dataGroup = dat_idl, projType = "cylin", custom = TRUE))


xy <- dat_idl[,c(1,2)]
dat_spdf <- sp::SpatialPointsDataFrame(
  data = dat_idl, coords = xy,
  proj4string =
    sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)

expect_silent(
  projectTracks(dataGroup = dat_spdf, projType = "azim", custom = TRUE))
expect_silent(
  projectTracks(dataGroup = dat_spdf, projType = "cylin", custom = TRUE))

