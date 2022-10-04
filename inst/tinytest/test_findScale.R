library(track2KBA)
library(lubridate)
library(sp)
library(dplyr)

dat <- data.frame(Longitude = c(1, 1.01, 1.02, 1.04, 1.05, 1), 
                  Latitude =  c(1, 1.01, 1.02, 1.03, 1.021, 1),
                  ID = rep("A", 6),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:5))
)

expect_error(findScale(dat), "data.frame", info = "data.frame")

colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])
trips <- tripSplit(dat, colony=colony, 
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = TRUE)

sumTrips <- tripSummary(trips, colony)
sumTrips$max_dist <- 19
trips <- projectTracks(trips, projType = "azim", custom = TRUE)

expect_message(findScale(trips), pattern = "sumTrips")
expect_message(findScale(trips), pattern = "res")

expect_true(is.na(findScale(trips)$med_max_dist))
expect_true(is.na(findScale(trips)$mag))

expect_false(is.na(findScale(trips, sumTrips=sumTrips)$med_max_dist))
expect_false(is.na(findScale(trips, sumTrips=sumTrips)$mag))

sumTrips_b <- sumTrips
sumTrips_b$max_dist <- 199

expect_false(findScale(trips, sumTrips=sumTrips_b)$med_max_dist == 
               findScale(trips, sumTrips=sumTrips)$med_max_dist,
             info = "frange differs")

expect_message(findScale(trips, res=1000000), pattern = "very large")

expect_silent(findScale(trips, peakMethod = "max"))
expect_silent(findScale(trips, peakMethod = "steep"))

## scales
sumTrips$max_dist <- 49
expect_silent(findScale(trips, sumTrips=sumTrips))
sumTrips$max_dist <- 75
expect_silent(findScale(trips, sumTrips=sumTrips))

## must be projected
xy <- dat[,c(1,2)]
dat_spdf <- sp::SpatialPointsDataFrame(
  data = dat, coords = xy,
  proj4string =
    sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)

expect_error(findScale(dat_spdf), pattern = "projected", 
             info = "projected data only")

trips2 <- trips[, -5]

expect_true(findScale(trips2)$step_length == findScale(trips)$step_length, 
             info = "step length by trip or only ID")


## data that cover massive distances
dat <- data.frame(Longitude = c(1, 1.1, 2, 2.1, 1.5, 1), 
                  Latitude =  c(1, 1.01, 1.02, 1.03, 1.021, 1),
                  ID = rep("A", 6),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:5))
)

trips <- tripSplit(dat, colony=colony, 
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = TRUE)

sumTrips <- tripSummary(trips, colony)
trips <- projectTracks(trips, projType = "azim", custom = TRUE)

expect_true(findScale(trips, sumTrips=sumTrips)$med_max_dist > 20)
