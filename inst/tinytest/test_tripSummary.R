library(track2KBA)
library(lubridate)

dat <- data.frame(Longitude = rep(c(1:10, 10:1), 2), 
                  Latitude =  rep(c(1:10, 10:1), 2),
                  ID = c(rep("A", 20), rep("B", 20)),
                  DateTime = as.character(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])

trips <- tripSplit(dat, colony=colony, 
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = FALSE)

## check basic use cases work w/out warnings/errors
expect_silent(tripSummary(trips, colony))
expect_silent(tripSummary(trips, colony, extraDist = TRUE))
expect_silent(tripSummary(trips@data, colony, extraDist = TRUE))

trips <- subset(trips, trips$tripID != "-1")
expect_silent(tripSummary(trips, colony))

## bad Lat/Long names don't work
colony_badname <- data.frame(longitude = dat$Longitude[1], 
                             latitude = dat$Latitude[1])

expect_error(tripSummary(trips, colony_badname, extraDist = TRUE))

## Test that nests option works
colony_nests <- data.frame(Longitude = c(1,2), 
                     Latitude = c(1,2),
                     ID = c("A", "B"))

trips <- tripSplit(dat, colony=colony_nests,
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = TRUE, nests = TRUE)

expect_true(nrow(tripSummary(trips, colony_nests, nests = TRUE))==3)


## warning message when some trips don't return to colony
dat <- data.frame(Longitude = rep(c(1:10)), 
                  Latitude =  rep(c(1:10)),
                  ID = rep("A", 20),
                  DateTime = as.character(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

trips <- tripSplit(dat, colony=colony,
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = TRUE)

expect_warning(tripSummary(trips, colony))
