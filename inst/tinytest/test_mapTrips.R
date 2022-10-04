library(track2KBA)

dat <- data.frame(Longitude = rep(c(1:10, 10:1), 2), 
                  Latitude =  rep(c(1:10, 10:1), 2),
                  ID = c(rep("A", 20), rep("B", 20)),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])

trips <- tripSplit(dat, colony=colony, 
                   innerBuff = 1, returnBuff = 1,
                   rmNonTrip=TRUE)

expect_error(mapTrips(trips))
expect_error(mapTrips(trips@data))
expect_message(mapTrips(trips, colony))
expect_silent(mapTrips(trips, colony, colorBy = "trip"))
expect_silent(mapTrips(trips, colony, IDs = 1))


## testing functionality along International Date Line
dat_idl <- data.frame(Longitude = rep(c(-170:-179, 179:170)), 
                  Latitude =  rep(c(1:10, 10:1)),
                  ID = c(rep("A", 20)),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

colony_idl <- data.frame(Longitude = dat_idl$Longitude[1], 
                         Latitude = dat_idl$Latitude[1])
trips_idl <- tripSplit(dat_idl, colony=colony_idl, 
                   innerBuff = 1, returnBuff = 1)

expect_silent(mapTrips(trips_idl, colony_idl))
