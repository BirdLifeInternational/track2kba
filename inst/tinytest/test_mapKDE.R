library(track2KBA)
library(lubridate)

dat <- data.frame(Longitude = c(1, 1.01, 1.02, 1.04, 1.05, 1.03, 1), 
                  Latitude =  c(1, 1.01, 1.02, 1.03, 1.021, 1.01, 1),
                  ID = rep("A", 7),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:6))
)
colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])
trips <- projectTracks(dat, projType = "azim", custom = TRUE)
KDEud <- estSpaceUse(trips, scale = 50, levelUD = 50)

expect_silent(mapKDE(KDEud))

KDEhr <- estSpaceUse(trips, scale = 50, levelUD = 50, polyOut = TRUE)

expect_error(mapKDE(KDEhr), pattern = "UDPolygons layer")
expect_silent(mapKDE(KDEhr$UDPolygons))
expect_silent(mapKDE(KDEhr$UDPolygons, colony = colony))
expect_silent(mapKDE(KDEhr$UDPolygons, show = FALSE))
