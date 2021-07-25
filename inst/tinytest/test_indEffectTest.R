library(track2KBA)
library(lubridate)

dat <- data.frame(Longitude = rep(c(seq(1, 1.1, length.out=10), 
                                    seq(1.2, 1, length.out=10)), 2), 
                  Latitude =  rep(c(seq(1, 1.1, length.out=10), 
                                    seq(1.2, 1, length.out=10)), 2),
                  ID = c(rep("A", 40)),
                  DateTime = as.character(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:39))
)

colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])
trips <- tripSplit(dat, colony=colony, 
                   innerBuff = 1, returnBuff = 1, duration = 0.5, 
                   rmNonTrip = FALSE)
trips <- projectTracks(trips, projType = "azim", custom = TRUE)

expect_error(indEffectTest(trips, tripID = "tripID", groupVar = "ID", 
                           method = "BA", scale = 50, iterations = 10), 
             pattern = "not enough")
