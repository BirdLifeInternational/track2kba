library(track2KBA)
library(lubridate)

dat <- data.frame(Longitude = c(1, 1.01, 1.02, 1.04, 1.05, 1.03, 1), 
                  Latitude =  c(1, 1.01, 1.02, 1.03, 1.021, 1.01, 1),
                  ID = rep("A", 7),
                  DateTime = as.character(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:6))
)

expect_error(estSpaceUse(dat, scale=50, levelUD = 50), pattern = "data.frame")

trips <- projectTracks(dat, projType = "azim", custom = TRUE)

expect_error(estSpaceUse(trips, levelUD = 50), pattern = "scale")
expect_message(estSpaceUse(trips, scale = 50, levelUD = 50), 
               pattern = "grid resolution")

if (ignore(expect_silent)(
  estSpaceUse(trips, scale=50, res=5, levelUD = 50, polyOut=TRUE)
  )) {
  expect_silent(
    estSpaceUse(trips, scale=50, res=5, levelUD = 50, polyOut=TRUE)
  )
} else {
  expect_warning(
    estSpaceUse(trips, scale=50, res=5, levelUD = 50, polyOut=TRUE),
    pattern = "CRS object"
  )
}
if (ignore(expect_silent)(
  estSpaceUse(trips, res=5, scale=50, levelUD = 50, polyOut=TRUE)
  )) {
  expect_silent(
    estSpaceUse(trips, res=5, scale=50, levelUD = 50, polyOut=TRUE)
  )
} else {
  expect_warning(
    estSpaceUse(trips, res=5, scale=50, levelUD = 50, polyOut=TRUE),
    pattern = "CRS object"
  )
}
expect_message(estSpaceUse(trips, res=5, scale=.5, levelUD = 50), 
               pattern = "very small")
expect_warning(estSpaceUse(trips, res=5, scale=.5, levelUD = 50, polyOut = TRUE)
               , pattern = "percent failed")
expect_error(estSpaceUse(trips, res=5, scale=50, levelUD = 101, polyOut = TRUE)
               , pattern = "levelUD")

## error when not enough points for any ID
dat2 <- dat[-c(1:3),]
trips2 <- projectTracks(dat2, projType = "azim", custom = TRUE)

expect_error(estSpaceUse(trips2, scale=50, res=5, levelUD = 50), 
             pattern = "No IDs")

dat3 <- dat2 
dat3$ID <- rep("B")
dat3 <- rbind(dat, dat3)
trips3 <- projectTracks(dat3, projType = "azim", custom = TRUE)

expect_message(estSpaceUse(trips3, scale=50, res=5, levelUD = 50), 
               pattern = "too few points")

## timestamp duplicates in data
dat4 <- rbind(dat2[1, ], dat2)
trips4 <- projectTracks(dat4, projType = "azim", custom = TRUE)

expect_message(estSpaceUse(trips4, scale=50, res=5, levelUD = 50), 
               pattern = "duplicated data")
