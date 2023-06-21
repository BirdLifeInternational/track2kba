library(track2KBA)
library(lubridate)
library(adehabitatHR)
library(raster)
library(sp)

dat <- read.csv("boobies_testdata.csv")
trips <- projectTracks(dat, projType = "azim", custom = TRUE)
KDE <- estSpaceUse(trips, scale = 50, levelUD = 50, res = 25)

expect_error(repAssess(trips, levelUD = 50), pattern = "KDE argument")
expect_message(repAssess(trips, KDE, levelUD = 50))
expect_error(repAssess(trips, KDE, levelUD = 101), pattern = "levelUD")
expect_silent(repAssess(trips, KDE, iteration = 10, levelUD = 50))

# weighted averaging option
expect_message(repAssess(trips, KDE, levelUD = 50, avgMethod = "weighted"))

# input KDE data types
KDEsp <- do.call(cbind, lapply(KDE, function(x) {
  as(x, "SpatialPixelsDataFrame")
} ))
names(KDEsp) <- make.names(unique(trips$ID))
expect_silent(repAssess(trips, KDEsp, levelUD = 50))

KDEraster <- raster::stack(KDEsp)
expect_silent(repAssess(trips, KDEraster, levelUD = 50))


one <- subset(trips, trips$ID %in% unique(trips$ID)[1])
one$ID <- rep("noKDEid")
trips <- rbind(trips, one)

expect_silent(repAssess(trips, KDE, levelUD = 50))
expect_silent(repAssess(trips, KDE, levelUD = 50, bootTable = TRUE))

trips3 <- subset(trips, trips$ID %in% unique(trips$ID)[1:4])
KDE3 <- estSpaceUse(trips3, scale = 50, levelUD = 50, res = 25)

expect_message(repAssess(trips3, KDE3, levelUD = 50), pattern = "unsuccessful")
