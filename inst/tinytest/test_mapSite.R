library(track2KBA)
library(rgdal)
library(sp)

dat <- read.csv("boobies_testdata.csv")
colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])
trips <- projectTracks(dat, projType = "azim", custom = TRUE)
KDE <- estSpaceUse(trips, scale = 50, levelUD = 50, res = 25)
sites <- findSite(KDE, represent = 100, levelUD = 50, popSize = 100)

expect_error(mapSite("sites"), pattern = "object")
expect_silent(mapSite(sites))
expect_silent(mapSite(sites, colony = colony))
expect_silent(mapSite(sites, colony = colony, show = FALSE))

sites_sf <- findSite(KDE, represent = 100, levelUD = 50, popSize = 100, 
                  polyOut = TRUE)

expect_silent(mapSite(sites_sf))
expect_silent(mapSite(sites_sf, colony = colony))
expect_silent(mapSite(sites_sf, colony = colony, show = FALSE))

sites2 <- findSite(KDE, represent = 100, levelUD = 50, polyOut = TRUE)

expect_silent(mapSite(sites2))
