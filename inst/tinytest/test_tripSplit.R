library(track2KBA)
library(sp)
library(sf)
library(rgdal)

dat <- data.frame(Longitude = rep(c(1:10, 10:1), 2), 
                  Latitude =  rep(c(1:10, 10:1), 2),
                  ID = c(rep("A", 20), rep("B", 20)),
                  DateTime = format(
                    ymd_hms("2021-01-01 00:00:00") + hours(0:19))
)

colony <- data.frame(Longitude = dat$Longitude[1], Latitude = dat$Latitude[1])

## test that data.frame input returns SPDF object 
expect_identical(class(
  tripSplit(dat, colony=colony, innerBuff = 1, returnBuff = 1))[[1]], 
  "SpatialPointsDataFrame")

## test that no NAs returned when DateTime is input as a ymd_HMS string
dat2 <- dat %>% mutate(DateTime = as.character(DateTime))

expect_true(
  sum(
    is.na(
      tripSplit(dat2, colony=colony, innerBuff = 1, returnBuff = 1)$DateTime)
    ) == 0)

## test that SPDFs and sf objects return error message 
xy <- dat[,c(1,2)]
dat_spdf <- sp::SpatialPointsDataFrame(
  data = dat, coords = xy,
  proj4string = 
    sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  )
dat_sf <- sf::st_as_sf(dat_spdf)

expect_error(tripSplit(dat_spdf, colony=colony, innerBuff = 1, returnBuff = 1))
expect_error(tripSplit(dat_sf, colony=colony, innerBuff = 1, returnBuff = 1))

## no duration specification still works based on buffs
expect_true(
  length(unique(tripSplit(dat, colony=colony, 
                          innerBuff = 1, returnBuff = 1,
                          rmNonTrip=TRUE)$tripID)) == 2,
  info = "When no duration specified"
)

## check that innerBuff works 
expect_true(
  all(unique( tripSplit(dat, colony=colony, 
                              innerBuff = 100000000, returnBuff = 1, 
                              rmNonTrip = FALSE)$tripID) == "-1"),
  info = "innerBuff check 1"
)

## function warns about timestamp duplicates in data

dat3 <- data.frame(Longitude = rep(c(1:2, 2:1)), 
                     Latitude =  rep(c(1:2, 2:1)),
                     ID = rep("A", 4),
                     DateTime = format(
                       ymd_hms("2021-01-01 00:00:00") + hours(c(0, 0:2)))
)

expect_message(
  tripSplit(dat3, colony=colony, innerBuff = 1, returnBuff = 1, duration = 0.5)
)

## check that returnBuff works 
dat2 <- dat[1:10, ]

expect_true(
  all(unique( tripSplit(dat2, colony=colony, 
                        innerBuff = 1, returnBuff = 1, 
                        rmNonTrip = TRUE)$Returns) == "No"),
  info = "returnBuff check 1"
)
expect_false(
  all(unique( tripSplit(dat2, colony=colony, 
                        innerBuff = 1, returnBuff = 2000000, 
                        rmNonTrip = TRUE)$Returns) == "No"),
  info = "returnBuff check 1"
)

## test verbose option
expect_silent(tripSplit(dat2, colony=colony, 
                        innerBuff = 1, returnBuff = 1, duration=0/5, 
                        verbose=FALSE))

## check that rmNonTrip option works 
expect_false(
  "-1" %in% unique( tripSplit(dat, colony=colony, 
                              innerBuff = 1, returnBuff = 1, 
                              rmNonTrip = TRUE)$tripID),
  info = "rmNonTrip check 1"
)
expect_true(
  "-1" %in% unique( tripSplit(dat, colony=colony, 
                              innerBuff = 1, returnBuff = 1, 
                              rmNonTrip = FALSE)$tripID),
  info = "rmNonTrip check 2"
)

## test colony specifications 
colony <- data.frame(Longitude = c(1,2), 
                     Latitude = c(1,2),
                     ID = c("A", "B"))
expect_true(
  length(unique(tripSplit(dat, colony=colony, 
                          innerBuff = 1, returnBuff = 1, 
                          nests = TRUE)$tripID)) == 4,
  info = "colony check 1"
)
expect_error(tripSplit(dat, colony=colony, 
                         innerBuff = 1, returnBuff = 1, duration = 0.5,
                       nests = FALSE),
             info = "colony check 2"
)

colnames(colony) <- c("longitude", "latitude")
expect_error(tripSplit(dat, colony=colony, innerBuff = 1, returnBuff = 1),
             info = "colony check 2"
)

