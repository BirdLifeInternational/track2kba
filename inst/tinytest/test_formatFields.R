library(track2KBA)
library(lubridate)

dat <- data.frame(x = rep(1:10, 2), y = rep(1:10, 2),
                  id = c(rep("A", 10), rep("B", 10)),
                  dt = as.character(ymd_hms("2021-01-01 00:00:00") + hours(0:9))
                  )

## test that misspecified DT column name gives error 
expect_error(formatFields(dat, 
                          fieldID = "id",
                          fieldDateTime = "t",
                          fieldLat = "y",
                          fieldLon = "x"))

## test that 'cleanDF' option works ## 
expect_true(!"dt" %in% colnames(formatFields(dat,
                                             fieldID = "id",
                                             fieldDateTime = "dt",
                                             fieldLat = "y",
                                             fieldLon = "x",
                                             formatDT = "ymd HMS",
                                             cleanDF = TRUE)
                                ) )
expect_false(!"dt" %in% colnames(formatFields(dat,
                                             fieldID = "id",
                                             fieldDateTime = "dt",
                                             fieldLat = "y",
                                             fieldLon = "x",
                                             formatDT = "ymd HMS",
                                             cleanDF = FALSE)
) )

## BirdLife STD standard format ##
dat2 <- data.frame(longitude = rep(1:10, 2), latitude = rep(1:10, 2),
                  track_id = c(rep("A", 10), rep("B", 10)),
                  date_gmt = rep("2021-01-01"),
                  time_gmt = do.call(rbind, strsplit(dat$dt, split=" "))[,2]
)

## test that 'formatBL' option works ## 
expect_true( "DateTime" %in% colnames(formatFields(dat2,
                                             formatBL = TRUE)) )            
expect_true( "ID" %in% colnames(formatFields(dat2,
                                                  formatBL = TRUE)) ) 
expect_true( "Latitude" %in% colnames(formatFields(dat2,
                                                  formatBL = TRUE)) )
expect_true( "Longitude" %in% colnames(formatFields(dat2,
                                                  formatBL = TRUE)) )
expect_true( "POSIXct" %in% class(formatFields(dat2, 
                                               formatBL = TRUE)$DateTime)  )