#' St. Helena Masked Boobies
#'
#' A GPS tracking data set of Masked Boobies during incubation and chick-rearing at St. Helena Island. 
#' Formatted following BirdLife International's Seabird Tracking Database standard \url{www.seabirdtracking.org}. Data from Oppel et al. 2015. 
#'
#' @format A data frame with 116355 obs. of  6 variables:
#' \describe{
#'   \item{track_id}{Unique identifier code for each bird}
#'   \item{date_gmt}{Character vector representing date (Greenwich Mean Time)}
#'   \item{time_gmt}{Character vector representing time (Greenwich Mean Time)}
#'   \item{longitude}{Longitudinal position of bird}
#'   \item{latitude}{Latitudinal position of bird}
#'   \item{lon_colony}{Longitudinal position of breeding colony}
#'   \item{lat_colony}{Latitudinal position of breeding colony}
#'   ...
#' }
#' @source \url{https://link.springer.com/article/10.1007/s00265-015-1903-3}
"boobies"