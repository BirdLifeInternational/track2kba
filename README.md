
<!-- README.md is generated from README.Rmd. Please edit that file -->
track2KBA
=========

This package is comprised of functions that facilitate the identification of areas of importance for biodiversity, such as Key Biodiversity Areas (KBAs) or Ecologically or Biologically Significant Areas (EBSAs), based on individual tracking data. Key functions include utilities to identify and summarise individual foraging trips, estimate utilization distributions, and overlay distributions to identify important aggregation areas. Utility functions are also included to download Movebank data, format data sets, as well functions to assess sample representativeness and space use independence.

Installation
------------

You can download the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("steffenoppel/track2iba", auth_key=ASK MARTIN FOR THIS!, dependencies=TRUE)
```

Example
-------

The `formatFields` function allows you to specify which columns in your data set correspond to the necessary fields for track2KBA analysis and re-format them. These are: a datetime field, latitude and longitude fields, and an ID field (i.e. individual animal, track, or trip).

In this example, we use a publicly available GPS dataset of Brown Pelicans, published in Geary et al. 2018, and stored in Movebank's data repository (<https://www.datarepository.movebank.org/>).

``` r
library(track2KBA)

data(boobies)

tracks <- formatFields(boobies, 
  field_ID = "track_id", 
  field_Date = "date_gmt", 
  field_Time = "time_gmt",
  field_Lon = "longitude", 
  field_Lat = "latitude"
  )
#> Warning in formatFields(boobies, field_ID = "track_id", field_Date =
#> "date_gmt", : No format was supplied for the input Date and Time fields,
#> a default format ('ymd_HMS') was attempted when combining the fields. If
#> an error is produced, see help page ('?lubridate::parse_date_time') for
#> information on date formats.

## basic example code
```

If your data come from a species which makes trips out from a centrally-located place, such as a nest in the case of a bird, or a beach colony in the case of a pinniped, you can use `tripSplit` to split up the data into discrete trips.

In order to do this, you must identify the location of the central place (e.g. nest or colony).

``` r
library(dplyr)

colony <- tracks %>% 
  summarise(
    Longitude = first(Longitude), 
    Latitude = first(Latitude))
```

Our *colony* dataframe tells us where trips originate from. Then we need to set some parameters to decide what constitutes a trip. To do that we should use our understanding of the movement ecology of the study species; Brown Pelicans belong to a coastal, nearshore species, and do not travel great distances on breeding season foraging trips. So in this case we set *InnerBuff* to 1 km, and *Duration* to 1 hour. *ReturnBuff* can be used to catch incomplete trips, where the animal began returning, but perhaps due to device failure the full trip wasn't captured. For short-ranging species with data from many trips this may be set to the same as *InnerBuff*.

Setting *rmColLocs* to TRUE will remove those points falling within the *InnerBuff*.

``` r
trips <- tripSplit(
  tracks = tracks, 
  colony, 
  InnerBuff = 3, 
  ReturnBuff = 10, 
  Duration = 1, 
  plotit = T, 
  rmColLocs = T)
#> [1] "track 693041 does not return to the colony"
#> [1] "track 693434 does not return to the colony"
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Then we can summarize the trip movements, using `tripSummary`.

``` r
tripSum <- tripSummary(Trips = trips, Colony = colony)
#> Warning in tripSummary(Trips = trips, Colony = colony): Some trips did not
#> return to the specified return buffer around the colony. The return date
#> given for these trips refers to the last location of the trip, and NOT the
#> actual return time to the colony.

tripSum
#> # A tibble: 215 x 10
#> # Groups:   ID [41]
#>    ID    trip_id n_locs departure           return              duration
#>    <chr> <chr>    <dbl> <dttm>              <dttm>                 <dbl>
#>  1 69302 693021     274 2012-07-22 07:52:11 2012-07-22 16:11:03     8.31
#>  2 69302 693022     124 2012-07-23 12:26:22 2012-07-23 15:54:05     3.46
#>  3 69302 693023     138 2012-07-25 08:30:53 2012-07-25 12:22:53     3.87
#>  4 69304 693041    1268 2013-08-22 11:50:41 2013-08-24 23:03:50    NA   
#>  5 69305 693051      71 2013-08-22 13:08:15 2013-08-22 15:10:59     2.05
#>  6 69306 693061      37 2014-01-06 16:28:42 2014-01-06 17:32:11     1.06
#>  7 69306 693062      83 2014-01-07 14:48:24 2014-01-07 17:10:21     2.37
#>  8 69306 693063     129 2014-01-08 14:25:11 2014-01-08 17:55:22     3.50
#>  9 69306 693064      50 2014-01-08 18:08:32 2014-01-08 19:30:40     1.37
#> 10 69306 693065     155 2014-01-09 14:47:04 2014-01-09 19:21:32     4.57
#> # ... with 205 more rows, and 4 more variables: total_dist <dbl>,
#> #   max_dist <dbl>, direction <dbl>, complete <chr>
```

Now that we have an idea how the animals are moving, we can start with the process of estimating their space use, and potential sites of aggregation!

`findScale` provides us with options for setting the all-important smoothing parameter in the Kernel Density Estimation.

If we know our animal uses area-restricted search to locate prey, then we can set the `ARSscale=T`. This will use First Passage Time analysis to identify the spatial scale at which area-restricted search is occuring.

``` r
Hvals <- findScale(trips,
  ARSscale = T,
  Colony = Colony,
  Trips_summary = tripSum)
#> Warning in findScale(trips, ARSscale = T, Colony = Colony, Trips_summary
#> = tripSum): No grid resolution ('Res') was specified, or the specified
#> resolution was >99 km and therefore ignored. Movement scale in the data was
#> compared to a 500-cell grid with cell size of 0.701 km squared.
#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used

#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used
#> [1] "No peak was found for: ID 69309"
#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used
#> [1] "No peak was found for: ID 69314"
#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used
#> [1] "No peak was found for: ID 69328"
#> [1] "No peak was found for: ID 69330"
#> [1] "No peak was found for: ID 69332"
#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used

#> Warning in if ((FirstPeak == Scales[length(Scales) - 1]) & (FirstPeak == :
#> the condition has length > 1 and only the first element will be used

Hvals
#>   med_max_dist  mag scaled_mag href ARSscale
#> 1        24.73 3.21        7.7 7.51    20.76
```

The other values calculated relate to the number of points in the data (`href`) and to the average foraging range (`med_max_dist`) estimated from the trips present in the data (`mag` and `scaled_mag`).

Then, we select a smoothing parameter value, based on our understanding of the species movement ecology, as well as our understanding of the management context within which these movements occur.

Using this smoothing value, we can run Kernel Density Estimation for each individual, with `estSpaceUse`. We need to specify the isopleth at which level we want to use utilisation distributions - this is by default set to 50, as the 50% utilisation distribution (where an animal spends about half of its time) is commonly used to define an animal's 'core range' (Lascelles et al. 2016).

``` r
KDEs <- estSpaceUse(
  DataGroup = trips, 
  Scale = Hvals$mag, 
  UDLev = 50, 
  polyOut = T)
#> Warning in estSpaceUse(DataGroup = trips, Scale = Hvals$mag, UDLev = 50, : No grid resolution ('Res') was specified, or the specified resolution was >99 km and therefore ignored.
#>                   Space use was calculated on a 500-cell grid, with cells of 0.727 square km
```

<img src="man/figures/README-estSpaceUse-1.png" width="100%" />

This gives us an estimate of the core areas in which each individual spends time while on foraging trips. At this step we should verify that the smoothing parameter value we selected is producing reasonable space use estimates, given what we know about our study animals.

The next step is to estimate to what degree this tracked sample is representative of the larger population. That is, how well does the variation in space use of these tracked individuals encapsulate variation in the wider population? To do this, we use the `repAssess` function. This function repeatedly samples a number of individual home ranges, averages them together, and quantifies how many locations from the unselected individuals fall within this combined home range area. This process is repeated for each sample size, and iterated a chosen number of time, from 1 to the 1 less than the total number of individuals in the study sample.

To speed up this procedure, we can supply the results of `estSpaceUse`, which will be randomly sampled and recombined in each iteration. We choose the number of times we want to re-sample at each sample size by setting the `Iteration` argument.

``` r
KDEs
#> $KDE.Surface
#> ********** Utilization distribution of several Animals ************
#> 
#> Type: probability density
#> Smoothing parameter estimated with a  specified smoothing parameter
#> This object is a list with one component per animal.
#> Each component is an object of class estUD
#> See estUD-class for more information
#> 
#> 
#> $UDPolygons
#> Simple feature collection with 41 features and 2 fields
#> geometry type:  MULTIPOLYGON
#> dimension:      XY
#> bbox:           xmin: -6.411805 ymin: -16.87262 xmax: -4.103352 ymax: -13.63024
#> epsg (SRID):    4326
#> proj4string:    +proj=longlat +datum=WGS84 +no_defs
#> First 10 features:
#>          id     area                       geometry
#> 69302 69302 477.2619 MULTIPOLYGON (((-5.691109 -...
#> 69304 69304 532.2372 MULTIPOLYGON (((-4.472789 -...
#> 69305 69305 305.7677 MULTIPOLYGON (((-5.704768 -...
#> 69306 69306 435.9407 MULTIPOLYGON (((-6.267572 -...
#> 69307 69307 424.6913 MULTIPOLYGON (((-5.920609 -...
#> 69308 69308 187.0009 MULTIPOLYGON (((-5.799668 -...
#> 69309 69309 187.6590 MULTIPOLYGON (((-5.806514 -...
#> 69310 69310 721.4156 MULTIPOLYGON (((-6.34918 -1...
#> 69311 69311 402.9683 MULTIPOLYGON (((-5.92123 -1...
#> 69312 69312 781.6377 MULTIPOLYGON (((-6.031113 -...
repr <- repAssess(trips, listKDE = KDEs$KDE.Surface, Scale = Hvals$mag, Iteration = 1, BootTable = F)
#> [1] "nls (non linear regression) successful, asymptote estimated for bootstrap sample."
```

The output is a table, with the estimated percentage of representativeness given in the `out` column.

The relationship between sample size and the inclusion of un-tested animals' space use areas is visualized via this automatic output plot, which is saved to the working directoty (i.e. `getwd()`). By quantifying this relationship across a range of different sample sizes, we can estimate how close we are to an information asymptote. Put another way, we estimate how much new space use information would be added by including more animals in the sample. In the case of this Brown Pelican dataset, we estimate that ~97% of the space used by this population is captured by the sample of 29 individuals. Highly representative!

<img src="man/figures/boobies_representativeness_h3_tripBuff3_N41.png" width="100%" height="10%" />

Using the space use estimates of each individual, we can now calculate where they overlap in space. Then, by including the representativeness value, we can estimate the proportion of the larger population using a given area and check whether this proportion meets threshold of importance at the population level. Here, if we have population size estimates, we can include this value to output a number of individuals aggregating in a given space instead of proportions, which can then use to compare against importance criteria (i.e KBA, EBSA criteria).

If you desire polygon output, instead of a raster surface, you can indicate this using the `polyOut` argument. This aggregates all cells with the same estimated number of individuals into to sinlgle polygons. If you stipulate `plotit=T` a plot will be produced (as below) showing all the polygons which meet the threshold proportion of birds using the area. If instead you just want the raster density distribution surface, simply use `plotit=F`.

``` r
KBAs <- findKBA(
  KDE.Surface = KDEs,
  Represent = repr$out,
  UDLev = 50,
  Col.size=500,
  polyOut = T,
  plotit = T)
```

<img src="man/figures/README-findKBA-1.png" width="100%" />

Since the output is in Simple Features (sf) visualizing these data is simple.

To visualize either the polygons or the density surface, use `plot(KBAs)`.
