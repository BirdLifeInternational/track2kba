# Patch change
This is a patch to deal with a DateTime related error in tripSplit(). In this version I have:

* Fixed errors related to midnight timestamps in POSIXct objects introduced to tripSplit() (i.e. tripSplit now accepts DateTime objects w/ midnight times).
* Attempted to remove maptools/rgdal/rgeos dependencies

# Version 1.1.0

## Test environments
* Windows (on local and GitHub Actions), R release
* Ubuntu 20.04 (on GitHub Actions), R release
* Ubuntu 20.04 (on GitHub Actions), devel
* macOS (on GitHub Actions), R release

## R CMD check results

0 errors | 0 warnings | 0 notes

There are currently no downstream dependencies for this package.