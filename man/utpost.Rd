% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{utpost}
\alias{utpost}
\title{Post office locations in Utah}
\format{A data frame with 253 rows and 7 variables:

 \describe{
 \item{CITY}{The character name of the post office city.}
 \item{COUNTY}{The character name of the county.}
 \item{STATE}{The two-character state abbrevation.}
 \item{LONGITUDE}{Longitude coordinate position.}
 \item{LATITUDE}{Latitude coordinate position.}
 \item{ELEVATION}{Elevation of the measurement location (meters).}
 \item{HUC}{The United States Geological Survey 8 digit
 hydrologic unit code.}
 }}
\usage{
utpost
}
\description{
These 253 post office locations provide the ideal set
of prediction locations to test the snow load functions.
These locations are also the same locations comprising the
city tables provided in the 2018 Utah Snow Load Study.
}
\references{
\insertRef{Bean2018-report}{snowload}
}
\keyword{datasets}
