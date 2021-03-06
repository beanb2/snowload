% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism.R
\name{prism}
\alias{prism}
\title{Title: PRISM}
\usage{
prism(formula = log(RESPONSE) ~ ELEVATION, locations, newdata,
  distImp = 0.8, minRad = 10, wdistance = 2, welevRange = list(lwr
  = 200, upr = 2500), welevation = c("ELEVATION", 1), wbasin = c(NA,
  3), weco = c(NA, 3), bound = FALSE, weights = FALSE, ...)
}
\arguments{
\item{locations}{The measurement locations, of class SpatialPointsDataFrame,
containing all the necessary colunmns to run the PRISM algorithm.}

\item{newdata}{The prediction locations, typically of class
SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
This object will be appended and returned by the prism function. Note that
the projection of this object must match the projection of the locations
variable.}

\item{distImp}{A weight controlling the distance importance factor relative
to the elevation importance factor. Must be a number between 0 and 1.}

\item{minRad}{The minimum radius of influence, defined in the units of
the dataset projection.}

\item{wdistance}{An exponent controlling the severity of the distance
weight.}

\item{welevRange}{The range of elevations considered in the elevation
weight. Differences less than the minimum receive weights of 1 while
differences greater than the maximum receive weights of 0.}

\item{welevation}{A vector of length 2. The first argument specifies the
name of the column where elevation information can be obtained. The second
argument specifies the exponent that controls the severity of the elevation
weight.}

\item{wbasin}{A vector of length 2. The first argument specifies the name of
column where water basin information can be obtained. The second argument
specifies the exponent that controls the severity of the water basin
weight. If the first argument is NA, then this weight is ignored.}

\item{weco}{Same structure as wbasin, but specifies weights based
on the EPA's ecoregions.}

\item{bound}{If true, predictions are not allowed to extrapolate beyond the
observed range of values in the input data. Thus predictions are capped on
both ends by the most extreme observed values. If false, predictions are
allowed to extrapolate.}

\item{weights}{If true, return the slope, intercept, and weights for each
of the predictions. This argument is only relevant if one wishes to
plot the PRISM predictions.}

\item{variables}{The variables used in the PRISM regression. The first variable
is the reponse and all remaining variables are treated as explanatory. This
is a replacement for a formula call as the use of lsfit greatly speeds the
computations.}

\item{score}{Additional arguments to}
}
\value{
If weights = FALSE, a numeric vector of predictions with length equal
  to the number of rows in newdata. Else, a list with the first component
  containing the prediction, intercept, and slope, and the second component
  containing a matrix of weights.
}
\description{
A function that implements the PRISM algorithm as developed by Daly et.al.
The PRISM model is weighted least squares regression primarily based on
similarities in location and elevation to the prediction location.
Note that a unique linear model is fit at each and every predictio location.
}
