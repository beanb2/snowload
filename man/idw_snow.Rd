% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otherMethods.R
\name{idw_snow}
\alias{idw_snow}
\title{Inverse Distance Weighting, the Idaho way}
\usage{
idw_snow(formula = RESPONSE ~ ELEVATION, locations, newdata,
  tlayer = 1220, power = c(2, 6), NGSL = TRUE,
  bound_output = FALSE, corMethod = "spearman", print = FALSE, ...)
}
\arguments{
\item{formula}{A formula that is passed to the idw function in the
sp package.}

\item{locations}{The measurement locations to be used for interpolation.
An object of class SpatialPointsDataFrame}

\item{newdata}{The prediction locations. An objet of class
SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.}

\item{tlayer}{The elevation (in meters) on which to split the data. The
Idaho snow load report splits the state into two layers using
4000 ft (1219.2 meters) as a cutoff. This variable can also be set
to NA, at which point an automatic separating elevation method is
applied (see details).}

\item{power}{A vector of length 2 indicating the exponent in the inverse
distance weighting for the lower and upper elevations respectively.}

\item{NGSL}{If true, used normalized ground snow loads, which are simply
snow loads divided by elevation. This is the way elevation is accounted
for in the 2015 Idaho study.}

\item{bound_output}{If true, the final predictions are capped at the range
of observed data.}

\item{corMethod}{(only relevant if tlayer = NA). The method by which
correlations are determined when deciding on a separating elevation.
Valid options include "pearson", "spearman" and "kendall"}

\item{print}{If true, gstat package messages are printed to the screen, as well as the
separating elevation.}
}
\value{
An appended version of the newdata class containing a column named
  "idw_snow".
}
\description{
A function that implements inverse distance weighting using
normalized ground snow loads.
}
\details{
One of the great challenges of this method is the proper selection
of a separating elevation.
}
