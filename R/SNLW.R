# Function to create old utah map given an input DEM
#' The 1992 Utah Snow Load Equations
#'
#' This function uses county-specific parameters stored internally in the
#' snowload package to make maps of the 1992 ground snow load equations
#' (referred to as SNLW).
#'
#' @param dem A spatial data frame containing the desired prediction locations.
#'   Note that non missing predictions are only provided for locations within the
#'   state of Utah.
#' @param elevVar The name of the variable in DEM corresponding to elevation.
#'   The elevation should be specified in meters or feet.
#' @param feet A logical indicating whether elevation is measured in meters
#'   or feet.
#'
#' @return A SpatialPointsDataFrame version of dem that includes a snlw
#'   prediction column.
#'
#' @export
snlwf = function(dem, elevVar = c("ELEVATION"), feet = FALSE){
  # Use some data internal to the snowload package. The first dataset is the
  # county shapefile for Utah and the second is a list of the coefficients
  # used for each county curve.

  if(!is.element(class(dem)[1], c("SpatialPointsDataFrame",
                                   "SpatialGridDataFrame",
                                   "SpatialPixelsDataFrame"))){
    stop("Prediction locations must be a spatial class from the sp package")
  }

  # Transform the input data to the same coordinate extent at the ut county
  # shapefile.
  dem <- as(dem, "SpatialPointsDataFrame")
  dem <- sp::spTransform(dem, sp::proj4string(utcounty))

  dem$COUNTYNBR <- as.numeric(as.character(sp::over(dem, utcounty)$COUNTYNBR))

  dem$COUNTYNBR[is.na(dem$COUNTYNBR)] <- nrow(snlw) + 1
  tdf <- data.frame(COUNTY = "None", P_0 = 0, S = 0, A_0 = 0)
  snlw <- rbind(snlw, tdf)

  dem$A_0 <- snlw$A_0[dem$COUNTYNBR]
  dem$P_0 <- snlw$P_0[dem$COUNTYNBR]
  dem$S <- snlw$S[dem$COUNTYNBR]

  # Note that the function assumes that the input data are in meters.
  adj <- 1/(1000*ifelse(feet, 1, .3048))
  dem$snlw <- sqrt(dem$P_0^2 + dem$S^2*((dem@data[, elevVar]*adj) - dem$A_0)^2)
  dem$snlw[dem@data[, elevVar]*adj <= dem$A_0] <-
    dem$P_0[dem@data[, elevVar]*adj <= dem$A_0]

  # Convert from psf to kpa
  dem$snlw <- dem$snlw*.04788

  # Remove the variables we no longer need now that predictions are
  # complete.
  dem$COUNTYNBR <- NULL
  dem$A_0 <- NULL
  dem$P_0 <- NULL
  dem$S <- NULL

  dem$snlw[dem$snlw == 0] <- NA

  return(dem)
}
