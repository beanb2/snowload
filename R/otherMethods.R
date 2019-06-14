#' Inverse Distance Weighting, the Idaho way
#'
#' A function that implements inverse distance weighting using
#' normalized ground snow loads.
#'
#' @param formula A formula that is passed to the idw function in the
#'   sp package.
#' @param locations The measurement locations to be used for interpolation.
#'   An object of class SpatialPointsDataFrame
#' @param newdata The prediction locations. An objet of class
#'   SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
#' @param tlayer The elevation (in meters) on which to split the data. The
#'   Idaho snow load report splits the state into two layers using
#'   4000 ft (1219.2 meters) as a cutoff.
#' @param power A vector of length 2 indicating the exponent in the inverse
#'   distance weighting for the lower and upper elevations respectively.
#' @param NGSL If true, used normalized ground snow loads, which are simply
#'   snow loads divided by elevation. This is the way elevation is accounted
#'   for in the 2015 Idaho study.
#' @param bound_output If true, the final predictions are capped at the range
#'   of observed data.
#'
#' @return An appended version of the newdata class containing a column named
#'   "idw_snow".
#'
#' @export
idw_snow = function(variables = c("RESPONSE", "ELEVATION"),
                    locations, newdata, tlayer = 1220, power = c(2, 6),
                    NGSL = TRUE, bound_output = FALSE, ...){
  # Run generic tests of model inputs
  input_check(locations, newdata)

  # Use the Normalized ground snow loads
  if(NGSL){
    # Prevent weird things from happening when elevation
    # is very close to zero or negative.
    locations@data[as.vector(locations[, variables[2]]@data < 1), variables[2]] <- 1
    newdata@data[as.vector(newdata[, variables[2]]@data < 1), variables[2]] <- 1

    locations@data[, variables[1]] <-
      locations@data[, variables[1]]/locations@data[, variables[2]]
  }

  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method)
  locations.low = locations[as.vector(locations[, variables[2]]@data < tlayer), ]
  locations.high = locations[as.vector(locations[, variables[2]]@data >= tlayer),]

  # Create a formula from the provided variables.
  tform <- as.formula(paste(variables[1], "1", sep = "~"))
  # Obtain the bi-layer results
  if(nrow(locations.low) > 0){
    temp2 = gstat::idw(tform, locations = locations.low, newdata = newdata,
                       idp = power[1], ...)$var1.pred
  }else{
    temp2 = NULL
    # print("No stations exist in the lower layer, using only one layer with exponent c2...",
    # call. = FALSE)
  }
  if(nrow(locations.high) > 0){
    temp6 = gstat::idw(tform, locations = locations.high, newdata = newdata,
                       idp = power[2], ...)$var1.pred
  }else{
    temp6 = NULL
  }

  # Store results by layer, assuming both layers are defined.
  results = vector("numeric", nrow(newdata))
  if(is.null(temp2)){
    results = temp6
  }else if(is.null(temp6)){
    results = temp2
  }else{
    results[newdata@data[, variables[2]] < tlayer] = temp2[newdata@data[, variables[2]] < tlayer]
    results[newdata@data[, variables[2]] >= tlayer] = temp6[newdata@data[, variables[2]] >= tlayer]
  }

  # Return ground snow load predictions as opposed to simply the NGSL.
  if(NGSL){results = results*newdata@data[, variables[2]]}

  # If requested restrict the response values to the highest valued snow load in the dataset.
  if(bound_output){
    bound_out <- c(min(locations@data[, variables[1]], na.rm = TRUE),
                   max(locations@data[, variables[1]], na.rm = TRUE))

    results[results < bound_out[1]] <- bound_out[1]
    results[results > bound_out[2]] <- bound_out[2]
  }

  newdata$idw_snow <- results
  return(newdata)
} # End the function


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


#' Linear Triangulation Interpolation
#'
#' A function that implements a linear triangulation interpolation.
#' The function is a wrapper to appropriate functions in the
#' akima package.
#'
#' @param variables A vector of length two specifying the column names
#'   of the response and elevation variables respectively.
#' @param locations The measurement locations to be used for interpolation.
#'   An object of class SpatialPointsDataFrame
#' @param newdata The prediction locations. An objet of class
#'   SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
#' @param density The density of the grid created by the interpolation.
#'   The first argument specifies the x density and the second the y density.
#' @param NGSL If true, normalized ground snow loads (which are snow loads
#'   divided by elevation) are used in prediction.
#' @param bound_ouput If TRUE, the final ouputs are restricted to the observed
#'   outputs at the measurement locations.
#'
#' @return An appended version of newdata that includes a column named
#'   "tri_snow" with the predictions.
#'
#' @export
tri_snow = function(variables = c("RESPONSE", "ELEVATION"),
                    locations, newdata, density = c(500, 500),
                    NGSL = TRUE, bound_output = FALSE, ...){
  # Run generic tests on model inputs (see checks.R)
  input_check(locations, newdata)

  # Determine if modeing raw snow loads or the NGSL
  if(NGSL){
    locations$NGSL = locations@data[,variables[1]]/locations@data[,variables[2]]
    # See Akima package for details behind this function
    temppredict = akima::interp(x = locations, z = "NGSL", linear = TRUE,
                                nx = density[1], ny = density[2])
  }else{
    # See Akima package for details behind this function
    temppredict = akima::interp(x = locations, z = variables[1], linear = TRUE,
                                nx = density[1], ny = density[2])
  }

  # One consideration is to make sure that the rows and columns line up
  # properly when creating the SpatialPixelsDataFrame. The determination
  # was made by taking a map with a known pattern and making sure that it
  # was plotted properly.
  temppredict <- data.frame(expand.grid(x = temppredict$x, y = temppredict$y),
                            z = as.vector(temppredict$z))

  # Reformat the output to be a spatial pixels data frame (should be this automatically)
  sp::coordinates(temppredict) <- c("x", "y")
  sp::proj4string(temppredict) <- sp::proj4string(locations)
  sp::gridded(temppredict) <- TRUE

  # Convert predictions to raster and extract raster information for the newdata DEM
  temppredict = raster::raster(temppredict)
  preds = raster::extract(temppredict, newdata)

  if(NGSL){
    preds = preds*newdata@data[, variables[2]]
  }

  if(bound_output){
    bound_out <- c(min(locations@data[, variables[1]], na.rm = TRUE),
                   max(locations@data[, variables[1]], na.rm = TRUE))

    preds[preds < bound_out[1]] <- bound_out[1]
    preds[preds > bound_out[2]] <- bound_out[2]
  }

  newdata$tri_snow <- preds

  return(newdata)
}


#' Linear regression for ground snow load prediction.
#'
#' This function is nothing more than a wrapper to the standard linear
#' models with inputs designed to look and feel similar to other
#' functions in the snowload package.
#'
#' @param formula A forumla object passed directly to the linear model.
#' @param locations An object of class SpatialPointsDataFrame. Specifies
#'   the measurement locations.
#' @param newdata A spatial data object. Specifies the prediction
#'   locations.
#' @param bound_elevation If TRUE, the elevations of the newdata locations
#'   are restricted to the range of observed data in the locations object.
#' @param bound_ouput If TRUE, the final ouputs are restricted to the observed
#'   outputs at the measurement locations.
#'
#' @return An appended newdata object that includes a "lm_snow" column
#'   with the model predictions. Note that if bound_elevation = TRUE then
#'   the elevations in this object may be altered from the original.
#'
#' @export
lm_snow = function(formula = log(RESPONSE) ~ ELEVATION,
                  locations, newdata,
                  bound_elevation = FALSE, bound_output = FALSE){
  # Run generic tests on model inputs (see checks.R)
  input_check(locations, newdata)

  vars <- all.vars(formula)

  # Restrct elevation trend estimates to be no greater than the trend value
  # of the highest elevation station location.
  if(bound_elevation){
    bound_el <- c(min(locations@data[, vars[2]], na.rm = TRUE),
                  max(locations@data[, vars[2]], na.rm = TRUE))

    newdata@data[,vars[2]][newdata@data[, vars[2]] > bound_el[2]] <-
      bound_el[2]
    newdata@data[,vars[2]][newdata@data[, vars[2]] < bound_el[1]] <-
      bound_el[1]
  }

  model <- lm(formula, data = as.data.frame(locations))

  preds = predict(model, as.data.frame(newdata))

  if(bound_output){
    temp_out <- lm(formula, data = as.data.frame(locations))
    temp_out <- temp_out$residuals + temp_out$fitted.values

    bound_out <- c(min(temp_out, na.rm = TRUE),
                   max(temp_out, na.rm = TRUE))

    preds[preds < bound_out[1]] <- bound_out[1]
    preds[preds > bound_out[2]] <- bound_out[2]
  }

  newdata$lm_snow <- preds

  return(newdata)
}


