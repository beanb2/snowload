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
#'   4000 ft (1219.2 meters) as a cutoff. This variable can also be set
#'   to NA, at which point an automatic separating elevation method is
#'   applied (see details).
#' @param power A vector of length 2 indicating the exponent in the inverse
#'   distance weighting for the lower and upper elevations respectively.
#' @param NGSL If true, used normalized ground snow loads, which are simply
#'   snow loads divided by elevation. This is the way elevation is accounted
#'   for in the 2015 Idaho study.
#' @param bound_output If true, the final predictions are capped at the range
#'   of observed data.
#' @param corMethod (only relevant if tlayer = NA). The method by which
#'   correlations are determined when deciding on a separating elevation.
#'   Valid options include "pearson", "spearman" and "kendall"
#' @param print If true, gstat package messages are printed to the screen, as well as the
#'   separating elevation.
#'
#' @return An appended version of the newdata class containing a column named
#'   "idw_snow".
#'
#' @details
#' One of the great challenges of this method is the proper selection
#' of a separating elevation.
#'
#' @export
idw_snow = function(formula = RESPONSE ~ ELEVATION,
                    locations, newdata, tlayer = 1220, power = c(2, 6),
                    NGSL = TRUE, bound_output = FALSE,
                    corMethod = "spearman", print = FALSE, ...){
  # Run generic tests of model inputs
  input_check(locations, newdata)

  modeldf <- stats::model.frame(formula, locations)
  vars <- all.vars(formula)

  # Use the Normalized ground snow loads
  if(NGSL){
    # NGSL can only be performed if a single explanatory variable
    # is provided (representing elevation)
    if(ncol(modeldf) != 2){
      stop("Exactly one response variable and one explanatory variable
           (corresponding to elevation) must be supplied when NGSL = TRUE")
    }

    # Prevent weird things from happening when elevation
    # is very close to zero or negative.
    modeldf[modeldf[, 2] < 1, 2] <- 1

    newdata$newelev <- newdata@data[, vars[2]]
    newdata$newelev[newdata$newelev < 1] <- 1

    locations$NGSL <- modeldf[, 1]/modeldf[, 2]
    formula = NGSL ~ 1

    # If no value for tlayer is supplied, automatically fit one that minimizes
    # a weighted average of absolute correlations between the upper and lower
    # layers.
    if(is.na(tlayer)){
      range = quantile(locations@data[, vars[2]], seq(0.25, 0.75, .1),
                       na.rm = TRUE)

      # Append the range to include all stations in either layer.
      range <- c(range, 1e16)

      # Cycle through the range and determine which separating elevation
      # is best based on percentages of 10 in the interquartile range.
      results <- matrix(0, ncol = 2, nrow = length(range))
      for(i in 1:length(range)){
        results[i, 1] <- range[i]
        results[i, 2] <- checkCor(locations, range[i], vars, corMethod)
      }

      # Set the separating elevation to this value.
      tlayer <- results[which.min(results[, 2]), 1]

      remove(results)
      if(print){print(paste("Separating elevation:", tlayer))}
    }
  }else{
    if(is.na(tlayer)){
      warning("Automatic elevation separating is only relevant when NGSL = TRUE.
              Placing all stations in the lower layer...")
      tlayer = 0
    }
  }

  # If print is requested, set debug.level = 1
  debug.level <- ifelse(print, 1, 0)

  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method)
  locations.low = locations[as.vector(locations[, vars[2]]@data < tlayer), ]
  locations.high = locations[as.vector(locations[, vars[2]]@data >= tlayer),]

  # Obtain the bi-layer results
  if(nrow(locations.low) > 0){
    temp2 = gstat::idw(formula, locations = locations.low, newdata = newdata,
                       idp = power[1], debug.level = debug.level, ...)$var1.pred
  }else{
    temp2 = NULL
    # print("No stations exist in the lower layer, using only one layer with exponent c2...",
    # call. = FALSE)
  }
  if(nrow(locations.high) > 0){
    temp6 = gstat::idw(formula, locations = locations.high, newdata = newdata,
                       idp = power[2], debug.level = debug.level, ...)$var1.pred
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
    results[newdata@data[, vars[2]] < tlayer] = temp2[newdata@data[, vars[2]] < tlayer]
    results[newdata@data[, vars[2]] >= tlayer] = temp6[newdata@data[, vars[2]] >= tlayer]
  }

  # Return ground snow load predictions as opposed to simply the NGSL.
  if(NGSL){results = results*newdata$newelev}

  # If requested restrict the response values to the highest valued snow load in the dataset.
  if(bound_output){
    bound_out <- c(min(modeldf[, 1], na.rm = TRUE),
                   max(modeldf[, 2], na.rm = TRUE))

    results[results < bound_out[1]] <- bound_out[1]
    results[results > bound_out[2]] <- bound_out[2]
  }

  return(results)
} # End the function

#
#
#
# Source code for the checkCor function that is called by idw_snow.
# tdata - (usually locations) to check for correlations
# vars - (passed from parent functions) the names of the variables
#   specified in the formula
# method - measure of correlation, either pearson, spearman,
#  or kendall
checkCor <- function(tdata, tlayer, vars, method){
  high <- tdata[tdata@data[, vars[2]] > tlayer, ]
  low <- tdata[tdata@data[, vars[2]] <= tlayer, ]

  lowcor <- stats::cor(low@data[, vars[2]], low$NGSL, method = method,
                       use = "na.or.complete")
  highcor <- stats::cor(high@data[, vars[2]], high$NGSL, method = method,
                        use = "na.or.complete")

  # If the correlation yields a missing value (due to a lack of data),
  # return a 0.
  if(is.na(lowcor)){lowcor <- 0}
  if(is.na(highcor)){highcor <- 0}

  # Return a weighted average of the absolute correlations based on elevation
  return( (nrow(high)*abs(highcor) + nrow(low)*abs(lowcor)) / nrow(tdata) )
}
#
#
#


# Function to create old utah map given an input DEM
#' The 1992 Utah Snow Load Equations
#'
#' This function uses county-specific parameters stored internally in the
#' snowload package to make maps of the 1992 ground snow load equations
#' (referred to as SNLW).
#'
#' @param A forumla of the form ~ELEVATION. The function assumes that the last
#'   given variable corresponds to the elevation column and all other arguments
#'   are ignored.
#' @param newdata A spatial data object containing the desired prediction locations.
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
snlwf = function(formula = ~ELEVATION, newdata, feet = FALSE){
  # Use some data internal to the snowload package. The first dataset is the
  # county shapefile for Utah and the second is a list of the coefficients
  # used for each county curve.

  if(!is.element(class(newdata)[1], c("SpatialPointsDataFrame",
                                      "SpatialGridDataFrame",
                                      "SpatialPixelsDataFrame"))){
    stop("Prediction locations must be a spatial class from the sp package")
  }

  # Transform the input data to the same coordinate extent at the ut county
  # shapefile.
  newdata <- as(newdata, "SpatialPointsDataFrame")
  newdata <- sp::spTransform(newdata, sp::proj4string(utcounty))

  newdata$COUNTYNBR <- as.numeric(as.character(sp::over(newdata, utcounty)$COUNTYNBR))

  newdata$COUNTYNBR[is.na(newdata$COUNTYNBR)] <- nrow(snlw) + 1
  tdf <- data.frame(COUNTY = "None", P_0 = 0, S = 0, A_0 = 0)
  snlw <- rbind(snlw, tdf)

  newdata$A_0 <- snlw$A_0[newdata$COUNTYNBR]
  newdata$P_0 <- snlw$P_0[newdata$COUNTYNBR]
  newdata$S <- snlw$S[newdata$COUNTYNBR]

  # Extract the column associated with elevation
  modeldf <- stats::model.frame(formula, newdata)
  modeldf <- modeldf[,ncol(modeldf)]

  # Note that the function assumes by default that the input data are in meters.
  adj <- 1/(1000*ifelse(feet, 1, .3048))

  snlw <- sqrt(newdata$P_0^2 + newdata$S^2*((modeldf*adj) - newdata$A_0)^2)

  snlw[modeldf*adj <= newdata$A_0] <- newdata$P_0[modeldf*adj <= newdata$A_0]

  # Convert from psf to kpa
  snlw <- snlw*.04788

  # 0 valued predictions are out of state, set these to na
  snlw[snlw == 0] <- NA

  return(snlw)
}


#' Linear Triangulation Interpolation
#'
#' A function that implements a linear triangulation interpolation.
#' The function is a wrapper to appropriate functions in the
#' akima package.
#'
#' @param A forumla argument specifying the response and explanatory variables.
#'   This formula should only be of the form RESPONSE ~ ELEVATION. No other
#'   formula types are allowed.
#' @param locations The measurement locations to be used for interpolation.
#'   An object of class SpatialPointsDataFrame
#' @param newdata The prediction locations. An objet of class
#'   SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
#' @param density The density of the grid created by the interpolation.
#'   The first argument specifies the x density and the second the y density.
#' @param bound_ouput If TRUE, the final ouputs are restricted to the observed
#'   outputs at the measurement locations.
#'
#' @return A numeric vector of predictions with length equal to the number
#'   of rows in newdata.
#'
#' @details Note that this function does not use the formula call directly.
#'   Rather, if an explanatory variable is present, it is assumed to be
#'   the variable name for elevation and a normalized response, which is
#'   simply response/elevation, is used.
#'
#' @export
tri_snow = function(formula = RESPONSE ~ ELEVATION,
                    locations, newdata, density = c(500, 500),
                    bound_output = FALSE, ...){
  # Run generic tests on model inputs (see checks.R)
  input_check(locations, newdata)

  # Obtain the desired model inputs.
  modeldf <- stats::model.frame(formula, locations)
  vars <- all.vars(formula)

  if(ncol(modeldf) > 2){
    stop("This function only allows one explanatory variable, which is elevation.")
  }

  # Determine if modeing raw snow loads or the NGSL
  if(ncol(modeldf) > 1){
    locations$NGSL = modeldf[, 1]/modeldf[, 2]
    # See Akima package for details behind this function
    temppredict = akima::interp(x = locations, z = "NGSL", linear = TRUE,
                                nx = density[1], ny = density[2])
  }else{
    # See Akima package for details behind this function
    temppredict = akima::interp(x = locations, z = vars[1], linear = TRUE,
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

  if(ncol(modeldf) > 1){
    preds = preds*newdata@data[, vars[length(vars)]]
  }

  if(bound_output){
    bound_out <- c(min(modeldf[, 1], na.rm = TRUE),
                   max(modeldf[, 1], na.rm = TRUE))

    preds[preds < bound_out[1]] <- bound_out[1]
    preds[preds > bound_out[2]] <- bound_out[2]
  }

  return(preds)
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
#' @return A numeric vector of predictions with length equal to the number
#'   of rows in newdata.
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

  model <- lm(formula, data = locations)

  preds = predict(model, newdata)

  if(bound_output){
    temp_out <- model$residuals + model$fitted.values

    bound_out <- c(min(temp_out, na.rm = TRUE),
                   max(temp_out, na.rm = TRUE))

    preds[preds < bound_out[1]] <- bound_out[1]
    preds[preds > bound_out[2]] <- bound_out[2]
  }

  return(preds)
}


