#' Title: PRISM
#'
#' A function that implements the PRISM algorithm as developed by Daly et.al.
#' The PRISM model is weighted least squares regression primarily based on
#' similarities in location and elevation to the prediction location.
#' Note that a unique linear model is fit at each and every predictio location.
#'
#' @param variables The variables used in the PRISM regression. The first variable
#'   is the reponse and all remaining variables are treated as explanatory. This
#'   is a replacement for a formula call as the use of lsfit greatly speeds the
#'   computations.
#' @param locations The measurement locations, of class SpatialPointsDataFrame,
#'   containing all the necessary colunmns to run the PRISM algorithm.
#' @param newdata The prediction locations, typically of class
#'   SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
#'   This object will be appended and returned by the prism function. Note that
#'   the projection of this object must match the projection of the locations
#'   variable.
#' @param distImp A weight controlling the distance importance factor relative
#'   to the elevation importance factor. Must be a number between 0 and 1.
#' @param minRad The minimum radius of influence, defined in the units of
#'   the dataset projection.
#' @param wdistance An exponent controlling the severity of the distance
#'   weight.
#' @param welevRange The range of elevations considered in the elevation
#'   weight. Differences less than the minimum receive weights of 1 while
#'   differences greater than the maximum receive weights of 0.
#' @param welevation A vector of length 2. The first argument specifies the
#'   name of the column where elevation information can be obtained. The second
#'   argument specifies the exponent that controls the severity of the elevation
#'   weight.
#' @param wbasin A vector of length 2. The first argument specifies the name of
#'   column where water basin information can be obtained. The second argument
#'   specifies the exponent that controls the severity of the water basin
#'   weight. If the first argument is NA, then this weight is ignored.
#' @param weco Same structure as wbasin, but specifies weights based
#'   on the EPA's ecoregions.
#' @param bound If true, predictions are not allowed to extrapolate beyond the
#'   observed range of values in the input data. Thus predictions are capped on
#'   both ends by the most extreme observed values. If false, predictions are
#'   allowed to extrapolate.
#' @param weights If true, return the slope, intercept, and weights for each
#'   of the predictions. This argument is only relevant if one wishes to
#'   plot the PRISM predictions.
#' @param score Additional arguments to
#'
#' @return A numeric vector of predictions with length equal to the number
#'   of rows in newdata.
#'
#' @export
prism = function(formula = log(RESPONSE) ~ ELEVATION,
                 locations,
                 newdata,
                 distImp = 0.8,
                 minRad = 10,
                 wdistance = 2,
                 welevRange = list(lwr = 200, upr = 2500),
                 welevation = c("ELEVATION", 1),
                 wbasin = c(NA, 3),
                 weco = c(NA, 3),
                 bound = FALSE,
                 weights = FALSE,
                 ...){

  # Run generic tests of model inputs
  input_check(locations, newdata)

  # Ensure that the arguments for the weights are properly specified
  tl <- c(length(welevation)-1, length(wbasin)-1, length(weco)-1,
          length(welevRange[[1]]), length(welevRange[[2]]))

  # Throw a warning if not enough weight arguments are supplied.
  if(any(tl < 1)){
    stop("Each weight argument must have one numeric weight supplied.")
  }

  # If more than one argument is supplied for any of the parameters,
  # call the tuning function.
  if(any(tl > 1)){
    tune_results <- prism_tune(formula = formula,
                               locations = locations,
                               distImp = distImp, minRad = minRad,
                               wdistance = wdistance, welevRange = welevRange,
                               welevation = welevation, wbasin = wbasin,
                               weco = weco, bound = FALSE, ...)

    tune_results <- as.vector(tune_results[1, ])
  }else{
    tune_results <- c(distImp, minRad, wdistance, welevRange[[1]],
                      welevRange[[2]], as.numeric(welevation[2]),
                      as.numeric(wbasin[2]), as.numeric(weco[2]))
  }
  # We store the final paramters (tune or not) in the tune_results variable.

  # Make sure that the required columns are in the dataset.

  # Ensure station data meet criteria.
  if(!(welevation[1] %in% colnames(locations@data) &
       welevation[1] %in% colnames(newdata@data))){
    stop("Specified elevation column not found in the data.")
  }

  if(!is.na(wbasin[1])){
    if(!(wbasin[1] %in% colnames(locations@data) &
         wbasin[1] %in% colnames(newdata@data))){
      stop("Specified basin column not found in the data.")
    }
  }

  if(!is.na(weco[1])){
    if(!(weco[1] %in% colnames(locations@data) &
         weco[1] %in% colnames(newdata@data))){
      stop("Specified eco-region column not found in the data.")
    }
  }

  # Find the distance and elevation weights (required)
  wd <- getDistWeight(locations, newdata, tune_results[2],
                      tune_results[3])
  wz <- getElevWeight(unlist(locations[, welevation[1]]@data),
                      unlist(newdata[, welevation[1]]@data),
                      tune_results[6], c(tune_results[4],
                                         tune_results[5]))

  wc <- getClusterWeight(locations, tune_results[2], tune_results[4])

  # Get the basin weight (if requested)
  if(!is.na(wbasin[1])){
    wb <- getHUCWeight(unlist(locations[, wbasin[1]]@data),
                       unlist(newdata[, wbasin[1]]@data),
                       tune_results[7])
  }else{
    wb <- 1
  }

  # Get the eco weight (if requested)
  if(!is.na(weco[1])){
    wf <- getEcoWeight(unlist(locations[, weco[1]]@data),
                       unlist(newdata[, weco[1]]@data),
                       tune_results[8])
  }else{
    wf <- 1
  }

  # Apply weighting equation from Daley et.al. 2008
  weightV <- wc*sqrt(tune_results[1]*wd^2 + (1-tune_results[1])*wz^2)*wb*wf

  # Scale the weights to sum to 1
  for(i in 1:ncol(weightV)){
    weightV[, i] <- weightV[, i]/sum(weightV[, i])
  }

  # Remove stuff we no longer need to save on memory
  remove(wc, wd, wz, wb, wf)
  ##################################

  # Preallocate for predictions
  modeldf <- stats::model.frame(formula, locations)
  vars = all.vars(formula)
  explanatory <- as.vector(modeldf[, -1])
  response <-  as.matrix(modeldf[, 1])

  tempReg <- matrix(0, nrow = length(newdata), ncol = ncol(modeldf))

  predictP <- as.matrix(newdata@data[, vars[2:length(vars)]])
  predictP <- cbind(1, predictP)

  for(i in 1:nrow(tempReg)){
    # Run linear regression on individual gridcell,
    # with the respective weighting vector
    if(any(is.na(predictP[i, ]))){
      tempReg[i, ] <- NA
    }else{
      tempReg[i, ] <- unname(lsfit(explanatory, response,
                                   wt = weightV[, i])$coef)
    }
  }

  predictP <- apply(predictP*tempReg, 1, sum)

  # If requested, do not let predictions extrapolate beyond readings
  if(bound){
    # Fit an unweighted linear model and extract the bounds of the data
    bound_out <- c(min(response, na.rm = TRUE), max(response, na.rm = TRUE))

    predictP[predictP < bound_out[1]] <- bound_out[1]
    predictP[predictP > bound_out[2]] <- bound_out[2]
  }

  # If weights are requested (for visualization) then return a list.
  # Otherwise, return the appended spatial data frame.
  if(weights){
    tdf <- data.frame(prism = predictP,
                      slope = tempReg[, 2],
                      intercept = tempReg[, 1])
    return(list(tdf, weightV))
  }else{
    return(predictP)
  }
}


# This file contains all model tuning functions.

#' Function to tune PRISM parameters using cross validation.
#'
#' This function conducts a grid search over the given list of
#' parameters and selects parameters that minimize the cross validated
#' error according to the provided scoring method.
#'
#'
#' @export
prism_tune <- function(formula = log(RESPONSE) ~ ELEVATION,
                       locations,
                       distImp = 0.8,
                       minRad = 10,
                       wdistance = 2,
                       welevRange = list(lwr = 200, upr = 2500),
                       welevation = c("ELEVATION", 1),
                       wbasin = c(NA, 3),
                       weco = c(NA, 3),
                       bound = FALSE,
                       tune_folds = 10,
                       tune_score = "MAE"){

  # Determine the number of combinations in the search.
  tlength <- length(distImp)*length(minRad)*length(wdistance)*
    length(welevRange[[1]])*length(welevRange[[2]])*
    (length(welevation)-1)*(length(wbasin)-1)*
    (length(weco)-1)

  # Create a matrix to store the results.
  tdf <- matrix(0, nrow = tlength, ncol = 9)
  colnames(tdf) <- c("distImp", "minRad", "wdistance", "lwr", "upr",
                     "welevation", "wbasin", "weco", "score")
  count <- 1
  for(w1 in as.numeric(distImp)){
    for(w2 in as.numeric(minRad)){
      for(w3 in as.numeric(wdistance)){
        for(w4 in as.numeric(welevRange[[1]])){
          for(w5 in as.numeric(welevRange[[2]])){
            for(w6 in as.numeric(welevation[-1])){
              for(w7 in as.numeric(wbasin[-1])){
                for(w8 in as.numeric(weco[-1])){
                  acc <- crv_pred(locations, formula,
                                  fun = "prism",
                                  distImp = w1, minRad = w2, wdistance = w3,
                                  welevRange = list(lwr = w4, upr = w5),
                                  welevation = c(welevation[1], w6),
                                  wbasin = c(wbasin[1], w7),
                                  weco = c(weco[1], w8),
                                  folds = tune_folds,
                                  score = tune_score)

                  if(length(acc) > 1){
                    stop("Model score length is greater than 1, make sure you specified a
                    proper scoring method.")
                  }

                  # Store the results.
                  tdf[count, ] <- c(w1, w2, w3, w4, w5, w6, w7, w8, acc)

                  # Advance
                  count <- count + 1
                } # end w8
              } # end w7
            } # end w6
          } # end w5
        } # end w4
      } # end w3
    } # end w2
  } # end w1

  tdf <- tdf[order(tdf[, 9]), ]
  return(tdf)
}


#
#
# These functions are used internally by the prism() function and
# will not be exported.
#
#

## Distance Weight Function ##
#=============================================================================
# Function determines the distance weight (W_d) in the PRISM algorithm.
# Distance is calculated as either geographical distance or euclidean distance
# depending on the preference of the user.
# gridC - Locations of predictions
# stationC - Locations of stations
# rm - minimum radius of influence
# a - distance weighting exponent.
#=============================================================================
getDistWeight = function(locations, newdata, rm, a){

  # Calculate distance with function from "fields" package.
  # Note that we transpose the results so weights are found in the columns,
  # not the rows.
  # Avoid singularity issues in getDistWeight by rounding.
  if(regexpr(sp::proj4string(locations), pattern = "+proj=longlat") > 1){
    dist <- round(fields::rdist.earth(locations@coords, newdata@coords, miles = FALSE),
                  digits = 4)
  }else{
    dist <- round(fields::rdist(locations@coords, newdata@coords), digits = 4)
  }

  # Set all stations within the minimum radius of influence to a weight of 1
  distWeight = dist - rm
  distWeight[distWeight <= 0] = 1

  # For all other stations, compute the weighting scheme.
  distWeight = distWeight^a
  distWeight = 1 / distWeight
  distWeight[distWeight > 1] = 1 # Set all values greater than 1 to the max

  # Normalize each weighting vector to sum to 1
  #distWeight = scale(distWeight, center = FALSE, scale = colSums(distWeight))
  for(i in 1:ncol(distWeight)){
    distWeight[, i] <- distWeight[, i]/sum(distWeight[, i])
  }

  return(distWeight)
}
#=============================================================================

## Elevation Weight Function ##
#=============================================================================
# Function determines the elevation weight (W_z) in the PRISM algorithm.
# It is assumed that elevation is measured in meters.
# gridElv - Elevations of predictions
# stationElv - Elevations of stations
# zm - minimum elevation precision (within this value = 1)
# zx - maximum elevation precision (outside this value = 0)
# b - elevation weighting exponent.
#=============================================================================
getElevWeight = function(stationElev, gridElev, b, elevRange){

  # Store elevation differences between each gridcell and the stations
  # in the columns of a matrix
  tempDiff = abs(outer(stationElev, gridElev, "-"))

  # Determine which elevation differences lie between the maximum and minimum
  # elevation of influence range
  tempDiffMin = tempDiff - elevRange[1]
  tempDiffMax = tempDiff - elevRange[2]

  # Apply weighting criterion to tempDiff based on min/max values
  tempDiff[tempDiffMin <= 0] = elevRange[1]
  tempDiff = tempDiff^b
  tempDiff = 1/tempDiff
  tempDiff[tempDiffMax >= 0] = 0

  # Scale elevation weights to sum to unity
  for(i in 1:ncol(tempDiff)){
    tempDiff[, i] <- tempDiff[, i]/sum(tempDiff[, i])
  }

  # In the event that ALL weights are equal to 0, scale will
  # return NAN values for that column.
  # We handle this by setting all NAN values equal to 0, meaning that the
  # elevation difference will have NO effect at this point
  # (and the weight will NOT sum to unity).
  tempDiff[is.nan(tempDiff)] = 0

  return(tempDiff)
}
#=============================================================================


## Cluster Weight Function ##
#=============================================================================
# The function down-weights stations that are clustered together to ensure
# that clusters of stations don't receive dissproportionate influence in the
# predictions.
# Inputs:
# locations - a spatialpointsdataframe of the measurement locations
# r - the minimum radius of influence, used in the distance cluster factor
# p - elevation precision: used in the elevation precision factor
#=============================================================================
getClusterWeight = function(locations, r, p){
  if(regexpr(sp::proj4string(locations), pattern = "+proj=longlat") > 1){
    distDiff = fields::rdist.earth(locations@coords, miles = FALSE)
  }else{
    distDiff = fields::rdist(locations@coords)
  }

  h = (.2*r - distDiff) / (.2*r)
  # Any stations outside the 20% of the radius of influence set to 0
  h[h < 0] = 0
  h[h > 1] = 1
  # Don't want to consider distance from a station to itself
  h = h - diag(diag(h))

  # Calculate s_ij matrix
  elevDiff = abs(outer(locations$ELEVATION,
                       locations$ELEVATION, "-")) - p
  elevDiff[elevDiff < 0] = 0
  # Calculate v_i
  v = (p - elevDiff)/p
  v[v < 0] = 0
  v[v > 1] = 1
  # Don't want to consider the elevation difference between a station
  # and itself
  v = v - diag(diag(v))

  Sc = h*v
  # +1 avoids dividing by 0, or getting values bigger than one for
  # near zero factors
  Sc = rowSums(Sc) + 1
  W_c  = 1 / Sc
  W_c = W_c / sum(W_c)
  return(as.vector(W_c))
}
#=============================================================================

## Eco-Region Weights ##
#=============================================================================
# The eco-region weight was a specific adaptation made for the Washington
# snow load report. This weight is intended to replace the water basin weight.
# Inputs:
# gridEco - a vector of strings with the eco-regions of the prediction
#           locations.
# stationEco - a vector of strings with the eco-regions of the measurement
#              locations.
# d - a user-defined exponent to control the severity of the applied weights.
#=============================================================================
getEcoWeight = function(stationEco, gridEco, d){

  rIII <- as.character(gridEco)
  sIII <- as.character(stationEco)

  # Check for missing values.
  na1 <- sum(is.na(rIII))
  na2 <- sum(is.na(sIII))

  if(na1 > 0){
    warning(paste("Ecoregions missing for", na1, "prediction locations..."))
    rIII[is.na(rIII)] <- " "
  }
  if(na2 > 0){
    warning(paste("Ecoregions missing for", na1, "station locations..."))
    sIII[is.na(sIII)] <- " "
  }

  # Determine hierachy of equality
  rII <- gsub(x = rIII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")
  sII <- gsub(x = sIII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")

  rI <- gsub(x = rII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")
  sI <- gsub(x = sII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")

  v1 <- outer(sI, rI, "==")
  v2 <- outer(sII, rII, "==")
  v3 <- outer(sIII, rIII, "==")

  # Stations receive weights proportional to the number of watersheds
  # they share.
  # C controls how drastic the weighting scheme will be
  tempweight = ((v1 + v2 + v3 + 1)/4)^d

  # stations in the row, cells in the columns
  for(i in 1:ncol(tempweight)){
    tempweight[, i] = tempweight[, i]/sum(tempweight[, i])
  }
  return(tempweight)
  #return(scale(tempweight, center = FALSE, scale = colSums(tempweight)))
}


## Water Basin Weights ##
#=============================================================================
# The water-basin weight was created in our original adapation of PRISM. It
# is intended to act as a surrogate for topographical facet. The program uses
# HUC codes 2 through 8 to define weights for a given region. Originally,
# the program used HUC12 designations.
# Inputs:
# basinDEM - the HUC codes of the prediction locations.
# stationBasin - the HUC codes of the measurement locations.
# d - a user-defined exponent to control the severity of the applied weights.
#=============================================================================
getHUCWeight = function(stationBasin, basinDEM, d){

  if(sum(is.na(stationBasin)) > 0){
    stop("Each station must have a specified HUC")
  }

  huc8.basin = floor(as.numeric(as.character(as.vector(basinDEM))))
  huc8.station = floor(as.numeric(as.character(stationBasin)))
  huc6.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 100)
  huc6.station = floor(as.numeric(as.character(stationBasin)) / 100)
  huc4.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 10000)
  huc4.station = floor(as.numeric(as.character(stationBasin)) / 10000)
  huc2.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 1000000)
  huc2.station = floor(as.numeric(as.character(stationBasin)) / 1000000)

  # Store elevation differences between each gridcell and the stations
  # in the columns of a matrix
  sameBasin8 = outer(huc8.station, huc8.basin, "==")
  sameBasin6 = outer(huc6.station, huc6.basin, "==")
  sameBasin4 = outer(huc4.station, huc4.basin, "==")
  sameBasin2 = outer(huc2.station, huc2.basin, "==")

  # Stations receive weights proportional to the number of
  # watersheds they share.
  # C controls how drastic the weighting scheme will be
  tempweight = ((sameBasin8+sameBasin6+sameBasin4+sameBasin2 + 1)/5)^d

  # stations in the row, cells in the columns
  for(i in 1:ncol(tempweight)){
    tempweight[, i] = tempweight[, i]/sum(tempweight[, i])
  }
  return(tempweight)
}

