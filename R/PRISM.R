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
#' @transform If "none" (the default). Regular least squares regression is
#'   applied. If "log" or "sqrt" then the appropriate transformations are
#'   applied. If not one of the options mentioned above, the system will
#'   throw an error.
#'
#' @return An appended version of the prediction locations spatial object
#'   with the added "prism" column. If weights = TRUE, then a list is returned
#'   with the second item being the matrix of PRISM weights.
#'
#' @export
prism = function(variables = c("RESPONSE", "ELEVATION"),
                 locations,
                 newdata,
                 distImp = 0.8,
                 minRad = 10,
                 wdistance = 2,
                 welevRange = c(200, 2500),
                 welevation = c("ELEVATION", 1),
                 wbasin = c(NA, 3),
                 weco = c(NA, 3),
                 bound = TRUE,
                 weights = FALSE,
                 transform = "none"){

  ### READ IN RASTER DATA ###
  # newdata MUST be a raster
  if(class(locations)[1] != "SpatialPointsDataFrame"){
    stop("Measurement locations must be of class \"SpatialPointsDataFrame\"")
  }

  if(!is.element(class(newdata)[1], c("SpatialPointsDataFrame",
                                      "SpatialGridDataFrame",
                                      "SpatialPixelsDataFrame"))){
    stop("Prediction locations must be a spatial class from the sp package")
  }

  if(sp::proj4string(locations) != sp::proj4string(newdata)){
    stop("Projections of inputs and outputs must be identical")
  }

  # Ensure that the arguments for the weights are properly specified
  tl <- c(length(welevation), length(wbasin), length(weco))
  if(min(tl) < 2 | max(tl) > 2){
    stop("The \'welevation\', \'wbasin\', and \'weco\' arguments must all be
         of length 2.")
  }

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

  # Make sure that any specified transformations are accepted.
  if(!is.element(transform, c("none", "log", "sqrt"))){
    stop("The only valid options for transform are \"log\",
    \"sqrt\", and \"none\"")
  }

  # Find the distance and elevation weights (required)
  wd <- getDistWeight(locations, newdata, minRad, wdistance)
  wz <- getElevWeight(unlist(locations[, welevation[1]]@data),
                      unlist(newdata[, welevation[1]]@data),
                      as.numeric(welevation[2]), welevRange)

  wc <- getClusterWeight(locations, minRad, welevRange[1])

  # Get the basin weight (if requested)
  if(!is.na(wbasin[1])){
    wb <- getHUCWeight(unlist(locations[, wbasin[1]]@data),
                       unlist(newdata[, wbasin[1]]@data),
                       as.numeric(wbasin[2]))
  }else{
    wb <- 1
  }

  # Get the eco weight (if requested)
  if(!is.na(weco[1])){
    wf <- getEcoWeight(unlist(locations[, weco[1]]@data),
                       unlist(newdata[, weco[1]]@data),
                       as.numeric(weco[2]))
  }else{
    wf <- 1
  }

  # Apply weighting equation from Daley et.al. 2008
  weightV <- wc*sqrt(distImp*wd^2 + (1-distImp)*wz^2)*wb*wf

  # Scale the weights to sum to 1
  for(i in 1:ncol(weightV)){
    weightV[, i] <- weightV[, i]/sum(weightV[, i])
  }

  # Remove stuff we no longer need to save on memory
  remove(wc, wd, wz, wb, wf)
  ##################################

  # Preallocate for predictions
  response <- unlist(locations[, variables[1]]@data)

  if(transform == "log"){response <- log(response + 1)}
  if(transform == "sqrt"){response <- sqrt(response)}

  explanatory <- as.matrix(locations[, variables[2:length(variables)]]@data)
  tempReg <- matrix(0, nrow = length(newdata), ncol = length(variables))

  predictP <- as.matrix(newdata[, variables[2:length(variables)]]@data)
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
    bounds <- max(response, na.rm = TRUE)

    predictP[predictP < 0] <- 0
    predictP[predictP > bounds] <- bounds
  }

  if(transform == "log"){predictP <- exp(predictP) - 1}
  if(transform == "sqrt"){predictP <- predictP^2}

  # If weights are requested (for visualization) then return a list.
  # Otherwise, return the appended spatial data frame.
  if(weights){
    newdata$prism <- predictP
    newdata$slope <- tempReg[, 2]
    newdata$intercept <- tempReg[, 1]
    return(list(newdata, weightV))
  }else{
    newdata$prism <- predictP
    return(newdata)
  }
}
