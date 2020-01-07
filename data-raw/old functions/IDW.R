#' Inverse Distance Weighting, the Idaho way
#' Confirmed this matched the old IVD2 on 10-25-2018
#'
#' A function that implements inverse distance weighting using
#' normalized ground snow loads.
#'
#' @param formula A formula that is passed to the idw function in the
#'   sp package.
#' @param train The measurement locations to be used for interpolation.
#'   An object of class SpatialPointsDataFrame
#' @param test The prediction locations. An objet of class
#'   SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.
#'
#' @return An appended version of the test class containing a column named
#'   IDW.
#'
#' @export
IDW_snow = function(variables = c("RESPONSE", "ELEVATION"),
                    train, test, tlayer = 1220, c1 = 2, c2 = 6,
                    NGSL = TRUE, maxout = TRUE, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    stop("Train must be an object of class SpatialPointsDataFrame")
  }
  if(!is.element(class(test)[1], c("SpatialPointsDataFrame",
                                      "SpatialGridDataFrame",
                                      "SpatialPixelsDataFrame"))){
    stop("Prediction locations must be a spatial class from the sp package")
  }
  # Make sure we have identical projections.
  if(sp::proj4string(train) != sp::proj4string(test)){
    stop("Projections for train and test must be identical")
  }

  if(maxout){
    max.out <- max(train@data[, variables[1]])
  }

  # Use the Normalized ground snow loads
  if(NGSL){
    # Prevent weird things from happening when elevation
    # is very close to zero or negative.
    train@data[as.vector(train[, variables[2]]@data < 1), variables[2]] <- 1
    test@data[as.vector(test[, variables[2]]@data < 1), variables[2]] <- 1

    train@data[, variables[1]] <-
      train@data[, variables[1]]/train@data[, variables[2]]
  }

  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method)
  train.low = train[as.vector(train[, variables[2]]@data < tlayer), ]
  train.high = train[as.vector(train[, variables[2]]@data >= tlayer),]

  # Create a formula from the provided variables.
  tform <- as.formula(paste(variables[1], "1", sep = "~"))
  # Obtain the bi-layer results
  if(nrow(train.low) > 0){
    temp2 = gstat::idw(tform, locations = train.low, newdata = test,
                       idp = c1, ...)$var1.pred
  }else{
    temp2 = NULL
    # print("No stations exist in the lower layer, using only one layer with exponent c2...",
    # call. = FALSE)
  }
  if(nrow(train.high) > 0){
    temp6 = gstat::idw(tform, locations = train.high, newdata = test,
                       idp = c2, ...)$var1.pred
  }else{
    temp6 = NULL
  }

  # Store results by layer, assuming both layers are defined.
  results = vector("numeric", nrow(test))
  if(is.null(temp2)){
    results = temp6
  }else if(is.null(temp6)){
    results = temp2
  }else{
    results[test@data[, variables[2]] < tlayer] = temp2[test@data[, variables[2]] < tlayer]
    results[test@data[, variables[2]] >= tlayer] = temp6[test@data[, variables[2]] >= tlayer]
  }

  # Return ground snow load predictions as opposed to simply the NGSL.
  if(NGSL){results = results*test@data[, variables[2]]}

  # If requested restrict the response values to the highest valued snow load in the dataset.
  if(maxout){
    results[results > max.out] <- max.out
  }

  return(results)
} # End the function

