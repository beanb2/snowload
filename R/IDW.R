#' Inverse Distance Weighting, the Idaho way
#' Confirmed this matched the old IVD2 on 10-25-2018
#' 
#' A function that implements inverse distance weighting using 
#' normalized ground snow loads. 
#'
#' @param 
#' 
#' @return An object with reponse variable predictions for 
#' each gridcell, and row and column names representing the latitude 
#' and longitude coordinates respectively 
#'
#' @examples
IDW = function(train, test, tlayer = 1220, c1 = 2, c2 = 6, 
               maxout = TRUE, NGSL = TRUE, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    sp::coordinates(train) = c("LONGITUDE", "LATITUDE")
    sp::proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }
  
  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train = sp::spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test = raster::rasterToPoints(test, spatial = TRUE)
    test$ELEVATION = test$layer
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      sp::coordinates(test) = c("LONGITUDE", "LATITUDE")
      sp::proj4string(test) = sp::proj4string(train)
    }
  }
  
  if(maxout){
    max.out <- max(train$RESPONSE)
  }
  
  # Use the Normalized ground snow loads
  if(NGSL){
    # Prevent weird things from happening when elevation 
    # is very close to zero or negative.
    train$ELEVATION[train$ELEVATION < 1] <- 1
    test$ELEVATION[test$ELEVATION < 1] <- 1
    
    train$RESPONSE <- train$RESPONSE/train$ELEVATION
  }
  
  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method) 
  train.low = train[train$ELEVATION < tlayer,]
  train.high = train[train$ELEVATION >= tlayer,]  
  
  
  # Obtain the bi-layer results
  if(nrow(train.low) > 0){
    temp2 = gstat::idw(RESPONSE~1, locations = train.low, newdata = test, 
                idp = c1, ...)$var1.pred
  }else{
    temp2 = NULL
   # print("No stations exist in the lower layer, using only one layer with exponent c2...", 
   # call. = FALSE)
  }
  if(nrow(train.high) > 0){
    temp6 = gstat::idw(RESPONSE~1, locations = train.high, newdata = test, 
                idp = c2, ...)$var1.pred
  }else{
    temp6 = NULL
    # print("No stations exist in the upper layer, using only one layer with exponent c1...", 
    # call. = FALSE)
  }
  
  # Store results by layer, assuming both layers are defined. 
  results = vector("numeric", nrow(test))
  if(is.null(temp2)){
    results = temp6
  }else if(is.null(temp6)){
    results = temp2
  }else{
    results[test$ELEVATION < tlayer] = temp2[test$ELEVATION < tlayer]
    results[test$ELEVATION >= tlayer] = temp6[test$ELEVATION >= tlayer]
  }
  
  # Return ground snow load predictions as opposed to simply the NGSL. 
  if(NGSL){results = results*test$ELEVATION}
  
  # If requested restrict the response values to the highest valued snow load in the dataset. 
  if(maxout){
    results[results > max.out] <- max.out
    }
  
  return(results)
} # End the function

