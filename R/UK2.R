#' Title: Universal Kriging - Original
#' 
#' A function that implements an ecoregion-based approach of universal kriging
#'
#' @param 
#' 
#' @return An object with reponse variable predictions for 
#' each gridcell, and row and column names representing the latitude 
#' and longitude coordinates respectively 
#'
#' @examples
UK2 = function(train, test, maxel = TRUE, maxout = TRUE, ...){
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
  
  # Restrict elevation values to max elevation observed in the data set. This helps
  # control for the effects of extrapolation. 
  if(maxel){
    max.el <- max(train$ELEVATION)
    test$ELEVATION[test$ELEVATION > max.el] <- max.el
  }
  
  # Fit Universal Kriging
  krigetest = gstat::krige(log(RESPONSE + 1) ~ ELEVATION, locations = train, 
                           newdata = test, ...)
  
  if(maxout){max.out <- max(train$RESPONSE)}
  
  # Add kriging coefficients and exponentiate
  if(israster){
    krigetest$var1.pred = exp(krigetest$var1.pred) - 1
    
    # Adjust high values if requested
    if(maxout){
      krigetest$var1.pred[krigetest$var1.pred > max.out] <- max.out
    }
    
    return(raster::rasterFromXYZ(krigetest))
  }else{
    # Estimate for witheld stations
    g.est2 = exp(krigetest$var1.pred) - 1
    
    # Adjust high values if requested
    if(maxout){
      g.est2[g.est2 > max.out] <- max.out
    }
    
    return(g.est2)
  }
  
} # End the function

