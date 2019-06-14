# Each of the functions are designed to receive similar inputs. This function
# ensures consistent input of the measurement and prediction locations
# across all functions. Its only design is informative error checking.
input_check <- function(locations, newdata){
  if(class(locations)[1] != "SpatialPointsDataFrame"){
    stop("locations must be an object of class SpatialPointsDataFrame")
  }
  if(!is.element(class(newdata)[1], c("SpatialPointsDataFrame",
                                      "SpatialGridDataFrame",
                                      "SpatialPixelsDataFrame"))){
    stop("Prediction locations must be a spatial class from the sp package")
  }
  # Make sure we have identical projections.
  if(sp::proj4string(locations) != sp::proj4string(newdata)){
    stop("Projections for locations and newdata must be identical")
  }
  return(NULL)
}
