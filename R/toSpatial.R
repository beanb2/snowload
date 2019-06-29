#' Convert snowload datasets to spatial data frames
#'
#' This function is a simple wrapper to sp package functions to create spatial
#'   data frame objects with a common projection. The function assumes that
#'   the variable name conventions of the datasets internal to the snowload
#'   package. It is not intended to be used on datasets outside of the
#'   snowload package.
#'
#' @param tdata A string naming one of the datasets internal to the snowload
#'   package.
#'
#' @return A spatial version with the same coordinate extent as the
#'   utdem SpatialPixelsDataFrame. (Also internal to the snowload package.)
#' @export
toSpatial <- function(tdata){
  data(list = tdata)

  # Use get function to make assignment
  tempd <- get(tdata)

  sp::coordinates(tempd) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(tempd) <- sp::proj4string(utdem)

  return(tempd)
}
