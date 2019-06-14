#' Regression-kriging functions for snow load estimation.
#'
#' The function combines simple kriging with varying local means (SKLM) and
#' universal kriging (UK) (more precisely, kriging with an external drift).
#' The only difference between the two regression-kriging estimators is that
#' UK leverages spatial correlations in the regression fitting.
#'
#' @param formula The formula used to account for elevation in the predictions.
#'   Typically of the form given as the default option.
#' @param locations The measurement locations (used to train the model).
#'   Of class SpatialPointsDataFrame from the sp package.
#' @param newdata The prediction locations. Must be one of the accepted spatial
#'   classes from the sp package.
#' @param model A variogram model created with the vgm function
#'   from the gstat package.
#' @param bound_elevation If true, restrict the trend predictions using the
#'   range of measurement location elevations.
#' @param bound_output If true, restrict the final predictions to the range
#'   of observed data at measurement locations.
#' @param sklm If true, use simple kriging with varying local means. If false,
#'   use "universal kriging". The only difference is that sklm assumes no
#'   correlation between locations when fitting the linear model for
#'   elevation.
#'
#' @return An appended version of newdata which includes a column
#'   named "rkriging" containing the predictions. Note that if
#'   bound_elevation-TRUE then the elevations in this object might be
#'   different than the original object.
#'
#' @export
rkriging = function(formula = log(RESPONSE) ~ ELEVATION, locations, newdata,
                    model, bound_elevation = FALSE, bound_output = FALSE,
                    sklm = FALSE, ...){
  # Generic check for inputs.
  input_check(locations, newdata)

  # Extract the variable names in the formula.
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

  if(sklm){
    # Fit Simple Kriging to Residuals
    gfit = lm(formula, data = as.data.frame(locations))
    locations$RESIDUALS = gfit$residuals
    krigetest = gstat::krige(RESIDUALS~1, locations = locations,
                      newdata = newdata, model = model, beta = 0, ...)
    krigetest = predict(gfit, newdata) + krigetest$var1.pred
  }else{
    # Fit Universal Kriging (or kriging with an external drift)
    krigetest = gstat::krige(formula, locations = locations,
                      newdata = newdata, model = model, ...)$var1.pred
  }

  # If requested, bound predictions by observations in the data
  if(bound_output){
    temp_out <- lm(formula, data = as.data.frame(locations))
    temp_out <- temp_out$residuals + temp_out$fitted.values
    bound_out <- c(min(temp_out, na.rm = TRUE),
                   max(temp_out, na.rm = TRUE))

    krigetest[krigetest < bound_out[1]] <- bound_out[1]
    krigetest[krigetest > bound_out[2]] <- bound_out[2]
  }

  newdata$rkriging <- krigetest

  return(newdata)
}
