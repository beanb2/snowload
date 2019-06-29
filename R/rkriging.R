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
#' @param weights If TRUE, bundle the weights and linear model coefficients
#'   as a list. Should be used only for plotting illustrations.
#'
#' @return A numeric vector of predictions with length equal to the number
#'   of rows in newdata.
#'
#' @export
rkriging = function(formula = log(RESPONSE) ~ ELEVATION, locations, newdata,
                    model, bound_elevation = FALSE, bound_output = FALSE,
                    sklm = FALSE, weights = FALSE, ...){
  # Generic check for inputs.
  input_check(locations, newdata)

  # Extract the variable names in the formula.
  vars <- all.vars(formula)

  if(any(vars == "RESIDUALS")){
    warning("Overwriting the RESIDUALS column in the locations data frame
    during prediction, consider renaming RESIDUALS column prior to
            function input.")
  }

  # Restrct elevation trend estimates to be no greater than the trend value
  # of the highest elevation station location.
  if(bound_elevation){
    if(length(vars) > 2){
      warning("More than one explanatory variable supplied. Assuming the
              first variable corresponds to elevation.")
    }
    bound_el <- c(min(locations@data[, vars[2]], na.rm = TRUE),
                  max(locations@data[, vars[2]], na.rm = TRUE))

    newdata@data[,vars[2]][newdata@data[, vars[2]] > bound_el[2]] <-
      bound_el[2]
    newdata@data[,vars[2]][newdata@data[, vars[2]] < bound_el[1]] <-
      bound_el[1]
  }

  # If the specified model has any NA values, call fit.variogram to
  # define the model.
  if(any(c(is.na(model$model), is.na(model$psill), is.na(model$range)))){
    g <- gstat::gstat(NULL, "vario", formula, locations)
    g_vario <- gstat::variogram(g)
    model <- gstat::fit.variogram(g_vario, model)
  }

  if(sklm){
    # Fit Simple Kriging to Residuals
    gfit = lm(formula, data = as.data.frame(locations))
    locations$RESIDUALS <- gfit$residuals
    krigetest = gstat::krige(RESIDUALS~1, locations = locations,
                             newdata = newdata, model = model, beta = 0, ...)
    krigetest = predict(gfit, newdata) + krigetest$var1.pred
  }else{
    # Fit Universal Kriging (or kriging with an external drift)
    krigetest = gstat::krige(formula, locations = locations,
                             newdata = newdata, model = model, ...)$var1.pred
  }

  # If requested, bound predictions by observations in the data
  # model.frame evaluates the provided expression for us.
  # The response is given in the first column.
  if(bound_output){
    temp_out <- stats::model.frame(formula, data = locations)[, 1]
    bound_out <- c(min(temp_out, na.rm = TRUE),
                   max(temp_out, na.rm = TRUE))

    krigetest[krigetest < bound_out[1]] <- bound_out[1]
    krigetest[krigetest > bound_out[2]] <- bound_out[2]
  }

  # If weights are requested, calculate the kriging weights "by hand"
  #===========================================================================
  if(weights){
    # Determine the response variable
    modeldf <- stats::model.frame(formula, locations)

    # If UK is requested, we need to determine the linear model coefficients
    if(!sklm){
      # Create a bogus data drame to infer the slope and intercept.
      # Help on matrix approach:
      # - https://stackoverflow.com/questions/32712301/create-empty-data-frame-with-column-names-by-assigning-a-string-vector
      fdf <-  as.data.frame(
        matrix(c(rep(mean(locations@coords[, 1], na.rm = TRUE), 2),
                 rep(mean(locations@coords[, 2], na.rm = TRUE), 2),
                 c(0, 1)), ncol = 3, byrow = FALSE)
      )
      # Ensure the final column has the same explanatory variable name.
      colnames(fdf) <- c("LONGITUDE", "LATITUDE", vars[length(vars)])
      sp::coordinates(fdf) <- c("LONGITUDE", "LATITUDE")
      sp::proj4string(fdf) <- sp::proj4string(newdata)

      testCoef <- gstat::gstat(formula = formula,
                               locations = locations,
                               model = model)

      coef <- predict(testCoef, fdf, BLUE = TRUE)

      params <- c(coef$var1.pred[1], coef$var1.pred[2] - coef$var1.pred[1])

      if(ncol(modeldf) > 2){
        stop("Weights = TRUE only possible if there is one and only one
             predictor variable.")
      }

      resid <- modeldf[, 1] - (params[1] + params[2]*modeldf[, 2])
      locations$RESIDUALS <- resid
    }else{
      # Retain the linear model coefficients if sklm = TRUE
      params <- unname(gfit$coefficients)
    }

    # Extract the distances to pass to the variogram.
    distM <- sp::spDists(locations)
    distL <- sp::spDists(locations, newdata)

    # Ensure diagonal exactly equal to 0 to avoid precision errors
    # when a nugget effect is present.
    for(i in 1:nrow(distM)){distM[i, i] <- 0}

    # Create the semi-variance matrices used in the linear model.
    # (SKLM uses the covariance specification)
    vgmM <- gstat::variogramLine(model, dist_vector = distM,
                                 covariance = sklm)
    vgmL <- gstat::variogramLine(model, dist_vector = distL,
                                 covariance = sklm)

    # -1 in the final column since semi-variances are used. If covariances
    # are used, this value should be +1.
    # However, SKLM is not subject to the same unbiasedness constraints.
    if(sklm){
      ls <- vgmM
    }else{
      ls <- cbind(vgmM, rep(-1, nrow(vgmM)))
      # Ensures the unbiaseness of the estimator (lagrange multiplier)
      ls <- rbind(ls, c(rep(1, ncol(vgmM)), 0))
    }

    tpred <- vector("numeric", nrow(newdata))
    mweights <- matrix(0, nrow(locations), ncol = nrow(newdata))
    for(i in 1:nrow(newdata)){
      # SKLM requires no lagrange multiplier, but SKLM = FALSE does
      if(sklm){
        rs <- vgmL[, i]
      }else{
        rs <- c(vgmL[, i], 1)
      }

      # Solve the kriging system
      tweights <- solve(ls, as.matrix(rs, ncol = 1))

      # Remove the lagrange multiplier estimate when SKLM = FALSE
      if(sklm){
        mweights[, i] <- tweights
      }else{
        mweights[, i] <- tweights[-length(tweights)]
      }

      tpred[i] <- sum(as.vector(locations$RESIDUALS) * as.vector(mweights[, i]))
    }

    finalPred <- params[1] +
      params[2]*newdata@data[, vars[length(vars)]] +
      tpred

    return(list(finalPred, tpred, mweights, params))
    #=========================================================================
  }else{
    return(unname(krigetest))
  }
}
