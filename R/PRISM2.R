#' Title: PRISM - Log Response (water basins)
#' 
#' A function that implements the PRISM algorithm as developed by Daly et.al. 
#'
#' @param eGrid A class "Raster" object representing a DEM.
#'  Must have a valid lat/lon extent defined in the 
#' layer for meaningful distance calculations 
#' @param stationD A dataframe with the following 
#' variables defined (case specific):
#' RESPONSE, ELEVATION, LATITUDE, LONGITUDE, HUC (Optional)
#' @param rm (Distance Weighting) - Minimum radus of influence 
#' (Typically 7-10km for precipitation data)
#' @param a (Distance Weighting) Distance weighting exponent.
#' @param b (Elevation Weighting) Elevation weighting exponent.
#' @param zm Minimum elevation difference of influence (100 - 300m)
#' @param zx Maximum elevation difference of influence (500 - 2500m)
#' @param Fd Importance factor for distance weighting must be 
#' between 0 and 1 (Fz = 1 - Fd)
#' @param zx Maximum elevation difference of influence (500 - 2500m)
#' @param d a coeffcient in the basin weighting
#' @param noNegative if true, all negative PRISM predictions are set to 0
#' @param maxout if true, PRISM predictions are not allowed to get beyond 
#' the maximum value of the response variable
#'
#' @return An object with reponse variable predictions for 
#' each gridcell, and row and column names representing the latitude 
#' and longitude coordinates respectively 
#'
#' @examples
PRISM2 = function(stationD, DEM, rm = 10, a = 2, b = 1, d = 1, 
                 zm = 200, zx = 2500, Fd = .8,
                 islog = TRUE, noNegative = TRUE, maxout = TRUE,
                 weights = TRUE){
  
  ### READ IN RASTER DATA ###
  # DEM MUST be a raster
  if(class(DEM)[1] != "SpatialPointsDataFrame" | 
     class(stationD)[1] != "SpatialPointsDataFrame"){
    stop("Both inputs must be of class \"SpatialPointsDataFrame\"")
  } 
  
  if(raster::projection(stationD) != raster::projection(DEM)){
    stop("Projections of inputs and outputs must be identical")
  }
  
  # Ensure station data meet criteria. 
  if(is.null(stationD$RESPONSE) | is.null(stationD$ELEVATION) | 
     is.null(stationD$HUC)){
    stop("The following variables MUST be included in StationD: 
         RESPONSE, ELEVATION, HUC")
  }
  
  if(is.null(DEM$ELEVATION) | is.null(DEM$HUC)){
    stop("The following variables MUST be included in DEM: 
          ELEVATION, HUC")
  }
  
  # Make copies of relevant variables. 
  expl <- stationD$ELEVATION
  if(islog){
    resp <- log(stationD$RESPONSE + 1)
  }else{
    resp <- stationD$RESPONSE
  }
  elev <- DEM$ELEVATION
  
  predCoord <- DEM@coords
  ### DETERMINE WEIGHTING SCHEME ###
  Wd = getDistWeight2.2(predCoord, stationD@coords, 
                      rm, a) 
  Wz = getElevWeight2.2(elev, expl, b, zm, zx)
  Wf = getHUCWeight(DEM$HUC, stationD$HUC, d)
  Wc = getClusterWeight.2(stationD, rm, zm)

  
  # Apply weighting equation from Daley et.al. 2008
  weightV = Wc*sqrt(Fd*Wd^2 + (1-Fd)*Wz^2)*Wf
  # weightV = base::scale(weightV, center=FALSE, scale = colSums(weightV))
  for(i in 1:ncol(weightV)){
    weightV[, i] <- weightV[, i]/sum(weightV[, i])
  }
  
  # Remove stuff we no longer need to save on memory
  remove(Wc,Wd,Wz,Wf)
  ##################################
  
  # Preallocate for predictions
  predictP = vector("numeric", length = length(DEM))
  slope = vector("numeric", length = length(DEM))
  intercept = vector("numeric", length = length(DEM))
  for(i in 1:length(predictP)){
    if(is.na(elev[i])){
      predictP[i] = NA
      slope[i] = NA
      intercept[i] = NA
    }else{
      # Run linear regression on individual gridcell, 
      # with the respective weighting vector
      tempReg = lsfit(as.matrix(expl, ncol = 1), resp, weightV[,i])
      # Predict based on centroid elevation of current grid (retransform data)
      slope[i] = tempReg$coefficients[2]
      intercept[i] = tempReg$coefficients[1]
      predictP[i] = tempReg$coefficients[1] + tempReg$coefficients[2]*elev[i]
    }
  }
  
  # Do not let predictions extrapolate beyond readings 
  if(maxout){predictP[predictP > max(resp)] <- max(resp)}
  
  # Retransform Results if log transform is being used
  if(islog){predictP = exp(predictP) - 1}
  
  # Set negative predicts to 0 if necessary. 
  if(noNegative){predictP[predictP < 0] <- 0}
  
  finalPredicts <- data.frame(LONGITUDE = predCoord[, 1], 
                              LATITUDE = predCoord[, 2],
                              predictions = predictP, 
                              intercept = intercept,
                              slope = slope)
  
  sp::coordinates(finalPredicts) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(finalPredicts) <- sp::proj4string(DEM)
  
  if(weights){
    return(list(finalPredicts, weightV))
  }else{
    return(finalPredicts)
  }
}


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
getDistWeight2.2 = function(gridC, stationC, rm, a){
  
  # Calculate distance with function from "fields" package.
  # Note that we transpose the results so weights are found in the columns, 
  # not the rows.
  # Avoid singularity issues in getDistWeight by rounding.
  dist = t(round(fields::rdist.earth(gridC, stationC, miles = FALSE), 
                 digits=4)) 
  
  # Set all stations within the minimum radius of influence to a weight of 1
  distWeight = dist - rm
  distWeight[distWeight <= 0] = 1
  
  # For all other stations, compute the weighting scheme. 
  distWeight = distWeight^a
  distWeight = 1 / distWeight
  distWeight[distWeight > 1] = 1 # Set all values greater than 1 to the max
  
  # Normalize each weighting vector to sum to 1
  # distWeight = scale(distWeight, center = FALSE, scale = colSums(distWeight))
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
## Find Elevation Weights ##
getElevWeight2.2 = function(gridElv, stationElv, b, zm, zx){
  
  # Store elevation differences between each gridcell and the stations 
  # in the columns of a matrix
  tempDiff = t(abs(outer(gridElv, stationElv, "-")))
  
  # Determine which elevation differences lie between the maximum and minimum
  # elevation of influence range
  tempDiffMin = tempDiff - zm
  tempDiffMax = tempDiff - zx
  
  # Apply weighting criterion to tempDiff based on min/max values
  tempDiff[tempDiffMin <= 0] = zm
  tempDiff = tempDiff^b
  tempDiff = 1/tempDiff
  tempDiff[tempDiffMax >= 0] = 0
  
  # Scale elevation weights to sum to unity
  # tempDiff = scale(tempDiff, center = FALSE, scale = colSums(tempDiff))
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

## Find Cluster Weights ##
getClusterWeight.2 = function(stationD, r, p){
  distDiff = fields::rdist.earth(stationD@coords, miles = FALSE)
  h = (.2*r - distDiff) / (.2*r)
  # Any stations outside the 20% of the radius of influence set to 0
  h[h < 0] = 0 
  h[h > 1] = 1
  # Don't want to consider distance from a station to itself
  h = h - diag(diag(h)) 
  
  # Calculate s_ij matrix
  elevDiff = abs(outer(stationD$ELEVATION, 
                       stationD$ELEVATION, "-")) - p
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

getHUCWeight = function(basinDEM, stationBasin, d){
  
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
  
  # Store elevation differences between each gridcell and the stations in the columns of a matrix
  sameBasin8 = t(outer(huc8.basin, huc8.station, "=="))
  sameBasin6 = t(outer(huc6.basin, huc6.station, "=="))
  sameBasin4 = t(outer(huc4.basin, huc4.station, "=="))
  sameBasin2 = t(outer(huc2.basin, huc2.station, "=="))
  
  # Stations receive weights proportional to the number of watersheds they share. 
  # C controls how drastic the weighting scheme will be
  tempweight = ((sameBasin8+sameBasin6+sameBasin4+sameBasin2 + 1)/5)^d
  
  # stations in the row, cells in the columns
  for(i in 1:ncol(tempweight)){
    tempweight[, i] = tempweight[, i]/sum(tempweight[, i])
  }
  return(tempweight)
 # return(scale(tempweight, center=FALSE, scale = colSums(tempweight))) # stations in the row, cells in the columns
}