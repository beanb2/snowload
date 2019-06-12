# This script contains the various weight functions for PRISM. These functions
# are used internally by the prism() function and will not be exported.

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
