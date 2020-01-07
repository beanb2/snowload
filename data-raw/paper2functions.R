#load("C:/Users/Brennan Bean/TNMDownloadManager/elevationDEM.R")

#' Title: PRISM - Log Response
#'
#' A function that implements the PRISM algorithm as developed by Daly et.al.
#'
#' @param eGrid A class "Raster" object representing a DEM. Must have a valid lat/lon extent defined in the
#' layer for meaningful distance calculations
#' @param stationD A dataframe with the following variables defined (case specific):
#' RESPONSE, ELEVATION, LATITUDE, LONGITUDE, HUC (Optional)
#' @param rm (Distance Weighting) - Minimum radus of influence (Typically 7-10km for precipitation data)
#' @param a (Distance Weighting) Distance weighting exponent.
#' @param b (Elevation Weighting) Elevation weighting exponent.
#' @param zm Minimum elevation difference of influence (100 - 300m)
#' @param zx Maximum elevation difference of influence (500 - 2500m)
#' @param Fd Importance factor for distance weighting must be between 0 and 1 (Fz = 1 - Fd)
#' @param zx Maximum elevation difference of influence (500 - 2500m)
#' @param cluster Defaults to -1 indicating that no cluster weighting will be used. Used if non-negative elevation precision
#' is specified (typically around 50m)
#' @param basin a rasterlayer or vector representing the HUC12 water basins of the DEM
#' @param c a coeffcient in the basin weighting
#' @param noNegative if true, all negative PRISM predictions are set to 0
#' @param maxout if true, PRISM predictions are not allowed to get beyond the maximum value of the response variable
#'
#' @return An object with reponse variable predictions for each gridcell, and
#' row and column names representing the latitude and longitude coordinates respectively
#'
#' @examples
prism2 = function(stationD, DEM, rm = 10, a = 2, b = 1,
                  zm = 200, zx = 2500, Fd = .8, cluster = -1,
                  basin = NULL, d = 1, feet = FALSE, islog = TRUE, NGSL = FALSE, noNegative = TRUE, maxout = TRUE){
  require(fields)
  require(raster)
  require(sp)
  ### READ IN RASTER DATA ###
  # DEM MUST be a raster
  if(class(DEM) != "RasterLayer" & class(DEM) != "data.frame"){
    stop("DEM must be of class \"RasterLayer\" or \"data.frame\"")}

  # Option to handle either a raster (map) or data.frame (cross validation)
  if(class(DEM) == "RasterLayer"){
    # Extract unique coordinate extent for RasterLayer
    gridC = sp::coordinates(DEM)
    lonCen = unique(gridC[,1])
    latCen = unique(gridC[,2])

    tproj <- projection(DEM)

    # Extract extent of raster layer
    temp.ext = extent(DEM)

    # Convert raster layer to vector (fills a matrix BY ROW)
    DEM = as.vector(DEM)
    israster = TRUE
  }else{
    if(is.null(DEM$ELEVATION) | is.null(DEM$LONGITUDE) | is.null(DEM$LATITUDE)){
      stop("ELEVATION, LATITUDE, and LONGITUDE, must be defined in order to proceed")
    }
    gridC = data.frame(LONGITUDE = DEM$LONGITUDE,
                       LATITUDE = DEM$LATITUDE)
    # Convert feet to meters if necessary
    if(feet){
      DEM = DEM$ELEVATION * .3048
    }else{
      DEM = DEM$ELEVATION
    }
    israster = FALSE
  }

  ### READ IN STATION DATA ###
  if(is.null(stationD$RESPONSE) | is.null(stationD$ELEVATION) | is.null(stationD$LONGITUDE)
     | is.null(stationD$LATITUDE)){
    stop("The following variables MUST be included in StationD: RESPONSE, ELEVATION, LONGITUDE, LATITUDE")
  }

  # If Station elevations are in feet, convert to meters to match the DEM
  if(feet){
    explanatory = stationD$ELEVATION * .3048 # Explanatory variable (elevation values in meters)
  }else{
    explanatory = stationD$ELEVATION
  }

  # Set snow loads as the response variable
  response = stationD$RESPONSE

  # Determine if modeling NGSL or raw snow loads
  if(NGSL){
    response = response/explanatory
  }

  # Determine if modeling log of values (either NGSL or raw snow loads)
  if(islog){
    response = log(response+1)
  }

  stationC = data.frame(LONGITUDE = stationD$LONGITUDE, LATITUDE = stationD$LATITUDE) # Lon/Lat coordiantes stored as dataframe

  ### DETERMINE WEIGHTING SCHEME ###

  ## FIND DISTANCE WEIGHTS ##
  getDistWeight2 = function(gridC, stationC, rm, a){

    # Calculate distance with function from "fields" package
    # Note that we transpose the results so weights are found in the columns, not the rows
    dist = t(round(fields::rdist.earth(as.matrix(gridC), as.matrix(stationC), miles = FALSE), digits=4)) # Avoid singularity issues in getDistWeight by rounding

    # Set all stations within the minimum radius of influence to a weight of 1
    distWeight = dist - rm
    distWeight[distWeight <= 0] = 1

    # For all other stations, compute the weighting scheme.
    distWeight = distWeight^a
    distWeight = 1 / distWeight
    distWeight[distWeight > 1] = 1 # Set all values greater than 1 to the max

    # Normalize each weighting vector to sum to 1
    distWeight = scale(distWeight, center = FALSE, scale = colSums(distWeight))

    return(distWeight)
  }

  ## Find Elevation Weights ##
  getElevWeight2 = function(gridElv, stationElv, b, zm, zx){

    # Store elevation differences between each gridcell and the stations in the columns of a matrix
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
    tempDiff = scale(tempDiff, center = FALSE, scale = colSums(tempDiff))

    # In the event that ALL weights are equal to 0, scale will return NAN values for that column
    # we handle this by setting all NAN values equal to 0, meaning that the elevation difference
    # will have NO effect at this point (and the weight will NOT sum to unity)
    tempDiff[is.nan(tempDiff)] = 0

    return(tempDiff)
  }

  ## Find Cluster Weights ##
  getClusterWeight = function(lonlat, stationElv, r, p){
    distDiff = rdist.earth(as.matrix(lonlat), miles = FALSE)
    h = (.2*r - distDiff) / (.2*r)
    h[h < 0] = 0 # Any stations outside the 20% of the radius of influence set to 0
    h[h > 1] = 1
    h = h - diag(diag(h)) # Don't want to consider distance from a station to itself

    # Calculate s_ij matrix
    elevDiff = abs(outer(stationElv, stationElv, "-")) - p
    elevDiff[elevDiff < 0] = 0
    # Calculate v_i
    v = (p - elevDiff)/p
    v[v < 0] = 0
    v[v > 1] = 1
    v = v - diag(diag(v)) # Don't want to consider the elevation difference between a station and itself

    Sc = h*v
    Sc = rowSums(Sc) + 1 # +1 avoids dividing by 0, or getting values bigger than one for near zero factors
    W_c  = 1 / Sc
    W_c = W_c / sum(W_c)
    return(as.vector(W_c))
  }

  ## Find Aspect Weights ##
  getaspect = function(basinDEM, stationBasin, d){

    if(sum(is.na(stationBasin)) >0){
      stop("Each station must have a specified HUC")
    }

    #huc12.basin = as.numeric(as.character(as.vector(basinDEM)))
    #huc12.station = as.numeric(as.character(stationBasin))
    #huc10.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 100)
    #huc10.station = floor(as.numeric(as.character(stationBasin)) / 100)
    huc8.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 10000)
    huc8.station = floor(as.numeric(as.character(stationBasin)) / 10000)
    huc6.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 1000000)
    huc6.station = floor(as.numeric(as.character(stationBasin)) / 1000000)
    huc4.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 100000000)
    huc4.station = floor(as.numeric(as.character(stationBasin)) / 100000000)
    huc2.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 10000000000)
    huc2.station = floor(as.numeric(as.character(stationBasin)) / 10000000000)

    # Store elevation differences between each gridcell and the stations in the columns of a matrix
    #sameBasin12 = t(outer(huc12.basin, huc12.station, "=="))
    #sameBasin10 = t(outer(huc10.basin, huc10.station, "=="))
    sameBasin8 = t(outer(huc8.basin, huc8.station, "=="))
    sameBasin6 = t(outer(huc6.basin, huc6.station, "=="))
    sameBasin4 = t(outer(huc4.basin, huc4.station, "=="))
    sameBasin2 = t(outer(huc2.basin, huc2.station, "=="))

    # Stations receive weights proportional to the number of watersheds they share.
    # C controls how drastic the weighting scheme will be
    tempweight = ((sameBasin8+sameBasin6+sameBasin4+sameBasin2 + 1)/5)^d

    # tempweight[is.na(tempweight)] <- 1 # used for Montana report to deal with NA values.

    return(scale(tempweight, center=FALSE, scale = colSums(tempweight))) # stations in the row, cells in the columns
  }

  Wd = getDistWeight2(gridC, stationC, rm, a) # Get Distance Weights
  Wz = getElevWeight2(DEM, explanatory, b, zm, zx) # Get Elevation Weights

  # Get cluster weight only if user specifies it
  if(length(cluster) > 1){
    stop("Cluster input must be single numeric value")
  }

  if(cluster > 0){
    # Use same minimum radius of influence used in distance weighting
    Wc = getClusterWeight(stationC, explanatory, rm, cluster)
  }else{
    Wc = 1 # Apply no cluster weighting unless the user specifies
  }

  # Get basin weight only if user specifies it
  if(!is.null(basin)){
    if(is.null(stationD$HUC)){
      stop("\'HUC\' must be a variable in StationD data frame when \'basin\' is defined")
    }
    #print(paste("Fitting aspect weighting with d =", d, "..."))
    Wf = getaspect(basin, stationD$HUC, d)
  }else{
    Wf = 1 # Apply no cluster weighting unless the user specifies
  }

  # Apply weighting equation from Daley et.al. 2008
  weightV = Wc*sqrt(Fd*Wd^2 + (1-Fd)*Wz^2)*Wf
  weightV = scale(weightV, center=FALSE, scale = colSums(weightV))

  # Remove stuff we no longer need to save on memory
  remove(Wc,Wd,Wz,Wf)
  ##################################

  # Preallocate for predictions
  predictP = vector("numeric", length = length(DEM))
  #slope = vector("numeric", length = length(DEM))
  #intercept = vector("numeric", length = length(DEM))
  for(i in 1:length(predictP)){
    if(is.na(DEM[i])){
      predictP[i] = NA
      #slope[i] = NA
      #intercept[i] = NA
    }else{
      # Run linear regression on individual gridcell, with the respective weighting vector
      tempReg = lsfit(as.matrix(explanatory,ncol = 1), response, weightV[,i])
      # Predict based on centroid elevation of current grid (retransform data)
      #intercept[i] = tempReg$coefficients[1]
      #slope[i] = tempReg$coefficients[2]
      predictP[i] = tempReg$coefficients[1] + tempReg$coefficients[2]*DEM[i]

    }
  }

  # Do not let predictions extrapolate beyond readings
  if(maxout){predictP[predictP > max(response)] <- max(response)}

  # Retransform Results if log transform is being used
  if(islog){predictP = exp(predictP) - 1}

  # Extract snow loads from NGSL if necessary
  if(NGSL){predictP = predictP*DEM}

  # Set negative predicts to 0 if necessary.
  if(noNegative){predictP[predictP < 0] <- 0}

  # Return results as a raster if the input DEM was a raster, if not return the vector
  if(israster){
    # Return results in same form as input, note that raster fills by ROW not Column
    predictP = matrix(predictP, nrow = length(latCen), ncol = length(lonCen), byrow = TRUE)
    #intercept = matrix(intercept, nrow = length(latCen), ncol = length(lonCen), byrow = TRUE)
    #slope = matrix(slope, nrow = length(latCen), ncol = length(lonCen), byrow = TRUE)

    # Extent defined at beginning of function.
    predictP = raster::raster(predictP, temp.ext@xmin, temp.ext@xmax, temp.ext@ymin, temp.ext@ymax)
    #intercept = raster::raster(intercept, temp.ext@xmin, temp.ext@xmax, temp.ext@ymin, temp.ext@ymax)
    #slope = raster::raster(slope, temp.ext@xmin, temp.ext@xmax, temp.ext@ymin, temp.ext@ymax)

    # Project final Raster into same projection as DEM
    raster::projection(predictP) <- tproj
    #raster::projection(intercept) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    #raster::projection(slope) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    #return(list(predictP, slope, intercept))
  }

  return(predictP)
}


### Inverse Distance Weighting ###


# Function employing the gstat package ivd functions
# This function replicated the idaho method of inverse distance weighting with one key difference:
# here distance is geographic distance based on latitude longitude coordinates, while Idaho uses
# a Tranverse Mercator projection and then computes euclidean distance. (i believe)
# layer - elevation to split the DEM into two separate layers
# c1 - lower elevation IDW exponent
# c2 - upper elevation IDW exponent
IVD2 = function(stationD, DEM, layer = 1219.2, c1 = 2, c2 = 6, NGSL = TRUE,
                maxout = FALSE, debug.level = 0){
  #require(fields)
  #require(raster)
  #require(sp)
  #require(gstat)
  ### READ IN RASTER DATA ###

  # Option to handle either a raster (map) or data.frame (cross validation)
  if(class(DEM) == "RasterLayer"){
    # Convert to spatialpointsdataframe
    israster = TRUE
    DEM = rasterToPoints(DEM, spatial = TRUE)
    names(DEM) = "ELEVATION"

    # If the raster doesn't have a predefined projection. Given it geographic coordinate reference.
    if(is.na(projection(DEM))){
      projection(DEM) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    }
  }else{
    if(is.null(DEM$ELEVATION) | is.null(DEM$LONGITUDE) | is.null(DEM$LATITUDE)){
      stop("ELEVATION, LATITUDE, and LONGITUDE, must be defined in order to proceed")
    }
    israster = FALSE
    coordinates(DEM) = c("LONGITUDE", "LATITUDE")
    # Ensure that our data frame has a geographic coordinate reference system if not already a raster.
    proj4string(DEM) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Project into geographic coordinates
  }

  # Ensure station data in acceptable format
  if(class(stationD) != "data.frame" & class(stationD) != "SpatialPointsDataFrame"){
    stop("STATIOND must be of class \"data.frame\" or \"SpatialPointsDataFrame\"")
  }

  # Make sure variables of interest are specified
  if(is.null(stationD$RESPONSE) | is.null(stationD$ELEVATION) | is.null(stationD$LONGITUDE)
     | is.null(stationD$LATITUDE)){
    stop("The following variables MUST be included in StationD: RESPONSE, ELEVATION, LONGITUDE, LATITUDE")
  }

  # Convert to spatial points data frame if not already done
  if(class(stationD) == "data.frame"){
    sp::coordinates(stationD) = c("LONGITUDE", "LATITUDE")
  }

  # Define projections to be the same if stationD does not have a projection pre-defined.
  if(is.na(proj4string(stationD))){
    proj4string(stationD) = proj4string(DEM)
  }

  # If both projections are defined yet different. Transform to DEM projection
  if(proj4string(stationD) != proj4string(DEM)){
    stationD = spTransform(stationD, CRS(proj4string(DEM)))
  }

  # Store elevation information separately
  dem.elevation = DEM$ELEVATION

  # If requested, store the maximum ground snow load of the data set
  # on which we can restrict values
  if(maxout){max.out <- max(stationD$RESPONSE)}

  # Prior to splitting define the NGSL variable as the response if the user requests
  if(NGSL){stationD$RESPONSE = (stationD$RESPONSE/stationD$ELEVATION)}

  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method)
  stationD.low = stationD[stationD$ELEVATION < layer,]
  stationD.high = stationD[stationD$ELEVATION >= layer,]


  # Obtain the bi-layer results
  if(nrow(stationD.low) > 0){
    temp2 = idw(RESPONSE~1, locations = stationD.low, newdata = DEM, idp = c1, debug.level = debug.level)$var1.pred
  }else{
    temp2 = NULL
    warning("No stations exist in the lower layer, using only one layer with exponent c2...", call. = FALSE)
  }
  if(nrow(stationD.high) > 0){
    temp6 = idw(RESPONSE~1, locations = stationD.high, newdata = DEM, idp = c2, debug.level = debug.level)$var1.pred
  }else{
    temp6 = NULL
    warning("No stations exist in the upper layer, using only one layer with exponent c1...", call. = FALSE)
  }

  # Store results by layer, assuming both layers are defined.
  results = vector("numeric", length(dem.elevation))
  if(is.null(temp2)){
    results = temp6
  }else if(is.null(temp6)){
    results = temp2
  }else{
    results[dem.elevation < layer] = temp2[dem.elevation < layer]
    results[dem.elevation >= layer] = temp6[dem.elevation >= layer]
  }

  # Return ground snow load predictions as opposed to simply the NGSL.
  if(NGSL){results = results*dem.elevation}

  # If requested restrict the response values to the highest valued snow load in the dataset.
  if(maxout){results[results > max.out] <- max.out}

  # Return results as a raster if the input DEM was a raster, if not return the vector
  if(israster){
    DEM$ELEVATION = NULL
    DEM$RESULTS = results
    DEM = rasterFromXYZ(DEM)
    return(DEM)
  }

  return(results)
}



### CREATE OLD UTAH MAP ###


# Function to create old utah map given an input DEM
snlwf = function(dem, path1 = NULL, path2 = NULL){
  if(is.null(path1)){
    path1 = "data-raw/Counties"
  }
  if(is.null(path2)){
    path2 = "data-raw/Utahsnowlaw.csv"
  }
  # Create Spatial Points Data Frame of the coordinates of our DEM
  if(!(class(dem) %in% c('Raster', 'RasterLayer'))){stop("Input object must be a raster.")}
  pnts = rasterToPoints(dem)
  colnames(pnts) = c("LONGITUDE", "LATITUDE", "ELEVATION")
  pnts = data.frame(pnts)
  coordinates(pnts) = c("LONGITUDE", "LATITUDE")
  proj4string(pnts) = projection(dem)
  utah = readOGR(dsn=path1, layer="Counties")
  # Reproject to the same coordinate extent as our DEM
  utah <- spTransform(utah, CRS(projection(dem)))
  points.in.counties = over(pnts, utah)
  points.in.counties$COUNTYNBR = as.character(points.in.counties$COUNTYNBR)
  points.in.counties$COUNTYNBR = as.numeric(points.in.counties$COUNTYNBR)
  # Set county number ofpoints outside of ALL counties to 0
  points.in.counties$COUNTYNBR[is.na(points.in.counties$COUNTYNBR)] = 0

  ### Read in Old Utah Snow Code Data ###
  snowlaw = read.csv(path2)
  snowcodepic = vector(mode = "numeric", length = nrow(points.in.counties))
  elev = as.vector(dem)
  elev = elev / (1000 * .3048) # Convert meters to THOUSANDS of FEET

  for(i in 1:nrow(points.in.counties)){
    tempc = points.in.counties$COUNTYNBR[i]
    # Set all values outside of Utah equal to 0
    if(tempc == 0){
      tempest = 0
    }else{
      # Set snowload equal to base snowload if elevation in question lies below the base elevation
      if(elev[i] - snowlaw$A_0[tempc] <= 0){
        tempest = snowlaw$P_0[tempc]
      }else{
        # Otherwise, compute the equation as defined by the Utah snow code for that county
        tempest = (snowlaw$P_0[tempc]^2 + snowlaw$S[tempc]^2*(elev[i] - snowlaw$A_0[tempc])^2)^.5
      }
    }
    snowcodepic[i] = tempest
  }

  snowcodepic = matrix(snowcodepic, nrow = length(unique(coordinates(dem)[,2])),
                       ncol = length(unique(coordinates(dem)[,1])), byrow = TRUE)

  snowcodepic = raster(snowcodepic, min(coordinates(dem)[,1]),
                       max(coordinates(dem)[,1]),
                       min(coordinates(dem)[,2]),
                       max(coordinates(dem)[,2]))
  extent(snowcodepic) = extent(dem)

  # Project final Raster into WGS coordinate plane
  projection(snowcodepic) <- projection(dem)

  return(snowcodepic)
}



### KRIGING ###


# Simple Kriging with varying local means
kriging = function(train, test, maxel = TRUE, maxout = TRUE, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    coordinates(train) = c("LONGITUDE", "LATITUDE")
    proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train = spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test = rasterToPoints(test, spatial = TRUE)
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) = c("LONGITUDE", "LATITUDE")
      proj4string(test) = proj4string(train)
    }
  }

  # Restrct elevation trend estimates to be no greater than the trend value
  # of the highest elevation station location.
  if(maxel){
    max.el <- max(train$ELEVATION)
    test$ELEVATION[test$ELEVATION > max.el] <- max.el
  }

  # Fit Simple Kriging to Residuals
  gfit = lm(log(RESPONSE)~ELEVATION, data = train)
  train$RESIDUALS = gfit$residuals
  krigetest = krige(RESIDUALS~1, locations = train, newdata = test, beta = 0, ...)

  # Store max from training set if requested for high valued restrictions
  if(maxout){max.out <- max(train$RESPONSE)}


  # Add kriging coefficients and exponentiate
  if(israster){
    g.est = gfit$coefficients[1] + gfit$coefficients[2]*test$ELEVATION
    krigetest$var1.pred = exp(g.est + krigetest$var1.pred)

    # Restrict high values if requested
    if(maxout){
      krigetest$var1.pred[krigetest$var1.pred > max.out] <- max.out
    }

    return(rasterFromXYZ(krigetest))
  }else{
    # Estimate for witheld stations
    g.est = gfit$coefficients[1] + gfit$coefficients[2]*test$ELEVATION
    g.est2 = exp(g.est + krigetest$var1.pred)

    # Adjust high values if requested
    if(maxout){
      g.est2[g.est2 > max.out] <- max.out
    }

    return(g.est2)
  }

} # End the function

# Universal Kriging
kriging2 = function(train, test, maxel = TRUE, maxout = TRUE, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    coordinates(train) = c("LONGITUDE", "LATITUDE")
    proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train = spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test = rasterToPoints(test, spatial = TRUE)
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) = c("LONGITUDE", "LATITUDE")
      proj4string(test) = proj4string(train)
    }
  }

  # Restrict elevation values to max elevation observed in the data set. This helps
  # control for the effects of extrapolation.
  if(maxel){
    max.el <- max(train$ELEVATION)
    test$ELEVATION[test$ELEVATION > max.el] <- max.el
  }

  # Fit Universal Kriging
  krigetest = krige(log(RESPONSE)~ELEVATION, locations = train, newdata = test, ...)

  if(maxout){max.out <- max(train$RESPONSE)}

  # Add kriging coefficients and exponentiate
  if(israster){
    krigetest$var1.pred = exp(krigetest$var1.pred)

    # Adjust high values if requested
    if(maxout){
      krigetest$var1.pred[krigetest$var1.pred > max.out] <- max.out
    }

    return(rasterFromXYZ(krigetest))
  }else{
    # Estimate for witheld stations
    g.est2 = exp(krigetest$var1.pred)

    # Adjust high values if requested
    if(maxout){
      g.est2[g.est2 > max.out] <- max.out
    }

    return(g.est2)
  }

} # End the function



### Akima's Linear Triangulation Interpolation


#(Denser grids needed for cross validation)
TRI = function(train, test, NGSL = FALSE, maxout = FALSE, density = 1000){

  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    train = data.frame(train)
    coordinates(train) = c("LONGITUDE", "LATITUDE")
    proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    # Project training set to be identical to that of the raster
    train = spTransform(train, CRS(projection(test))) # Reproject to CRS of raster

    # Preserve the matrix dimensions and extent of the original raster later.
    dims = dim(as.matrix(test))
    ext = extent(test)

    # Convert raster layer to spatial points.
    test = rasterToPoints(test, spatial = TRUE)
  }else{
    # If not a raster, ensure that the test set is converted to a spatial points data frame
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) = c("LONGITUDE", "LATITUDE")
      proj4string(test) = proj4string(train)
    }
  }

  # Determine if modeing raw snow loads or the NGSL
  if(NGSL){
    train$NGSL = train$RESPONSE/train$ELEVATION
    temppredict = akima::interp(train, z = "NGSL", nx = density, ny = density) # See Akima package for details behind this function
  }else{
    temppredict = akima::interp(train, z = "RESPONSE", nx = density, ny = density)
  }

  # Convert predictions to raster and extract raster information for the test DEM
  temppredict = raster(temppredict)
  preds = raster::extract(temppredict, test)

  if(israster){
    # Convert NGSL to raw snow loads using the DEM elevations
    if(NGSL){preds = preds*test$layer}
    # Return results in same form as input, note that raster fills by ROW not Column
    preds = matrix(preds, nrow = dims[1], ncol = dims[2], byrow = TRUE)

    # Extent defined at beginning of function.
    preds = raster::raster(preds, ext@xmin, ext@xmax, ext@ymin, ext@ymax,
                           crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }else{
    # Convert NGSL to raw snow loads using the DEM elevations
    if(NGSL){preds = preds*test$ELEVATION}
  }

  if(maxout){
    max.out <- max(train$RESPONSE)
    preds[preds > max.out] <- max.out
  }

  return(preds)
}



### REGULAR LINEAR REGRESSION


# Not well documented because it is not a primary method
lmPred = function(train, test, maxel = TRUE, maxout = TRUE, islog = TRUE){
  if(islog){
    tpred = lm(log(RESPONSE) ~ ELEVATION, data = train)
  }else{
    tpred = lm(RESPONSE ~ ELEVATION, data = train)
  }

  # Restrict elevation to the highest valued input station
  if(maxel){
  max.el <- max(train$ELEVATION)
  test$ELEVATION[test$ELEVATION > max.el] <- max.el
  }

  preds = predict(tpred, test)

  # Exponentiate if using log transform
  if(islog){preds = exp(preds)}

  # Restrict predictions if requested.
  if(maxout){
    max.out <- max(train$RESPONSE)
    preds[preds > max.out] <- max.out
  }

  return(preds)
}
#===========================================================================================================================
# CROSS VALIDATION


### FUNCTION TO TEST 9 PRISM weight dimensions ###


gridTest = function(k, tdata, rm, a, b,
                    zm, zx, Fd,
                    basin = NULL, d = 1, ...){

  values = matrix(nrow = length(rm)*length(a)*length(b)*length(zm)*length(zx)*length(d)*length(Fd),
                  ncol = 10)

  #print(dim(values))
  count = 1
  for(i in rm){
    for(j in a){
      for(l in b){
        for(m in zm){
          for(n in zx){
            for(o in d){
              for(q in Fd){
   #             for(w in cluster){
                  temp = crs_val(k = k, method = "PRISM",
                                 tdata = tdata, rm = i, a = j, b = l,
                                 zm = m, zx = n, cluster = m,
                                 Fd = q, d = o, ...)
                  abserr <- mean(abs(temp$ERROR))
                  mederr <- median(abs(temp$ERROR))
                  relerr <- mean(abs(tdata$RESPONSE - temp$PREDICT)/tdata$RESPONSE)
                  values[count,] = c(i, j, l, m, n, o, q,
                                     abserr, mederr, relerr)
                  #print(count)
                  count = count + 1
  #              }
              }
            }
          }
        }
      }
    }
  }
  colnames(values) = c("rm", "a", "b", "zm", "zx", "d", "Fd", "ABSERR", "MEDERR", "RELERR")
  return(values)
}



### Universal Cross Validation Function ###



# Method is a character vector that corresponds to one of the following:
# PRISM
# Kriging
# IVD
# TRI
# LinReg
crs_val = function(k, tdata, method, maxel=FALSE, maxout = FALSE, ...){

  # Make sure user specfies a method that we have defined
  if(!is.element(method, c("PRISM", "Kriging", "IVD", "IVD3", "TRI", "LinReg", "Kriging2", "Kriging0", "Kriging.loess"))){
    error("Please select one of the following methods: PRISM, Kriging, Kriging2, IVD, IVD3, TRI, or LinReg")
  }
  # Split data into k groups
  groups = rep(1:k, length = nrow(tdata))
  # Randomize group assignment
  groups = sample(groups)

  # Assign each station to a group
  tdata$GROUPS = groups
  tdata$PREDICT = 0

  # Do we want restrictions on the extrapolated predictions?
  # We will store these results.
  maxel.2 <- maxel
  maxout.2 <- maxout

  for(i in 1:k){
    test = tdata[tdata$GROUPS == i,]
    train = tdata[tdata$GROUPS != i,]

    # Predict accoring to the appropriate method
    if(method == "Kriging"){tdata$PREDICT[tdata$GROUPS == i] = kriging(train, test,
                                                                       maxel = maxel.2,
                                                                       maxout = maxout.2, ...)}

    if(method == "Kriging2"){tdata$PREDICT[tdata$GROUPS == i] = kriging2(train, test,
                                                                         maxel = maxel.2,
                                                                         maxout = maxout.2, ...)}

    if(method == "Kriging0"){tdata$PREDICT[tdata$GROUPS == i] = kriging.0(train, test, layer.low = 1275, ...)}

    if(method == "Kriging.loess"){tdata$PREDICT[tdata$GROUPS == i] = kriging.loess(train, test, layer.low = 1275, ...)}

    if(method == "PRISM"){tdata$PREDICT[tdata$GROUPS == i] = prism2(train, test, basin = test$HUC,
                                                                    maxout = maxout.2, ...)}

    if(method == "IVD"){tdata$PREDICT[tdata$GROUPS == i] = IVD2(train, test,
                                                                maxout = maxout.2, ...)}

    if(method == "IVD3"){tdata$PREDICT[tdata$GROUPS == i] = IVD3(train, test, ...)}

    if(method == "LinReg"){tdata$PREDICT[tdata$GROUPS == i] = lmPred(train, test, maxel = maxel.2,
                                                                     maxout = maxout.2, ...)}

    if(method == "TRI"){tdata$PREDICT[tdata$GROUPS == i] = TRI(train, test,
                                                               maxout = maxout.2, ...)}
  }

  # Calculate error (same way we would calculate a residual)
  tdata$ERR = tdata$RESPONSE - tdata$PREDICT

  return(data.frame(STATION_NAME = tdata$STATION_NAME, PREDICT = tdata$PREDICT, ERROR = tdata$ERR))
}



### SNOWLOAD "CROSS VALIDATION"



# Predict at Snow sites given by old utah snow report
#snowlaw = read.csv("F:/Snow Project/Data/Current Utah Snow Loads/Utahsnowlaw.csv")
snlwCV2 = function(tempdata, snowlaw, state = utah){
  # Extract elevation values
  if(is.null(tempdata$COUNTYNBR)){
    coordinates(tempdata) = c("LONGITUDE", "LATITUDE")
    proj4string(tempdata) = projection(state)
    tempdata$COUNTYNBR = over(tempdata, state)$COUNTYNBR
  }

  # Convert from meters to thousands of feet
  elev = (tempdata$ELEVATION*3.28084)/1000

  # Determine county of each point.
  tempc = tempdata$COUNTYNBR

  snowcodepic.st = vector(mode = "numeric", length = nrow(tempdata))
  for(i in 1:nrow(tempdata)){
    if(is.na(tempc[i])){
      snowcodepic.st[i] = NA
    }else{
      # Set snowload equal to base snowload if elevation in question lies below the base elevation
      if(elev[i] - snowlaw$A_0[tempc[i]] <= 0){
        snowcodepic.st[i] = snowlaw$P_0[tempc[i]] * .04788
      }else{
        # Otherwise, compute the equation as defined by the Utah snow code for that county
        snowcodepic.st[i] = ((snowlaw$P_0[tempc[i]]^2 + snowlaw$S[tempc[i]]^2*(elev[i] - snowlaw$A_0[tempc[i]])^2)^.5)*0.04788
      }
    }
  }

  # Convert to KPA from PSF
  ERR = tempdata$RESPONSE - snowcodepic.st
  tfinal = data.frame(PREDICT = snowcodepic.st, ERROR = ERR)
  return(tfinal)
}




#===========================================================================================================================
# Graphics

# Function to plot snow load maps given an input state
plotSL = function(trast, state, tbreaks = c(0,30,60,90,120,150,200,300, 600) * 0.0478803, tcol = brewer.pal(8,"Blues"),
                  plotmap = TRUE, ylegend = TRUE, labind = c(2,4,6,8) , ...){
  test = mask(trast, state)
  plot(test, breaks = tbreaks, col = tcol, legend = FALSE, ...)
  if(ylegend){
    plot(test, legend.only=TRUE,  col = tcol,
         legend.width=1, legend.shrink=0.75, breaks = tbreaks,
         axis.args=list(at=labind, labels=round(labind,1)), legend.args=list(text='kPa', side=3, line=1, cex=1),box = FALSE)
  }
  if(plotmap){plot(state, border = "grey", add=TRUE)}
}



#============================================================================================================================
# Extra function that tries non-linear NGSL

# c1 - lower elevation IDW exponent
# c2 - upper elevation IDW exponent
# layer - splits the data into two regions (layer = 0 performs regular ivd using c2)
# NGSL - exponent of NGSL adjustement (NGSL = 0 performs raw inverse distance weighting)
IVD3 = function(stationD, DEM, layer = 1219.2, c1 = 2, c2 = 6, NGSL = 1, debug.level = 0){
  #require(fields)
  #require(raster)
  #require(sp)
  #require(gstat)
  ### READ IN RASTER DATA ###

  # Option to handle either a raster (map) or data.frame (cross validation)
  if(class(DEM) == "RasterLayer"){
    # Convert to spatialpointsdataframe
    israster = TRUE
    DEM = rasterToPoints(DEM, spatial = TRUE)
    names(DEM) = "ELEVATION"

    # If the raster doesn't have a predefined projection. Given it geographic coordinate reference.
    if(is.na(projection(DEM))){
      projection(DEM) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    }
  }else{
    if(is.null(DEM$ELEVATION) | is.null(DEM$LONGITUDE) | is.null(DEM$LATITUDE)){
      stop("ELEVATION, LATITUDE, and LONGITUDE, must be defined in order to proceed")
    }
    israster = FALSE
    coordinates(DEM) = c("LONGITUDE", "LATITUDE")
    # Ensure that our data frame has a geographic coordinate reference system if not already a raster.
    proj4string(DEM) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Project into geographic coordinates
  }

  # Ensure station data in acceptable format
  if(class(stationD) != "data.frame" & class(stationD) != "SpatialPointsDataFrame"){
    stop("DEM must be of class \"data.frame\" or \"SpatialPointsDataFrame\"")
  }

  # Make sure variables of interest are specified
  if(is.null(stationD$RESPONSE) | is.null(stationD$ELEVATION) | is.null(stationD$LONGITUDE)
     | is.null(stationD$LATITUDE)){
    stop("The following variables MUST be included in StationD: RESPONSE, ELEVATION, LONGITUDE, LATITUDE")
  }

  # Convert to spatial points data frame if not already done
  if(class(stationD) == "data.frame"){
    sp::coordinates(stationD) = c("LONGITUDE", "LATITUDE")
  }

  # Define projections to be the same if stationD does not have a projection pre-defined.
  if(is.na(proj4string(stationD))){
    proj4string(stationD) = proj4string(DEM)
  }

  # If both projections are defined yet different. Transform to DEM projection
  if(proj4string(stationD) != proj4string(DEM)){
    stationD = spTransform(stationD, CRS(proj4string(DEM)))
  }

  # Store elevation information separately
  dem.elevation = DEM$ELEVATION

  # Prior to splitting define the NGSL variable as the response if the user requests
  stationD$RESPONSE = stationD$RESPONSE/(stationD$ELEVATION^NGSL)

  # Split stations into two layers (above and below 4000 ft (1219.2m) in the idaho method)
  stationD.low = stationD[stationD$ELEVATION < layer,]
  stationD.high = stationD[stationD$ELEVATION >= layer,]


  # Obtain the bi-layer results
  if(nrow(stationD.low) > 0){
    temp2 = idw(RESPONSE~1, locations = stationD.low, newdata = DEM, idp = c1, debug.level = debug.level)$var1.pred
  }else{
    temp2 = NULL
    warning("No stations exist in the lower layer, using only one layer with exponent c2...", call. = FALSE)
  }
  if(nrow(stationD.high) > 0){
    temp6 = idw(RESPONSE~1, locations = stationD.high, newdata = DEM, idp = c2, debug.level = debug.level)$var1.pred
  }else{
    temp6 = NULL
    warning("No stations exist in the upper layer, using only one layer with exponent c1...", call. = FALSE)
  }

  # Store results by layer, assuming both layers are defined.
  results = vector("numeric", length(dem.elevation))
  if(is.null(temp2)){
    results = temp6
  }else if(is.null(temp6)){
    results = temp2
  }else{
    results[dem.elevation < layer] = temp2[dem.elevation < layer]
    results[dem.elevation >= layer] = temp6[dem.elevation >= layer]
  }

  # Return ground snow load predictions as opposed to simply the NGSL.
  results = results*(dem.elevation^NGSL)

  # Return results as a raster if the input DEM was a raster, if not return the vector
  if(israster){
    DEM$ELEVATION = NULL
    DEM$RESULTS = results
    DEM = rasterFromXYZ(DEM)
    return(DEM)
  }

  return(results)
}


# Simple Kriging with varying local means - including a conservative adjustment using the kriging prediction variance
kriging_adj = function(train, test, sd = 2, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    coordinates(train) = c("LONGITUDE", "LATITUDE")
    proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train = spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test = rasterToPoints(test, spatial = TRUE)
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) = c("LONGITUDE", "LATITUDE")
      proj4string(test) = proj4string(train)
    }
  }


  # Fit Simple Kriging to Residuals
  gfit = lm(log(RESPONSE)~ELEVATION, data = train)
  train$RESIDUALS = gfit$residuals
  krigetest = krige(RESIDUALS~1, locations = train, newdata = test, beta = 0, ...)



  # Add kriging coefficients and exponentiate
  if(israster){
    g.est = gfit$coefficients[1] + gfit$coefficients[2]*test$layer
    krigetest$var1.pred = exp(g.est + krigetest$var1.pred + sd*sqrt(krigetest$var1.var))
    return(rasterFromXYZ(krigetest))
  }else{
    # Estimate for witheld stations
    g.est = gfit$coefficients[1] + gfit$coefficients[2]*test$ELEVATION
    g.est2 = exp(g.est + krigetest$var1.pred)
    return(g.est2)
  }

} # End the function



# Kriging for April 1st SWE data (handles 0 values)
# layer.low - layer below which all SWE values are predicted as 0.
# loess = TRUE - use a loess smoothing curve to fit rather than a linear model
kriging.0 = function(train, test, layer.low = 1275, maxout = TRUE, maxel = TRUE, ...){
  if(class(train)[1] != "SpatialPointsDataFrame"){
    # Create spatial points data frame2
    coordinates(train) = c("LONGITUDE", "LATITUDE")
    proj4string(train) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster = is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train = spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test = rasterToPoints(test, spatial = TRUE)
    test$ELEVATION = test$layer
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) = c("LONGITUDE", "LATITUDE")
      proj4string(test) = proj4string(train)
    }
  }

  # Restrict elevation values to max elevation observed in the data set. This helps
  # control for the effects of extrapolation.
  if(maxel){
    max.el <- max(train$ELEVATION)
    test$ELEVATION[test$ELEVATION > max.el] <- max.el
  }

  # Fit Universal Kriging
  krigetest = krige(log(WESD + 1)~ELEVATION, locations = train, newdata = test, ...)

  if(maxout){max.out <- max(train$WESD)}

  # Add kriging coefficients and exponentiate
  if(israster){
    krigetest$var1.pred = exp(krigetest$var1.pred)

    # Adjust high values if requested
    if(maxout){
      krigetest$var1.pred[krigetest$var1.pred > max.out] <- max.out
    }

    # Set values in the lower layer equal to 0
    krigetest$var1.pred[test$layer < layer.low] <- 0

    return(rasterFromXYZ(krigetest))
  }else{
    # Estimate for witheld stations
    g.est2 = exp(krigetest$var1.pred)

    # Adjust high values if requested
    if(maxout){
      g.est2[g.est2 > max.out] <- max.out
    }

    g.est2[test$ELEVATION < layer.low] <- 0

    return(g.est2)
  }

} # End the function

# Kriging for April 1st SWE data (handles 0 values)
# layer.low - layer below which all SWE values are predicted as 0.
# loess = TRUE - use a loess smoothing curve to fit rather than a linear model
kriging.loess = function(train, test, layer.low = 1275, maxout = TRUE, ...){
  require(dplyr) # This package is necessary for manipulations

  # Create spatial points data frame
  coordinates(train) <- c("LONGITUDE", "LATITUDE")
  proj4string(train) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

  # If prediction layer is a raster, set projection of training data equal to that of the raster
  israster <- is.element(class(test)[1], c("Raster", "RasterLayer"))
  if(israster){
    train <- spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
    test <- rasterToPoints(test, spatial = TRUE)
    test$ELEVATION <- test$layer
  }else{
    if(class(test)[1] != "SpatialPointsDataFrame"){
      coordinates(test) <- c("LONGITUDE", "LATITUDE")
      proj4string(test) <- proj4string(train)
    }
  }

  # Set all elevations greater than the maximum observed elevation
  # equal to the maximum observed elevation.
  maxel <- max(train$ELEVATION)
  test$ELEVATION[test$ELEVATION > maxel] <- maxel

  # Cap predictions if requested
  if(maxout){max.out <- max(train$WESD)}

  # Extract coefficients.
  fit3 <- loess(log(WESD+1) ~ ELEVATION, data = train,
                control = loess.control(surface = "direct"))

  train$RESIDUAL <- fit3$residuals
  # bw <- npregbw(formula=log(WESD+1)~ELEVATION, data = as.data.frame(train))
  # fit3 <- npreg(bw, newdata = as.data.frame(train), residuals = TRUE)
  # train$RESIDUAL <- fit3$resid


  # Fit Universal Kriging
  krigetest = krige(RESIDUAL~1, locations = train, newdata = test, ...)

  # Add kriging coefficients and exponentiate
  if(israster){
    g.est = predict(fit3, newdata = as.data.frame(test))
    krigetest$var1.pred = exp(g.est + krigetest$var1.pred) - 1

    # Set values in the lower layer equal to 0
    krigetest$var1.pred[test$layer < layer.low] <- 0

    # No negative SWE values so we will set them to 0.
    krigetest$var1.pred[krigetest$var1.pred < 0] <- 0

    # Adjust high values if requested
    if(maxout){
      krigetest$var1.pred[krigetest$var1.pred > max.out] <- max.out
    }

    return(rasterFromXYZ(krigetest))
  }else{
    g.est = predict(fit3, newdata = as.data.frame(test))
    g.est2 = exp(g.est + krigetest$var1.pred) - 1

    # Set values in the lower layer equal to 0
    g.est2[test$ELEVATION < layer.low] <- 0

    g.est2[g.est2 < 0] <- 0

    # Adjust high values if requested
    if(maxout){
      g.est2[g.est2 > max.out] <- max.out
    }

    return(g.est2)
  }

} # End the function




### NEW PRISM FUNCTION FOR PICTURES
prism2.2 = function(stationD, DEM, rm = 10, a = 2, b = 1,
                    zm = 200, zx = 2500, Fd = .8, cluster = -1,
                    basin = NULL, d = 1, feet = FALSE, islog = TRUE, NGSL = FALSE){
  require(fields)
  require(raster)
  require(sp)
  ### READ IN RASTER DATA ###
  # DEM MUST be a raster
  if(class(DEM) != "RasterLayer" & class(DEM) != "data.frame"){
    stop("DEM must be of class \"RasterLayer\" or \"data.frame\"")}

  # Option to handle either a raster (map) or data.frame (cross validation)
  if(class(DEM) == "RasterLayer"){
    # Extract unique coordinate extent for RasterLayer
    gridC = sp::coordinates(DEM)
    lonCen = unique(gridC[,1])
    latCen = unique(gridC[,2])

    # Extract extent of raster layer
    temp.ext = extent(DEM)

    # Convert raster layer to vector (fills a matrix BY ROW)
    DEM = as.vector(DEM)
    israster = TRUE
  }else{
    if(is.null(DEM$ELEVATION) | is.null(DEM$LONGITUDE) | is.null(DEM$LATITUDE)){
      stop("ELEVATION, LATITUDE, and LONGITUDE, must be defined in order to proceed")
    }
    gridC = data.frame(LONGITUDE = DEM$LONGITUDE,
                       LATITUDE = DEM$LATITUDE)
    # Convert feet to meters if necessary
    if(feet){
      DEM = DEM$ELEVATION * .3048
    }else{
      DEM = DEM$ELEVATION
    }
    israster = FALSE
  }

  ### READ IN STATION DATA ###
  if(is.null(stationD$RESPONSE) | is.null(stationD$ELEVATION) | is.null(stationD$LONGITUDE)
     | is.null(stationD$LATITUDE)){
    stop("The following variables MUST be included in StationD: RESPONSE, ELEVATION, LONGITUDE, LATITUDE")
  }

  # If Station elevations are in feet, convert to meters to match the DEM
  if(feet){
    explanatory = stationD$ELEVATION * .3048 # Explanatory variable (elevation values in meters)
  }else{
    explanatory = stationD$ELEVATION
  }

  # Set snow loads as the response variable
  response = stationD$RESPONSE

  # Determine if modeling NGSL or raw snow loads
  if(NGSL){
    response = response/explanatory
  }

  # Determine if modeling log of values (either NGSL or raw snow loads)
  if(islog){
    response = log(response)
  }

  stationC = data.frame(LONGITUDE = stationD$LONGITUDE, LATITUDE = stationD$LATITUDE) # Lon/Lat coordiantes stored as dataframe

  ### DETERMINE WEIGHTING SCHEME ###

  ## FIND DISTANCE WEIGHTS ##
  getDistWeight2 = function(gridC, stationC, rm, a){

    # Calculate distance with function from "fields" package
    # Note that we transpose the results so weights are found in the columns, not the rows
    dist = t(round(fields::rdist.earth(as.matrix(gridC), as.matrix(stationC), miles = FALSE), digits=4)) # Avoid singularity issues in getDistWeight by rounding

    # Set all stations within the minimum radius of influence to a weight of 1
    distWeight = dist - rm
    distWeight[distWeight <= 0] = 1

    # For all other stations, compute the weighting scheme.
    distWeight = distWeight^a
    distWeight = 1 / distWeight
    distWeight[distWeight > 1] = 1 # Set all values greater than 1 to the max

    # Normalize each weighting vector to sum to 1
    distWeight = scale(distWeight, center = FALSE, scale = colSums(distWeight))

    return(distWeight)
  }

  ## Find Elevation Weights ##
  getElevWeight2 = function(gridElv, stationElv, b, zm, zx){

    # Store elevation differences between each gridcell and the stations in the columns of a matrix
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
    tempDiff = scale(tempDiff, center = FALSE, scale = colSums(tempDiff))

    # In the event that ALL weights are equal to 0, scale will return NAN values for that column
    # we handle this by setting all NAN values equal to 0, meaning that the elevation difference
    # will have NO effect at this point (and the weight will NOT sum to unity)
    tempDiff[is.nan(tempDiff)] = 0

    return(tempDiff)
  }

  ## Find Cluster Weights ##
  getClusterWeight = function(lonlat, stationElv, r, p){
    distDiff = rdist.earth(as.matrix(lonlat), miles = FALSE)
    h = (.2*r - distDiff) / (.2*r)
    h[h < 0] = 0 # Any stations outside the 20% of the radius of influence set to 0
    h[h > 1] = 1
    h = h - diag(diag(h)) # Don't want to consider distance from a station to itself

    # Calculate s_ij matrix
    elevDiff = abs(outer(stationElv, stationElv, "-")) - p
    elevDiff[elevDiff < 0] = 0
    # Calculate v_i
    v = (p - elevDiff)/p
    v[v < 0] = 0
    v[v > 1] = 1
    v = v - diag(diag(v)) # Don't want to consider the elevation difference between a station and itself

    Sc = h*v
    Sc = rowSums(Sc) + 1 # +1 avoids dividing by 0, or getting values bigger than one for near zero factors
    W_c  = 1 / Sc
    W_c = W_c / sum(W_c)
    return(as.vector(W_c))
  }

  ## Find Aspect Weights ##
  getaspect = function(basinDEM, stationBasin, d){

    #huc12.basin = as.numeric(as.character(as.vector(basinDEM)))
    #huc12.station = as.numeric(as.character(stationBasin))
    #huc10.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 100)
    #huc10.station = floor(as.numeric(as.character(stationBasin)) / 100)
    huc8.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 10000)
    huc8.station = floor(as.numeric(as.character(stationBasin)) / 10000)
    huc6.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 1000000)
    huc6.station = floor(as.numeric(as.character(stationBasin)) / 1000000)
    huc4.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 100000000)
    huc4.station = floor(as.numeric(as.character(stationBasin)) / 100000000)
    huc2.basin = floor(as.numeric(as.character(as.vector(basinDEM))) / 10000000000)
    huc2.station = floor(as.numeric(as.character(stationBasin)) / 10000000000)

    # Store elevation differences between each gridcell and the stations in the columns of a matrix
    #sameBasin12 = t(outer(huc12.basin, huc12.station, "=="))
    #sameBasin10 = t(outer(huc10.basin, huc10.station, "=="))
    sameBasin8 = t(outer(huc8.basin, huc8.station, "=="))
    sameBasin6 = t(outer(huc6.basin, huc6.station, "=="))
    sameBasin4 = t(outer(huc4.basin, huc4.station, "=="))
    sameBasin2 = t(outer(huc2.basin, huc2.station, "=="))

    # Stations receive weights proportional to the number of watersheds they share.
    # C controls how drastic the weighting scheme will be
    tempweight = ((sameBasin8+sameBasin6+sameBasin4+sameBasin2 + 1)/5)^d

    return(scale(tempweight, center=FALSE, scale = colSums(tempweight))) # stations in the row, cells in the columns
  }

  Wd = getDistWeight2(gridC, stationC, rm, a) # Get Distance Weights
  Wz = getElevWeight2(DEM, explanatory, b, zm, zx) # Get Elevation Weights

  # Get cluster weight only if user specifies it
  if(length(cluster) > 1){
    stop("Cluster input must be single numeric value")
  }

  if(cluster > 0){
    # Use same minimum radius of influence used in distance weighting
    Wc = getClusterWeight(stationC, explanatory, rm, cluster)
  }else{
    Wc = 1 # Apply no cluster weighting unless the user specifies
  }

  # Get basin weight only if user specifies it
  if(!is.null(basin)){
    if(is.null(stationD$HUC)){
      stop("\'HUC\' must be a variable in StationD data frame when \'basin\' is defined")
    }
    #print(paste("Fitting aspect weighting with d =", d, "..."))
    Wf = getaspect(basin, stationD$HUC, d)
  }else{
    Wf = 1 # Apply no cluster weighting unless the user specifies
  }

  # Apply weighting equation from Daley et.al. 2008
  weightV = Wc*sqrt(Fd*Wd^2 + (1-Fd)*Wz^2)*Wf
  weightV = scale(weightV, center=FALSE, scale = colSums(weightV))

  # Remove stuff we no longer need to save on memory
  ##################################

  # Preallocate for predictions
  predictP = vector("numeric", length = length(DEM))
  slope = vector("numeric", length = length(DEM))
  intercept = vector("numeric", length = length(DEM))
  for(i in 1:length(predictP)){
    if(is.na(DEM[i])){
      predictP[i] = NA
      #slope[i] = NA
      #intercept[i] = NA
    }else{
      # Run linear regression on individual gridcell, with the respective weighting vector
      tempReg = lsfit(as.matrix(explanatory,ncol = 1), response, weightV[,i])
      # Predict based on centroid elevation of current grid (retransform data)
      intercept[i] = tempReg$coefficients[1]
      slope[i] = tempReg$coefficients[2]
      predictP[i] = tempReg$coefficients[1] + tempReg$coefficients[2]*DEM[i]

    }
  }

  # Retransform Results if log transform is being used
  if(islog){predictP = exp(predictP)}

  # Extract snow loads from NGSL if necessary
  if(NGSL){predictP = predictP*DEM}

  return(list(predictP,intercept,slope,weightV,gridC, DEM))
}




# # OLD KRIGING FUNCTION WHERE 0 VALUES WERE NOT USED IN LINEAR MODEL FIT
#
# # Kriging for April 1st SWE data (handles 0 values)
# # layer.low - layer below which all SWE values are predicted as 0.
# # loess = TRUE - use a loess smoothing curve to fit rather than a linear model
# kriging.0 = function(train, test, layer.low = 1275, ...){
#   require(dplyr) # This package is necessary for manipulations
#
#   # Filter 0 values prior to making spatial data frames
#   # train.0 <- train %>% filter(RESPONSE > 0)
#
#   # Create spatial points data frame
#   coordinates(train) <- c("LONGITUDE", "LATITUDE")
#   proj4string(train) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#
#   # coordinates(train.0) <- c("LONGITUDE", "LATITUDE")
#   # proj4string(train.0) <- proj4string(train)
#
#   # Create spatial data frame to extact coefficients from GLS in gstat
#   testD <- data.frame(ELEVATION = c(0,1), LATITUDE = c(39.5, 39.6),
#                       LONGITUDE = c(-111,-111.1))
#   coordinates(testD) <- c("LONGITUDE", "LATITUDE")
#   proj4string(testD) <- proj4string(train)
#
#   # If prediction layer is a raster, set projection of training data equal to that of the raster
#   israster <- is.element(class(test)[1], c("Raster", "RasterLayer"))
#   if(israster){
#     train <- spTransform(train, CRS(projection(test))) # Reproject to CRS of raster
#     test <- rasterToPoints(test, spatial = TRUE)
#     test$ELEVATION <- test$layer
#   }else{
#     if(class(test)[1] != "SpatialPointsDataFrame"){
#       coordinates(test) <- c("LONGITUDE", "LATITUDE")
#       proj4string(test) <- proj4string(train)
#     }
#   }
#
#   # # Extract coefficients.
#   # fit3 <- gstat(id = "tkrige",
#   #               formula = log(WESD+1) ~ ELEVATION,
#   #               locations = train.0, ...)
#   #
#   # # Predict slope and intercept
#   # coef <- predict(fit3, newdata = testD, BLUE = TRUE)
#   # t.int <- coef$tkrige.pred[1]
#   # t.slope <- coef$tkrige.pred[2]-coef$tkrige.pred[1]
#   #
#   # # Obtain residuals for ALL values (including 0 valued stations)
#   # train$PREDICT <- t.int + t.slope*train$ELEVATION
#   # train$RESIDUAL <- log(train$RESPONSE + 1) - train$PREDICT
#   #
#   # # Fit Universal Kriging
#   # krigetest = krige(RESIDUAL~1, locations = train, newdata = test, ...)
#
#   krigetest = krige(log(WESD+1) ~ ELEVATION, locations = train,
#                     newdata = test, ...)
#
#   # Add kriging coefficients and exponentiate
#   if(israster){
#     g.est = t.int + t.slope*test$ELEVATION
#     krigetest$var1.pred = exp(g.est + krigetest$var1.pred) - 1
#
#     # Set values in the lower layer equal to 0
#     krigetest$var1.pred[test$layer < layer.low] <- 0
#
#     # No negative SWE values so we will set them to 0.
#     krigetest$var1.pred[krigetest$var1.pred < 0] <- 0
#
#     return(rasterFromXYZ(krigetest))
#   }else{
#     g.est = t.int + t.slope*test$ELEVATION
#     g.est2 = exp(g.est + krigetest$var1.pred) - 1
#
#     # Set values in the lower layer equal to 0
#     g.est2[test$ELEVATION < layer.low] <- 0
#
#     g.est2[g.est2 < 0] <- 0
#
#     return(g.est2)
#   }
#
# } # End the function
#
#
#
#
