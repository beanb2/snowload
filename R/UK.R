#' Title: Universal Kriging - Ecoregions
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
UK = function(train, test, maxel = TRUE, maxout = TRUE, ...){
  ### READ IN RASTER DATA ###
  # test MUST be a raster
  if(class(test)[1] != "SpatialPointsDataFrame" | 
     class(train)[1] != "SpatialPointsDataFrame"){
    stop("Both inputs must be of class \"SpatialPointsDataFrame\"")
  } 
  
  if(raster::projection(train) != raster::projection(test)){
    stop("Projections of inputs and outputs must be identical")
  }
  
  # Ensure station data meet criteria. 
  if(is.null(train$RESPONSE) | is.null(train$ELEVATION) | 
     is.null(train$REGION)){
    stop("The following variables MUST be included in train: 
         RESPONSE, ELEVATION, REGION")
  }
  
  if(is.null(test$ELEVATION) | is.null(test$REGION)){
    stop("The following variables MUST be included in test: 
          ELEVATION, REGION")
  }
  
  
  tpredict <- vector("numeric", nrow(test))
  for(i in 1:nrow(test)){
    # Subset all stations in the same level III ecoregion as the area of interest or within 50km
    isRegion <- as.character(train$REGION) == as.character(test$REGION[i])
    tdist <- fields::rdist.earth(train@coords, 
                                 data.frame(x = test@coords[i, 1], 
                                            y = test@coords[i, 2]), 
                                 miles = FALSE)
    isClose <- tdist <= 50
    # If there are 50 of these stations, suspend the search, 
    # if not, keep looking...
    if(sum(as.logical(isRegion | isClose)) < 50){
      # Move up one ecoregion and check again
      rII <- gsub(train$REGION, pattern = "[[:punct:]][[:digit:]]+$", 
                  replacement = "")
      sII <- gsub(test$REGION[i], pattern = "[[:punct:]][[:digit:]]+$", 
                  replacement = "")
      nextRegion <- rII == sII
      if(sum(as.logical(isRegion | isClose | nextRegion)) < 50){
        # If necessary, move up one ecoregion and check again
        rI <- gsub(rII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")
        sI <- gsub(sII, pattern = "[[:punct:]][[:digit:]]+$", replacement = "")
        nextnextRegion <- rI == sI
        if(sum(as.logical(isRegion | isClose | nextRegion | nextnextRegion)) < 50){
          # If there are still not enough stations, find the closest 
          # remaining stations necessary to achieve a sample size of 50. 
          tord <- order(tdist, decreasing = TRUE)[!(isRegion | isClose | 
                                                      nextRegion | nextnextRegion)]
          tord <- tord[1:(50-sum(as.logical(isRegion | isClose | nextRegion | 
                                   nextnextRegion)))]
          finalRegion <- rep(FALSE, length(rI))
          finalRegion[tord] <- TRUE
          
          train.sub <- train[as.logical(isRegion | isClose | nextRegion | 
                               nextnextRegion | finalRegion), ]
        }else{
          train.sub <- train[as.logical(isRegion | isClose | nextRegion | 
                               nextnextRegion), ]
        }
      }else{
        train.sub <- train[as.logical(isRegion | isClose | nextRegion), ]
      }
    }else{
      train.sub <- train[as.logical(isRegion | isClose), ]
    }

    # Fit Universal Kriging
    krigetest = gstat::krige(log(RESPONSE + 1) ~ ELEVATION, 
                             locations = train.sub, newdata = test[i, ], ...)
    tpredict[i] <- exp(krigetest$var1.pred) - 1
  } # End the for-loop
  
  # If requested, cap the predictions. 
  if(maxout){
    tpredict[tpredict > max(train$RESPONSE)] <- max(train$RESPONSE)
  }
  
  test$predict <- tpredict
  
  return(test)
} # End the function