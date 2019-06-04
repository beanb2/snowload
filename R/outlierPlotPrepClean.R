#=============================================================================
# Function to clean the downloaded snow load data. 
# Author: Brennan Bean
# Created 10-8-2018
# Last Updated: 10-8-2018
# Ensure that working directory is set to source file location prior to 
# running this script. 
#=============================================================================

#=============================================================================
# Helpful references. 
# Acceleration due to gravity:
# - https://physics.nist.gov/cgi-bin/cuu/Value?gn 
# Conditional Mutation
# - https://stackoverflow.com/questions/24459752/can-dplyr-package-be-used-for-conditional-mutating
#=============================================================================


#=============================================================================
# Function: Clean the data for a given state. 
# fileName - name of the .csv file that will be cleaned. 
# stations - R data frame with all station meta data
# year1    - the earliest calendar year to consider
# year2    - the most recent calendar year to consider
# remove.outlier - (true/false) should outliers be removed as part of analysis?
# method   - ("Sturm" or "RMCD") the method for converting snow weight from 
#             snow depth. 
# h        - maximum similarity score for the cluster analysis. 
# distAdj  - constant value that scales the dissimilarity score for geographic
#            distance. 
# elevAdj  - constant value that scales elevation dissimilarity
# cluster  - (true/false) should nearby stations be clustered together as part
#            of the analysis? 
#=============================================================================
dataPrep <- function(fileName, stations, year1 = 1967, year2 = 2018, 
                     remove.outlier = TRUE, method = "Sturm", h = 2,
                     distAdj = 4, elevAdj = 50, cluster = TRUE){
  
  # Define the pipeline operator:
  # - https://stackoverflow.com/questions/27386694/using-operator-from-dplyr-without-loading-dplyr-in-r
  `%>%` <- magrittr::`%>%`
  
  # Step 1 - Data Consolidation
  #=============================================================================
  # Note that SNWD needs conversion from mm to cm. 
  # and WESD needs conversion from mm/10 to cm.
  tempDF <- try(read.csv(fileName), TRUE)
  # Return error if data fails to read. 
  if(inherits(tempDF, "try-error")){
    print(paste("No viable file found at ", fileName, " returning NULL...",
                sep = ""))
    return(NULL)
  }
  
  tempDF <- tempDF %>%
    dplyr::mutate(YEAR = as.numeric(as.character(YEAR)),
                  MONTH = as.numeric(as.character(MONTH)),
                  DAY = as.numeric(as.character(DAY)),
                  VALUE = as.numeric(as.character(VALUE)),
                  # value conversion for SNWD and WESD are different. 
                  VALUE = dplyr::if_else(ELEMENT == "SNWD", VALUE*0.1, VALUE*.01),
                  # Define a water year from October to May (move October forward 1 year)
                  wYEAR = dplyr::if_else(MONTH < 9, YEAR, YEAR+1),
                  QFLAG = dplyr::if_else(is.na(QFLAG), " ", as.character(QFLAG)),
                  MFLAG = dplyr::if_else(is.na(MFLAG), " ", as.character(MFLAG)),
                  SFLAG = dplyr::if_else(is.na(SFLAG), " ", as.character(SFLAG))) %>% 
    dplyr::filter(wYEAR >= year1, wYEAR <= year2) %>% # Only keep relevant years
    dplyr::filter(QFLAG == " ") %>% # remove observations failing NDCD quality check
    dplyr::filter(MFLAG != "P") %>% # Filter missing presumed 0.
    dplyr::select(-QFLAG, -MFLAG, -SFLAG) 
  
  if(nrow(tempDF) < 1){
    print(paste("No viable records in ", fileName, " returning NULL...", sep = ""))
    return(NULL)
  }
  
  ### Add station cluster information if requested.  ######################
  stID <- as.character(unique(tempDF$ID))
  
  # Only cluster if there is more than one station in the set. 
  stations.sub <- stations[is.element(stations$ID, stID), ]
  if(cluster & length(stID) > 1){  
    tdist <- fields::rdist.earth(stations.sub[, c("LONGITUDE", "LATITUDE")])
    telev <- abs(outer(stations.sub$ELEVATION, stations.sub$ELEVATION, "-"))
    teco <- outer(stations.sub$NA_L3CODE, stations.sub$NA_L3CODE, "!=")
    # Prevent clusters from occuring across eco regions by setting the eco
    # score beyond the cutting value. 
    teco <- (h+1)*teco 
    
    # Produce "similarity matrix": 
    # - distance of 2 km produces a score of 1
    # - elevation difference of 50m produces a score of 1
    dist.clust <- stats::as.dist((tdist/distAdj) + (telev/elevAdj) + teco)
    
    # Cluster all stations with the farthest neighbors in a cluster having a 
    # score of no more than 2.
    tclust <- stats::hclust(dist.clust, method = "complete")
    cu.id <- stats::cutree(tclust, h = h)
    
    stations.sub$cuID <- cu.id
    stations.sub$cuID <- paste(stations.sub$STATE, cu.id, sep = "-")
  }else{ # Return the original ID variables if clustering is not desired. 
    stations.sub$cuID <- stations.sub$ID
  }
  
  # Now summarize the results based on this information. 
  stations.final <- stations.sub %>% 
    dplyr::filter(ELEVATION > -100) %>% # Remove missing elevation information
    dplyr::group_by(cuID) %>%
    dplyr::summarize(NAME = NAME[1],
                     STATE = STATE[1],
                     LATITUDE = mean(LATITUDE, na.rm = TRUE),
                     LONGITUDE = mean(LONGITUDE, na.rm = TRUE),
                     ELEVATION = mean(ELEVATION, na.rm = TRUE),
                     REGIONIII = unique(NA_L3CODE),
                     REGION_name = unique(NA_L3NAME),
                     climate = unique(climate),
                     numST = dplyr::n())
  
  #############################################################################
  
  tempDF <- dplyr::left_join(tempDF, stations.sub, by = "ID")
  
  # Define the date using the lubridate package and determine the 
  # day of the year (DOY)
  tempDF$DATE <- paste(tempDF$YEAR, tempDF$MONTH, tempDF$DAY, sep = "-")
  tempDF$DATE <- lubridate::as_date(tempDF$DATE)
  tempDF$DOY <- lubridate::yday(tempDF$DATE)
  tempDF$DOYA <- tempDF$DOY
  tempDF$DOYA[tempDF$DOYA > 250] = tempDF$DOYA[tempDF$DOYA > 250] - 366
  # Must avoid 0 value that occurs for leap years.
  tempDF$DOYA[tempDF$DOYA == 0] <- -1
  
  params <- base::data.frame(climate = c("Alpine", "Maritime", "Prairie", "Tundra",
                                         "Taiga", "Ephemeral"),
                             pmax = c(0.5975,	0.5979,	0.594, 0.363, 0.217, 
                                      (0.5975 + 0.5979 + 0.594)/3),
                             po = c(0.2237, 0.2578, 0.2332,	0.2425,	0.217, 
                                    (0.2237 + 0.2578 + 0.2332)/3),
                             k1 = c(0.0012, 0.001, 0.0016, 0.0029, 0, 
                                    (0.0012 + 0.001 + 0.0016)/3),
                             k2 = c(0.0038,	0.0038, 0.0031, 0.0049, 0, 
                                    (0.0038 + 0.0038 + 0.0031)/3))
  
  tempDF.snwd <- tempDF %>% 
    dplyr::filter(ELEMENT == "SNWD") %>%
    dplyr::select(ID, cuID, DATE, wYEAR, DOYA, VALUE, climate) %>%
    dplyr::left_join(., params, by = "climate") %>%
    dplyr::mutate(WEIGHT1 = VALUE*((pmax - po)*(1 - exp((-k1*VALUE) - (k2*DOYA))) + po)*0.09806665,
                  WEIGHT2 = dplyr::if_else(VALUE/2.54 < 22, 
                                           .9 * .04788 * VALUE / 2.54,
                                           ((2.36*VALUE/2.54) - 31.9)*.04788),
                  SNWD = VALUE) %>%
    dplyr::select(-pmax, -po, -k1, -k2, -VALUE)
  
  tempDF.wesd <- tempDF %>% 
    dplyr::filter(ELEMENT == "WESD") %>%
    dplyr::select(ID, cuID, DATE, wYEAR, DOYA, VALUE, climate) %>%
    dplyr::left_join(., params, by = "climate") %>%
    dplyr::mutate(WEIGHT3 = VALUE*0.09806665,
                  WESD = VALUE) %>%
    dplyr::select(-pmax, -po, -k1, -k2, -VALUE, -climate)
  
  tempDF <- dplyr::full_join(tempDF.snwd, tempDF.wesd, by = c("ID", "cuID", "DATE", "DOYA", "wYEAR"))
  
  # Now determine a final weight for each observation, giving preference to direct 
  # calculations of WESD. 
  tempDF$WEIGHT <- tempDF$WEIGHT3
  tempDF$direct <- 1 # Indicator that measurement was calculated directly from WESD
  if(method == "Sturm"){ # By default, use Sturm's method
    tempDF$WEIGHT[is.na(tempDF$WEIGHT)] <- tempDF$WEIGHT1[is.na(tempDF$WEIGHT)]
    tempDF$direct[is.na(tempDF$WEIGHT)] <- 0
    tempDF$WEIGHT1[is.na(tempDF$WEIGHT1)] <- -1 # Replace missing with negative
    tempDF$WEIGHT[is.na(tempDF$WEIGHT)] <- -1 # Replace missing with negative
    tempDF$WEIGHT[tempDF$WEIGHT <= 0 & tempDF$WEIGHT1 > 0] <- 
      tempDF$WEIGHT1[tempDF$WEIGHT <= 0 & tempDF$WEIGHT1 > 0]
    tempDF$direct[tempDF$WEIGHT <= 0 & tempDF$WEIGHT1 > 0] <- 0
  }else{ # Use RMCD if specified. 
    tempDF$WEIGHT[is.na(tempDF$WEIGHT)] <- tempDF$WEIGHT2[is.na(tempDF$WEIGHT)]
    tempDF$direct[is.na(tempDF$WEIGHT)] <- 0
    tempDF$WEIGHT2[is.na(tempDF$WEIGHT2)] <- -1 # Replace missing with negative
    tempDF$WEIGHT[is.na(tempDF$WEIGHT)] <- -1 # Replace missing with negative
    tempDF$WEIGHT[tempDF$WEIGHT <= 0 & tempDF$WEIGHT2 > 0] <- 
      tempDF$WEIGHT2[tempDF$WEIGHT <= 0 & tempDF$WEIGHT2 > 0]
    tempDF$direct[tempDF$WEIGHT <= 0 & tempDF$WEIGHT2 > 0] <- 0
  }
  
  # Retain only relevant columns. Add back month column. 
  tempDF <- tempDF %>%
    dplyr::select(ID, cuID, DATE, wYEAR, WEIGHT, direct) %>%
    dplyr::filter(WEIGHT >= 0) %>%
    dplyr::mutate(MONTH = lubridate::month(DATE))
  #=============================================================================
  
  # Step 2 - Outlier Detection
  #=============================================================================
  # Collect the initial maximums, these will be used to search 
  # for high outliers. 
  # First block gives preference to observations directly obtained from WESD
  # when determining the maximum value within clusters. 
  tempMax <- tempDF %>% 
    dplyr::group_by(cuID, DATE) %>% 
    dplyr::arrange(dplyr::desc(direct), dplyr::desc(WEIGHT)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cuID, wYEAR) %>%
    dplyr::summarize(maxW = max(WEIGHT, na.rm = TRUE),
                     nday = n(),
                     nmonth = length(unique(MONTH)),
                     wmonth = MONTH[which.max(WEIGHT)[1]],
                     ID = ID[which.max(WEIGHT)][1])
  
  # For each station, determine the median and 10th percentile 
  # of the data points. 
  tempQuant <- tempMax %>%
    dplyr::group_by(cuID) %>%
    dplyr::summarize(q50 = median(maxW)) %>%
    dplyr::ungroup()
  
  # Filter all maximums in water years that lack at least 30 records 
  # spanning over 5 different months 
  # (UNLESS the maximum is strictly above the median).
  tempMax <- dplyr::left_join(tempMax, tempQuant, by = "cuID") %>%
    dplyr::filter((nday >= 30 & nmonth >= 5) | maxW > q50)
  
  # Determine the bounds for which we will flag "outliers". 
  # Because we assume the maximum values are right skewed, we take the log 
  # of the data for this calculation. 
  tempBounds <- tempMax %>%
    dplyr::group_by(cuID) %>%
    dplyr::summarize(q25 = stats::quantile(log(maxW + 1), 0.25),
                     q75 = stats::quantile(log(maxW + 1), 0.75),
                     maxl = max(log(maxW + 1)),
                     # 0.7871 is the equivalent to 25 psf after transformation 
                     IQR2 = dplyr::if_else(3*(q75-q25) + q75 > 0.7871, 
                                           3*(q75-q25) + q75, 0.7871)) %>%
    dplyr::mutate(CUTOFF = exp(IQR2) - 1) %>%
    dplyr::select(cuID, CUTOFF, IQR2) %>%
    dplyr::ungroup()
  
  outliers <- dplyr::left_join(tempDF, tempBounds, by = "cuID") %>%
    dplyr::filter(WEIGHT > CUTOFF) 
  
  # If a station has 5 or more observations above the cutoff in a given year
  # remove the outlier flag from these values. 
  outSum <- outliers %>% 
    dplyr::group_by(cuID, wYEAR) %>% 
    dplyr::tally() %>%
    dplyr::filter(n < 5) %>%
    dplyr::ungroup() %>%
    dplyr::select(cuID, wYEAR, n)
  
  # Only retain "outliers" that have less than 5 values above the cutoff. 
  outliers <- dplyr::inner_join(outliers, outSum, by = c("cuID", "wYEAR")) %>%
    dplyr::select(cuID, ID, DATE, WEIGHT, CUTOFF)
  #=============================================================================
  
  return(list(final.load, outliers, tempMax))
}
#=============================================================================
#=============================================================================



