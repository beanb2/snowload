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
    dplyr::filter(MONTH < 6 | MONTH > 9) %>% # only keep relevant months
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
  if(remove.outlier){
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
    
    tempDF <- dplyr::anti_join(tempDF, outliers, by = c("cuID", "ID", "DATE"))
  }
  #=============================================================================
  
  # Step 3 - Distributional estimates. 
  #=============================================================================
  # With the outliers removed, redo the collection of maximums. 
  # (see earlier code in Step 2)
  # First block gives preference to observations directly obtained from WESD
  # when determining the maximum value within clusters. 
  tempMax <- tempDF %>% 
    dplyr::group_by(cuID, DATE) %>%
    dplyr::arrange(dplyr::desc(direct), dplyr::desc(WEIGHT)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cuID, wYEAR) %>% 
    dplyr::summarize(maxW = max(WEIGHT),
                     nday = length(unique(DATE)),
                     nmonth = length(unique(MONTH)),
                     wmonth = which.max(WEIGHT)[1])
  
  tempQuant <- tempMax %>%
    dplyr::group_by(cuID) %>%
    dplyr::summarize(q50 = stats::median(maxW),
                     q10 = stats::quantile(maxW, 0.1),
                     q25 = stats::quantile(maxW, 0.25)) %>%
    dplyr::ungroup()
  
  # Filter all observations at or below the 10th percentile.
  # Unless the 25th percentile is also 0, then keep them. 
  tempMax <- dplyr::left_join(tempMax, tempQuant, by = "cuID") %>%
    dplyr::filter((nday >= 30 & nmonth >= 5) | maxW > q50) # nday >= 30 & 
  
  stScreen <- tempMax %>% dplyr::group_by(cuID) %>% dplyr::tally()
  
  tempMax <- dplyr::left_join(tempMax, stScreen, by = "cuID") %>%
    dplyr::filter(n > 11) %>%
    dplyr::filter(maxW > q10 | q25 == 0)
  
  final.load <- tempMax %>%
    dplyr::group_by(cuID) %>%
    dplyr::summarize(medMax = median(maxW),
                     mMax = max(maxW),
                     DLL = mle50(maxW),
                     DLL2 = mle50(maxW, pcnt = 0.99),
                     DLG = mle50(maxW, dbn = "gamma"),
                     DLG2 = mle50(maxW, pcnt = 0.99, dbn = "gamma"),
                     DLI = mle50(maxW, dbn = "gumbel"),
                     DLI2 = mle50(maxW, pcnt = 0.99, dbn = "gumbel"),
                     DLII = mle50(maxW, dbn = "frechet"),
                     DLII2 = mle50(maxW, pcnt = 0.99, dbn = "frechet"),
                     DLIII = mle50(maxW, dbn = "revweibull"),
                     DLIII2 = mle50(maxW, pcnt = 0.99, dbn = "revweibull"),
                     DLW = mle50(maxW, dbn = "weibull"),
                     DLW2 = mle50(maxW, pcnt = 0.99, dbn = "weibull"),
                     DLN = mle50(maxW, dbn = "gamma"),
                     DLN2 = mle50(maxW, pcnt = 0.99, dbn = "gamma"),
                     #tskew = e1071::skewness(maxW),
                     numMax = n())
  #=============================================================================
  
  # Combine information and return. 
  final.load <- dplyr::inner_join(stations.final, final.load, by = "cuID") %>%
    dplyr::mutate(REGIONII = gsub(REGIONIII, 
                                  pattern = "[[:punct:]][[:digit:]]+$", 
                                  replacement = ""),
                  REGIONI = gsub(REGIONII, 
                                 pattern = "[[:punct:]][[:digit:]]+$", 
                                 replacement = "")) %>%
    dplyr::select(cuID, NAME, STATE, LATITUDE, LONGITUDE, ELEVATION, 
                  REGIONI, REGIONII, REGIONIII, REGION_name, 
                  climate, numST, numMax, medMax, mMax, DLL, DLL2,
                  DLG, DLG2, DLI, DLI2, DLII, DLII2, DLIII, DLIII2,
                  DLW, DLW2, DLN, DLN2)
  
  if(nrow(outliers) > 1){
    outliers$relevant = 0
    outliers$relevant[is.element(outliers$cuID, final.load$cuID)] <- 1
  }
  
  return(list(final.load, outliers, tempMax))
}
#=============================================================================
#=============================================================================


mle50 = function(x, pcnt = 0.98, dbn = "lnorm", method = "mle", fplot = FALSE){
  if(!is.element(dbn, c("lnorm", "gamma", "norm", "weibull",
                        "frechet", "gumbel", "revweibull"))){
    stop("Only viable dbn options are lnorm, norm, weibull, 
         frechet, gumbel, revweibull, and gamma")
  }
  # Ensure no na values in the vector (there shouldn't be any)
  tempvec = as.vector(na.omit(x))
  
  # If a station has NEVER had snow, then simply return 0
  if(max(tempvec) == 0){return(0)}
  
  # Determine the proportion of zero-valued maximums (no snow in a year).
  tempvec.sub <- tempvec[tempvec > 0]
  prop.zero <- 1 - (length(tempvec.sub)/length(tempvec))
  
  # If there are not enough values for the distribution fitting 
  # (at least 5) then simply return the maximum max. 
  if(length(tempvec.sub) < 5){return(max(tempvec))}
  
  # If there are enough non-zero maxiums, commence with the
  # distribution fitting process. 
  
  # First, perform appropriate transformations 
  # (see "Continuous Univariate Distributions Vol 2", 
  #  Johnson, Kotz, and Balakrishnan)
  if(is.element(dbn, c("weibull", "frechet", 
                       "gumbel", "revweibull"))){
    dbn2 = "weibull"
    if(dbn == "frechet"){
      tempvec.sub <- 1/tempvec.sub
    }
    if(dbn == "gumbel"){
      tmean <- mean(tempvec.sub)
      tempvec.sub <- 1/exp(tempvec.sub - tmean)
    }
    if(dbn == "revweibull"){
      tmax <- max(tempvec.sub) + 1
      tempvec.sub <- tmax - tempvec.sub
    }
  }else{
    dbn2 = dbn
  }
  
  tempdist = try(fitdistrplus::fitdist(tempvec.sub, distr = dbn2,
                                       method = method), TRUE)
  
  if(inherits(tempdist, "try-error")){
    print(paste("No viable MLE estimators for", dbn, "distribution. Returning NA...",
                sep = ""))
    return(NA)
  }
  
  # 50 year estimate is the 98th percentile of the above fitted distribution
  tpcnt <- (pcnt - prop.zero)/(1-prop.zero)
  
  if(tpcnt <= 0){
    return(0)
  }
  
  if(is.element(dbn, c("frechet", "gumbel", "revweibull"))){
    tpcnt = 1-tpcnt
  }
  
  # If a station has no snow more often than the desired percentile,
  # return 0. 
  # if not, return the adjusted recurrence value.
  if(dbn2 == "lnorm"){
    recurrencevalue = stats::qlnorm(tpcnt, meanlog = tempdist$estimate[1],
                                    sdlog = tempdist$estimate[2])
  }
  if(dbn2 == "gamma"){
    recurrencevalue = stats::qgamma(tpcnt, shape = tempdist$estimate[1],
                                    rate = tempdist$estimate[2])
  }
  if(dbn2 == "norm"){
    recurrencevalue = stats::qnorm(tpcnt, mean = tempdist$estimate[1],
                                   sd = tempdist$estimate[2])
  }
  if(dbn2 == "weibull"){
    recurrencevalue = stats::qweibull(tpcnt, shape = tempdist$estimate[1],
                                      scale = tempdist$estimate[2])
  }
  
  
  if(fplot){
    plot(tempdist)
    # Title idea:
    # - https://stat.ethz.ch/pipermail/r-help/2008-August/170110.html
    title(paste("Distribution:", dbn, sep = " "), outer=TRUE)
  }
  
  if(dbn == "frechet"){
    recurrencevalue = 1/recurrencevalue
  }
  if(dbn == "gumbel"){
    recurrencevalue = -log(recurrencevalue) + tmean
  }
  if(dbn == "revweibull"){
    recurrencevalue = tmax - recurrencevalue
  }
  
  # I confirmed these results were equivalent to the qzmlnorm results 
  # (for the lognormal distribution only)
  # in the EnvStats package on 10-9-2018.
  return(recurrencevalue)
}


# # Script to bug check the cleaning function.
# fileName =  "../Data/Raw/BCsnow2.csv"
# stations = read.csv("../Data/stationList.csv")
# year1 = 1967
# year2 = 3000
# remove.outlier = TRUE
# method = "Sturm"
# h = 2
# distAdj = 4
# elevAdj = 50
# cluster = TRUE
