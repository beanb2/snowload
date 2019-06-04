#=============================================================================
# Functions used for web scraping of snow records. 
# Author: Brennan Bean
# Created 10-8-2018
# Last Updated: 10-8-2018
# Ensure that working directory is set to source file location prior to 
# running this script. 
# Run before:  - dataCleanFunctions.R
# Run after:   - getMetaData.R
#=============================================================================

#=============================================================================
# Function: Organize data downloaded from the web. 
# This function is designed to take NOAA weather data organized in .dly format,
# extract snow depth and WESD information, and return a neatly organized 
# data frame. Each function call takes an incredible amount of time to run.
# URL - the full URL of the web query. 
# elem - the desired climate variable. 
#=============================================================================
getStationData <- function(URL, elem = c("SNWD", "WESD")){
  # % Define the piping command (avoid loading dplyr specifically)
  # https://stackoverflow.com/questions/27386694/using-operator-from-dplyr-without-loading-dplyr-in-r
  `%>%` <- magrittr::`%>%`
  
  # Test to see if the given URL contains a file
  #vi.url <- RCurl::url.exists(URL)
  
  #if(!vi.url){
  #  warning(paste("File ", URL, " does not exist."))
  #  return(NULL)
  #}
  
  # Read the table from the FTP file. Convert to character vector.
  dly <- try(read.table(URL, sep = "\t"), silent = TRUE)
  if(inherits(dly, "try-error")){
    print(paste("No viable file found at ", URL, sep = ""))
    return(NULL)
  }
  dly <- as.character(dly$V1)
  
  # Create a list that can hold all the entries of dly
  tdata <- vector("list", length = length(dly))
  # Create a separate count variable as we will not keep all variables.
  count = 1
  for(i in 1:length(tdata)){
    tentry <- dly[i]
    
    # Extract the variable name and only proceed if the variable is desired. 
    telem <- substring(tentry, 18, 21)
    
    # Create a data frame for each month of data
    if(is.element(telem, elem)){
      tdf <- data.frame(
        ID = rep(substring(tentry, 1, 11), 31),
        YEAR = rep(substring(tentry, 12, 15), 31),
        MONTH = rep(substring(tentry, 16, 17), 31),
        DAY = 1:31,
        ELEMENT = rep(substring(tentry, 18, 21), 31),
        VALUE = rep("-9999", 31),
        MFLAG = rep(" ", 31),
        QFLAG = rep(" ", 31),
        SFLAG = rep(" ", 31),
        # Avoid creating factor varibles from the character strings
        stringsAsFactors = FALSE 
      )
      # Iteration of k taken from documentation:
      # https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
      k = 22
      for(j in 1:31){
        tdf$VALUE[j] <- substring(tentry, k, k+4)
        tdf$MFLAG[j] <- substring(tentry, k+5, k+5)
        tdf$QFLAG[j] <- substring(tentry, k+6, k+6)
        tdf$SFLAG[j] <- substring(tentry, k+7, k+7)
        k = k+8
      }
      tdata[[count]] <- tdf
      count = count + 1
    }
  }
  
  tdata2 <- data.table::rbindlist(tdata) 
  
  if(nrow(tdata2) > 0){
    tdata2 <- tdata2 %>% 
      dplyr::mutate(YEAR = as.numeric(as.character(YEAR)),
                    MONTH = as.numeric(as.character(MONTH)),
                    DAY = as.numeric(as.character(DAY)),
                    VALUE = as.numeric(as.character(VALUE))) %>%
      dplyr::filter(VALUE != -9999)
    
    if(nrow(tdata2) > 0){
      return(tdata2)
    }else{
      return(NULL)
    }
  }else{
    return(NULL)
  }
}
#=============================================================================


#=============================================================================
# Function: download snow data for a particular state. 
# Uses the getStationData function in a loop to organize results by state. 
# state - the abbreviated state name as provided in the station meta-data. 
# stations - the station meta data run previously. 
# mainpath - The main URL to download station data from the web. 
# ext - the desired file extension. 
#=============================================================================
getState <- function(state, stations, mainpath, ext = ".dly", ...){
  # Define the pipeline operator:
  # - https://stackoverflow.com/questions/27386694/using-operator-from-dplyr-without-loading-dplyr-in-r
  `%>%` <- magrittr::`%>%`
  
  # Subset the stations according to the given state. Blank will give return
  # all stations without a state designation. 
  stations.sub <- stations[stations$STATE == state, ]
  st.ID <- stations.sub$ID
  data.list <- vector("list", length(st.ID))
  for(i in 1:length(st.ID)){
    data.list[[i]] <- getStationData(URL = paste(mainpath, st.ID[i], ext, sep = ""), ...)
  }
  
  tdata2 <- data.table::rbindlist(data.list) 
  
  if(nrow(tdata2) > 0){
    tdata2 <- tdata2 %>% 
      dplyr::mutate(YEAR = as.numeric(as.character(YEAR)),
                    MONTH = as.numeric(as.character(MONTH)),
                    DAY = as.numeric(as.character(DAY)),
                    VALUE = as.numeric(as.character(VALUE)))
    
    if(nrow(tdata2) > 0){
      return(tdata2)
    }else{
      return(NULL)
    }
  }else{
    return(NULL)
  }
}
#=============================================================================

