# A script to download station data for NWS first order stations.
# Taken from: https://www.ncdc.noaa.gov/homr/
library(tidyverse)
spots <- c(1, 7, 18, 39, 50, 54, 60, 69, 76, 88, 139, 200,
           226, 229, 232, 243, 253, 264, 271, 278, 280)
spotsdiff <- diff(spots)

# Read in the first order station list.
stationList <- read.fwf("data-raw/lcd-stations.txt",
                         widths = spotsdiff)

sum(is.element(as.character(stationList$V13), c("AZ ", "CA ", "CO ", "ID ",
                                "MT ", "NM ", "NV ", "OR ", "UT ",
                                "WA ", "WY ")))

sum(is.element(as.character(stationList$V13), c("UT ")))

# Read in eco regions.

stations <- data.frame(ID = substring(stationList, 1, 12),
                       LATITUDE = substring(stationList, 13, 21),
                       LONGITUDE = substring(stationList, 22, 30),
                       ELEVATION = substring(stationList, 31, 38),
                       STATE = substring(stationList, 39, 41),
                       NAME = substring(stationList, 42, 72),
                       NETWORK1 = substring(stationList, 73, 76),
                       NETWORK2 = substring(stationList, 77, 80),
                       ALTID = substring(stationList, 81),
                       stringsAsFactors = FALSE) %>%
  dplyr::mutate(ID = gsub(ID, pattern = "[[:space:]]", replacement = ""),
                LATITUDE = gsub(LATITUDE, pattern = "[[:space:]]", replacement = ""),
                LONGITUDE = gsub(LONGITUDE, pattern = "[[:space:]]", replacement = ""),
                ELEVATION = gsub(ELEVATION, pattern = "[[:space:]]", replacement = ""),
                STATE = gsub(STATE, pattern = "[[:space:]]", replacement = ""),
                NAME = gsub(NAME, pattern = "[[:space:]]+$", replacement = ""),
                NAME = gsub(NAME, pattern = "^[[:space:]]+", replacement = ""),
                NETWORK1 = gsub(NETWORK1, pattern = "[[:space:]]", replacement = ""),
                NETWORK2 = gsub(NETWORK2, pattern = "[[:space:]]", replacement = ""),
                ALTID = gsub(ALTID, pattern = "[[:space:]]", replacement = "")) %>%
  dplyr::mutate(LATITUDE = as.numeric(LATITUDE),
                LONGITUDE = as.numeric(LONGITUDE),
                ELEVATION = as.numeric(ELEVATION)) %>%
  dplyr::filter(STATE != "") # Only want stations in U.S. and Canada

