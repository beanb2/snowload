
### 2015 Idaho Dataset
#=============================================================================
# Code for Idaho data - reading it into R and some exploratory data analysis
# Read in station data as well as the metadata, merge results.
idaho_ngsl = read.csv("data-raw/Idaho-NGSL.csv")
idaho_ngsl = dplyr::select(idaho_ngsl, STATION, YRS, RESPONSE, NGSL)

idaho_meta = read.csv("data-raw/Idaho-NGSL-metadata.csv")

idaho_final = dplyr::inner_join(idaho_meta,idaho_ngsl, by = "STATION")

idaho_final$LONGITUDE <- round(idaho_final$LONGITUDE, 3)
idaho_final$LATITUDE <- round(idaho_final$LATITUDE, 3)

# Assign stations to appropriate water basins
# Read in water basin data assign to various HUCS
watersheds = rgdal::readOGR(dsn = "data-raw",
                        layer = "naWatershed")

sp::coordinates(idaho_final) = c("LONGITUDE", "LATITUDE")
sp::proj4string(idaho_final) = sp::proj4string(watersheds)

idaho_final$HUC = sp::over(idaho_final, watersheds)$HUC8

# Check for identical station locations
# We will jitter points by .001 degree latitude or longitude
# for the second of every pair of duplicates.
ident_locs <- sp::zerodist(idaho_final, zero = 0.0)[, 2]

id2015 <- as.data.frame(idaho_final)

# For places with identical locations, we add an arbitarily
# small amount to the coordinates to make the locations unique.
set.seed(123)
id2015$LONGITUDE[ident_locs] <- id2015$LONGITUDE[ident_locs] +
  sample(c(-.001, .001), size = length(ident_locs), replace = TRUE)
id2015$LATITUDE[ident_locs] <- id2015$LATITUDE[ident_locs] +
  sample(c(-.001, .001), size = length(ident_locs), replace = TRUE)

# Convert elevation to meters
id2015$ELEVATION <- round(id2015$ELEVATION*.3048)
# Convert load to kpa
id2015$RESPONSE <- round(id2015$RESPONSE*0.04788, 3)
id2015 <- dplyr::select(id2015, STATION, STATION.NAME, STATE,
                        LONGITUDE, LATITUDE, ELEVATION, YRS, RESPONSE, HUC)
colnames(id2015) <- c("STATION", "STATION_NAME", "STATE", "LONGITUDE",
                      "LATITUDE", "ELEVATION", "YRS", "yr50", "HUC")

# usethis::use_data(id2015, overwrite = TRUE)
#=============================================================================

### 1992 Utah Dataset
#=============================================================================
ut1992 <- read.csv("data-raw/oldutahwithlocations.csv")
ut1992$STATE <- "UT"

# Round to three decimal places
ut1992$LONGITUDE <- round(ut1992$LONGITUDE, 3)
ut1992$LATITUDE <- round(ut1992$LATITUDE, 3)

sp::coordinates(ut1992) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(ut1992) <- sp::proj4string(watersheds)
ut1992$HUC <- sp::over(ut1992, watersheds)$HUC8

ident_locs <- sp::zerodist(ut1992, zero = 0.0)[, 2]

ut1992 <- as.data.frame(ut1992)

# For places with identical locations, we add an arbitarily
# small amount to the coordinates to make the locations unique.
ut1992$LONGITUDE[ident_locs] <- ut1992$LONGITUDE[ident_locs] +
  sample(c(-.001, .001), size = length(ident_locs), replace = TRUE)
ut1992$LATITUDE[ident_locs] <- ut1992$LATITUDE[ident_locs] +
  sample(c(-.001, .001), size = length(ident_locs), replace = TRUE)

# Convert elevation to meters
ut1992$ELEVATION <- round(ut1992$ELEVATION*.3048)
# Convert load to kpa
ut1992$yr50 <- round(ut1992$yr50*0.04788, 3)
ut1992$maxobs <- round(ut1992$maxobs*0.04788, 3)
ut1992$approx[is.na(ut1992$approx)] <- 0
ut1992 <- dplyr::select(ut1992, STATION, STATION_NAME, STATE, LONGITUDE,
                        LATITUDE, ELEVATION, approx, maxobs, yr50, HUC)

# usethis::use_data(ut1992, overwrite = TRUE)
#=============================================================================

### 2017 Utah Dataset
#=============================================================================
ut2017 <- read.csv("data-raw/finalUtah.csv")

ut2017$STATE <- stringr::str_extract_all(ut2017$STATION_NAME,
                                     "[[:space:]]+(UT|WY|ID|CO|NM|NV|AZ)[[:space:]]+US")
ut2017$STATE <- gsub(pattern = "[[:space:]]|US", replacement = "",
                     ut2017$STATE)
ut2017$STATION_NAME <- gsub(ut2017$STATION_NAME,
                            pattern = "[[:space:]]+(UT|WY|ID|CO|NM|NV|AZ)[[:space:]]+US",
                            replacement = "")

# Manually fill in blank states
# From the file (snotelextra2.txt)
ut2017$STATE[ut2017$STATE == "character(0)"] <- c("CO", "ID", "WY", "WY", "ID", "ID")

sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(ut2017) <- sp::proj4string(watersheds)
ut2017$HUC <- sp::over(ut2017, watersheds)$HUC8

ut2017 <- as.data.frame(ut2017)
ut2017$LONGITUDE <- round(ut2017$LONGITUDE, 3)
ut2017$LATITUDE <- round(ut2017$LATITUDE, 3)
# Convert elevation to meters
ut2017$ELEVATION <- round(ut2017$ELEVATION*.3048)
# Convert load to kpa
ut2017$yr50 <- round(ut2017$yr50*0.04788, 3)
ut2017$maxobs <- round(ut2017$maxobs*0.04788, 3)
ut2017 <- dplyr::select(ut2017, STATION, STATION_NAME, STATE, LONGITUDE,
                        LATITUDE, ELEVATION, n, maxobs, yr50, HUC)
colnames(ut2017) <- c("STATION", "STATION_NAME", "STATE", "LONGITUDE",
                      "LATITUDE", "ELEVATION", "YRS", "maxobs", "yr50", "HUC")



usethis::use_data(ut2017, overwrite = TRUE)
#=============================================================================

#=============================================================================
`%>%` <- magrittr::`%>%`
postOffice <- read.csv("data-raw/utpostOffices.csv")
postOffice <- postOffice %>%
  dplyr::mutate(CITY = gsub(x = NAME, pattern = "[[:space:]]+Post Office.*$",
                     replacement = ""),
         LON = gsub(x = LON, pattern = "W", replacement = ""),
         LAT = gsub(x = LAT, pattern = "N", replacement = "")) %>%
  dplyr::mutate(CITY = gsub(x = CITY, pattern = "Saint", replacement = "St.")) %>%
  dplyr::mutate(
    LON = as.numeric(LON),
    LAT = as.numeric(LAT),
    LONGITUDE = floor(LON/10000) +
      (floor((LON %% 10000)/100)/60) +
      ((LON %% 100)/3600),
    LATITUDE = floor(LAT/10000) +
      (floor((LAT %% 10000)/100)/60) +
      ((LAT %% 100)/3600),
    LONGITUDE = -LONGITUDE, # West degrees are negative
    ELEVATION = ELEVATION * 0.3048) %>%
  dplyr::select(CITY, COUNTY, STATE, LONGITUDE, LATITUDE, ELEVATION)

sp::coordinates(postOffice) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(postOffice) <- sp::proj4string(watersheds)

postOffice$HUC <- as.numeric(as.character(sp::over(postOffice, watersheds)$HUC8))

utpost <- as.data.frame(postOffice)

usethis::use_data(utpost, overwrite = TRUE)
#=============================================================================

head(id2015)
tail(id2015)
head(ut1992)
tail(ut1992)
head(ut2017)
tail(ut2017)

### TEST LOCATIONS






