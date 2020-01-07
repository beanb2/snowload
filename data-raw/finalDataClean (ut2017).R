library(plyr); library(dplyr) # Load Plyr first to not cause problems
library(data.table)
library(fitdistrplus)
#library(evd)
#library(EnvStats)
library(RColorBrewer)
library(rgdal)
library(sp)
library(fields)
require(tidyr)
require(moments)
#library(extRemes)

# All file paths in the document are relative to the location of where the code is located on the thumb drive. 
# Read in elevation DEM (see "unzipandstitch.R" code in data folder for how this was done)
load("../../R Objects/elevationDEM2.R")
### READ IN SHAPEFILES ###
watersheds = readOGR(dsn = "../../Shapefiles and DEMs/ShapeFiles/extractHUC/hydrologic_units",
                     layer = "wbdhu12_a_extract")
watersheds <- spTransform(watersheds, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
                                          +datum=WGS84 +no_defs +towgs84=0,0,0"))

# Read In Utah Counties Shapefiles as Obtained from gis.utah.gov
map <- readOGR(dsn="../../Shapefiles and DEMs/ShapeFiles/Counties_shp/Counties", layer="Counties")
# Reproject to the same coordinate extent as our DEM
map.2 <- spTransform(map, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
                              +datum=WGS84 +no_defs +towgs84=0,0,0"))

# Read In Utah STATE Shapefile as Obtained from gis.utah.gov
utah <- readOGR(dsn="../../Shapefiles and DEMs/ShapeFiles/Utah_shp/Utah", layer="Utah")
# Reproject to the same coordinate extent as our DEM
utah <- spTransform(utah, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
                              +datum=WGS84 +no_defs +towgs84=0,0,0"))



### COOP STATION DATA - ORIGINAL ###
### File Read ###
tfiles = list.files("../../Data/Snow Project Stations/CoopNew/", pattern = ".csv")
tfiles2 = list.files("../../Data/Snow Project Stations/CoopNew2", pattern = ".csv")
tfiles3 = list.files("../../Data/Snow Project Stations/COOPSimon/StationData", pattern = ".csv")

coopl.original <- lapply(paste("../../Data/Snow Project Stations/CoopNew/",
                               tfiles, sep = ""), fread, 
                         colClasses = rep("Character", 16))

coopl.new <- lapply(paste("../../Data/Snow Project Stations/CoopNew2/",
                          tfiles2, sep = ""), fread)

coopl.new2 <- lapply(paste("../../Data/Snow Project Stations/COOPSimon/StationData/",
                          tfiles3, sep = ""), fread)

coopl = rbindlist(coopl.original, fill = TRUE)
colnames(coopl) <- c("STATION","STATION_NAME","ELEVATION","LATITUDE",	
                     "LONGITUDE","DATE","SNWD","Measurement.Flag",
                     "Quality.Flag","Source.Flag",
                     "Time.of.Observation","WESD","Measurement.Flag.1",
                     "Quality.Flag.1","Source.Flag.1", "Time.of.Observation.1")

# All values are read in as character values. 
coopl2 = rbindlist(coopl.new, fill = TRUE)
coopl3 = rbindlist(coopl.new2, fill = TRUE)

## coopl (original) formatting steps:
# Includes the following:
# 1 - remove "GHCND" tag on station numbers
# 2 - convert from inches to centimeters
coopl <- coopl %>%
  mutate(STATION = gsub(x = STATION, pattern = "GHCND:", replacement = "")) %>%
  mutate(SNWD = as.numeric(SNWD),
         WESD = as.numeric(WESD),
         ELEVATION = as.numeric(ELEVATION),
         LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE))


## coopl2 (new) formatting steps:
# 1 - separate SNWD_ATTRIBUTES and WESD_ATTRIBUTE columns 
#     to match format of original
# 2 - reorder columns to match order of original
# 3 - remove "-" symbols in date variable and convert to numeric
# 4 - replace NA values with the respective missing value designations 
#     from the original file
# 5 - convert from mm to cm
coopl2 <- coopl2 %>% 
  separate(col = SNWD_ATTRIBUTES, 
           into = c("Measurement.Flag", "Quality.Flag", 
                    "Source.Flag", "Time.of.Observation"), 
           fill = "right",
           sep = ",") %>% 
  separate(col = WESD_ATTRIBUTES,
           into = c("Measurement.Flag.1", "Quality.Flag.1", 
                    "Source.Flag.1", "Time.of.Observation.1"), 
           fill = "right",
           sep = ",") %>%
  dplyr::select(STATION, NAME, ELEVATION, LATITUDE, LONGITUDE, DATE, SNWD, 
                Measurement.Flag, Quality.Flag, Source.Flag, Time.of.Observation, 
                WESD, Measurement.Flag.1, Quality.Flag.1, Source.Flag.1, 
                Time.of.Observation.1) %>% 
  mutate(DATE = gsub(x = DATE, pattern = "-", replacement = "")) %>%
  replace_na(list(Measurement.Flag = "",
                  Quality.Flag = "",
                  Source.Flag = "", 
                  Time.of.Observation = "9999",
                  Measurement.Flag.1 = "",
                  Quality.Flag.1 = "",
                  Source.Flag.1 = "",
                  Time.of.Observation.1 = "9999")) %>% 
  mutate(SNWD = as.numeric(SNWD),
         WESD = as.numeric(WESD))

# Do again for coopl3
coopl3 <- coopl3 %>% 
  separate(col = SNWD_ATTRIBUTES, 
           into = c("Measurement.Flag", "Quality.Flag", 
                    "Source.Flag", "Time.of.Observation"), 
           fill = "right",
           sep = ",") %>% 
  separate(col = WESD_ATTRIBUTES,
           into = c("Measurement.Flag.1", "Quality.Flag.1", 
                    "Source.Flag.1", "Time.of.Observation.1"), 
           fill = "right",
           sep = ",") %>%
  dplyr::select(STATION, NAME, ELEVATION, LATITUDE, LONGITUDE, DATE, SNWD, 
                Measurement.Flag, Quality.Flag, Source.Flag, Time.of.Observation, 
                WESD, Measurement.Flag.1, Quality.Flag.1, Source.Flag.1, 
                Time.of.Observation.1) %>% 
  mutate(DATE = gsub(x = DATE, pattern = "-", replacement = "")) %>%
  replace_na(list(Measurement.Flag = "",
                  Quality.Flag = "",
                  Source.Flag = "", 
                  Time.of.Observation = "9999",
                  Measurement.Flag.1 = "",
                  Quality.Flag.1 = "",
                  Source.Flag.1 = "",
                  Time.of.Observation.1 = "9999")) %>% 
  mutate(SNWD = as.numeric(SNWD)*0.0393701, # Convert from mm to inches
         WESD = as.numeric(WESD)*0.0393701)

# Make sure the column names match
colnames(coopl2) <- colnames(coopl)
colnames(coopl3) <- colnames(coopl)


coopl.combined <- bind_rows(coopl, coopl2, coopl3)

### NOW ADD THE SNOTEL DATA
sntl = read.csv("../../Data/Snow Project Stations/SNOTEL/snotelextra2.txt", 
                header = TRUE, comment.char = "#")
colnames(sntl) = c("DATE", "STATION", "STATION_NAME", "ELEVATION", "LATITUDE",
                   "LONGITUDE", "WESD", "SNWD", "HUC")

sntl <- sntl %>% mutate(DATE = gsub(x=DATE, pattern = "-", "", DATE))

# Create and mutate a few columns to get forms to match
sntl <- sntl %>% 
  mutate(ELEVATION = ELEVATION*0.3048,
         Measurement.Flag = "",
         Measurement.Flag.1 = "",
         Quality.Flag = "",
         Quality.Flag.1 = "",
         Source.Flag = "",
         Source.Flag.1 = "",
         Time.of.Observation = 9999,
         Time.of.Observation.1 = 9999,
         STATION = as.character(STATION),
         ELEVATION = as.numeric(as.character(ELEVATION)),
         LATITUDE = as.numeric(as.character(LATITUDE)),
         LONGITUDE = as.numeric(as.character(LONGITUDE)),
         Time.of.Observation = as.character(Time.of.Observation),
         Time.of.Observation.1 = as.character(Time.of.Observation.1),
         WESD = WESD*0.0393701, # convert from mm to inches
         SNWD = SNWD*0.393701) # convert from cm to inches

sntl <- sntl %>% dplyr::select(STATION, STATION_NAME, ELEVATION, LATITUDE, 
                               LONGITUDE, DATE, SNWD, Measurement.Flag, 
                               Quality.Flag, Source.Flag, Time.of.Observation,
                               WESD, Measurement.Flag.1, Quality.Flag.1,
                               Source.Flag.1, Time.of.Observation.1) %>%
  filter(STATION_NAME != "Mancos") %>% 
  filter(STATION_NAME != "Sharkstooth") # Filter out duplicated stations (now that we are using expanded dataset)

coopl.combined <- bind_rows(coopl.combined, sntl)

# Convert to tbl
coopl = tbl_df(coopl.combined) # 9629175

# Remove other verions
remove(coopl.new, coopl.original, coopl.combined, coopl2) 
remove(coopl3, sntl, coopl.new2)
#coopl.test = coopl

# sum(is.na(coopl$STATION)) # No station names are missing (good)


# We don't care what time of day observations were made, eliminate these variables
coopl = dplyr::select(coopl, -Time.of.Observation, -Time.of.Observation.1)

# Convert elevation, latitude, and longitude to numeric 
coopl = mutate(coopl, ELEVATION = as.numeric(ELEVATION),
               LATITUDE = as.numeric(LATITUDE),
               LONGITUDE = as.numeric(LONGITUDE))

# Replace WESD and SNWD NA values with -9999
coopl = coopl %>% replace_na(list(SNWD = -9999, WESD = -9999)) %>%
  filter(SNWD >= 0 | WESD >= 0) %>%
  mutate(STATION_NAME = toupper(STATION_NAME),
         STATION_NAME = gsub("[[:punct:]]", "", STATION_NAME)) # 7045305

coopl = coopl %>% 
  mutate(DATE = as.numeric(DATE),
         YEAR = floor(DATE/10000), 
         MONTH = floor(DATE / 100) %% 100) %>%
  filter(MONTH < 6 | MONTH > 10)  %>%
  filter(YEAR >= 1969) # 4138747 observations

# Code to determine unique station numbers
stnum <- dplyr::select(coopl, STATION_NAME) %>% unique(.) # 1350 unique station numbers

### Station Metdata Preparation

# Extract unique station meta data
tempstations <- dplyr::select(coopl, STATION_NAME, STATION, LATITUDE, LONGITUDE, ELEVATION) %>% 
  unique(.)


# Previously, we looked at station names and separated stations that had differences of more 
# than 100m or 10km. It seems like in most cases these are simply lapses in the meta data. 

# Rather we will separate by station NUMBER and we will take the median of each meta data value, 
# no longer choosing to split the data into two groups like we did previously. 
# Based on what we see, we trust that fishiness will be resolved or filtered out at a later point. 

# We are operating under the assumption that major deviances in location/elevation for the 
# same station number are errors in reporting, not true differences in the location of 
# measurements. 

# For each of these stations, average the latitude, longitudes and elevations
# so that there is only ONE station location and ONE elevation for each location of interest
combined_stations = tempstations %>% dplyr::select(-STATION_NAME) %>% group_by(STATION) %>% 
  na.omit() %>% summarise_all(funs(median))

# For differing station numbers at this point, assign it to the first station in terms of alphanumeric ordering
station_numbers = tempstations %>% dplyr::select(STATION_NAME, STATION) %>% group_by(STATION) %>% 
  slice(1)

combined_stations.2 = inner_join(combined_stations, station_numbers, by = "STATION")

# Look at stations "left out" of averaging
lostations = dplyr::anti_join(tempstations, combined_stations.2, by = "STATION")

# Two stations need meta data information which we obtain from
# - https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt
lostations[1,3:5] = c(39,-110.5,1239.9) # Green River, UT   
lostations[2,3:5] = c(40.1500, -110.0667, 1623.1) # Myton UT
lostations[3,3:5] = c(36.15,-109.53,1690.1) # Canyon de Chelly, AZ

lostations = dplyr::select(lostations, STATION_NAME, LATITUDE, LONGITUDE, ELEVATION, STATION)

finalStations = bind_rows(combined_stations.2, lostations)

# We now deal with the fact that there are still more station numbers than names. 
# We want to see how signficantly these stations differ. 

elevDiff = finalStations %>% na.omit(.) %>% group_by(STATION_NAME) %>%
  summarise(max(abs(outer(ELEVATION, ELEVATION, "-")))) %>% dplyr::arrange(STATION_NAME)
maxdist = finalStations %>% na.omit(.) %>% group_by(STATION_NAME) %>%
  summarise(max(rdist.earth(matrix(c(LONGITUDE,LATITUDE),ncol = 2),
                            matrix(c(LONGITUDE,LATITUDE),ncol = 2), miles = FALSE))) %>%
  dplyr::arrange(STATION_NAME)
differences = left_join(elevDiff, maxdist)
colnames(differences) = c("STATION_NAME", "elevdiff", "distdiff")
goodStations = filter(differences, elevdiff < 100 & distdiff < 5)
goodStations.2 = finalStations %>% filter(is.element(STATION_NAME, goodStations$STATION_NAME)) %>% dplyr::arrange(STATION_NAME)
badStations = filter(differences, !(elevdiff < 100 & distdiff < 5))
badStations.2 = finalStations %>% filter(is.element(STATION_NAME, badStations$STATION_NAME)) %>% dplyr::arrange(STATION_NAME)

# Rename the three stations that should be different.
# Both in the final station list and in the original files 
coopl$STATION_NAME[coopl$STATION == "USC00051020"] =  "BROWNS PARK REFUGE CO US - 2"
coopl$STATION_NAME[coopl$STATION == "USS0012M13S"] =  "CASTLE VALLEY UT US - 2"
coopl$STATION_NAME[coopl$STATION == "USS0011J42S"] =  "SNOWBIRD UT US - 2"

finalStations$STATION_NAME[finalStations$STATION == "USC00051020"] =  "BROWNS PARK REFUGE CO US - 2"
finalStations$STATION_NAME[finalStations$STATION == "USS0012M13S"] =  "CASTLE VALLEY UT US - 2"
finalStations$STATION_NAME[finalStations$STATION == "USS0011J42S"] =  "SNOWBIRD UT US - 2"


# Now, take the final station list, group by STATION_NAME instead of number, and take the median
# meta data values. We will then group everything following on station name instead of station 
# number. 

# We MAY run into a place where two equations will be measuring the same location at the same 
# time. In such caes, we will simply take the larger of the two measurements in the case of an 
# overlap. 

finalStations <- finalStations %>% arrange(STATION) %>% group_by(STATION_NAME) %>% slice(1)


# We now have both unique station names and station numbers. Lastly, we need to combine any
# station information that shares nearly identical location. 

pw.dist <- rdist.earth(matrix(c(finalStations$LONGITUDE,finalStations$LATITUDE),ncol = 2),
            matrix(c(finalStations$LONGITUDE,finalStations$LATITUDE),ncol = 2), 
            miles = FALSE)
pw.elev <- abs(outer(finalStations$ELEVATION, finalStations$ELEVATION, "-"))

# Add arbitrary number to diagonal to make sure this is not counted as
# the minimum distance. 
pw.dist <- pw.dist + diag(nrow(pw.dist))*1000
pw.elev <- pw.elev + diag(nrow(pw.elev))*1000

sim.ind <- which(pw.dist < .01, arr.ind=TRUE)
sim.vec <- which(pw.dist < .01)

station.zero = data.frame(station1 = finalStations$STATION_NAME[sim.ind[,1]],
                    station2 = finalStations$STATION_NAME[sim.ind[,2]],
                    distance = pw.dist[sim.vec],
                    elevation = pw.elev[sim.vec])

station.zero <- station.zero %>% group_by(elevation) %>% slice(1) # 1286 stations

# Two sets of locations have distances between them of 0 with elevation differences less than 50 meters. 
# We are going to combine each set into one station. 
coopl$STATION_NAME[coopl$STATION_NAME == "PARACHUTE CO US"] <- "GRAND VALLEY CO US - combined"
coopl$STATION_NAME[coopl$STATION_NAME == "GRAND VALLEY CO US"] <- "GRAND VALLEY CO US - combined"

coopl$STATION_NAME[coopl$STATION_NAME == "SNOW BASIN UT US"] <- "SNOW BASIN UT US - combined"
coopl$STATION_NAME[coopl$STATION_NAME == "HUNTSVILLE SNOW BSN UT US"] <- "SNOW BASIN UT US - combined"

finalStations$STATION_NAME[finalStations$STATION_NAME == "PARACHUTE CO US"] <- "GRAND VALLEY CO US - combined"
finalStations$STATION_NAME[finalStations$STATION_NAME == "GRAND VALLEY CO US"] <- "GRAND VALLEY CO US - combined"

finalStations$STATION_NAME[finalStations$STATION_NAME == "SNOW BASIN UT US"] <- "SNOW BASIN UT US - combined"
finalStations$STATION_NAME[finalStations$STATION_NAME == "HUNTSVILLE SNOW BSN UT US"] <- "SNOW BASIN UT US - combined"

finalStations <- finalStations %>% ungroup() %>% arrange(STATION) %>% group_by(STATION_NAME) %>% slice(1) # 1284 stations

### DATA FILTERS ###

# We want to exclude stations that are more than 100km outside of the state of Utah.
finalStations.sp= finalStations
coordinates(finalStations.sp) = c("LONGITUDE", "LATITUDE")
proj4string(finalStations.sp) = "+init=epsg:4326 +proj=longlat +ellps=WGS84
                              +datum=WGS84 +no_defs +towgs84=0,0,0"
finalStations.sp$COUNTYNBR = over(finalStations.sp, map.2)$COUNTYNBR

toborder = round(rdist.earth(data.frame(Longitude = finalStations.sp$LONGITUDE, Latitude = finalStations.sp$LATITUDE),
                             utah@polygons[[2]]@Polygons[[1]]@coords, miles = FALSE))
toborder = apply(toborder, 1, min)
finalStations.sp$toborder = toborder <= 100

# Get list of stations within 100km of Utah OR IN UTAH
relevantstations = finalStations[!is.na(finalStations.sp$COUNTYNBR) | finalStations.sp$toborder, ] # 966 stations

# check station number and names at this point
coopl = filter(coopl, is.element(STATION_NAME, relevantstations$STATION_NAME)) # 3170837 observations

relevantstations2 = coopl %>% dplyr::select(STATION_NAME, YEAR) %>% group_by(STATION_NAME) %>% unique(.) %>%
  count(.) %>% filter(n >= 10) # 506 stations

relevantstations_final = inner_join(relevantstations, relevantstations2, by = 'STATION_NAME')
relevantstations_final = dplyr::arrange(relevantstations_final, n) # 473 stations for now 

# Separate stations based on elevation (the value given represents the third quantile of station elevations)
stations.low  = relevantstations_final %>% dplyr::select(STATION_NAME, ELEVATION) %>% filter(ELEVATION < 2113.6)
stations.high  = relevantstations_final %>% dplyr::select(STATION_NAME, ELEVATION) %>% filter(ELEVATION >= 2113.6)

coopl = filter(coopl, is.element(STATION_NAME, relevantstations_final$STATION_NAME)) # 2992976

# NA's are all in lat/lon and elevation (quality control variables do not count)
# Extract Date information that we can use when calculating bulk snow density in Sturm's equation
# Use this to recode days according to snow season
# Day of year adjusted has November 1st = Day -62 and May 1st = 120 and we exclude leap years
temp = as.POSIXlt(as.character(coopl$DATE), format = "%Y%m%d")
coopl = mutate(coopl, DOY = temp$yday+1, DOYA = DOY)
coopl$DOYA[coopl$DOYA >= 305] =  coopl$DOYA[coopl$DOYA >= 305] - 366
coopl$DOYA[coopl$DOYA == 0] = -1 # No leap years in calculation - set extra day to december 31st

# Adjust years to reflect water year (starting in October, ending in June) 
# rather than calendar year. We do this by advancing the year of negative 
# values of the snow season by 1. 
coopl$YEAR[coopl$DOYA < 0] <- coopl$YEAR[coopl$DOYA < 0] + 1

# Remove all readings from the unfinished water year of 2018, or from before prior to November 1, 1969.
coopl <- coopl %>% filter(YEAR < 2018, YEAR > 1969) #2956525

# Check for moments when both WESD and SNWD are defined but disagree
sum(coopl$SNWD == 0 & coopl$WESD > 0)
sum(coopl$SNWD > 0 & coopl$WESD == 0)

coopl.1 = filter(coopl, WESD > 0 | (WESD == 0 & SNWD <= 0)) # 883382
coopl.2 = filter(coopl, !(WESD > 0 | (WESD == 0 & SNWD <= 0))) # 2073143 


# Quality Control for coop2.1
# Define quality control flags we wish to look out for
flags = c("D", "G", "I", "K", "L", "M", "N", "O", "R", "S", "T", "W", "X", "Z")

# Remove observations flagged for quality control
coopl.1 = filter(coopl.1, !is.element(Quality.Flag.1, flags))  # 883381 (removed 1)
coopl.2 = filter(coopl.2, !is.element(Quality.Flag, flags)) # 2066769 (removed 6374)


# Remove "Missing Values Presumed 0"
# We dont need to remove "P" flags in coopl.1 as there aren't any, this is the reason the variable gets read in 
# as a logical (T - should stand for "trace" and not "true" but either way we know what it means)
coopl.2 = filter(coopl.2, Measurement.Flag != "P") # 1271236

# Use the maximum density of water (english units) in pressure calculations
coopl.1 = mutate(coopl.1, PRESSURE = WESD * (62.426/12))

# When SWE is not defined, use Sturm's equation OR the RCMD to estimate
coopl.2.low = coopl.2 %>%  filter(is.element(STATION_NAME, stations.low$STATION_NAME)) %>%
  mutate(p_d =  .3608*(1 - exp(-.0016*(SNWD*2.54) - .0031*DOYA)) + .2332)
coopl.2.high = coopl.2 %>%  filter(is.element(STATION_NAME, stations.high$STATION_NAME)) %>%
  mutate(p_d =  .3738*(1 - exp(-.0012*(SNWD*2.54) - .0038*DOYA)) + .2237)

coopl.2 = bind_rows(coopl.2.low, coopl.2.high)

coopl.2 <- coopl.2 %>% 
  mutate(PRESSURE = SNWD*p_d*(62.426/12),
         PRESSURE2 = (2.36*SNWD - 31.9)) %>% 
  dplyr::select(-p_d)

coopl.2$PRESSURE2[coopl.2$SNWD < 22] = .9 * coopl.2$SNWD[coopl.2$SNWD < 22]

coopl = bind_rows(coopl.1, coopl.2) # 2154617

# Determine number of unique stations left at this point. 
# check <- coopl %>% dplyr::select(STATION_NAME) %>% unique(.)

# Rename the coopl.1 to use as a copy for the next few data manipulations
coopl.1 = coopl
coopl.1 = coopl.1 %>% dplyr::select(STATION_NAME, YEAR, DOY, PRESSURE) %>% dplyr::arrange(STATION_NAME, YEAR, DOY)
# Determine pressure changes of more than 30 psf over consecutive measurements.
diffs1 = abs(c(0,diff(coopl.1$PRESSURE)))
diffs2 = abs(c(diff(coopl.1$PRESSURE),0))
# Measurements must be within 10 days of each other to be flagged as a suspcious jump in pressure
diffs.d1 = abs(c(0,diff(coopl.1$DOY)))
diffs.d2 = abs(c(diff(coopl.1$DOY),0))
# Measurements needs to be from the same year to be compared
diffs.y1 = abs(c(0,diff(coopl.1$YEAR)))
diffs.y2 = abs(c(diff(coopl.1$YEAR),0))
# Measurements also need to be from the same station
samestation = duplicated(coopl.1$STATION_NAME, fromLast = TRUE)

tester = data.frame(diffs1 > 30,diffs2 > 30,diffs.d1 < 10,diffs.d2 < 10,diffs.y1 < 1,diffs.y2 < 1, samestation)

flagged = apply(tester,1,prod)

# Suspect results
suspect = c(which(as.logical(flagged)),which(as.logical(flagged))-1,which(as.logical(flagged))-2, 
            which(as.logical(flagged))+1, which(as.logical(flagged))+2)
suspect = suspect[order(suspect)]

suspect = coopl.1[suspect,]

# I am not sure why but none of these readings look abnormal except for the ones in indian canyon. 
# UPDATE 4-22-2017 - I needed to look at the ABSOLUTE value of differences and had not done so until now. This has been updated. 
# We still remove these readings at Indian Canyon, but also remove some additional points as given below. 

coopl = coopl %>% filter(!(STATION_NAME == "INDIAN CANYON UT US" & YEAR == 1979 & is.element(DOY, 12:20))) %>%
  filter(!(STATION_NAME == "BEAVER DAMS UT US" & YEAR == 1991 & DOY == 330)) %>%
  filter(!(STATION_NAME == "MACK 5 NW CO US" & YEAR == 2007 & DOY == 340)) %>%
  filter(!(STATION_NAME == "PRESTON ID US" & YEAR == 1982 & DOY == 21)) %>% 
  filter(!(STATION_NAME == "SEELEY CREEK UT US" & YEAR == 1979 & DOY == 30)) %>%
  filter(!(STATION_NAME == "SILVER LAKE BRIGHTON UT US" & YEAR == 1982 & DOY == 89)) %>%
  filter(!(STATION_NAME == "ELY YELLAND FIELD AIRPORT NV US" & YEAR == 1990 & DOY == 48)) %>%
  filter(!(STATION_NAME == "COTTONWOOD HEIGHTS 15 SE UT US" & YEAR == 2010 & DOY == 18)) %>% 
  filter(!(STATION_NAME == "COTTONWOOD HEIGHTS 15 SE UT US" & YEAR == 2010 & DOY == 119)) %>%
  filter(!(STATION_NAME == "SILVER LAKE BRIGHTON UT US" & YEAR == 1982 & DOY == 89)) %>%
  filter(!(STATION_NAME == "LOGAN 17 ESE UT US" & YEAR == 2009 & DOY == 55)) %>%
  filter(!(STATION_NAME == "LOGAN 17 ESE UT US" & YEAR == 2009 & DOY == 69)) %>%
  filter(!(STATION_NAME == "LOGAN 17 ESE UT US" & YEAR == 2010 & DOY == 357)) # 2154601

# Make NA values for pressure2 (THE RCMD estimate) equal to -9999
coopl$PRESSURE2[is.na(coopl$PRESSURE2)] = -9999 
coopl$PRESSURENEW = coopl$PRESSURE2
coopl$PRESSURENEW[coopl$PRESSURE2 < 0] = coopl$PRESSURE[coopl$PRESSURE2 < 0]

# Use a filtered set of variables for the remainder of the analysis
coopl.2 = dplyr::select(coopl, STATION, STATION_NAME, LATITUDE, LONGITUDE, ELEVATION, YEAR, MONTH, DOY, PRESSURE, PRESSURENEW)

coopl.low = coopl.2 %>% filter(is.element(STATION_NAME, stations.low$STATION_NAME))
# Adjust the high elevation stations to reflect the "alpine" parameters. 
coopl.high = coopl.2 %>% filter(is.element(STATION_NAME, stations.high$STATION_NAME)) 

maxes.low = coopl.low %>% group_by(STATION_NAME, YEAR) %>% summarize(maxp = max(PRESSURE)) # 10672
maxes.low2 = coopl.low %>% group_by(STATION_NAME, YEAR) %>% summarize(maxp = max(PRESSURENEW))
coverage.low = coopl.low %>% 
  group_by(STATION_NAME, YEAR) %>% 
  summarize(COUNT = sum(is.element(unique(MONTH), c(12,1,2,3))))

maxes.high = coopl.high %>% group_by(STATION_NAME, YEAR) %>% summarize(maxp = max(PRESSURE)) # 4668
maxes.high2 = coopl.high %>% group_by(STATION_NAME, YEAR) %>% summarize(maxp = max(PRESSURENEW))
coverage.high = coopl.high %>% 
 group_by(STATION_NAME, YEAR) %>% 
  summarize(COUNT = sum(is.element(unique(MONTH), c(2,3,4,5))))

# Determine the median snow load at each station, 
# keep all maximums above the median regardless of coverage
medians.low = maxes.low %>% ungroup() %>% 
  group_by(STATION_NAME) %>% summarize(meds = median(maxp))
medians.low2 = maxes.low2 %>% ungroup() %>% 
  group_by(STATION_NAME) %>% summarize(meds = median(maxp))

medians.high = maxes.high %>% ungroup() %>% 
  group_by(STATION_NAME) %>% summarize(meds = median(maxp))
medians.high2 = maxes.high2 %>% ungroup() %>% 
  group_by(STATION_NAME) %>% summarize(meds = median(maxp))

maxes.low <- left_join(maxes.low, coverage.low, 
                       by = c("STATION_NAME", "YEAR")) %>% 
  left_join(., medians.low, 
            by = "STATION_NAME") %>% 
  replace_na(list(COUNT = 0)) %>% filter(COUNT > 3 | maxp > meds) %>%
  dplyr::select(-COUNT) # 5979

maxes.high <- left_join(maxes.high, coverage.high, 
                        by = c("STATION_NAME", "YEAR")) %>% 
  left_join(., medians.high, 
            by = "STATION_NAME") %>%
  replace_na(list(COUNT = 0)) %>% filter(COUNT > 3 | maxp > meds) %>%
  dplyr::select(-COUNT) # 3981

maxes.low2 <- left_join(maxes.low2, coverage.low, 
                        by = c("STATION_NAME", "YEAR")) %>% 
  left_join(., medians.low2, 
            by = "STATION_NAME") %>%
  replace_na(list(COUNT = 0)) %>% filter(COUNT > 3 | maxp > meds) %>%
  dplyr::select(-COUNT) 

maxes.high2 <- left_join(maxes.high2, coverage.high, 
                         by = c("STATION_NAME", "YEAR")) %>% 
  left_join(., medians.high2, 
            by = "STATION_NAME") %>%
  replace_na(list(COUNT = 0)) %>% filter(COUNT > 3 | maxp > meds) %>%
  dplyr::select(-COUNT)

maxes = rbind(maxes.low, maxes.high)
maxes2 = rbind(maxes.low2, maxes.high2)

# ENSURE THAT THERE ARE AT LEAST 12 MAXIMUMS PER STATION (REQUIRED)
numMax = maxes %>% group_by(STATION_NAME) %>% count(.) %>% filter(n >= 12)

maxes.2 = maxes %>% filter(is.element(STATION_NAME, numMax$STATION_NAME))
maxes2.2 = maxes2 %>% filter(is.element(STATION_NAME, numMax$STATION_NAME)) 

maxes.med <- maxes.2 %>% ungroup() %>% group_by(STATION_NAME) %>% summarize(q10 = quantile(maxp, .1))
maxes2.med <- maxes2.2 %>% ungroup() %>% group_by(STATION_NAME) %>% summarize(q10 = quantile(maxp, .1))

maxes.2 <- left_join(maxes.2, maxes.med, by = "STATION_NAME") %>%
filter(maxp > q10)

maxes2.2 <- left_join(maxes2.2, maxes2.med, by = "STATION_NAME") %>%
 filter(maxp > q10)

# Now filter out zero valued maximums as they will not be used in the lognormal distribution fitting
# check to make sure that we have at least 5 non zero years prior to fit. 
maxes.2 = maxes.2 %>% filter(maxp > 0) 
maxes2.2 = maxes2.2 %>% filter(maxp > 0) 

numMax.2 = maxes.2 %>% group_by(STATION_NAME) %>% count(.) %>% filter(n >= 5)

maxes.final = maxes.2 %>% filter(is.element(STATION_NAME, numMax.2$STATION_NAME)) # 10776
maxes2.final = maxes2.2 %>% filter(is.element(STATION_NAME, numMax.2$STATION_NAME)) 

# Check to make sure that ALL values are strictly positive
min(maxes.final$maxp)
min(maxes2.final$maxp)


### FIND 98 PERCENTILE FOR COOP STATIONS SWE ESTIMATES
# Function to compute 50 year estimate via lognormal distribution for a vector of values. 
require(fitdistrplus)
mle50 = function(x, pcnt = 0.98, yplot = FALSE){
  tempvec = as.vector(na.omit(x))
  if(min(tempvec <= 0)){stop("All values must be strictly positive for lognormal distribution fitting.")}
  # Fit the lognormal distribution via maximum likihood estimation
  tempdist = fitdist(tempvec, distr = "lnorm", method = "mle") 
  # 50 year estimate is the 98th percentile of the above fitted distribution
  recurrencevalue = qlnorm(pcnt, meanlog = tempdist$estimate[1],
                           sdlog = tempdist$estimate[2])
  # Automatically Generate Goodness of fit plots if desired. 
  if(yplot){
    plot(tempdist)
  }
  
  return(recurrencevalue)
}

mle50_est = maxes.final %>% group_by(STATION_NAME) %>% 
  summarize(yr50 = mle50(maxp), 
            skew = skewness(maxp)) %>% 
  dplyr::arrange(STATION_NAME)
mle50_est2 = maxes2.final %>% group_by(STATION_NAME) %>% 
  summarize(yr50.2 = mle50(maxp)) %>% 
  dplyr::arrange(STATION_NAME)

mle100_est = maxes.final %>% group_by(STATION_NAME) %>% 
  summarize(yr100 = mle50(maxp, 0.99)) %>% 
  dplyr::arrange(STATION_NAME)
mle100_est2 = maxes2.final %>% group_by(STATION_NAME) %>% 
  summarize(yr100.2 = mle50(maxp, 0.99)) %>% 
  dplyr::arrange(STATION_NAME)

maxobs = maxes.final %>% group_by(STATION_NAME) %>% summarize(maxobs = max(maxp)) %>% dplyr::arrange(STATION_NAME)
maxobs2 = maxes2.final %>% group_by(STATION_NAME) %>% summarize(maxobs2 = max(maxp)) %>% dplyr::arrange(STATION_NAME)

# Now its time to throw everything together...
# Filter from the station metadata all stations that passed all imposed filter (should be 403)
# Then add the absolute maximum pressure, recording years, and both 50 yr estimates to the data frame
coop_final = relevantstations_final %>% filter(is.element(STATION_NAME, numMax.2$STATION_NAME)) %>% 
  dplyr::arrange(STATION_NAME) %>% dplyr::select(-n)

# Combine all relevant results together into one common data frame
coop_final = left_join(coop_final, numMax) # retain maximums passing coverage filter
coop_final = left_join(coop_final, numMax.2, by = "STATION_NAME") # retain maximums used in lognormal distribution fitting
coop_final = left_join(coop_final, maxobs)
coop_final = left_join(coop_final, maxobs2)
coop_final = left_join(coop_final, mle50_est)
coop_final = left_join(coop_final, mle50_est2)
coop_final = left_join(coop_final, mle100_est)
coop_final = left_join(coop_final, mle100_est2)

newUtah.final <- coop_final

newUtah.final$RELDIFF <- (newUtah.final$maxobs - newUtah.final$yr50) / newUtah.final$maxobs
newUtah.final$RELDIFF2 <- (newUtah.final$maxobs - newUtah.final$yr100) / newUtah.final$maxobs

newUtah.final$yr50[newUtah.final$RELDIFF < -0.5] <- newUtah.final$maxobs[newUtah.final$RELDIFF < -0.5] * 1.5
newUtah.final$yr100[newUtah.final$RELDIFF2 < -1] <- newUtah.final$maxobs[newUtah.final$RELDIFF2 < -1] * 2

newUtah.final$RELDIFF <- (newUtah.final$maxobs - newUtah.final$yr50) / newUtah.final$maxobs

newUtah.final.sp <- newUtah.final
coordinates(newUtah.final.sp) <- c("LONGITUDE", "LATITUDE")
projection(newUtah.final.sp) <- projection(watersheds)
hucs <- over(newUtah.final.sp, watersheds)
newUtah.final$HUC <- as.character(hucs$HUC12)

save(newUtah.final, file = "../../R Objects/newUtah2.R")

# Checked this with the previous newUtah file and things are very similar. 
# No more than 10 psf differnt and most stations are identical

# Determine stations in Utah as well as SNOTEL stations
newUtah.final <- newUtah.final %>% 
  mutate(SNOTEL = regexpr(pattern = "S$", STATION) > 0 | regexpr(pattern = "^[[:digit:]]", STATION) > 0,
         UT = regexpr(pattern = " UT ", STATION_NAME) > 0)

ut.sub <- newUtah.final %>% filter(UT)
nout.sub <-newUtah.final %>% filter(!UT)

# Outlier detection process
st.name <- newUtah.final$STATION_NAME[abs(newUtah.final$RELDIFF) > 0.5]
st.scale.low <- vector("numeric", length = length(st.name))
st.scale.high <- vector("numeric", length = length(st.name))
count = 1
for(i in st.name){
  temp <- final_maxes %>% filter(STATION_NAME == i)
  obs <- temp$maxp[temp$maxp > 0]
  par(mfrow = c(2,2))
  boxplot(log(obs), main = i)
  boxplot(obs, main = i)
  hist(obs, main = paste(newUtah.final$RELDIFF[newUtah.final$STATION_NAME == i]))
  hist(log(obs), main = paste(newUtah.final$RELDIFF[newUtah.final$STATION_NAME == i]))
  obs.scale <- (obs - mean(obs)) / sd(obs)
  st.scale.low[count] <- min(obs.scale)
  st.scale.high[count] <- max(obs.scale)
  
  count = count + 1
}

newUtah.final$st.low <- st.scale.low
newUtah.final$st.high <- st.scale.high



### COMPARISONS TO THE ORIGINAL FORM
load("../../R Objects/cleanedDataSets.R")

# Make the station names match the form in the new utah dataset
ut.final <- utahdata.final2 %>% 
  mutate(STATION_NAME = gsub(x=STATION_NAME, pattern = "[[:punct:]]", replacement = ""), 
         STATION_NAME = toupper(STATION_NAME))

# See which stations were left out in the new version
added <- anti_join(newUtah.final, ut.final, by = "STATION_NAME")

# Now see which stations were added
left.out <- anti_join(ut.final, newUtah.final, by = "STATION_NAME")


compare <- left_join(newUtah.final, ut.final, by = "STATION_NAME") %>% dplyr::select(STATION_NAME, maxobs.x, maxobs.y, yr50.x, yr50.y)

compare <- compare %>% mutate(maxRelDiff = ((maxobs.x - maxobs.y)/maxobs.y),
                              yr50RelDiff = ((yr50.x - yr50.y)/yr50.y))




# Code to prepare final version for publication

load("../../R Objects/newUtah2.R")

newUtah.sp <- newUtah.final
coordinates(newUtah.sp) <- c("LONGITUDE", "LATITUDE")
projection(newUtah.sp) <- projection(utah)
newUtah.final$COUNTY <- over(newUtah.sp, utah)$COUNTYNBR

utahFinal <- newUtah.final %>% ungroup(.) %>%
  mutate(LATITUDE = round(LATITUDE,3),
         LONGITUDE = round(LONGITUDE,3),
         ELEVATION = round(ELEVATION*3.28084),
         n = n.x,
         maxobs = round(maxobs),
         yr50 = round(yr50),
         yr100 = round(yr100), 
         STATE = str_extract(STATION_NAME, " [ACINUW][TODYMZV] "),
         STATE = gsub(STATE, pattern = " ", replacement = ""),
         STATION_NAME = gsub(STATION_NAME, pattern = "[[:alpha:]][[:alpha:]][[:space:]]US$", replacement = ""),
         STATION_NAME = gsub(STATION_NAME, pattern = "HEADQUARTERS", replacement = "HQS"),
         STATION_NAME = gsub(STATION_NAME, pattern = "INTERNATIONAL", replacement = "INTL"),
         STATION_NAME = gsub(STATION_NAME, pattern = "AIRPORT", replacement = "ARPT"),
         STATION_NAME = gsub(STATION_NAME, pattern = "- combined", replacement = ""),
         STATION_NAME = gsub(STATION_NAME, pattern = "POWERHOU", replacement = "PWRHS"),
         STATION_NAME = gsub(STATION_NAME, pattern = "- [[:digit:]]", replacement = ""),
         STATION_NAME = gsub(STATION_NAME, pattern = "NUMBER", replacement = "")
         ) %>%
  dplyr::select(STATION, STATION_NAME, STATE, COUNTY, LATITUDE, LONGITUDE, ELEVATION, n, maxobs, yr50, yr100)

utahFinal$STATE[is.na(utahFinal$STATE)] <- c("CO", "ID", "WY", "WY", "ID", "ID")

utahFinal.ut <- filter(utahFinal, STATE == "UT") %>% dplyr::select(-STATE) %>%
  arrange(COUNTY, STATION_NAME)
utahFinal.other <- filter(utahFinal, STATE != "UT") %>% dplyr::select(-COUNTY) %>%
  arrange(STATE, STATION_NAME)

write.csv(utahFinal.ut, file = "../../utFinal.csv", row.names = FALSE, quote = FALSE)
write.csv(utahFinal.other, file = "../../otherFinal.csv", row.names = FALSE, quote = FALSE)



# Now format post office table
newUtah.final$RESPONSE <- newUtah.final$yr50
newUtah.final <- data.frame(newUtah.final)
postOffice <- read.csv("../../Data/postOffices.csv")
postOffice <- postOffice %>% 
  mutate(CITY = gsub(x = NAME, pattern = " Post Office.*$", replacement = ""),
         LON = gsub(x = LON, pattern = "W", replacement = ""),
         LAT = gsub(x = LAT, pattern = "N", replacement = ""),
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
  dplyr::select(-LON, -LAT)


watersheds = readOGR(dsn = "../../Shapefiles and DEMs/ShapeFiles/extractHUC/hydrologic_units",
                     layer = "wbdhu12_a_extract")
watersheds <- spTransform(watersheds, CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
                                          +datum=WGS84 +no_defs +towgs84=0,0,0"))
postOffice.sp <- postOffice
coordinates(postOffice.sp) <- c("LONGITUDE", "LATITUDE")
projection(postOffice.sp) <- projection(watersheds)

postOffice$HUC <- over(postOffice.sp, watersheds)$HUC12

test.prism <- prism2(newUtah.final, postOffice, a = 2, b = 1, d = 3, Fd = .8, zx = 1500, zm = 100,
                     rm = 50, cluster = 100, basin = postOffice$HUC, feet = FALSE)
test.uk <- kriging2(newUtah.final, postOffice, 
                    model = vgm(psill = .17, model = "Sph", range = 198, nugget = .07))
test.idw <- IVD2(newUtah.final, postOffice, NGSL = TRUE)

postOffice$PRISM <- test.prism
postOffice$UK <- test.uk
postOffice$IDW <- test.idw

po.final <- postOffice %>% arrange(COUNTY) %>%
  mutate(COUNTY = as.numeric(as.factor(COUNTY)),
         NAME = gsub(NAME, pattern = "Post Office$", replacement = ""),
         NAME = gsub(NAME, pattern = "Post Office", replacement = "-"),
         NAME = gsub(NAME, pattern = "- \\(", replacement = " \\("),
         NAME = gsub(NAME, pattern = "- -", replacement = "-"),
         NAME = gsub(NAME, pattern = "Utah State University", replacement = "USU"),
         NAME = gsub(NAME, pattern = "[[:space:]]+$", replacement = ""),
         NAME = gsub(NAME, pattern = "[[:space:]]+", replacement = " "),
         ELEVATION = ELEVATION*3.28084,
         LATITUDE = round(LATITUDE,3),
         LONGITUDE = round(LONGITUDE,3),
         ELEVATION = round(ELEVATION),
         PRISM = round(PRISM),
         UK = round(UK)) %>%
  dplyr::select(NAME, COUNTY, LATITUDE, LONGITUDE, ELEVATION, PRISM, UK)

write.csv(po.final, file = "../../postoffices.csv", row.names = FALSE, quote = FALSE)

