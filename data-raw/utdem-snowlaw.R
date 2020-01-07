# Run datasetPrep prior to this file.

# This is an aggregated raster from the National Map
# (originally on a 30m by 30m resolution)
# We wish to save an aggregated version of this in our package.
utdem <- raster::raster("data-raw/UTDEMC")
# Aggregate for speed
utdem2 <- raster::aggregate(utdem, fact=8, fun=mean)

utdem <- as(utdem, "SpatialPixelsDataFrame")
utdem2 <- as(utdem2, "SpatialPixelsDataFrame")

watersheds <- rgdal::readOGR(dsn = "data-raw", layer = "naWatershed")
watersheds <- sp::spTransform(watersheds, utdem@proj4string)

utdem$HUC <- over(utdem, watersheds)$HUC8
utdem2$HUC <- over(utdem2, watersheds)$HUC8
utdem$ELEVATION <- utdem$layer
utdem$layer <- NULL

utdem2$ELEVATION <- utdem2$layer
utdem2$layer <- NULL

usethis::use_data(utdem, overwrite = TRUE)


# Read in the Utah county shapefile.
counties <- rgdal::readOGR(dsn = "data-raw/Counties", layer = "Counties")

counties <- sp::spTransform(counties,
                            sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

counties_smooth <- rgeos::gSimplify(counties, tol = .00005,
                                    topologyPreserve = TRUE)
counties_smooth <- as(counties_smooth, "SpatialPolygonsDataFrame")
counties_smooth$COUNTYNBR <- counties$COUNTYNBR
counties_smooth$dummy <- NULL
# counties_smooth <- counties

plot(counties_smooth)
object.size(counties)
object.size(counties_smooth)

utcounty <- counties_smooth

# Now read in the coefficients for the Utah snowload and save it also.
snlw <- read.csv("data-raw/Utahsnowlaw.csv")
snlw$COUNTY <- as.character(snlw$COUNTY)
snlw$P_0 <- as.numeric(as.character(snlw$P_0))
snlw$S <- as.numeric(as.character(snlw$S))
snlw$A_0 <- as.numeric(as.character(snlw$A_0))

# Script to extract predictions for the old functions and ensure that the new
# functions return an equivalent result.

# Load some necessary libraries
library(sp)
library(rgdal)
library(gstat)
library(raster)
library(rgeos)

source("data-raw/paper2functions.R")

data_prism <- ut2017
# Get HUC 12 designations
data_prism$HUC <- as.numeric(as.character(data_prism$HUC))*10000
data_prism$RESPONSE <- data_prism$yr50

dem <- raster(utdem2["ELEVATION"])

prism_predict <- prism2(as.data.frame(data_prism), dem, cluster = 200,
                        basin = as.numeric(as.character(utdem2$HUC))*10000,
                        d = 3, noNegative = FALSE, maxout = FALSE)
# as.vector(prism_predict)

idw_predict <- IVD2(data_prism, dem)
# as.vector(idw_predict)

snlw_predict <- snlwf(dem)
# as.vector(snlw_predict)

sklm_predict <- kriging(as.data.frame(data_prism), dem,
                        maxel = FALSE, maxout = FALSE,
                        model = vgm(psill = .21, model = "Sph",
                                    range = 200, nugget = .06))
# as.vector(sklm_predict)

uk_predict <- kriging2(as.data.frame(data_prism), dem,
                       maxel = FALSE, maxout = FALSE,
                       model = vgm(psill = .21, model = "Sph",
                                   range = 200, nugget = .06))
# as.vector(uk_predict)

tri_predict <- TRI(as.data.frame(data_prism), dem,
                   density = 100)
# as.vector(tri_predict)

lm_predict <- lmPred(as.data.frame(data_prism),
                     as.data.frame(dem), maxout = FALSE)
# as.vector(lm_predict)

test_df <- data.frame(prism = as.vector(prism_predict),
                      idw = as.vector(idw_predict),
                      snlw = as.vector(snlw_predict)*.04788,
                      sklm = as.vector(sklm_predict),
                      uk = as.vector(uk_predict),
                      tri = as.vector(tri_predict),
                      lm = as.vector(lm_predict)
                      )

# A few post processing steps to ensure consistent output
test_df$snlw[test_df$snlw == 0] <- NA

usethis::use_data(utcounty, snlw, utdem2, test_df,
                  internal = TRUE, overwrite = TRUE)




