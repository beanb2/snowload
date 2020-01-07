library(intkrige)
data(utsnow)
utsnow$RESPONSE <- utsnow$pointDL

utdem <- raster::raster("data-raw/UTDEMC")
utdem2 <- utdem
utdem <- as(utdem, "SpatialPixelsDataFrame")
coordinates(utsnow) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(utsnow) <- sp::proj4string(utdem)

watersheds <- rgdal::readOGR(dsn = "data-raw", layer = "naWatershed")
watersheds <- sp::spTransform(watersheds, utdem@proj4string)

utdem$HUC <- over(utdem, watersheds)$HUC8
utdem$ELEVATION <- utdem$layer
utsnow$HUC <- over(utsnow, watersheds)$HUC8

test1 <- IDW(utsnow, utdem)
test2 <- IDW_snow(c("RESPONSE", "ELEVATION"), utsnow, utdem)

all.equal(test1, test2)



tester <- snlwf2(utdem, "layer")
sp::gridded(tester) <- TRUE
testRast <- raster::raster(tester["snlw"])
plot(testRast)

x <- as.vector(tester$snlw)
x[is.na(x)] <- 0

tester2 <- snlwf(utdem2, path1 = "data-raw/Counties", path2 = "data-raw/Utahsnowlaw.csv")

all.equal(x, as.vector(tester2))




