# Check to see how post office predictions in Utah compare to those given in the report.
library(snowload)

# Assign stations to appropriate water basins
# Read in water basin data assign to various HUCS
watersheds = rgdal::readOGR(dsn = "data-raw",
                            layer = "naWatershed")

data(ut2017)
sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(ut2017) <- sp::proj4string(watersheds)

# Now Read in the post office information
postOffice <- read.csv("data-raw/postofficepred.csv")
postOffice$ELEVATION <- postOffice$ELEVATION*0.3048
sp::coordinates(postOffice) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(postOffice) <- sp::proj4string(watersheds)

postOffice$HUC <- sp::over(postOffice, watersheds)$HUC8

mont <- postOffice[postOffice$NAME == "Monticello", ]

temp <- prism(log(yr50) ~ ELEVATION, locations = ut2017,
              newdata = mont, distImp = 0.8, minRad = 50,
              wdistance = 2, welevRange = list(lwr = 100, upr = 1500),
              welevation = c("ELEVATION", 1), wbasin = c("HUC", 3))
temp2 <- rkriging(log(yr50) ~ ELEVATION, locations = ut2017, newdata = mont,
                  model = gstat::vgm("Sph"))

exp(temp)/.04788
exp(temp2)/.04788


mont2 <- mont
mont2$ELEVATION <- 2078.736

temp.2 <- prism(log(yr50) ~ ELEVATION, locations = ut2017,
              newdata = mont2, distImp = 0.8, minRad = 50,
              wdistance = 2, welevRange = list(lwr = 100, upr = 1500),
              welevation = c("ELEVATION", 1), wbasin = c("HUC", 3))
temp2.2 <- rkriging(log(yr50) ~ ELEVATION, locations = ut2017, newdata = mont2,
                  model = gstat::vgm("Sph"))

exp(temp.2)/.04788
exp(temp2.2)/.04788






postOffice$PRISM2 <- exp(temp)/.04788

View(as.data.frame(postOffice))
