Sys.time()
library(snowload)
source("otherScripts/cross_validation_results.R")

data(ut2017)
tdata <- ut2017
sp::coordinates(tdata) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(tdata) <- sp::proj4string(utdem)

# Specify cores and number of trials.
cores = 4
trials = 4

# Create an register the clusters.
# `%dopar%` <- foreach::`%dopar%`
# https://stackoverflow.com/questions/8358098/how-to-set-seed-for-random-simulations-with-foreach-and-domc-packages/41109556
`%dorng%` <- doRNG::`%dorng%`
`%dopar%` <- foreach::`%dopar%`
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

set.seed(123)
errors1 <- try(foreach::foreach(i = 1:trials,
                                .packages = c("snowload",
                                              "sp",
                                              "gstat",
                                              "stats")) %dorng%
                 {snowload_trials(tdata)}
)

if(inherits(errors1, "try-error")){
  parallel::stopCluster(cl) # Close the parallel connection.
  print("Error in parallel execution. Closing connection...")
}

set.seed(123)
errors2 <- try(foreach::foreach(i = 1:trials,
                                .packages = c("snowload",
                                              "sp",
                                              "gstat",
                                              "stats")) %dorng%
                 {snowload_trials(tdata)}
)

parallel::stopCluster(cl) # Close the parallel connection.

# save(errors1, file = "data-raw/cvResults1.R")
Sys.time()

# identical(errors1, errors2)
# # Test script to see how large I really need to make the grid for tri_snow.
# dens <- seq(100, 900, 200)
# check <- check2 <- vector("numeric", length = length(dens))
# for(i in 1:length(dens)){
#   print(i)
#   check[i] <- crv_pred(tdata, formula = yr50 ~ 1, fun = "tri_snow",
#                        density = c(dens[i], dens[i]), score = "MAE")
#   check2[i] <- crv_pred(tdata, formula = yr50 ~ ELEVATION, fun = "tri_snow",
#                     density = c(dens[i], dens[i]), score = "MAE")
# }

data(id2015)
tdata <- id2015
sp::coordinates(tdata) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(tdata) <- sp::proj4string(utdem)

fun = "idw_snow"
preds <- try(crv_pred(tdata, formula = log(yr50) ~ ELEVATION,
                      fun = fun,
                      NGSL = TRUE, tlayer = 2500,
                      debug.level = 0), silent = TRUE)

# mean(abs(preds - id2015$yr50))
mean(abs(exp(preds) - id2015$yr50))

check <- function(tdata, layer, ...){
  high <- tdata[tdata$ELEVATION > layer, ]
  low <- tdata[tdata$ELEVATION <= layer, ]

  lowcor <- cor(low$ELEVATION, low$yr50/low$ELEVATION, use = "na.or.complete", ...)
  highcor <- cor(high$ELEVATION, high$yr50/high$ELEVATION, use = "na.or.complete", ...)

  # If the correlation yields a missing value (due to a lack of data),
  # return a 0.
  if(is.na(lowcor)){lowcor <- 0}
  if(is.na(highcor)){highcor <- 0}

  # Return a weighted average of the absolute correlations based on elevation
  return( (nrow(high)*abs(highcor) + nrow(low)*abs(lowcor)) / nrow(tdata) )
}

