library(snowload)
source("cross_validation_results.R")

data(ut1992)
tdata <- ut1992
sp::coordinates(tdata) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(tdata) <- sp::proj4string(utdem)

# Specify cores and number of trials.
cores = 25
trials = 100

# Create an register the clusters.
`%dopar%` <- foreach::`%dopar%`
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

errors2 <- try(foreach::foreach(i = 1:trials,
                                .packages = c("snowload",
                                              "sp",
                                              "gstat",
                                              "stats")) %dopar%
                 snowload_trials(tdata)
)

if(inherits(errors2, "try-error")){
  parallel::stopCluster(cl) # Close the parallel connection.
  print("Error in parallel execution. Closing connection...")
}

parallel::stopCluster(cl) # Close the parallel connection.

save(errors2, "./data/cvResults2.R")
