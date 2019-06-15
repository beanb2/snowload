library(snowload)
source("cross_validation_results.R")

data(ut2017)
tdata <- ut2017
sp::coordinates(tdata) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(tdata) <- sp::proj4string(utdem)

# Specify cores and number of trials.
cores = 2
trials = 2

# Create an register the clusters.
`%dopar%` <- foreach::`%dopar%`
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

errors1 <- try(foreach::foreach(i = 1:trials,
                                .packages = c("snowload",
                                              "sp",
                                              "gstat",
                                              "stats")) %dopar%
                 snowload_trials(tdata)
)

if(inherits(errors1, "try-error")){
  parallel::stopCluster(cl) # Close the parallel connection.
  print("Error in parallel execution. Closing connection...")
}

parallel::stopCluster(cl) # Close the parallel connection.

save(errors1, file = "cvResults1.R")
