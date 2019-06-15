library(snowload)
source("cross_validation_results.R")

data(id2015)
tdata <- id2015
sp::coordinates(tdata) <- c("LONGITUDE", "LATITUDE")
sp::proj4string(tdata) <- sp::proj4string(utdem)

# Specify cores and number of trials.
cores = 25
trials = 100

# Create an register the clusters.
`%dopar%` <- foreach::`%dopar%`
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

errors3 <- try(foreach::foreach(i = 1:trials,
                                .packages = c("snowload",
                                              "sp",
                                              "gstat",
                                              "stats")) %dopar%
                 snowload_trials(ut2017)
)

if(inherits(errors3, "try-error")){
  parallel::stopCluster(cl) # Close the parallel connection.
  print("Error in parallel execution. Closing connection...")
}

parallel::stopCluster(cl) # Close the parallel connection.

save(errors3, "./data/cvResults3.R")
