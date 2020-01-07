#' Title: Cross Validation Function Driver
#' 
#' A function that implements an ecoregion-based approach of universal kriging
#'
#' @param 
#' 
#' @return An object with reponse variable predictions for 
#' each gridcell, and row and column names representing the latitude 
#' and longitude coordinates respectively 
#'
#' @examples
crv = function(k, tdata, method, maxel = FALSE, maxout = FALSE, 
               d = 1:10, rm = c(30, 50, 75, 100), 
               zx = c(500, 1000, 1500, 2000, 2500),
               zm = c(100, 200, 300), print = FALSE, RMSE = FALSE,
               lc = c(0, 500, 1000, 1500, 2000, 2500, 10000),
               ...){
  
  # Make sure user specfies a method that we have defined
  if(!is.element(method, c("PRISM", "PRISM2", "UK", "UK2",
                           "UK2.1", "UK3", "IDW"))){
    error("Please select one of the following methods: 
          PRISM, PRISM2, UK, UK2, UK3, IDW.")
  }
  # Split data into k groups
  groups = rep(1:k, length = nrow(tdata))
  # Randomize group assignment
  groups = sample(groups)
  
  # Assign each station to a group
  tdata$GROUPS = groups
  tdata$PREDICT = 0
  
  # Do we want restrictions on the extrapolated predictions?
  # We will store these results. 
  maxel.2 <- maxel
  maxout.2 <- maxout
  
  tvalues <- vector("list", k)
  for(i in 1:k){
    if(print){
      print(Sys.time())
      print(paste("Round", i, "of", k, "...", sep = " "))
    }
    test = tdata[tdata$GROUPS == i,]
    train = tdata[tdata$GROUPS != i,]
    
    # Predict accoring to the appropriate method
    ### UK TRIALS ###
    if(method == "UK"){
      tvgm <- gstat::gstat(NULL, id = "tkrige", 
                           formula = log(RESPONSE + 1) ~ ELEVATION, 
                           locations = train)
      
      tvalues[[i]] <- gstat::fit.variogram(gstat::variogram(tvgm), gstat::vgm("Sph"))
      tdata$PREDICT[tdata$GROUPS == i] = UK(train, test, 
                                            maxel = maxel.2, 
                                            maxout = maxout.2, model = tvalues[[i]], ...)$predict
    }
    ###                                                                        ### 
    
    ### UK2 TRIALS ###
    if(method == "UK2"){
      tvgm <- gstat::gstat(NULL, id = "tkrige", 
                           formula = log(RESPONSE + 1) ~ ELEVATION, 
                           locations = train)
      
      tvalues[[i]] <- gstat::fit.variogram(gstat::variogram(tvgm), gstat::vgm("Sph"))
      tdata$PREDICT[tdata$GROUPS == i] = UK2(train, test, 
                                             maxel = maxel.2, 
                                             maxout = maxout.2, model = tvalues[[i]], ...)
    }
    ###                                                                        ### 
    
    ### UK2 TRIALS ###
    if(method == "UK2.1"){
      tvgm <- gstat::gstat(NULL, id = "tkrige", 
                           formula = log(RESPONSE + 1) ~ ELEVATION, 
                           locations = train)
      
      tvalues[[i]] <- gstat::fit.variogram(gstat::variogram(tvgm), gstat::vgm("Sph"))
      tdata$PREDICT[tdata$GROUPS == i] = UK2(train, test, 
                                             maxel = maxel.2, 
                                             maxout = maxout.2, model = tvalues[[i]], 
                                             maxdist = tvalues[[i]][2, 3], ...)
    }
    ###                                                                        ### 
    
    ### UK3 TRIALS ###
    if(method == "UK3"){
      tvgm <- gstat::gstat(NULL, id = "tkrige", 
                           formula = log(RESPONSE + 1) ~ ELEVATION, 
                           locations = train)
      
      tvalues[[i]] <- gstat::fit.variogram(gstat::variogram(tvgm), gstat::vgm("Sph"))
      tdata$PREDICT[tdata$GROUPS == i] = UK3(train, test, 
                                             maxel = maxel.2, 
                                             maxout = maxout.2, model = tvalues[[i]], ...)
    }
    ###                                                                        ### 
    
    ### PRISM TRIALS ### 
    if(method == "PRISM"){
      # Rough tune using only 9/10 of the data
      tlength <- length(d)*length(rm)*length(zx)*length(zm)
      if(tlength > 1){
        values <- matrix(0, nrow = tlength, ncol = 5)
        count = 1
        for(q in d){
          for(j in rm){
            for(m in zx){
              for(l in zm){
                # Split data into k groups
                groups2 = rep(1:k, length = nrow(train))
                # Randomize group assignment
                groups2 = sample(groups2)
                
                # Assign each station to a group
                train$GROUPS2 = groups2
                train$PREDICT2 = 0
                for(z in 1:k){
                  test.sub = train[train$GROUPS2 == z,]
                  train.sub = train[train$GROUPS2 != z,]
                  
                  train$PREDICT2[train$GROUPS2 == z] = 
                    PRISM(train.sub, test.sub, a = 2, b = 1, d = q, 
                          Fd = .8, zx = m, zm = l, rm = j)[[1]]$predictions
                } # end cross validation for loop
                # By default, minimize RMSE, but minimize MAE if specified. 
                if(RMSE){
                  rmse <- sqrt(mean((train$PREDICT2 - train$RESPONSE)^2))
                }else{
                  rmse <- mean(abs(train$PREDICT2 - train$RESPONSE))
                }
                values[count,] = c(q, j, m, l, rmse)
                count = count + 1
              } # end minimum elevation loop
            } # end maximum elevation loop
          } # end radius of influence loop
        } # end eco region loop
        colnames(values) = c("d", "rm", "zx", "zm", "rmse")
        values <- values[base::order(values[, 5]), ]
        values <- values[1, ]
      }else{
        values <- c(d, rm, zx, zm, NA)
      }
      tvalues[[i]] <- values
      tdata$PREDICT[tdata$GROUPS == i] = 
        PRISM(train, test, a = 2, b = 1, d = values[1], 
              Fd = .8, zx = values[3], 
              zm = values[4], rm = values[2])[[1]]$predictions
    }
    ###                                                                        ### 
    
    ### PRISM TRIALS ### 
    if(method == "PRISM2"){
      # Rough tune using only 9/10 of the data
      tlength <- length(d)*length(rm)*length(zx)*length(zm)
      if(tlength > 1){
        values <- matrix(0, nrow = length(d)*length(rm)*length(zx)*length(zm), 
                         ncol = 5)
        count = 1
        for(q in d){
          for(j in rm){
            for(m in zx){
              for(l in zm){
                # Split data into k groups
                groups2 = rep(1:k, length = nrow(train))
                # Randomize group assignment
                groups2 = sample(groups2)
                
                # Assign each station to a group
                train$GROUPS2 = groups2
                train$PREDICT2 = 0
                for(z in 1:k){
                  test.sub = train[train$GROUPS2 == z,]
                  train.sub = train[train$GROUPS2 != z,]
                  
                  train$PREDICT2[train$GROUPS2 == z] = 
                    PRISM2(train.sub, test.sub, a = 2, b = 1, d = q, 
                           Fd = .8, zx = m, zm = l, rm = j)[[1]]$predictions
                } # end cross validation for loop
                if(RMSE){
                  rmse <- sqrt(mean((train$PREDICT2 - train$RESPONSE)^2))
                }else{
                  rmse <- mean(abs(train$PREDICT2 - train$RESPONSE))
                }
                values[count,] = c(q, j, m, l, rmse)
                count = count + 1
              } # end minimum elevation loop
            } # end maximum elevation loop
          } # end radius of influence loop
        } # end eco region loop
        colnames(values) = c("d", "rm", "zx", "zm", "rmse")
        values <- values[base::order(values[, 5]), ]
        values <- values[1, ]
      }else{
        values <- c(d, rm, zx, zm, NA)
      }
      tvalues[[i]] <- values
      tdata$PREDICT[tdata$GROUPS == i] = 
        PRISM2(train, test, a = 2, b = 1, d = values[1], 
               Fd = .8, zx = values[3], 
               zm = values[4], rm = values[2])[[1]]$predictions
    }
    ###                                                                        ### 
    
    ### IDW TRIALS ### 
    if(method == "IDW"){
      values <- matrix(0, nrow = length(lc), ncol = 2)
      if(length(lc) > 1){
        count = 1
        for(q in lc){
          # Split data into k groups
          groups2 = rep(1:k, length = nrow(train))
          # Randomize group assignment
          groups2 = sample(groups2)
          
          # Assign each station to a group
          train$GROUPS2 = groups2
          train$PREDICT2 = 0
          for(z in 1:k){
            test.sub = train[train$GROUPS2 == z,]
            train.sub = train[train$GROUPS2 != z,]
            
            train$PREDICT2[train$GROUPS2 == z] = IDW(train.sub, test.sub, tlayer = q,
                                                     maxout = maxout.2, ...)
          } # end cross validation for loop
          if(RMSE){
            rmse <- sqrt(mean((train$PREDICT2 - train$RESPONSE)^2))
          }else{
            rmse <- mean(abs(train$PREDICT2 - train$RESPONSE))
          }
          values[count,] = c(q, rmse)
          count = count + 1
        } # end the layer loop
        colnames(values) = c("tlayer", "rmse")
        values <- values[base::order(values[, 2]), ]
        values <- values[1, ]
      }else{
        values <- lc
      }
      tvalues[[i]] <- values
      tdata$PREDICT[tdata$GROUPS == i] = 
        IDW(train, test, tlayer = values[1],
            maxout = maxout.2, ...)
      
    }
    ###                                                                        ### 
    
  } # End cross validation for-loop
  
  # Calculate error (same way we would calculate a residual)
  tdata$ERR = tdata$RESPONSE - tdata$PREDICT 
  
  return(list(data.frame(cuID = tdata$cuID, PREDICT = tdata$PREDICT, 
                         ERROR = tdata$ERR), tvalues))
}




