# Function to compute the weights of the universal kriging 
# for purposes of visualization.

getKW <- function(train, test, tmodel){
  fdf <- data.frame(LONGITUDE = rep(mean(train$LONGITUDE), 2), 
                    LATITUDE = rep(mean(train$LATITUDE), 2),
                    ELEVATION = c(0, 1))
  sp::coordinates(fdf) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(fdf) <- raster::projection(train)
  
  test2 <- gstat::gstat(formula = log(RESPONSE+1)~ELEVATION, 
                        locations = train, 
                        model = tmodel)
  
  coef <- predict(test2, fdf, BLUE = TRUE)
  
  resid <- log(train$RESPONSE + 1) - 
    (coef$var1.pred[1] + (coef$var1.pred[2] - 
                            coef$var1.pred[1])*train$ELEVATION)
  
  train$RESIDUALS <- resid
  
  # Recreate kriging of the residuals.
  distM <- fields::rdist.earth(train@coords, train@coords, miles = FALSE)
  distL <- fields::rdist.earth(train@coords, test@coords, miles = FALSE)
  
  vgmM <- gstat::variogramLine(tmodel, dist_vector = distM)
  vgmM <- vgmM - diag(diag(vgmM))
  vgmL <- gstat::variogramLine(tmodel, dist_vector = distL)
  
  ls <- cbind(vgmM, c(rep(-1, nrow(vgmM))))
  ls <- rbind(ls, c(rep(1, ncol(vgmM)), 0))
  
  tpred <- vector("numeric", nrow(test))
  weights <- matrix(0, nrow(train), ncol = nrow(test))
  for(i in 1:nrow(test)){
    rs <- c(vgmL[, i], 1)
    
    tweights <- solve(ls, as.matrix(rs, ncol = 1))
    weights[, i] <- tweights[-length(tweights)]
    
    tpred[i] <- sum(as.vector(train$RESIDUALS) * as.vector(weights[, i]))
  }
  finalPred <- (coef$var1.pred[1] + (coef$var1.pred[2] - 
                                       coef$var1.pred[1])*test$ELEVATION) + tpred
  return(list(finalPred, tpred, weights, as.vector(coef$var1.pred)))
}