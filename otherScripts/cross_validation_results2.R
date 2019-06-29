# Write a function that will return the cross validated error for one of the three
# snowload datasets we are comparing to.
snowload_trials2 <- function(tdata, message = FALSE){
  # Pre-allocate a three dimensional array to store the errors.
  errors <- matrix(0, nrow = nrow(tdata), ncol = 7)
  ### PRISM predictions
  # Notice that we tune the predictions along the way for
  # some parameters.
  fun <- "prism"
  preds <- try(crv_pred(tdata, formula = log(yr50 + 1) ~ ELEVATION,
                    fun = fun,
                    message = message,
                    distImp = 0.8,
                    minRad = c(30, 50, 75, 100),
                    wdistance = 2,
                    welevRange = list(lwr = c(100, 200, 300),
                                      upr = c(500, 1000, 1500, 2000, 2500)),
                    welevation = c("ELEVATION", 1),
                    wbasin = c("HUC", 1:5),
                    bound = FALSE,
                    weights = FALSE,
                    tune_folds = 10,
                    tune_score = "MAE"), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }


  # Record the prism CV errors
  errors[, 1] <- tdata$yr50 - (exp(preds) - 1)

  ### SKLM predictions
  fun <- "rkriging"
  preds <- try(crv_pred(tdata, formula = log(yr50 + 1) ~ ELEVATION,
                        fun = fun,
                        model = gstat::vgm("Sph"), sklm = TRUE,
                        debug.level = 0), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  # record the sklm cv errors
  errors[, 2] <- tdata$yr50 - (exp(preds) - 1)

  ### UK predictions
  fun <- "rkriging"
  preds <- try(crv_pred(tdata, formula = log(yr50 + 1) ~ ELEVATION,
                    fun = fun,
                    model = gstat::vgm("Sph"), sklm = FALSE,
                    debug.level = 0), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }


  # record the sklm cv errors
  errors[, 3] <- tdata$yr50 - (exp(preds) - 1)

  ### IDW_snow preds (as implemented in Idaho)
  fun <- "idw_snow"
  # Notice no log transform when doing NGSL predictions
  preds <- try(crv_pred(tdata, formula = log(yr50 + 1) ~ ELEVATION,
                    fun = fun,
                    NGSL = TRUE, tlayer = NA,
                    print = FALSE), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  # record the idw_snow cv errors
  errors[, 4] <- tdata$yr50 - (exp(preds) - 1)

  ### tri_snow preds
  fun <- "tri_snow"
  # Notice no log transform when doing NGSL predictions
  preds <- try(crv_pred(tdata, formula = log(yr50+1) ~ ELEVATION,
                    fun = fun), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  # record the tri_snow cv errors
  errors[, 5] <- tdata$yr50 - (exp(preds) - 1)

  ### lr_snow preds
  fun <- "lm_snow"
  preds <- try(crv_pred(tdata, formula = log(yr50 + 1) ~ ELEVATION,
                    fun = fun), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  # record the lr_snow cv errors
  errors[, 6] <- tdata$yr50 - (exp(preds) - 1)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  ### snlwf preds
  fun <- "snlwf"
  preds <- try(crv_pred(tdata, formula = ~ELEVATION,
                    fun = fun), silent = TRUE)

  if(inherits(preds, "try-error")){
    print(paste("Method", fun, "failed to predict..."))
    preds <- rep(NA, nrow(tdata))
  }

  # record the snlwf cv errors
  errors[, 7] <- tdata$yr50 - preds

  return(errors)
}

