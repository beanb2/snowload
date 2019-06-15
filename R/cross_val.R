# This file contains a set of functions to perform cross validation.
# The main function is the wrapper function "crv" which provides
# a common cross validation framework for a set of allowable functions.

#' Perform cross validation for the snowload functions
#'
#' This function is a wrapper function to perform k-fold cross validation
#' on the various functions in the snowload package.
#'
#' @param tdata The dataset on which cross validation will be performed.
#' @param fun The name of the function on which to perform cross validation
#'   this must be a function from the snowload package.
#' @param folds The number of groups to separate the data when performing
#'   cross validation. The default is 10, the minimum is 2 and the maximum
#'   is the number of observations.
#' @param score A character string of "MAE" for mean absolute error,
#'   "RMSE" for root mean square error, "ME" for mean error, and
#'   "none" when you would like the function to simply return the
#'   cross validation predictions.
#' @param ... Additional arguments to the function specified in the fun
#'   argument.
#'
#' @return A list of two, with the first item being a single numeric value
#'   specifying the MAE or RMSE. The second is a vector with the full set of
#'   errors resulting from the cross validation.
#'
#' @export
crv_pred <- function(tdata, formula, fun, folds = 10,
                score = "None", message = FALSE, ...){
  # Split data into k groups
  groups = rep(1:folds, length = nrow(tdata))
  # Randomize group assignment
  groups = sample(groups)

  # Provide a caveat for the snlwf function
  if(fun == "snlwf"){
    preds <- snlwf(formula, tdata, ...)
  }else{
    preds <- vector("numeric", nrow(tdata))
    for(i in 1:folds){
      # If requested, print the time when each of the folds is being processed.
      if(message){
        print(Sys.time())
        print(paste("Round", i, "of", folds, "...", sep = " "))
      }

      # get() treats the character call as a function and makes the predictions
      # This form relies on the fact that the predictions are stored in a
      # column with the same name as the function.
      f <- get(fun)
      preds[groups == i] <- f(formula, tdata[groups != i,],
                              tdata[groups == i, ], ...)
    } # End the for-loop
  }

  # If a score is supplied, only return the score.
  if(score != "None"){
    tdf <- model.frame(formula, tdata)

    errors <- tdf[, 1] - preds

    if(score == "MAE"){
      return(mean(abs(errors)))
    }else if(score== "RMSE"){
      return(sqrt(mean(errors^2)))
    }else if(score== "ME"){
      return(mean(errors))
    }else{ # Throw an error for invalid scoring s.
      stop("Valid scoring methodss are \'MAE\', \'RMSE\' and \'ME\'")
    }

  }else{   # If a score is not supplied, return the predictions.
    return(preds)
  }
}

