#' Predict Function for Ranger
#' 
#' Returns prediction function for different modes of ranger.
#' 
#' @noRd
#' @keywords internal
#' @param treetype The value of `fit$treetype` in a fitted ranger model.
#' @param survival Cumulative hazards "chf" (default) or probabilities "prob" per time.
#' 
#' @returns A function with signature f(model, newdata, ...).
create_ranger_pred_fun <- function(treetype, survival = c("chf", "prob")) {
  survival <- match.arg(survival)
  
  if (treetype != "Survival") {
    pred_fun <- function(model, newdata, ...) {
      stats::predict(model, newdata, ...)$predictions
    }
    return(pred_fun)
  }
  
  if (survival == "prob") {
    survival <- "survival"
  }
  
  pred_fun <- function(model, newdata, ...) {
    pred <- stats::predict(model, newdata, ...)
    out <- pred[[survival]]
    colnames(out) <- paste0("t", pred$unique.death.times)
    return(out)
  }
  return(pred_fun)
}

