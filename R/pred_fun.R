#' Predict Function for Ranger
#' 
#' Internal function that prepares the predictions of different types of ranger models,
#' including survival models.
#' 
#' @noRd
#' @keywords internal
#' @param model Fitted ranger model.
#' @param newdata Data to predict on.
#' @param survival Cumulative hazards "chf" (default) or probabilities "prob" per time.
#' @param ... Additional arguments passed to ranger's predict function.
#' 
#' @returns A vector or matrix with predictions.
pred_ranger <- function(model, newdata, survival = c("chf", "prob"), ...) {
  survival <- match.arg(survival)
  
  pred <- stats::predict(model, newdata, ...)
  
  if (model$treetype == "Survival") {
    out <- if (survival == "chf") pred$chf else pred$survival
    colnames(out) <- paste0("t", pred$unique.death.times)
  } else {
    out <- pred$predictions
  }
  return(out)
}

