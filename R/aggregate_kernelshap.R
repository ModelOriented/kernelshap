#' Aggregate SHAP objects
#'
#' Useful for reporting SHAP for multiple models in a cross validation schema.
#'
#' @param x a `kernelshap` object or a list of them.
#' @param ... other `kernelshap` objects when `x` is not a list.
#'
#' @details
#' SHAP values, data, standar errors, number of iterations, convergence and predictions of
#' each object are concatenated in the result. The other components of the results are
#' averaged or unique.
#'
#' @return A `kernelshap` object aggregating all the inputs.
#' @export
aggregate_kernelshap <- function(x, ...) {
  if (is.kernelshap(x)) {
    x <- list(x, ...)
  }
  aggr_shap <- list(
    S = do.call(rbind, lapply(x, function(s) s$S)),
    X = do.call(rbind, lapply(x, function(s) s$X)),
    baseline = mean(sapply(x, function(s) s$baseline)),
    SE = do.call(rbind, lapply(x, function(s) s$SE)),
    n_iter = unname(do.call(c, lapply(x, function(s) s$n_iter))),
    converged = unname(do.call(c, lapply(x, function(s) s$converged))),
    m = unique(sapply(x, function(s) s$m)),
    m_exact = unique(sapply(x, function(s) s$m_exact)),
    prop_exact = mean(sapply(x, function(s) s$prop_exact)),
    exact = unique(sapply(x, function(s) s$exact)),
    txt = unique(sapply(x, function(s) s$txt)),
    predictions = do.call(rbind, lapply(x, function(s) s$predictions))
  )
  class(aggr_shap) <- "kernelshap"
  return(aggr_shap)
}
