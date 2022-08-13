#' Prints "kernelshap" Object
#'
#' @param x An object of class "kernelshap".
#' @param n Maximum number of rows of SHAP values, standard errors and feature values to print.
#' @param ... Further arguments passed from other methods.
#' @return Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' pred_fun <- function(X) stats::predict(fit, X)
#' s <- kernelshap(iris[1:3, -1], pred_fun = pred_fun, iris[-1])
#' s
#' @seealso \code{\link{kernelshap}}.
print.kernelshap <- function(x, n = 2L, ...) {
  S <- ks_extract(x, "S")
  n <- min(n, nrow(S))
  cat(
    "'kernelshap' object representing \n  - SHAP matrix of dimension",
    nrow(S), "x", ncol(S),
    "\n  - feature data.frame/matrix of dimension",  nrow(S), "x", ncol(S),
    "\n  - baseline value of", ks_extract(x, "baseline"),
    "\n  - average number of iterations of", mean(ks_extract(x, "n_iter"))
  )
  cat("\n\n")
  cat("SHAP values of first", n, "observations:\n")
  print(utils::head(S, n))
  cat("\n Corresponding standard errors:\n")
  print(utils::head(ks_extract(x, "SE"), n))
  cat("\n And the feature values:\n")
  print(utils::head(ks_extract(x, "X"), n))
  cat("\n")
  invisible(x)
}

#' Check for kernelshap
#'
#' Is object of class "kernelshap"?
#'
#' @param object An R object.
#' @return Returns \code{TRUE} if \code{object} has "\code{kernelshap}" among its classes, and \code{FALSE} otherwise.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' pred_fun <- function(X) stats::predict(fit, X)
#' s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[-1])
#' is.kernelshap(s)
#' is.kernelshap("a")
is.kernelshap <- function(object){
  inherits(object, "kernelshap")
}

#' Extractor Function
#'
#' Function to extract an element of a "kernelshap" object, e.g., the SHAP values "S".
#'
#' @param object Object to extract something.
#' @param what Element to extract. One of "S", "X", "baseline", "SE", "n_iter", or "converged".
#' @param ... Currently unused.
#' @return The corresponding object is returned.
#' @export
ks_extract <- function(object, ...){
  UseMethod("ks_extract")
}

#' @describeIn ks_extract Method for "kernelshap" object.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' pred_fun <- function(X) stats::predict(fit, X)
#' s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[-1])
#' ks_extract(s, what = "S")
ks_extract.kernelshap = function(object, what = c("S", "X", "baseline", "SE", "n_iter", "converged"), ...) {
  what <- match.arg(what)
  object[[what]]
}

#' @describeIn ks_extract No default method available.
#' @export
ks_extract.default = function(object, ...) {
  stop("ks_extract() only accepts 'kernelshap' objects.")
}
