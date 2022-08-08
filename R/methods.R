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
  S <- ks_shap_values(x)
  n <- min(n, nrow(S))
  cat(
    "'kernelshap' object representing \n  - SHAP matrix of dimension",
    nrow(S), "x", ncol(S),
    "\n  - feature data.frame/matrix of dimension",  nrow(S), "x", ncol(S),
    "\n  - baseline value of", ks_baseline(x)
  )
  cat("\n\n")
  cat("SHAP values of first", n, "observations:\n")
  print(utils::head(S, n))
  cat("\n Corresponding standard errors:\n")
  print(utils::head(ks_standard_errors(x), n))
  cat("\n And the feature values:\n")
  print(utils::head(ks_feature_values(x), n))
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

#' Extractor Functions
#'
#' Functions to extract SHAP values, feature values, standard errors etc. from a "kernelshap" object.
#'
#' @name extractors
#' @param object Object to extract something.
#' @param ... Currently unused.
#' @return The corresponding object is returned, i.e.,
#' \itemize{
#'   \item \code{ks_shap_values()} returns the matrix of SHAP values, 
#'   \item \code{ks_feature_values()} the \code{data.frame} of feature values, 
#'   \item \code{ks_baseline()} the numeric baseline value of the input, 
#'   \item \code{ks_standard_errors()} the matrix of standard errors of SHAP values, 
#'   \item \code{ks_converged()} returns the vector of convergence flags, and finally
#'   \item \code{ks_n_iter()} the number of iterations per row.
#' }
NULL

#' @rdname extractors
#' @export
ks_shap_values <- function(object, ...){
  UseMethod("ks_shap_values")
}

#' @rdname extractors
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' pred_fun <- function(X) stats::predict(fit, X)
#' s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[-1])
#' ks_shap_values(s)
ks_shap_values.kernelshap = function(object, ...) {
  object[["S"]]
}

#' @rdname extractors
#' @export
ks_shap_values.default = function(object, ...) {
  stop("No default method available.")
}

#' @rdname extractors
#' @export
ks_feature_values <- function(object, ...){
  UseMethod("ks_feature_values")
}

#' @rdname extractors
#' @export
ks_feature_values.kernelshap = function(object, ...) {
  object[["X"]]
}

#' @rdname extractors
#' @export
ks_feature_values.default = function(object, ...) {
  stop("No default method available.")
}

#' @rdname extractors
#' @export
ks_baseline <- function(object, ...){
  UseMethod("ks_baseline")
}

#' @rdname extractors
#' @export
ks_baseline.kernelshap = function(object, ...) {
  object[["baseline"]]
}

#' @rdname extractors
#' @export
ks_baseline.default = function(object, ...) {
  stop("No default method available.")
}

#' @rdname extractors
#' @export
ks_standard_errors <- function(object, ...){
  UseMethod("ks_standard_errors")
}

#' @rdname extractors
#' @export
ks_standard_errors.kernelshap = function(object, ...) {
  object[["SE"]]
}

#' @rdname extractors
#' @export
ks_standard_errors.default = function(object, ...) {
  stop("No default method available.")
}

#' @rdname extractors
#' @export
ks_n_iter <- function(object, ...){
  UseMethod("ks_n_iter")
}

#' @rdname extractors
#' @export
ks_n_iter.kernelshap = function(object, ...) {
  object[["n_iter"]]
}

#' @rdname extractors
#' @export
ks_n_iter.default = function(object, ...) {
  stop("No default method available.")
}

#' @rdname extractors
#' @export
ks_converged <- function(object, ...){
  UseMethod("ks_converged")
}

#' @rdname extractors
#' @export
ks_converged.kernelshap = function(object, ...) {
  object[["converged"]]
}

#' @rdname extractors
#' @export
ks_converged.default = function(object, ...) {
  stop("No default method available.")
}

