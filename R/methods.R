#' Prints "kernelshap" Object
#'
#' @param x An object of class "kernelshap".
#' @param compact Set to \code{TRUE} to hide printing the top n SHAP values, 
#' standard errors and feature values. 
#' @param n Maximum number of rows of SHAP values, standard errors and feature values to print.
#' @param ... Further arguments passed from other methods.
#' @return Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:3, -1], bg_X = iris[-1])
#' s
#' @seealso \code{\link{kernelshap}}.
print.kernelshap <- function(x, compact = FALSE, n = 2L, ...) {
  S <- ks_extract(x, "S")
  SE <- ks_extract(x, "SE")
  X <- ks_extract(x, "X")
  if (!is.list(S)) {
    n <- min(n, nrow(S))
    s_text <- paste("- SHAP matrix of dimension",  nrow(S), "x", ncol(S))
  } else {
    n <- min(n, nrow(S[[1L]]))
    s_text <- paste(
      "-", length(S), "SHAP matrices of dimension",  nrow(S[[1L]]), "x", ncol(S[[1L]])
    )
  }
  cat(
    "'kernelshap' object representing \n ", s_text,
    "\n  - feature data.frame/matrix of dimension",  nrow(X), "x", ncol(X),
    "\n  - baseline:", ks_extract(x, "baseline"),
    "\n  - average iterations:", mean(ks_extract(x, "n_iter")),
    "\n  - rows not converged:", sum(!ks_extract(x, "converged"))
  )
  cat("\n")
  if (!compact) {
    cat("\nSHAP values of first", n, "observations:\n")
    if (!is.list(S)) print(utils::head(S, n)) else print(lapply(S, utils::head, n))
    cat("\n Corresponding standard errors:\n")
    if (!is.list(S)) print(utils::head(SE, n)) else print(lapply(SE, utils::head, n))
    cat("\n And the feature values:\n")
    print(utils::head(X, n))
    cat("\n")
  }
  invisible(x)
}

#' Check for kernelshap
#'
#' Is object of class "kernelshap"?
#'
#' @param object An R object.
#' @return Returns \code{TRUE} if \code{object} has "\code{kernelshap}" among its classes, 
#' and \code{FALSE} otherwise.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:2, -1], bg_X = iris[-1])
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
#' @param what Element to extract. One of "S", "X", "baseline", "SE", "n_iter", 
#' "converged", "m", "m_exact", "prop_exact", or "txt".
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
#' s <- kernelshap(fit, iris[1:2, -1], bg_X = iris[-1])
#' ks_extract(s, what = "S")
ks_extract.kernelshap = function(object, 
                                 what = c("S", "X", "baseline", "SE", "n_iter", "converged",
                                          "m", "m_exact", "prop_exact", "txt"), 
                                 ...) {
  what <- match.arg(what)
  object[[what]]
}

#' @describeIn ks_extract No default method available.
#' @export
ks_extract.default = function(object, ...) {
  stop("ks_extract() only accepts 'kernelshap' objects.")
}
