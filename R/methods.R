#' Print Method
#' 
#' Prints the first two rows of the matrix (or matrices) of SHAP values. 
#'
#' @param x An object of class "kernelshap".
#' @param n Maximum number of rows of SHAP values to print.
#' @param ... Further arguments passed from other methods.
#' @return Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:3, -1], bg_X = iris[-1])
#' s
#' @seealso \code{\link{kernelshap}}.
print.kernelshap <- function(x, n = 2L, ...) {
  cat("SHAP values of first", n, "observations:\n")
  head_list(ks_extract(x, "S"), n = n)
  invisible(x)
}

#' Summary Method
#'
#' @param object An object of class "kernelshap".
#' @param compact Set to \code{TRUE} to hide printing the top n SHAP values, 
#' standard errors and feature values. 
#' @param n Maximum number of rows of SHAP values, standard errors and feature values to print.
#' @param ... Further arguments passed from other methods.
#' @return Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:3, -1], bg_X = iris[-1])
#' summary(s)
#' @seealso \code{\link{kernelshap}}.
summary.kernelshap <- function(object, compact = FALSE, n = 2L, ...) {
  cat(ks_extract(object, "txt"))

  S <- ks_extract(object, "S")
  if (!is.list(S)) {
    n <- min(n, nrow(S))
    cat(paste("\n  - SHAP matrix of dim",  nrow(S), "x", ncol(S)))
  } else {
    n <- min(n, nrow(S[[1L]]))
    cat(
      "\n  -", length(S), "SHAP matrices of dim", nrow(S[[1L]]), "x", ncol(S[[1L]])
    )
  }
  cat("\n  - baseline:", ks_extract(object, "baseline"))
  ex <- ks_extract(object, "exact")
  if (!ex) {
    cat(
      "\n  - average number of iterations:", mean(ks_extract(object, "n_iter")),
      "\n  - rows not converged:", sum(!ks_extract(object, "converged")),
      "\n  - proportion exact:", ks_extract(object, "prop_exact"),
      "\n  - m/iter:", ks_extract(object, "m")
    )
  }
  cat("\n  - m_exact:", ks_extract(object, "m_exact"))
  if (!compact) {
    cat("\n\nSHAP values of first", n, "observations:\n")
    head_list(S, n = n)
    if (!ex) {
      cat("\nCorresponding standard errors:\n")
      head_list(ks_extract(object, "SE"), n = n) 
    }
  }
  invisible(object)
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
#' "converged", "m", "m_exact", "prop_exact", "exact", or "txt".
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
                                          "m", "m_exact", "prop_exact", "exact", "txt"), 
                                 ...) {
  what <- match.arg(what)
  object[[what]]
}

#' @describeIn ks_extract No default method available.
#' @export
ks_extract.default = function(object, ...) {
  stop("ks_extract() only accepts 'kernelshap' objects.")
}
