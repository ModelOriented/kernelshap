#' Print Method
#' 
#' Prints the first two rows of the matrix (or matrices) of SHAP values. 
#'
#' @param x An object of class "kernelshap".
#' @param n Maximum number of rows of SHAP values to print.
#' @param ... Further arguments passed from other methods.
#' @returns Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:3, -1], bg_X = iris[-1])
#' s
#' @seealso [kernelshap()]
print.kernelshap <- function(x, n = 2L, ...) {
  cat("SHAP values of first observations:\n")
  print(head_list(getElement(x, "S"), n = n))
  invisible(x)
}

#' @describeIn print.kernelshap Print method for "permshap" object
#' @export
print.permshap <- function(x, n = 2L, ...) {
  print.kernelshap(x, n = n, ...)
}

#' Summary Method
#'
#' @param object An object of class "kernelshap".
#' @param compact Set to `TRUE` to hide printing the top n SHAP values, 
#'   standard errors and feature values. 
#' @param n Maximum number of rows of SHAP values, standard errors and feature values 
#'   to print.
#' @param ... Further arguments passed from other methods.
#' @returns Invisibly, the input is returned.
#' @export
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:3, -1], bg_X = iris[-1])
#' summary(s)
#' @seealso [kernelshap()]
summary.kernelshap <- function(object, compact = FALSE, n = 2L, ...) {
  cat(getElement(object, "txt"))

  S <- getElement(object, "S")
  if (!is.list(S)) {
    n <- min(n, nrow(S))
    cat(paste("\n  - SHAP matrix of dim",  nrow(S), "x", ncol(S)))
  } else {
    n <- min(n, nrow(S[[1L]]))
    cat(
      "\n  -", length(S), "SHAP matrices of dim", nrow(S[[1L]]), "x", ncol(S[[1L]])
    )
  }
  cat("\n  - baseline:", getElement(object, "baseline"))
  ex <- getElement(object, "exact")
  if (!ex) {
    cat(
      "\n  - average number of iterations:", mean(getElement(object, "n_iter")),
      "\n  - rows not converged:", sum(!getElement(object, "converged")),
      "\n  - proportion exact:", getElement(object, "prop_exact"),
      "\n  - m/iter:", getElement(object, "m")
    )
  }
  cat("\n  - m_exact:", getElement(object, "m_exact"))
  if (!compact) {
    cat("\n\nSHAP values of first observations:\n")
    print(head_list(S, n = n))
    if (!ex) {
      cat("\nCorresponding standard errors:\n")
      print(head_list(getElement(object, "SE"), n = n))
    }
  }
  invisible(object)
}

#' Check for kernelshap
#'
#' Is object of class "kernelshap"?
#'
#' @param object An R object.
#' @returns `TRUE` if `object` is of class "kernelshap", and `FALSE` otherwise.
#' @export
#' @examples
#' fit <- lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:2, -1], bg_X = iris[-1])
#' is.kernelshap(s)
#' is.kernelshap("a")
#' @seealso [kernelshap()]
is.kernelshap <- function(object){
  inherits(object, "kernelshap")
}

#' Check for permshap
#'
#' Is object of class "permshap"?
#'
#' @param object An R object.
#' @returns `TRUE` if `object` is of class "permshap", and `FALSE` otherwise.
#' @export
#' @examples
#' fit <- lm(Sepal.Length ~ ., data = iris)
#' s <- permshap(fit, iris[1:2, -1], bg_X = iris[-1])
#' is.permshap(s)
#' is.permshap("a")
#' @seealso [kernelshap()]
is.permshap <- function(object){
  inherits(object, "permshap")
}
