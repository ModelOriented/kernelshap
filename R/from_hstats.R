# These functions have originally been implemented in {hstats}

#' Fast Index Generation
#' 
#' For not too small m, much faster than `rep(seq_len(m), each = each)`.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param m Integer. See `each`.
#' @param each Integer. How many times should each value in `1:m` be repeated?
#' @returns Like `x`, but converted to matrix.
#' @examples
#' rep_each(10, 2)
#' rep(1:10, each = 2)  # Dito
rep_each <- function(m, each) {
  out <- .col(dim = c(each, m))
  dim(out) <- NULL
  out 
}

#' Fast OHE
#' 
#' Turns vector/factor into its One-Hot-Encoding.
#' Ingeniouly written by Mathias Ambuehl.
#' 
#' Working with integers instead of doubles would be slightly faster, but at the price
#' of potential integer overflows in subsequent calculations.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param x Object representing model predictions.
#' @returns Like `x`, but converted to matrix.
fdummy <- function(x) {
  x <- as.factor(x)
  lev <- levels(x)
  out <- matrix(0, nrow = length(x), ncol = length(lev))
  out[cbind(seq_along(x), as.integer(x))] <- 1
  colnames(out) <- lev
  out 
}
