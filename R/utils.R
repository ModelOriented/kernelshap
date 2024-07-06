#' Fast Row Subsetting
#' 
#' Internal function used to row-subset data.frames.
#' Brings a massive speed-up for data.frames. All other classes (tibble, data.table,
#' matrix) are subsetted in the usual way.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param x A matrix-like object.
#' @param i Logical or integer vector of rows to pick.
#' @returns Subsetted version of `x`.
rep_rows <- function(x, i) {
  if (!(all(class(x) == "data.frame"))) {
    return(x[i, , drop = FALSE])  # matrix, tibble, data.table, ...
  }
  # data.frame
  out <- lapply(x, function(z) if (length(dim(z)) != 2L) z[i] else z[i, , drop = FALSE])
  attr(out, "row.names") <- .set_row_names(length(i))
  class(out) <- "data.frame"
  out
}

#' Weighted Version of colMeans()
#' 
#' Internal function used to calculate column-wise weighted means.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param x A matrix-like object.
#' @param w Optional case weights.
#' @returns A (1 x ncol(x)) matrix of column means.
wcolMeans <- function(x, w = NULL, ...) {
  x <- as.matrix(x)
  out <- if (is.null(w)) colMeans(x) else colSums(x * w) / sum(w)
  t.default(out)
}

#' All on-off Vectors
#'
#' Internal function that creates matrix of all on-off vectors of length `p`.
#'
#' @noRd
#' @keywords internal
#'
#' @param p Number of features.
#' @param feature_names Feature names.
#' @param keep_extremes Should extremes be kept? Defaults to `FALSE` (for kernelshap).
#' @returns An integer matrix of all on-off vectors of length `p`.
exact_Z <- function(p, feature_names, keep_extremes = FALSE) {
  Z <- as.matrix(do.call(expand.grid, replicate(p, 0:1, simplify = FALSE)))
  colnames(Z) <- feature_names
  if (keep_extremes) Z else Z[2:(nrow(Z) - 1L), , drop = FALSE]
}

#' Masker
#'
#' Internal function. 
#' For each on-off vector (rows in `Z`), the (weighted) average prediction is returned.
#'
#' @noRd
#' @keywords internal
#'
#' @inheritParams kernelshap
#' @param X Row to be explained stacked m*n_bg times.
#' @param bg Background data stacked m times.
#' @param Z A (m x p) matrix with on-off values.
#' @param w A vector with case weights (of the same length as the unstacked
#'   background data).
#' @returns A (m x K) matrix with vz values.
get_vz <- function(X, bg, Z, object, pred_fun, w, ...) {
  m <- nrow(Z)
  not_Z <- !Z
  n_bg <- nrow(bg) / m   # because bg was replicated m times
  
  # Replicate not_Z, so that X, bg, not_Z are all of dimension (m*n_bg x p)
  g <- rep_each(m, each = n_bg)
  not_Z <- not_Z[g, , drop = FALSE]
  
  if (is.matrix(X)) {
    # Remember that columns of X and bg are perfectly aligned in this case
    X[not_Z] <- bg[not_Z]
  } else {
    for (v in colnames(Z)) {
      s <- not_Z[, v]
      X[[v]][s] <- bg[[v]][s]
    }
  }
  preds <- align_pred(pred_fun(object, X, ...))
  
  # Aggregate (distinguishing fast 1-dim case)
  if (ncol(preds) == 1L) {
    return(wrowmean_vector(preds, ngroups = m, w = w))
  }
  if (is.null(w)) {
    return(rowsum(preds, group = g, reorder = FALSE) / n_bg)
  }
  rowsum(preds * w, group = g, reorder = FALSE) / sum(w)
}

#' Combine Matrices
#'
#' Binds list of matrices along new first axis.
#'
#' @noRd
#' @keywords internal
#'
#' @param a List of n (p x K) matrices.
#' @returns A (n x p x K) array.
abind1 <- function(a) {
  out <- array(
    dim = c(length(a), dim(a[[1L]])), 
    dimnames = c(list(NULL), dimnames(a[[1L]]))
  )
  for (i in seq_along(a)) {
    out[i, , ] <- a[[i]]
  }
  out
}

#' Reorganize List
#'
#' Internal function that turns list of n (p x K) matrices into list of K (n x p) 
#' matrices. Reduce if K = 1.
#'
#' @noRd
#' @keywords internal
#'
#' @param alist List of n (p x K) matrices.
#' @returns List of K (n x p) matrices.
reorganize_list <- function(alist) {
  if (!is.list(alist)) {
    stop("alist must be a list")
  }
  out <- asplit(abind1(alist), MARGIN = 3L)
  if (length(out) == 1L) {
    return(as.matrix(out[[1L]]))
  }
  lapply(out, as.matrix)
}

#' Aligns Predictions
#'
#' Turns predictions into matrix.
#'
#' @noRd
#' @keywords internal
#'
#' @param x Object representing model predictions.
#' @returns Like `x`, but converted to matrix.
align_pred <- function(x) {
  if (is.data.frame(x) && ncol(x) == 1L) {
    x <- x[[1L]]
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("Predictions must be numeric!")
  }
  return(x)
}

#' Head of List Elements
#'
#' Internal function that returns the top n rows of each element in the input list.
#'
#' @noRd
#' @keywords internal
#'
#' @param x A list or a matrix-like.
#' @param n Number of rows to show.
#' @returns List of first rows of each element in the input.
head_list <- function(x, n = 6L) {
  if (!is.list(x)) utils::head(x, n) else lapply(x, utils::head, n)
}

# Summarize details about the chosen algorithm (exact, hybrid, sampling)
summarize_strategy <- function(p, exact, deg) {
  if (exact || trunc(p / 2) == deg) {
    txt <- "Exact Kernel SHAP values"
    if (!exact) {
      txt <- paste(txt, "by the hybrid approach")
    }
    return(txt)
  }
  if (deg == 0L) {
    return("Kernel SHAP values by iterative sampling")
  } 
  paste("Kernel SHAP values by the hybrid strategy of degree", deg)
}

# Case p = 1 returns exact Shapley values
case_p1 <- function(n, feature_names, v0, v1, X, verbose) {
  txt <- "Exact Shapley values (p = 1)"
  if (verbose) {
    message(txt)
  }
  S <- v1 - v0[rep(1L, n), , drop = FALSE]                        #  (n x K)
  SE <- matrix(numeric(n), dimnames = list(NULL, feature_names))  #  (n x 1)
  if (ncol(v1) > 1L) {
    SE <- replicate(ncol(v1), SE, simplify = FALSE)
    S <- lapply(
      asplit(S, MARGIN = 2L), function(M) 
        as.matrix(M, dimnames = list(NULL, feature_names))
    )
  } else {
    colnames(S) <- feature_names      
  }
  out <- list(
    S = S, 
    X = X, 
    baseline = as.vector(v0), 
    SE = SE, 
    n_iter = integer(n), 
    converged = rep(TRUE, n),
    m = 0L,
    m_exact = 0L,
    prop_exact = 1,
    exact = TRUE,
    txt = txt,
    predictions = v1,
    algorithm = "kernelshap"
  )
  class(out) <- "kernelshap"
  out
}

#' Fast Index Generation (from {hstats})
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

#' Grouped Means for Single-Column Matrices (adapted from {hstats})
#'
#' Grouped means for matrix with single column over fixed-length groups.
#'
#' @noRd
#' @keywords internal
#'
#' @param x Matrix with one column.
#' @param ngroups Number of subsequent, equals sized groups.
#' @param w Optional vector of case weights of length `NROW(x) / ngroups`.
#' @returns Matrix with one column.
wrowmean_vector <- function(x, ngroups = 1L, w = NULL) {
  if (ncol(x) != 1L) {
    stop("x must have a single column")
  }
  nm <- colnames(x)
  dim(x) <- c(length(x) %/% ngroups, ngroups)
  out <- if (is.null(w)) colMeans(x) else colSums(x * w) / sum(w)
  dim(out) <- c(ngroups, 1L)
  if (!is.null(nm)) {
    colnames(out) <- nm
  }
  out
}

#' Basic Input Checks
#' 
#' @noRd
#' @keywords internal
#' 
#' @inheritParams kernelshap
#' 
#' @returns TRUE or an error
basic_checks <- function(X, bg_X, feature_names, pred_fun) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.matrix(X) == is.matrix(bg_X),
    dim(X) >= 1L,
    dim(bg_X) >= 1L,
    !is.null(colnames(X)),
    !is.null(colnames(bg_X)),
    length(feature_names) >= 1L,
    all(feature_names %in% colnames(X)),
    all(feature_names %in% colnames(bg_X)),  # not necessary, but clearer
    all(colnames(X) %in% colnames(bg_X)),
    is.function(pred_fun),
    "If X is a matrix, feature_names must equal colnames(X)" = 
      !is.matrix(X) || identical(colnames(X), feature_names)
  )
  TRUE
}

#' Warning on Slow Computations
#' 
#' @noRd
#' @keywords internal
#' 
#' @param m Number of on-off vectors.
#' @param bg_n Number of rows in the background data.
#' 
#' @returns TRUE.
warning_burden <- function(m, bg_n) {
  warning("\nPredictions on large data sets with ", m, "x", bg_n,
          " observations are being done.\n",
          "Consider reducing the computational burden (e.g. use smaller X_bg)")
  TRUE
}

#' Prepare Case Weights
#' 
#' @noRd
#' @keywords internal
#' 
#' @param w Vector of case weights.
#' @param bg_n Number of rows in the background data.
#' 
#' @returns TRUE or an error
prep_w <- function(w, bg_n) {
  stopifnot(
    length(w) == bg_n, 
    all(w >= 0), 
    !all(w == 0)
  )
  if (!is.double(w)) as.double(w) else w
}

