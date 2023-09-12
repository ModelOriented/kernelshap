# Functions required only for handling (partly) exact cases

# Provides fixed input for the exact case:
# - Z: Matrix with all 2^p-2 on-off vectors z
# - w: Vector with row weights of Z ensuring that the distribution of sum(z) matches
#      the SHAP kernel distribution
# - A: Exact matrix A = Z'wZ
input_exact <- function(p, feature_names) {
  Z <- exact_Z(p, feature_names = feature_names)
  # Each Kernel weight(j) is divided by the number of vectors z having sum(z) = j
  w <- kernel_weights(p) / choose(p, 1:(p - 1L))
  list(Z = Z, w = w[rowSums(Z)], A = exact_A(p, feature_names = feature_names))
}

#' Exact Matrix A
#'
#' Internal function that calculates exact A. 
#' Notice the difference to the off-diagnonals in the Supplement of 
#' Covert and Lee (2021). Credits to David Watson for figuring out the correct formula,
#' see our discussions in https://github.com/ModelOriented/kernelshap/issues/22
#'
#' @noRd
#' @keywords internal
#'
#' @param p Number of features.
#' @param feature_names Feature names.
#' @returns A (p x p) matrix.
exact_A <- function(p, feature_names) {
  S <- 1:(p - 1L)
  c_pr <- S * (S - 1) / p / (p - 1)
  off_diag <- sum(kernel_weights(p) * c_pr)
  A <- matrix(
    off_diag, nrow = p, ncol = p, dimnames = list(feature_names, feature_names)
  )
  diag(A) <- 0.5
  A
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
#' @returns An integer ((2^p - 2) x p) matrix of all on-off vectors of length `p`.
exact_Z <- function(p, feature_names) {
  Z <- as.matrix(do.call(expand.grid, replicate(p, 0:1, simplify = FALSE)))
  colnames(Z) <- feature_names
  Z[2:(nrow(Z) - 1L), , drop = FALSE]
}

# List all length p vectors z with sum(z) in {k, p - k}
partly_exact_Z <- function(p, k, feature_names) {
  if (k < 1L) {
    stop("k must be at least 1")
  }
  if (p < 2L * k) {
    stop("p must be >=2*k")
  }
  if (k == 1L) {
    Z <- diag(p)
  } else {
    Z <- t(
      utils::combn(seq_len(p), k, FUN = function(z) {x <- numeric(p); x[z] <- 1; x})
    )
  }
  if (p != 2L * k) {
    Z <- rbind(Z, 1 - Z)
  }
  colnames(Z) <- feature_names
  Z
}

# Create Z, w, A for vectors z with sum(z) in {k, p-k} for k in {1, ..., deg}.
# The total weights do not sum to one, except in the special (exact) case deg=p-deg.
# (The remaining weight will be added via input_sampling(p, deg=deg)).
# Note that for a given k, the weights are constant.
input_partly_exact <- function(p, deg, feature_names) {
  if (deg < 1L) {
    stop("deg must be at least 1")
  }
  if (p < 2L * deg) {
    stop("p must be >=2*deg")
  }
  
  kw <- kernel_weights(p)
  Z <- w <- vector("list", deg)
  
  for (k in seq_len(deg)) {
    Z[[k]] <- partly_exact_Z(p, k = k, feature_names = feature_names)
    n <- nrow(Z[[k]])
    w_tot <- kw[k] * (2 - (p == 2L * k))
    w[[k]] <- rep(w_tot / n, n)
  }
  w <- unlist(w, recursive = FALSE, use.names = FALSE)
  Z <- do.call(rbind, Z)
  
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
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
    predictions = v1
  )
  class(out) <- "kernelshap"
  out
}