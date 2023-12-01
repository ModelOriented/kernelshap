#' Shapley Weights
#'
#' Weights used in Shapley's formula. Vectorized over `p` and/or `ell`.
#'
#' @noRd
#' @keywords internal
#'
#' @param p Number of features.
#' @param ell Size of subset (i.e., sum of on-off vector z).
#' @returns Shapley weights.
shapley_weights <- function(p, ell) {
  1 / choose(p, ell) / (p - ell)
}

#' SHAP values for one row
#'
#' Calculates permutation SHAP values for a single row.
#'
#' @noRd
#' @keywords internal
#'
#' @inheritParams permshap
#' @param v1 Prediction of `x`.
#' @param v0 Average prediction on background data.
#' @param x A single row to be explained.
#' @param precalc A list with precalculated values that are identical for all rows.
#' @return A (p x K) matrix of SHAP values.
permshap_one <- function(x, v1, object, pred_fun, bg_w, v0, precalc, ...) {
  Z <- precalc[["Z"]]                                                    # ((m_ex+2) x K)
  vz <- get_vz(                                                          # (m_ex x K)
    X = rep_rows(x, rep.int(1L, times = nrow(precalc[["bg_X_rep"]]))),   # (m_ex*n_bg x p)
    bg = precalc[["bg_X_rep"]],                                          # (m_ex*n_bg x p)
    Z = Z[2:(nrow(Z) - 1L), , drop = FALSE],                             # (m_ex x p)
    object = object,
    pred_fun = pred_fun,
    w = bg_w,
    ...
  )
  vz <- rbind(v0, vz, v1)  # we add the cheaply calculated v0 and v1
  rownames(vz) <- precalc[["Z_code"]]
  shapley_formula(Z, vz = vz)
}

#' Shapley's formula
#'
#' Evaluates Shapley's formula for each feature.
#'
#' @noRd
#' @keywords internal
#'
#' @param Z Matrix of on-off row vectors.
#' @param vz Named vector of vz values.
#' @returns SHAP values organized as (p x K) matrix.
shapley_formula <- function(Z, vz) {
  p <- ncol(Z)
  out <- matrix(nrow = p, ncol = ncol(vz), dimnames = list(colnames(Z), colnames(vz)))
  for (j in seq_len(p)) {
    s1 <- Z[, j] == 1L
    vz1 <- vz[s1, , drop = FALSE]
    L <- rowSums(Z[s1, -j, drop = FALSE])  # how many players are playing with j?
    s0 <- rownames(vz1)
    substr(s0, j, j) <- "0"
    vz0 <- vz[s0, , drop = FALSE]
    w <- shapley_weights(p, L)
    out[j, ] <- wcolMeans(vz1 - vz0, w = w)
  }
  out
}

#' Rowwise Paste
#'
#' Fast version of `apply(Z, 1L, FUN = paste0, collapse = "")`.
#'
#' @noRd
#' @keywords internal
#'
#' @param Z A (n x p) matrix.
#' @returns A length n vector.
rowpaste <- function(Z) {
  do.call(paste0, asplit(Z, 2L))
}
