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
#' @param precalc A list with pre-calculated values that are identical for all rows.
#' @return A list with (p x K) matrix of SHAP values, a (p x K) matrix of standard
#'   errors, number of iterations, and convergence status.
permshap_one <- function(
    x,
    v1,
    object,
    pred_fun,
    bg_w,
    v0,
    precalc,
    feature_names,
    exact,
    low_memory,
    tol,
    max_iter,
    ...) {
  bg <- precalc[["bg_X_rep"]]
  X <- rep_rows(x, rep.int(1L, times = nrow(bg)))

  if (exact) {
    Z <- precalc[["Z"]] # ((m_ex+2) x K)
    vz <- get_vz( # (m_ex x K)
      X = X,
      bg = bg,
      Z = Z[2L:(nrow(Z) - 1L), , drop = FALSE], # (m_ex x p)
      object = object,
      pred_fun = pred_fun,
      w = bg_w,
      ...
    )
    vz <- rbind(v0, vz, v1) # we add the cheaply calculated v0 and v1
    rownames(vz) <- precalc[["Z_code"]]
    beta <- shapley_formula(Z, vz = vz)
    out <- list(beta = beta, sigma = 0 * beta, n_iter = 1L, converged = TRUE)
    return(out)
  }

  # Approximate SHAP values
  K <- ncol(v0)
  p <- length(feature_names)

  beta_n <- matrix(
    data = 0, nrow = p, ncol = K, dimnames = list(feature_names, colnames(v0))
  )
  est_m <- list()
  converged <- FALSE
  n_iter <- 0L
  stride <- 2L * (p - 1L)

  while (!converged && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    chains <- balanced_chains(p)
    Z <- lapply(chains, sample_Z_from_chain, feature_names = feature_names)
    if (!low_memory) { # predictions for all chains at once
      Z <- do.call(rbind, Z)
      vz <- get_vz(
        X = X, bg = bg, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
      )
    } else { # predictions for each chain separately
      vz <- vector("list", length = p)
      for (j in seq_len(p)) {
        vz[[j]] <- get_vz(
          X = X,
          bg = bg,
          Z = Z[[j]],
          object = object,
          pred_fun = pred_fun,
          w = bg_w,
          ...
        )
      }
      vz <- do.call(rbind, vz)
    }
    for (j in seq_len(p)) {
      chain <- chains[[j]]
      vzj <- vz[(1L + (j - 1L) * stride):(j * stride), , drop = FALSE]
      vzj <- pad_vz(vzj, v0 = v0, v1 = v1)
      J <- order(chain)
      forward <- vzj[J, , drop = FALSE] - vzj[J + 1L, , drop = FALSE]
      backward <- vzj[p + J + 1L, , drop = FALSE] - vzj[p + J, , drop = FALSE]
      est_m[[length(est_m) + 1L]] <- delta <- (forward + backward) / (2 * p)
      beta_n <- beta_n + delta
    }
    sigma_n <- get_sigma(est_m, iter = n_iter)
    rownames(sigma_n) <- feature_names
    converged <- all(conv_crit(sigma_n, beta_n / n_iter) < tol)
  }
  out <- list(
    beta = beta_n / n_iter, sigma = sigma_n, n_iter = n_iter, converged = converged
  )
  return(out)
}

#' Shapley Weights
#'
#' Weights used in Shapley's formula. Vectorized over `p` and/or `ell`.
#' Required for exact permutation SHAP only.
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

#' Shapley's formula
#'
#' Evaluates Shapley's formula for each feature.
#' Required for exact permutation SHAP only.
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
    L <- rowSums(Z[s1, -j, drop = FALSE]) # how many players are playing with j?
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
#' Required for exact permutation SHAP only.
#'
#' @noRd
#' @keywords internal
#'
#' @param Z A (n x p) matrix.
#' @returns A length n vector.
rowpaste <- function(Z) {
  do.call(paste0, asplit(Z, 2L))
}

#' Balanced Chains
#'
#' Creates `p` permutations of the numbers `1:p` such that each number appears once
#' as the first element.
#' Only used in approximate permutation SHAP.
#'
#' @noRd
#' @keywords internal
#'
#' @param p An integer number of features.
#' @returns A list of length `p`. Each element is a permutation vector of `1:p`.
balanced_chains <- function(p) {
  out <- vector("list", length = p)
  J <- seq_len(p)
  for (j in J) {
    out[[j]] <- c(j, J[-j][sample.int(p - 1L)])
  }
  return(out)
}

#' Z Matrix for iterative permutation SHAP
#'
#' Creates a (2 * (p - 1) x p) on-off-matrix with antithetic rows.
#' The first and last rows (all `TRUE`) and the middle one (all `FALSE`) are skipped
#' (because we already know their values).
#' Only used in approximate permutation SHAP.
#'
#' @noRd
#' @keywords internal
#'
#' @param J A permutation vector of length `p`.
#' @param feature_names A character vector of feature names.
#' @returns A (2 * (p - 1) x p) on-off-matrix with antithetic rows.
sample_Z_from_chain <- function(J, feature_names) {
  m <- length(J) - 1L
  Z <- matrix(TRUE, nrow = m, ncol = m + 1L, dimnames = list(NULL, feature_names))
  for (i in seq_len(m)) {
    Z[i:m, J[i]] <- FALSE
  }
  return(rbind(Z, !Z))
}

#' Fills Gaps in vz
#'
#' Fills the first and last rows of `vz` with `v1`, and the middle one with `v0`.
#' Only used in approximate permutation SHAP.
#'
#' @noRd
#' @keywords internal
#'
#' @param vz A (2 * (p - 1) x K) matrix of vz values.
#' @param v0 A (1 x K) vector of average predictions on background data.
#' @param v1 A (1 x K) vector with the prediction for `x`.
#' @param feature_names A character vector of feature names.
#' @returns A (2 * (p - 1) x K) matrix with vz values.
pad_vz <- function(vz, v0, v1) {
  m <- nrow(vz) %/% 2
  out <- rbind(
    v1, vz[1L:m, , drop = FALSE], v0, vz[(m + 1L):(2L * m), , drop = FALSE], v1
  )
  dimnames(out) <- list(NULL, colnames(v1))
  return(out)
}
