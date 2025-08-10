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
#' @return A list with (p x K) matrix of SHAP values, and for the sampling version
#'   a (p x K) matrix of standard errors, number of iterations,
#'   and convergence status.
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
    bg_n,
    ...) {
  p <- length(feature_names)
  K <- ncol(v1)
  K_names <- colnames(v1)
  beta_n <- matrix(0, nrow = p, ncol = K, dimnames = list(feature_names, K_names))

  if (exact) {
    vz <- get_vz(
      x = x,
      bg_rep = precalc$bg_exact_rep,
      Z_rep = precalc$Z_exact_rep,
      object = object,
      pred_fun = pred_fun,
      w = bg_w,
      bg_n = bg_n,
      ...
    )
    vz <- rbind(v0, vz, v1)

    for (j in seq_len(p)) {
      pos <- precalc$positions[[j]]
      beta_n[j, ] <- wcolMeans(
        vz[pos$on, , drop = FALSE] - vz[pos$off, , drop = FALSE],
        w = precalc$shapley_w[pos$on]
      )
    }
    return(list(beta = beta_n))
  }

  # Sampling version
  est_m <- array(
    data = 0, dim = c(max_iter, p, K), dimnames = list(NULL, feature_names, K_names)
  )
  converged <- FALSE
  n_iter <- 0L
  stride <- 2L * (p - 3L)

  # Pre-calculate part of Z with rowsum 1 or p - 1
  vz_balanced <- get_vz( # (2p x K)
    x = x,
    bg_rep = precalc$bg_balanced_rep,
    Z_rep = precalc$Z_balanced_rep,
    object = object,
    pred_fun = pred_fun,
    w = bg_w,
    bg_n = bg_n,
    ...
  )

  # vzj has constant first, middle, and last row
  vzj <- init_vzj(p, v0 = v0, v1 = v1)

  # Important positions to be filled in vzj
  from_balanced <- c(2L, 2L + p, p, 2L * p)
  from_iter <- c(3L:(p - 1L), (p + 3L):(2L * p - 1L))

  bg_sampling_rep <- precalc$bg_sampling_rep
  g <- rep_each(nrow(bg_sampling_rep) %/% bg_n, each = bg_n)

  while (!converged && n_iter < max_iter) {
    chains <- balanced_chains(p)
    Z <- lapply(chains, sample_Z_from_chain, feature_names = feature_names)
    if (!low_memory) { # predictions for all chains at once
      Z <- do.call(rbind, Z)
      vz <- get_vz(
        x = x,
        bg_rep = bg_sampling_rep,
        Z_rep = Z[g, , drop = FALSE],
        object = object,
        pred_fun = pred_fun,
        w = bg_w,
        bg_n = bg_n,
        ...
      )
    } else { # predictions for each chain separately
      vz <- vector("list", length = p)
      for (j in seq_len(p)) {
        vz[[j]] <- get_vz(
          x = x,
          bg_rep = bg_sampling_rep,
          Z_rep = Z[[j]][g, , drop = FALSE],
          object = object,
          pred_fun = pred_fun,
          w = bg_w,
          bg_n = bg_n,
          ...
        )
      }
      vz <- do.call(rbind, vz)
    }

    for (j in seq_len(p)) {
      n_iter <- n_iter + 1L
      chain <- chains[[j]]

      # Fill vzj by pre-calculated masks
      vzj[from_balanced, ] <- vz_balanced[c(j, j + p, chain[p] + p, chain[p]), ]
      vzj[from_iter, ] <- vz[(1L + (j - 1L) * stride):(j * stride), , drop = FALSE]

      J <- order(chain)
      forward <- vzj[J, , drop = FALSE] - vzj[J + 1L, , drop = FALSE]
      backward <- vzj[p + J + 1L, , drop = FALSE] - vzj[p + J, , drop = FALSE]
      est_m[n_iter, , ] <- delta <- (forward + backward) / 2
      beta_n <- beta_n + delta
    }
    sigma_n <- get_sigma(est_m[seq_len(n_iter), , , drop = FALSE])
    converged <- check_convergence(beta = beta_n / n_iter, sigma = sigma_n, tol = tol)
  }
  rownames(sigma_n) <- feature_names
  out <- list(
    beta = beta_n / n_iter, sigma = sigma_n, n_iter = n_iter, converged = converged
  )
  return(out)
}

#' Shapley Weights
#'
#' Weights used in Shapley's formula. Vectorized over `p` and/or `ell`.
#' Required for exact permutation SHAP only. Note that the result is `Inf` when
#' ell < 0.
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

#' Balanced Chains
#'
#' Creates `p` permutations of the numbers `1:p` such that each number appears once
#' as the first element.
#' Only used in iterative permutation SHAP.
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
#' Creates a (2 * (p - 2) x p) on-off-matrix with antithetic scheme for the backward
#' pass. The first and last rows (all `TRUE`) and the middle one (all `FALSE`) are dropped
#' (because we already know their values). Also rows with rowsum 1 or p-1
#' are dropped. Only used in iterative permutation SHAP.
#'
#' @noRd
#' @keywords internal
#'
#' @param J A permutation vector of length `p`.
#' @param feature_names A character vector of feature names.
#' @returns A (2 * (p - 2) x p) on-off-matrix.
sample_Z_from_chain <- function(J, feature_names) {
  m <- length(J) - 1L
  if (m < 3L) {
    stop("J must have at least 3 elements")
  }
  Z <- matrix(TRUE, nrow = m, ncol = m + 1L, dimnames = list(NULL, feature_names))
  for (i in seq_len(m)) {
    Z[i:m, J[i]] <- FALSE
  }
  # Drop rows with rowsum 1 or p-1 (could be done more efficiently in advance
  Z <- Z[2L:(m - 1L), , drop = FALSE]
  return(rbind(Z, !Z))
}

#' Exact Z part for p balanced permutations
#'
#' Creates a (2p x p) on-off-matrix with the starting on-off vector for the forward and
#' backward pass. Rows j and j + p refer to the same starting feature.
#'
#' @noRd
#' @keywords internal
#'
#' @param p An integer number of features.
#' @param feature_names A character vector of feature names.
#' @returns A (2p x p) on-off-matrix with antithetic rows.
exact_Z_balanced <- function(p, feature_names) {
  Z <- diag(p)
  storage.mode(Z) <- "logical"
  colnames(Z) <- feature_names
  return(rbind(!Z, Z))
}

#' Exact Z part for p balanced permutations
#'
#' Creates a (2p x p) on-off-matrix with the starting on-off vector for the forward and
#' backward pass. Rows j and j + p refer to the same starting feature.
#'
#' @noRd
#' @keywords internal
#'
#' @param p An integer number of features.
#' @param v0 A (1 x K) matrix with the average prediction on background data.
#' @param v1 A (1 x K) matrix with the prediction on the row to be explained.
#' @returns A ((2p + 1) x K) matrix with initialized first, middle and last row.
init_vzj <- function(p, v0, v1) {
  vzj <- matrix(0, nrow = 2L * p + 1L, ncol = ncol(v1), dimnames = list(NULL, colnames(v1)))
  vzj[2L * p + 1L, ] <- vzj[1L, ] <- v1
  vzj[p + 1L, ] <- v0
  return(vzj)
}


#' On-Off Positions
#'
#' Creates a list with on- and off-positions of vz for each feature. This can
#' be calculated once for all rows to be explained.
#'
#' @noRd
#' @keywords internal
#'
#' @param mask Logical matrix. The output of `exact_Z(p, feature_names)`.
#' @returns List with p elements, each containing an `on` and `off` vector.
positions_for_exact <- function(mask) {
  p <- ncol(mask)
  codes <- seq_len(nrow(mask)) # Row index = binary code of the row

  positions <- vector("list", p)
  for (j in seq_len(p)) {
    on <- codes[mask[, j]]
    off <- on - 2^(p - j) # trick to turn "bit" off
    positions[[j]] <- list(on = on, off = off)
  }
  return(positions)
}
