# TODO
# 1. run multiple chains together
# 2. balanced permutations: generate max_iter * chains_per_iter chains
#  so that the first value is balanced within chains_per_iter. Then,
# produce corresponding list of Z. Theoretically, we can do this before entering
# the main function, so that all observations use the same chains.


balanced_chains <- function(p) {
  out <- vector("list", length = p)
  J <- seq_len(p)
  for (j in J) {
    out[[j]] <- c(j, J[-j][sample.int(p - 1L)])
  }
  return(out)
}

sample_Z_from_chain <- function(J, feature_names) {
  m <- length(J) - 1L
  Z <- matrix(1L, nrow = m, ncol = m + 1L, dimnames = list(NULL, feature_names))
  for (i in seq_len(m)) {
    Z[i:m, J[i]] <- 0L
  }
  return(rbind(Z, 1L - Z))
}

pad_vz <- function(vz, v0, v1) {
  m <- nrow(vz) %/% 2
  out <- rbind(
    v1, vz[1L:m, , drop = FALSE], v0, vz[(m + 1L):(2L * m), , drop = FALSE], v1
  )
  dimnames(out) <- list(NULL, colnames(v1))
  return(out)
}

library(ranger)
set.seed(1)
object <- lm(Sepal.Length ~ ., data = iris)
pred_fun <- predict

object <- ranger(Sepal.Length ~ ., data = iris)
pred_fun <- function(m, x) predict(m, x)$predictions

... <- NULL

bg_X <- iris[, -1]
bg_w <- NULL
x <- iris[1, -1, drop = F]

p <- ncol(x)
bg_n <- nrow(bg_X)
v1 <- pred_fun(object, x, ...)
v0 <- wcolMeans(pred_fun(object, bg_X, ...), w = bg_w)

tol <- 0.005
max_iter <- 10L
low_memory <- FALSE

bg_rep <- 2L * (p - 1L) * (if (low_memory) 1L else p)
precalc <- list(
  bg_X_rep = rep_rows(bg_X, rep.int(seq_len(bg_n), bg_rep))
)



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
#' @param chains_per_iter Number of chains to sample.
#' @return A (p x K) matrix of SHAP values.
permshap_sampling_one <- function(
    x, v1, object, pred_fun, bg_w, v0, precalc, low_memory, tol, max_iter, ...) {
  K <- ncol(v0)
  p <- ncol(x)
  feature_names <- colnames(x)
  bg <- precalc[["bg_X_rep"]]
  X <- rep_rows(x, rep.int(1L, times = nrow(bg)))

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
    converged <- all(conv_crit(sigma_n, beta_n / n_iter) < tol)
  }
  out <- list(
    beta = beta_n / n_iter,
    sigma = sigma_n,
    n_iter = n_iter,
    converged = converged
  )
  return(out)
}

permshap_sampling_one(x, v1, object, predict, w, v0, precalc, chains_per_iter = 100)


library(kernelshap)

permshap(object, X = iris[1, -1], bg_X = iris[, -1])

# Sepal.Width Petal.Length Petal.Width   Species
#   0.2195135    -1.955357   0.3149451 0.5823533

# Sepal.Width Petal.Length Petal.Width    Species
#   0.1734389   -0.5202173  -0.2416169 -0.1420209
