sample_Z_from_chain <- function(J, feature_names) {
  m <- length(J) - 1L
  Z <- matrix(1L, nrow = m, ncol = m + 1L, dimnames = list(NULL, feature_names))
  for (i in seq_len(m)) {
    Z[i:m, J[i]] <- 0L
  }
  return(rbind(Z, 1L - Z))
}

pad_vz <- function(vz, v0, v1) {
  m <- nrow(vz) / 2
  out <- rbind(
    v1, vz[1:m, , drop = FALSE], v0, vz[(m + 1):(2 * m), , drop = FALSE], v1
  )
  dimnames(out) <- list(NULL, colnames(v1))
  return(out)
}

object <- lm(Sepal.Length ~ ., data = iris)
bg_X <- iris[, -1]
bg_w <- NULL
x <- iris[1, -1, drop = F]

p <- ncol(x)
bg_n <- nrow(bg_X)
v1 <- predict(object, x)
v0 <- wcolMeans(predict(object, bg_X), w = w)
pred_fun <- predict
... <- NULL

precalc <- list(
  bg_X_rep = rep_rows(bg_X, rep.int(seq_len(bg_n), 2L * p - 2L))
)

chains_per_iter <- 100L / p
tol <- 0.005
max_iter <- 100L


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
    x, v1, object, pred_fun, bg_w, v0, precalc, chains_per_iter, tol, max_iter, ...) {
  K <- ncol(v0)
  p <- ncol(x)
  feature_names <- colnames(x)
  bg <- precalc[["bg_X_rep"]]
  X <- rep_rows(x, rep.int(1L, times = nrow(bg)))

  J <- seq_len(p)
  beta_n <- matrix(
    data = 0, nrow = p, ncol = K, dimnames = list(feature_names, colnames(v0))
  )
  est_m <- list()
  converged <- FALSE
  n_iter <- 0L

  while (!converged && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    chain <- sample(p)
    Z <- sample_Z_from_chain(chain, feature_names = feature_names)
    vz <- get_vz(
      X = X, bg = bg, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
    )
    vz <- pad_vz(vz, v0 = v0, v1 = v1)
    forward <- vz[J, , drop = FALSE] - vz[J + 1L, , drop = FALSE]
    backward <- vz[p + J + 1L, , drop = FALSE] - vz[p + J, , drop = FALSE]
    est_m[[n_iter]] <- (forward + backward)[order(chain), , drop = FALSE] / 2
    beta_n <- beta_n + est_m[[n_iter]]

    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      sigma_n <- get_sigma(est_m, iter = n_iter)
      converged <- all(conv_crit(sigma_n, beta_n / n_iter) < tol)
    }
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

# Sepal.Width Petal.Length Petal.Width   Species
#   0.2195135    -1.955357   0.3149451 0.5823533


library(kernelshap)
library(ranger)
set.seed(1)
fit <- ranger(Sepal.Length ~ ., data = iris)
ks <- kernelshap(fit, X = iris[1:10, -1], bg_X = iris[, -1], exact = F, hybrid_degree = 0)
