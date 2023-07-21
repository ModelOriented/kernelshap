# Functions required only for handling exact cases

# Provides fixed input for the exact case:
# - Z: Matrix with all 2^p-2 on-off vectors z
# - w: Vector with row weights of Z ensuring that the distribution of sum(z) matches
#      the SHAP kernel distribution
# - A: Exact matrix A = Z'wZ
input_exact <- function(p) {
  Z <- exact_Z(p)
  # Each Kernel weight(j) is divided by the number of vectors z having sum(z) = j
  w <- kernel_weights(p) / choose(p, 1:(p - 1L))
  list(Z = Z, w = w[rowSums(Z)])
}

# Creates (2^p-2) x p matrix with all on-off vectors z of length p
# Instead of calculating this object, we could evaluate it for different p <= p_max
# and store it as a list in the package.
exact_Z <- function(p) {
  Z <- as.matrix(do.call(expand.grid, replicate(p, 0:1, simplify = FALSE)))
  dimnames(Z) <- NULL
  Z[2:(nrow(Z) - 1L), , drop = FALSE]
}

# List all length p vectors z with sum(z) in {k, p - k}
partly_exact_Z <- function(p, k) {
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
  if (p == 2L * k) {
    return(Z)
  }
  return(rbind(Z, 1 - Z))
}

# Create Z, w, A for vectors z with sum(z) in {k, p-k} for k in {1, ..., deg}.
# The total weights do not sum to one, except in the special (exact) case deg=p-deg.
# (The remaining weight will be added via input_sampling(p, deg=deg)).
# Note that for a given k, the weights are constant.
input_partly_exact <- function(p, deg) {
  if (deg < 1L) {
    stop("deg must be at least 1")
  }
  if (p < 2L * deg) {
    stop("p must be >=2*deg")
  }

  kw <- kernel_weights(p)
  Z <- w <- vector("list", deg)

  for (k in seq_len(deg)) {
    Z[[k]] <- partly_exact_Z(p, k = k)
    n <- nrow(Z[[k]])
    w_tot <- kw[k] * (2 - (p == 2L * k))
    w[[k]] <- rep(w_tot / n, n)
  }
  w <- unlist(w, recursive = FALSE, use.names = FALSE)
  Z <- do.call(rbind, Z)

  list(Z = Z, w = w)
}
