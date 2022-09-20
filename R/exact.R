# Functions required only for handling exact cases

# Provides fixed input for the exact case:
# - Z: Matrix with all 2^p-2 on-off vectors z
# - w: Vector with row weights of Z ensuring that the distribution of sum(z) matches
#      the SHAP kernel distribution
# - A: Exact matrix A = Z'wZ
input_exact <- function(p) {
  Z <- exact_Z(p)
  w <- kernel_weights(p) / choose(p, 1:(p - 1L))
  list(Z = Z, w = w[rowSums(Z)], A = exact_A(p))
}

# Calculates exact A. Notice the difference in the off-diagnonals to Appendix 1 of 
# Covert and Lee (2021), which seems slightly off
exact_A <- function(p) {
  S <- 1:(p - 1L)
  c_pr <- S * (S - 1) / p / (p - 1)
  off_diag <- sum(kernel_weights(p) * c_pr)
  A <- matrix(off_diag, nrow = p, ncol = p)
  diag(A) <- 0.5
  A
}

# Creates (2^p-2) x p matrix with all on-off vectors z of length p
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

# Create Z, w, A for vectors with sum(z) in {k, p-k} for k in deg
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
  
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}

# Case p = 1 returns exact Shapley values
case_p1 <- function(n, nms, v0, v1, X, verbose) {
  txt <- "Exact Shapley values (p = 1)"
  if (verbose) {
    message(txt)
  }
  S <- v1 - v0[rep(1L, n), , drop = FALSE]
  SE <- matrix(numeric(n), dimnames = list(NULL, nms))
  if (ncol(v1) > 1L) {
    SE <- replicate(ncol(v1), SE, simplify = FALSE)
    S <- lapply(
      asplit(S, MARGIN = 2L), function(M) as.matrix(M, dimnames = list(NULL, nms))
    )
  } else {
    colnames(S) <- nms      
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
    txt = txt
  )
  class(out) <- "kernelshap"
  out
}