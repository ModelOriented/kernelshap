# Functions required only for handling exact cases

# Provides fixed input for the exact case:
# - Z: Matrix with all 2^p-2 on-off vectors z
# - w: Vector with row weights of Z ensuring that the distribution of sum(z) matches
#      the SHAP kernel distribution
# - A: Exact matrix A = Z'wZ
exact_input <- function(p) {
  Z <- enumerate_Z(p)
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
enumerate_Z <- function(p) {
  Z <- as.matrix(do.call(expand.grid, replicate(p, 0:1, simplify = FALSE)))
  dimnames(Z) <- NULL
  Z[2:(nrow(Z) - 1L), , drop = FALSE]
}

# Case p = 1 returns exact Shapley values
case_p1 <- function(n, nms, v0, v1, X) {
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
    m = 0L
  )
  class(out) <- "kernelshap"
  out
}