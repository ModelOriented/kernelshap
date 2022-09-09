# Functions required only for the sampling version of Kernel SHAP

# m permutations distributed according to Kernel SHAP weights -> (m x p) matrix
sample_Z <- function(m, p, paired) {
  if (p < 2L) {
    stop("Sampling impossible for p < 2")
  }
  if (paired) { # Note that m is an even number
    m <- m / 2
  }
  
  # First draw number of elements in S
  len_S <- sample(1:(p - 1L), m, replace = TRUE, prob = kernel_weights(p))
  
  # Then, conditional on that number, set random positions to 1
  # Can this be done without loop/vapply?
  Z <- t(
    vapply(
      len_S, 
      function(z) {out <- numeric(p); out[sample(1:p, z)] <- 1; out}, 
      FUN.VALUE = numeric(p)
    )
  )
  if (paired) rbind(Z, 1 - Z) else Z
}

# Calculate standard error from list of m estimates
get_sigma <- function(est, iter) {
  apply(abind1(est), 3L, FUN = function(Y) sqrt(diag(stats::cov(Y)) / iter))
}

# Convergence criterion
conv_crit <- function(sig, bet) {
  if (any(dim(sig) != dim(bet))) {
    stop("sig must have same dimension as bet")
  }
  apply(sig, 2L, FUN = max) / apply(bet, 2L, FUN = function(z) diff(range(z)))
}