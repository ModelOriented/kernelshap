# Functions required only for the sampling version of Kernel SHAP

# Draw m binary vectors z of length p with sum(z) distributed according to Kernel SHAP 
# weights -> (m x p) matrix. The argument S can be used to restrict the range of sum(z)
sample_Z <- function(m, p, S = 1:(p - 1L)) {
  # First draw s = sum(z)
  probs <- kernel_weights(p, S = S)
  len_S <- S[sample.int(length(S), m, replace = TRUE, prob = probs)]
  
  # Then, conditional on that number, set random positions of z to 1
  # Can this be done without loop/vapply?
  out <- vapply(
    len_S, 
    function(z) {out <- numeric(p); out[sample(1:p, z)] <- 1; out}, 
    FUN.VALUE = numeric(p)
  )
  t(out)
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

# Create Z, w, A for strategy "simple"
input_simple <- function(m, p) {
  Z <- sample_Z(m = m, p = p)
  list(Z = Z, w = 1 / m, A = crossprod(Z, Z / w))
}

# Create Z, w, A for strategy "paired"
input_paired <- function(m, p) {
  Z <- sample_Z(m = m / 2, p = p)
  Z <- rbind(Z, 1 - Z)
  list(Z = Z, w = 1 / m, A = crossprod(Z, Z / m))
}

# Create Z, w, A for strategy "hybrid"
input_hybrid <- function(m, p) {
  kw1 <- kernel_weights(p)[1L]

  # Enumerate all z with sum(z) = 1 and store their effective weights
  Z <- diag(p)
  w <- rep(kw1 / p, p)
  
  if (p >= 4L) {
    m_rest <- m / 2 - p
    Z <- rbind(Z, sample_Z(m_rest, p, S = 2:(p - 2L)))
    w <- c(w, rep(0.5 - kw1, times = m_rest) / m_rest)
  }
  Z <- rbind(Z, 1 - Z)
  w <- c(w, w)
  
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}
