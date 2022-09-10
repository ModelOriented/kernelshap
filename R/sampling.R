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

# List all z vectors with sum(z) = 2
all_pairs <- function(p) {
  t(utils::combn(seq_len(p), 2L, FUN = function(z) {x <- numeric(p); x[z] <- 1; x}))
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

# Create Z, w, A for strategy "simple" and "paired"
input_simple_paired <- function(p, m, paired) {
  Z <- sample_Z(m = if (paired) m / 2 else m, p = p)
  if (paired) {
    Z <- rbind(Z, 1 - Z)
  }
  w <- rep(1 / m, m)
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}

# Create Z, w, A for strategy "hybrid"
input_hybrid <- function(p, m, pairs = NULL) {
  deg <- hybrid_degree(p = p, m = m)
  kw <- kernel_weights(p)
  need_sampling <- FALSE

  if (deg == 0L) {
    stop("m is too small for the hybrid strategy, please make it larger")
  }

  # Enumerate all z with sum(z) = 1 and store their effective weights
  Z <- diag(p)
  w <- rep(kw[1L] / p, p)
  m_rest <- m / 2 - p
  
  # Can we also enumerate all z with sum(z) = 2 while keeping a rest proportional to
  # Kernel SHAP weights?
  n_pairs <- p * (p - 1L) / 2
  if (deg == 2L) {
    if (is.null(pairs)) {
      pairs <- all_pairs(p)
    }
    Z <- rbind(Z, pairs)
    w_pairs <- rep(kw[2L], n_pairs) / n_pairs
    w <- c(w, w_pairs)
    m_rest <- m_rest - n_pairs
    if (p >= 6L) {
      need_sampling <- TRUE
      K <- 3L
    }
  } else {
    need_sampling <- m_rest > 0L
    K <- 2L
  }
  if (need_sampling) {
    # Paired sampling for the rest
    Z_rest <- sample_Z(m = m_rest, p = p, S = K:(p - K))
    w_rest <- rep(0.5 - sum(kw[seq_len(K - 1L)]), times = m_rest) / m_rest
    Z <- rbind(Z, Z_rest)
    w <- c(w, w_rest)
  }
  
  # Mirror everything, except when p = 4 where the pairs are already mirrored
  if (p == 4L) {
    Z <- rbind(Z, 1 - diag(p))
    w <- c(w, rep(kw[1L] / p, p))
  } else {
    Z <- rbind(Z, 1 - Z)
    w <- c(w, w)
  }
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}
