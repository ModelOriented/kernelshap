# Functions required only for the sampling version of Kernel SHAP

# Functions required only for the sampling version of Kernel SHAP

# Draw m binary vectors z of length p with sum(z) distributed according to Kernel SHAP 
# weights -> (m x p) matrix. The argument S can be used to restrict the range of sum(z)
sample_Z <- function(p, m, S = 1:(p - 1L)) {
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

# Vectorized version of the vapply part, credits: Mathias Ambuehl
# f <- function(N, p = 4L) {
#   l <- length(N)
#   out <- rep(rep(0:1, l), as.vector(rbind(p-N, N)))
#   dim(out) <- c(p, l)
#   ord <- order(col(out), sample.int(l*p))
#   out[] <- out[ord]
#   t.default(out)
# }

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

# Create Z, w, A from sampling
input_sampling <- function(p, m, deg, paired) {
  S <- (deg + 1L):(p - deg - 1L)
  Z <- sample_Z(m = if (paired) m / 2 else m, p = p, S = S)
  if (paired) {
    Z <- rbind(Z, 1 - Z)
  }
  w_total <- if (deg == 0L) 1 else 1 - 2 * sum(kernel_weights(p)[seq_len(deg)])
  w <- w_total / m
  list(Z = Z, w = rep(w, m), A = crossprod(Z) * w)
}

