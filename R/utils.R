# Calculates regression coefficients
solver <- function(A, b, v1, v0) {
  p <- ncol(A)
  # Ainv <- solve(A)
  Ainv <- MASS::ginv(A)
  s <- (matrix(colSums(Ainv %*% b), nrow = 1L) - v1 + v0) / sum(Ainv)  # (1 x K)
  Ainv %*% (b - s[rep(1L, p), , drop = FALSE])                         # (p x K)
}

# m permutations distributed according to Kernel SHAP weights -> (p x m) matrix
make_Z <- function(m, p) {
  if (p <= 1L) {
    stop("p must be 2 or larger")
  }
  S <- 1:(p - 1)
  # First draw number of elements in S
  probs <- (p - 1) / (choose(p, S) * S * (p - S))
  len_S <- sample(S, m, replace = TRUE, prob = probs)
  
  # Then, conditional on that number, set random positions to 1
  # Can this be done without loop/vapply?
  vapply(
    len_S, 
    function(z) {out <- numeric(p); out[sample(1:p, z)] <- 1; out}, 
    FUN.VALUE = numeric(p)
  )
}

# Calculates all vz of an iteration and thus takes time
get_vz <- function(X, bg, Z, pred_fun, w) {
  not_Z <- t(!Z)
  n_Z <- nrow(not_Z)
  n_bg <- nrow(bg) / n_Z   # Remember that bg was replicated n_Z times
  
  # Replicate not_Z, so that X, bg, not_Z are all of dimension (n_Z*n_bg x p)
  g <- rep(seq_len(n_Z), each = n_bg)
  not_Z <- not_Z[g, , drop = FALSE]
  
  if (is.matrix(X)) {
    X[not_Z] <- bg[not_Z]
  } else {
    for (j in seq_len(ncol(bg))) {
      s <- not_Z[, j, drop = TRUE]
      X[[j]][s] <- bg[[j]][s]
    }
  }
  preds <- check_pred(pred_fun(X), n = nrow(X))
  
  # Aggregate
  if (is.null(w)) {
    return(rowsum(preds, group = g, reorder = FALSE) / n_bg)
  }
  rowsum(preds * rep(w, times = n_Z), group = g, reorder = FALSE) / sum(w)
}

# Weighted colMeans(). Always returns a (1 x ncol(x)) matrix
weighted_colMeans <- function(x, w = NULL, ...) {
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  if (is.null(w)) {
    out <- colMeans(x, ...)
  } else {
    out <- colSums(x * w, ...) / sum(w)  
  }
  matrix(out, nrow = 1L)
}

# Binds list of matrices along new first axis
abind1 <- function(a) {
  out <- array(dim = c(length(a), dim(a[[1L]])))
  for (i in seq_along(a)) {
    out[i, , ] <- a[[i]]
  }
  out
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

# Turn list of n (p x K) matrices into list of K (n x p) matrices. Reduce if K = 1.
reorganize_list <- function(alist, nms) {
  if (!is.list(alist)) {
    stop("alist must be a list")
  }
  out <- abind1(alist)
  dimnames(out) <- list(NULL, nms, NULL)
  out <- asplit(out, MARGIN = 3L)
  if (length(out) == 1L) {
    return(as.matrix(out[[1L]]))
  }
  lapply(out, as.matrix)
}

# Checks and reshapes predictions to (n x K) matrix
check_pred <- function(x, n) {
  if (!is.numeric(x)) {
    stop("Predictions must be numeric")
  }
  if (is.matrix(x) && nrow(x) == n) {
    return(x)
  }
  if (length(x) == n) {
    return(matrix(x, nrow = n))
  }
  stop("Predictions must be a length n vector or a matrix with n rows.")
}

