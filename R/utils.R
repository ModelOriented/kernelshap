# Calculates regression coefficients
solver <- function(A, b, v1, v0) {
  p <- ncol(A)
  Ainv <- solve(A)
  s <- (matrix(colSums(crossprod(Ainv, b)), nrow = 1L) - v1 + v0) / sum(Ainv)  # (1 x K)
  Ainv %*% (b - repn(s, p))                                                    # (p x K)
}

# Generates m permutations distributed according to Kernel SHAP weights
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
  out <- vapply(
    len_S, 
    function(z) {out <- numeric(p); out[sample(1:p, z)] <- 1; out}, 
    FUN.VALUE = numeric(p)
  )
  t(out)
}

# This function is by far the most expensive: data manipulation and prediction
# We try to work with matrices or data.tables. Note that stacking many evaluation
# data and then calling pred_fun once is much less expensive than calling pred_fun
# for each evaluation data separately.
get_vz <- function(X, bg, Z, pred_fun, w, use_dt) {
  if (inherits(X, "data.table") && use_dt) {
    modify_X <- function(not_z) {
      X_mod <- data.table::copy(X)
      for (j in which(not_z)) {
        data.table::set(X_mod, j = j, value = bg[[j]])
      }
      X_mod
    }
    data_list <- apply(!Z, 1L, modify_X, simplify = FALSE)
    pred_data <- data.table::rbindlist(data_list)
  } else {
    modify_X <- function(not_z) {
      X[, not_z] <- bg[, not_z, drop = FALSE]
      X
    }
    data_list <- apply(!Z, 1L, modify_X, simplify = FALSE)
    if (!use_dt) {
      pred_data <- do.call(rbind, data_list)
    } else {
      pred_data <- as.data.frame(data.table::rbindlist(data_list))
    }
  }
  preds <- check_pred(pred_fun(pred_data), n = nrow(pred_data))
  rowmean(preds, n_bg = nrow(bg), n_z = nrow(Z), w = w)
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

# Average prediction per row of Z (optimized for speed)
rowmean <- function(P, n_bg, n_z, w = NULL) {
  g <- rep(seq_len(n_z), each = n_bg)
  if (is.null(w)) {
    return(rowsum(P, group = g, reorder = FALSE) / n_bg)
  }
  rowsum(P * rep(w, times = n_z), group = g, reorder = FALSE) / sum(w)
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

# Replicate single row data.frame or matrix n times
repn <- function(x, n) {
  if ((!is.matrix(x) && !is.data.frame(x)) || nrow(x) != 1L) {
    stop("x must be a matrix or data.frame with 1 row")
  }
  x[rep(1L, times = n), , drop = FALSE]
}
