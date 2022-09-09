# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(object, X, bg_X, pred_fun, bg_w, v0, v1, 
                           sampling_strategy, m, precalc, tol, max_iter, ...) {
  p <- ncol(X)
  v0_ext <- v0[rep(1L, m), , drop = FALSE]                        #  (m x K)
  est_m = list()
  converged <- FALSE
  n_iter <- 0L
  Asum <- matrix(0, nrow = p, ncol = p)                           #  (p x p)
  bsum <- matrix(0, nrow = p, ncol = ncol(v0))                    #  (p x K)
  
  while(!isTRUE(converged) && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    
    # Get Z, w, A for any strategy
    input <- switch(
      sampling_strategy,
      hybrid = input_hybrid(m = m, p = p, pairs = precalc),
      simple = input_simple_paired(m = m, p = p, paired = FALSE),
      paired = input_simple_paired(m = m, p = p, paired = TRUE),
      exact = precalc
    )
    Z <- input[["Z"]]                                             #  (m x p)
    
    # Expensive step                                              #  (m x K)
    vz <- get_vz(
      X = X, bg = bg_X, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
    )
    
    # Least-squares with constraint that beta_1 + ... + beta_p = v_1 - v_0. 
    # The additional constraint beta_0 = v_0 is dealt via offset
    Atemp <- input[["A"]]                                         #  (p x p)
    btemp <- crossprod(Z, input[["w"]] * (vz - v0_ext))           #  (p x K)
    Asum <- Asum + Atemp                                          #  (p x p)
    bsum <- bsum + btemp                                          #  (p x K)
    
    # Constrained regression -> parameter matrix                  #  (p x K)
    est_m[[n_iter]] <- beta <- solver(Atemp, btemp, constraint = v1 - v0)
    
    if (sampling_strategy == "exact") {
      return(list(beta = beta, sigma = 0 * beta, n_iter = 1L, converged = TRUE))
    }

    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(Asum / n_iter, bsum / n_iter, constraint = v1 - v0)  # (p x K)
      sigma_n <- get_sigma(est_m, iter = n_iter)                            # (p x K)
      converged <- all(conv_crit(sigma_n, beta_n) < tol)
    }
  }
  list(beta = beta_n, sigma = sigma_n, n_iter = n_iter, converged = converged)
}

# Regression coefficients given sum(beta) = constraint
# A: (p x p), b: (p x k), constraint: (1 x K)
solver <- function(A, b, constraint) {
  p <- ncol(A)
  Ainv <- MASS::ginv(A)
  s <- (matrix(colSums(Ainv %*% b), nrow = 1L) - constraint) / sum(Ainv)  # (1 x K)
  Ainv %*% (b - s[rep(1L, p), , drop = FALSE])                            # (p x K)
}

# Calculates all vz of an iteration and thus takes time
get_vz <- function(X, bg, Z, object, pred_fun, w, ...) {
  m <- nrow(Z)
  not_Z <- !Z
  n_bg <- nrow(bg) / m   # because bg was replicated m times
  
  # Replicate not_Z, so that X, bg, not_Z are all of dimension (m*n_bg x p)
  g <- rep(seq_len(m), each = n_bg)
  not_Z <- not_Z[g, , drop = FALSE]
  
  if (is.matrix(X)) {
    X[not_Z] <- bg[not_Z]
  } else {
    for (j in seq_len(ncol(bg))) {
      s <- not_Z[, j, drop = TRUE]
      X[[j]][s] <- bg[[j]][s]
    }
  }
  preds <- check_pred(pred_fun(object, X, ...), n = nrow(X))
  
  # Aggregate
  if (is.null(w)) {
    return(rowsum(preds, group = g, reorder = FALSE) / n_bg)
  }
  rowsum(preds * rep(w, times = m), group = g, reorder = FALSE) / sum(w)
}

# Weighted colMeans(). Always returns a (1 x ncol(x)) matrix
weighted_colMeans <- function(x, w = NULL, ...) {
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  out <- if (is.null(w)) colMeans(x, ...) else colSums(x * w, ...) / sum(w)
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
  if (
    !is.vector(x) && 
    !is.matrix(x) && 
    !is.data.frame(x) && 
    !(is.array(x) && length(dim(x)) <= 2L)
  ) {
    stop("Predictions must be a vector, matrix, data.frame, or <=2D array")
  }
  if (is.data.frame(x) || is.array(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) {
    stop("Predictions must be numeric")
  }
  if (is.matrix(x) && nrow(x) == n) {
    return(x)
  }
  if (length(x) == n) {
    return(matrix(x, nrow = n))
  }
  stop("Predictions must be a length n vector or a matrix/data.frame/array with n rows.")
}

# Informative warning if background data is small or large
check_bg_size <- function(n) {
  if (n > 1000L) {
    warning("Your background data 'bg_X' is large, which will slow down the process. Consider using 100-200 rows.")
  }
  if (n < 20L) {
    warning("Your background data 'bg_X' is small, which might lead to imprecise SHAP values. Consider using 100-200 rows.")
  }
}

# Kernel weights normalized to a non-empty subset S of {1, ..., p-1}
kernel_weights <- function(p, S = seq_len(p - 1L)) {
  if (length(S) == 0L) {
    stop("S must be non-empty")
  }
  if (!all(S %in% 1:(p - 1L))) {
    stop("S must be subset of 1:(p-1)")
  }
  probs <- (p - 1L) / (choose(p, S) * S * (p - S))
  probs / sum(probs)
}

