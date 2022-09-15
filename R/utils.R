# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(x, v1, object, pred_fun, bg_w, exact, deg, 
                           paired, m, tol, max_iter, v0, precalc, ...) {
  p <- ncol(x)

  # Calculate A_exact and b_exact
  if (exact || deg >= 1L) {
    A_exact <- precalc[["A"]]                                     #  (p x p)
    bg_X_exact <- precalc[["bg_X_exact"]]                         #  (m_ex*n_bg x p)
    Z <- precalc[["Z"]]                                           #  (m_ex x p)
    m_exact <- nrow(Z)
    v0_m_exact <- v0[rep(1L, m_exact), , drop = FALSE]            #  (m_ex x K)
    
    # Most expensive part
    vz <- get_vz(                                                 #  (m_ex x K)
      X = x[rep(1L, times = nrow(bg_X_exact)), , drop = FALSE],   #  (m_ex*n_bg x p)
      bg = bg_X_exact,                                            #  (m_ex*n_bg x p)
      Z = Z,                                                      #  (m_ex x p)
      object = object, 
      pred_fun = pred_fun, 
      w = bg_w, 
      ...
    )
    # Note: w is correctly replicated along columns of (vz - v0_m_exact)
    b_exact <- crossprod(Z, precalc[["w"]] * (vz - v0_m_exact))   #  (p x K)
    
    # Some of the hybrid cases are exact as well
    if (exact || p %in% (0:1 + (2L * deg))) {
      beta <- solver(A_exact, b_exact, constraint = v1 - v0)      #  (p x K)
      return(list(beta = beta, sigma = 0 * beta, n_iter = 1L, converged = TRUE))  
    }
  } 
  
  # Iterative sampling part, always using A_exact and b_exact to fill up the weights
  bg_X_m <- precalc[["bg_X_m"]]                                   #  (m*n_bg x p)
  X <- x[rep(1L, times = nrow(bg_X_m)), , drop = FALSE]           #  (m*n_bg x p)
  v0_m <- v0[rep(1L, m), , drop = FALSE]                          #  (m x K)

  est_m = list()
  converged <- FALSE
  n_iter <- 0L
  A_sum <- matrix(0, nrow = p, ncol = p)                          #  (p x p)
  b_sum <- matrix(0, nrow = p, ncol = ncol(v0))                   #  (p x K)
  if (deg == 0L) {
    A_exact <- A_sum
    b_exact <- b_sum
  }
  
  while(!isTRUE(converged) && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    input <- input_sampling(p = p, m = m, deg = deg, paired = paired)
    Z <- input[["Z"]]
      
    # Expensive                                                              #  (m x K)
    vz <- get_vz(
      X = X, bg = bg_X_m, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
    )
    
    # The sum of weights of A_exact and input[["A"]] is 1, same for b
    A_temp <- A_exact + input[["A"]]                                         #  (p x p)
    b_temp <- b_exact + crossprod(Z, input[["w"]] * (vz - v0_m))             #  (p x K)
    A_sum <- A_sum + A_temp                                                  #  (p x p)
    b_sum <- b_sum + b_temp                                                  #  (p x K)
    
    # Least-squares with constraint that beta_1 + ... + beta_p = v_1 - v_0. 
    # The additional constraint beta_0 = v_0 is dealt via offset   
    est_m[[n_iter]] <- solver(A_temp, b_temp, constraint = v1 - v0)          #  (p x K)

    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(A_sum / n_iter, b_sum / n_iter, constraint = v1 - v0) #  (p x K)
      sigma_n <- get_sigma(est_m, iter = n_iter)                             #  (p x K)
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
  s <- (matrix(colSums(Ainv %*% b), nrow = 1L) - constraint) / sum(Ainv)     #  (1 x K)
  Ainv %*% (b - s[rep(1L, p), , drop = FALSE])                               #  (p x K)
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

# Given p and maximum m, determine hybrid degree (currently not used)
find_degree <- function(p, m_max) {
  if (p < 2L) {
    "p must be at least 2"
  }
  if (m_max < 2L * p) {
    return(list(degree = 0L, m_exact = 0L))
  }
  # Non-integers are rounded down
  S <- seq_len(p / 2)
  const <- (2L - (p == 2L * S))
  m_kum <- cumsum(const * choose(p, S))
  degree <- max(S[m_kum <= m_max])
  list(degree = degree, m_exact = m_kum[degree])
}

# Describe what is happening
summarize_strategy <- function(p, exact, deg, m_exact, m) {
  if (exact || p %in% (0:1 + (2L * deg))) {
    return(
      sprintf(
        "Exact Kernel SHAP values %s(m_exact = %s)", 
        if (exact) "" else "by the hybrid approach ",
        m_exact
      )
    )
  }
  if (deg == 0L) {
    return(sprintf("Kernel SHAP values by iterative sampling (m/iter = %s)", m))
  } 
  sprintf(
    "Kernel SHAP values by the iterative hybrid strategy of degree %s (m_exact = %s, m/iter = %s)", 
    deg, 
    m_exact, 
    m
  )
}

# Kernel weights normalized to a non-empty subset S of {1, ..., p-1}
kernel_weights <- function(p, S = seq_len(p - 1L)) {
  probs <- (p - 1L) / (choose(p, S) * S * (p - S))
  probs / sum(probs)
}
