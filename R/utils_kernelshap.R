# Kernel SHAP algorithm for a single row x
# If exact, a single call to predict() is necessary.
# If sampling is involved, we need at least two additional calls to predict().
kernelshap_one <- function(
    x,
    v1,
    object,
    pred_fun,
    feature_names,
    bg_w,
    exact,
    deg,
    m,
    tol,
    max_iter,
    v0,
    precalc,
    ...) {
  p <- length(feature_names)
  K <- ncol(v1)
  K_names <- colnames(v1)

  # Calculate A_exact and b_exact
  if (exact || deg >= 1L) {
    A_exact <- precalc[["A"]] #  (p x p)
    bg_X_exact <- precalc[["bg_X_exact"]] #  (m_ex*n_bg x p)
    Z <- precalc[["Z"]] #  (m_ex x p)
    m_exact <- nrow(Z)
    v0_m_exact <- v0[rep.int(1L, m_exact), , drop = FALSE] #  (m_ex x K)

    # Most expensive part
    vz <- get_vz(
      x = x, bg = bg_X_exact, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
    )
    # Note: w is correctly replicated along columns of (vz - v0_m_exact)
    b_exact <- crossprod(Z, precalc[["w"]] * (vz - v0_m_exact)) #  (p x K)

    # Some of the hybrid cases are exact as well
    if (exact || trunc(p / 2) == deg) {
      beta <- solver(A_exact, b_exact, constraint = v1 - v0) #  (p x K)
      return(list(beta = beta))
    }
  }

  # Iterative sampling part, always using A_exact and b_exact to fill up the weights
  bg_X_m <- precalc[["bg_X_m"]] #  (m*n_bg x p)
  v0_m <- v0[rep.int(1L, m), , drop = FALSE] #  (m x K)
  est_m <- array(
    data = 0, dim = c(max_iter, p, K), dimnames = list(NULL, feature_names, K_names)
  )
  converged <- FALSE
  n_iter <- 0L
  A_sum <- matrix( #  (p x p)
    data = 0, nrow = p, ncol = p, dimnames = list(feature_names, feature_names)
  )
  b_sum <- matrix(0, nrow = p, ncol = K, dimnames = list(feature_names, K_names)) # (p x K)
  if (deg == 0L) {
    A_exact <- A_sum
    b_exact <- b_sum
  }

  while (!converged && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    input <- input_sampling(p = p, m = m, deg = deg, feature_names = feature_names)
    Z <- input[["Z"]]

    # Expensive                                                              #  (m x K)
    vz <- get_vz(
      x = x, bg = bg_X_m, Z = Z, object = object, pred_fun = pred_fun, w = bg_w, ...
    )

    # The sum of weights of A_exact and input[["A"]] is 1, same for b
    A_temp <- A_exact + input[["A"]] #  (p x p)
    b_temp <- b_exact + crossprod(Z, input[["w"]] * (vz - v0_m)) #  (p x K)
    A_sum <- A_sum + A_temp #  (p x p)
    b_sum <- b_sum + b_temp #  (p x K)

    # Least-squares with constraint that beta_1 + ... + beta_p = v_1 - v_0.
    # The additional constraint beta_0 = v_0 is dealt via offset
    est_m[n_iter, , ] <- solver(A_temp, b_temp, constraint = v1 - v0) #  (p x K)

    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(A_sum / n_iter, b_sum / n_iter, constraint = v1 - v0) #  (p x K)
      sigma_n <- get_sigma(est_m[seq_len(n_iter), , , drop = FALSE]) #  (p x K)
      converged <- check_convergence(beta = beta_n, sigma = sigma_n, tol = tol)
    } else if (max_iter == 1L) { # if n_iter == 1 and max_iter == 1
      beta_n <- solver(A_sum, b_sum, constraint = v1 - v0) #  (p x K)
      converged <- FALSE
      sigma_n <- NA * beta_n #  (p x K)
    }
  }
  list(beta = beta_n, sigma = sigma_n, n_iter = n_iter, converged = converged)
}

# Regression coefficients given sum(beta) = constraint
# A: (p x p), b: (p x k), constraint: (1 x K)
# Full credits: https://github.com/iancovert/shapley-regression/blob/master/shapreg/shapley.py
solver <- function(A, b, constraint) {
  Ainv1 <- solve(A, matrix(1, nrow = nrow(A)))
  Ainvb <- solve(A, b)
  num <- rbind(colSums(Ainvb)) - constraint
  return(Ainvb - Ainv1 %*% num / sum(Ainv1))
}


# Draw m binary vectors z of length p with sum(z) distributed according
# to Kernel SHAP weights -> (m x p) matrix.
# The argument S can be used to restrict the range of sum(z).
sample_Z <- function(p, m, feature_names, S = 1:(p - 1L)) {
  probs <- kernel_weights(p, per_coalition_size = TRUE, S = S)
  N <- S[sample.int(length(S), m, replace = TRUE, prob = probs)]

  # Then, conditional on that number, set random positions of z to 1
  # Original, unvectorized code
  # out <- vapply(
  #   N,
  #   function(z) {out <- numeric(p); out[sample(1:p, z)] <- 1; out},
  #   FUN.VALUE = numeric(p)
  # )
  # t(out)

  # Vectorized by Mathias Ambuehl
  out <- rep(rep.int(c(FALSE, TRUE), m), as.vector(rbind(p - N, N)))
  dim(out) <- c(p, m)
  ord <- order(col(out), sample.int(m * p))
  out[] <- out[ord]
  rownames(out) <- feature_names
  t(out)
}

# Provides random input for paired SHAP sampling:
# - Z: Matrix with m on-off vectors z with sum(z) following Kernel weight distribution.
# - w: Vector (1/m, 1/m, ...) of length m (if pure sampling)
# - A: Matrix A = Z'wZ
# The weights are constant (Kernel weights have been used to draw the z vectors).
#
# If deg > 0, vectors z with sum(z) restricted to [deg+1, p-deg-1] are sampled.
# This case is used in combination with input_partly_hybrid(). Consequently, sum(w) < 1.
input_sampling <- function(p, m, deg, feature_names) {
  if (p < 2L * deg + 2L) {
    stop("p must be >=2*deg + 2")
  }
  S <- (deg + 1L):(p - deg - 1L)
  Z <- sample_Z(p = p, m = m / 2, feature_names = feature_names, S = S)
  Z <- rbind(Z, !Z)
  w <- if (deg == 0L) 1 else 1 - prop_exact(p, deg = deg)
  w <- w / m
  list(Z = Z, w = rep.int(w, m), A = crossprod(Z) * w)
}

# Functions required only for handling (partly) exact cases

# Provides fixed input for the exact case:
# - Z: Matrix with all 2^p-2 on-off vectors z
# - w: Vector with row weights of Z ensuring that the distribution of sum(z) matches
#      the SHAP kernel distribution
# - A: Exact matrix A = Z'wZ
input_exact <- function(p, feature_names) {
  Z <- exact_Z(p, feature_names = feature_names)
  Z <- Z[2L:(nrow(Z) - 1L), , drop = FALSE]
  kw <- kernel_weights(p, per_coalition_size = FALSE) # Kernel weights for all subsets
  w <- kw[rowSums(Z)] # Corresponding weight for each row in Z
  w <- w / sum(w)
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}

# List all length p vectors z with sum(z) in {k, p - k}
partly_exact_Z <- function(p, k, feature_names) {
  if (k < 1L) {
    stop("k must be at least 1")
  }
  if (p < 2L * k) {
    stop("p must be >=2*k")
  }
  if (k == 1L) {
    Z <- diag(p)
    storage.mode(Z) <- "logical"
  } else {
    Z <- t(
      utils::combn(seq_len(p), k, FUN = function(z) {
        x <- logical(p)
        x[z] <- TRUE
        x
      })
    )
  }
  if (p != 2L * k) {
    Z <- rbind(Z, !Z)
  }
  colnames(Z) <- feature_names
  Z
}

# Create Z, w, A for vectors z with sum(z) in {k, p-k} for k in {1, ..., deg}.
# The total weights do not sum to one, except in the special (exact) case deg=p-deg.
# (The remaining weight will be added via input_sampling(p, deg=deg)).
# Note that for a given k, the weights are constant.
input_partly_exact <- function(p, deg, feature_names) {
  if (deg < 1L) {
    stop("deg must be at least 1")
  }
  if (p < 2L * deg) {
    stop("p must be >=2*deg")
  }

  kw <- kernel_weights(p, per_coalition_size = FALSE)

  Z <- vector("list", deg)
  for (k in seq_len(deg)) {
    Z[[k]] <- partly_exact_Z(p, k = k, feature_names = feature_names)
  }
  Z <- do.call(rbind, Z)
  w <- kw[rowSums(Z)]
  w_target <- prop_exact(p, deg = deg) # How much of total weight to spend here
  w <- w / sum(w) * w_target
  list(Z = Z, w = w, A = crossprod(Z, w * Z))
}

# Kernel weight distribution
#
# `per_coalition_size = TRUE` is required, e.g., when one wants to sample random masks
# according to the Kernel SHAP distribution: Pick a coalition size as per
# these weights, then randomly place "on" positions. `FALSE` refer to weights
# if all masks has been calculated and one wants to calculate their weights based
# on the number of "on" positions.
kernel_weights <- function(p, per_coalition_size, S = seq_len(p - 1L)) {
  const <- if (per_coalition_size) 1 else choose(p, S)
  probs <- (p - 1) / (const * S * (p - S)) # could drop the numerator
  return(probs / sum(probs))
}

# How much Kernel SHAP weights do coalitions of size
# {1, ..., deg, ..., p-deg-1 ..., p-1} have?
prop_exact <- function(p, deg) {
  if (deg == 0) {
    return(0)
  }
  w <- kernel_weights(p, per_coalition_size = TRUE)
  w_total <- 2 * sum(w[seq_len(deg)]) - w[deg] * (p == 2 * deg)
  return(w_total)
}
