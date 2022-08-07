#' Kernel SHAP
#'
#' This function implements the Kernel SHAP Algorithm 1 of Covert and Lee (2021)
#' with pairwise sampling, see reference.
#' It is applied separately to each row in \code{X}. Due to its iterative nature,
#' standard errors of the resulting SHAP values can be provided. During each iteration,
#' \code{2m} feature subsets are evaluated, until the worst standard error is small enough.
#' 
#' @param X Matrix or data.frame containing the observations to be explained.
#' Should only contain features required in \code{pred_fun}, i.e., it should only contain
#' features.
#' @param pred_fun A function taking objects like \code{X} as input and providing numeric 
#' predictions.
#' @param bg_X Matrix or data.frame used as background dataset to calculate marginal 
#' expectations. Its column structure must be similar to \code{X}.
#' If too large (>120 rows). Use subsampling or some more sophisticated strategy. 
#' @param bg_w An optional vector of case weights for each row of \code{bg_X}.
#' @param m The number of feature subsets S to be evaluated during one iteration. 
#' By default \code{trunc(20 * sqrt(ncol(bg_X)))}. Since we use the pairwise strategy,
#' the actual number of evaluations is 2m.
#' @param tol Tolerance determining when to stop. The algorithm keeps sampling until
#' max(sigma_n) / diff(range(beta_n)) < tol, where sigma_n are the standard errors 
#' and beta_n are the SHAP values of a given observation.
#' @param verbose Set to \code{FALSE} to suppress messages and progress bar.
#' @param ... Currently unused.
#' @return An object of class "kernelshap" with the following components:
#' \itemize{
#'   \item \code{S}: Matrix with SHAP values.
#'   \item \code{X}: Same as parameter \code{X}.
#'   \item \code{baseline}: The average prediction on the background data.
#'   \item \code{SE}: Standard errors corresponding to \code{S}.
#' }
#' @export
#' @references Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(head(iris[-1], 1), function(X) stats::predict(fit, X), iris[-1])
#' s
kernelshap <- function(X, pred_fun, bg_X = X, bg_w = NULL, 
                       m = trunc(20 * sqrt(ncol(bg_X))), tol = 0.01, 
                       verbose = TRUE, ...) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.function(pred_fun),
    dim(X) >= 1L,
    dim(bg_X) >= 1L,
    colnames(X) == colnames(bg_X)
  )
  if (!is.null(bg_w)) {
    stopifnot(
      length(bg_w) == nrow(bg_X),
      all(bg_w >= 0),
      !all(bg_w == 0)
    )
  }
  if (verbose && nrow(bg_X) > 120) {
    message("Your background data is large. This will slow down the process.")
  }
  preds <- pred_fun(bg_X)
  stopifnot(
    "Predictions should be numeric" = is.numeric(preds), 
    length(preds) == nrow(bg_X)
  )
  # Baseline
  v0 <- weighted_mean(preds, bg_w)
  
  # Handle simple (exact) cases
  p <- ncol(X)
  n <- nrow(X)
  S <- SE <- matrix(0, nrow = n, ncol = p, dimnames = list(NULL, colnames(X)))

  if (p <= 1L) {
    if (p == 1L) {
      S <- as.matrix(preds - v0, ncol = 1, dimnames = list(NULL, colnames(X)))
    }
    out <- list(S = S, X = X, baseline = v0, SE = SE)
    class(out) <- "kernelshap"
    out
  }
  
  # KernelSHAP: loop over rows of X
  if (verbose && n >= 2L) {
    pb <- utils::txtProgressBar(1, n, style = 3)  
  }
  for (i in seq_len(n)) {
    res <- kernelshap_one(
      X[i, , drop = FALSE], 
      pred_fun = pred_fun, 
      bg_X = bg_X, 
      bg_w = bg_w, 
      v0 = v0,
      m = m,
      tol = tol
    )
    S[i, ] <- res$beta
    SE[i, ] <- res$sigma
    if (verbose && n >= 2L) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  out <- list(S = S, X = X, baseline = v0, SE = SE)
  class(out) <- "kernelshap"
  out
}

# Little helpers

# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(x, pred_fun, bg_X, bg_w, v0, m, tol) {
  v1 <- pred_fun(x)
  p <- ncol(x)
  X <- x[rep(1, nrow(bg_X)), ]
  group <- rep(1:m, each = nrow(bg_X))
  w <- if (!is.null(bg_w)) rep(bg_w, times = m)

  n_runs <- n <- counter <- 0L
  A <- Atemp <- matrix(0, nrow = p, ncol = p)
  b <- btemp <- numeric(p)
  est_m = list()

  converged <- FALSE

  while(!isTRUE(converged)) {
    Z <- make_Z(m, p)
    
    # Pairwise strategy
    v_z1 <- rowmean(
      pred_fun(modify_and_stack(X, bg = bg_X, Z = Z)), group = group, w = w
    )
    v_z2 <- rowmean(
      pred_fun(modify_and_stack(X, bg = bg_X, Z = 1 - Z)), group = group, w = w
    )

    # Maybe loop can be replaced by vectorized operation
    for (i in 1:m) { # i <- 1
      z <- Z[i, ]
      Asample <- (tcrossprod(z) + tcrossprod(1 - z)) / 2
      bsample <- (z * v_z1[i] + (1 - z) * v_z2[i] - v0) / 2

      # Welford's algorithm to iteratively calculate covariances etc
      n <- n + 1L
      A <- A + (Asample - A) / n
      b <- b + (bsample - b) / n
      counter <- counter + 1L
      Atemp = Atemp + (Asample - Atemp) / counter
      btemp = btemp + (bsample - btemp) / counter
    }

    est_m[[length(est_m) + 1L]] <- solver(Atemp, btemp, v1, v0)

    counter <- 0L
    Atemp <- matrix(0, nrow = p, ncol = p)
    btemp <- numeric(p)

    n_runs <- n_runs + 1L

    if (n_runs >= 2L) {
      beta_n <- solver(A, b, v1, v0)
      sigma_beta <- m * stats::cov(do.call(rbind, est_m))
      sigma_n <- sqrt(diag(sigma_beta) / n)
      converged <- max(sigma_n) / diff(range(beta_n)) < tol
    }
  }
  list(beta = beta_n, sigma = sigma_n, nruns = n_runs)
}

# Simple version without standard errors
# kernelshap_one <- function(x, pred_fun, bg_X, bg_w, v0, tol = 0) {
#   v1 <- pred_fun(x)
#   p <- ncol(x)
#   m <- trunc(500 * sqrt(p))
#   X <- x[rep(1, nrow(bg_X)), ]
#   group <- rep(1:m, each = nrow(bg_X))
#   w <- if (!is.null(bg_w)) rep(bg_w, times = m)
#   
#   Z <- make_Z(m, p)
#   vz <- rowmean(
#     pred_fun(modify_and_stack(X, bg = bg_X, Z = Z)), group = group, w = w
#   )
#   A <- crossprod(Z) / m
#   b <- crossprod(Z, (vz - v0)) / m
# 
#   list(beta = solver(A, b, v1, v0), sigma = 0)
# }

# Calculates regression coefficients
solver <- function(A, b, v1, v0) {
  Ainv <- solve(A)
  out <- Ainv %*% (b - (sum(crossprod(Ainv, b)) - v1 + v0) / sum(Ainv))
  as.numeric(out)
}

# Generates m permutations with unequal probability
make_Z <- function(m, p) {
  if (p <= 1L) {
    stop("p must be 2 or larger")
  }
  S <- 1:(p - 1)
  # First draw number of elements in S
  probs <- (p - 1) / (choose(p, S) * S * (p - S))
  len_S <- sample(S, m, replace = TRUE, prob = probs)
  
  # Then, conditional on that number, set random positions to 1
  # Can this be done without loop?
  out <- vapply(len_S, function(z) {
    out <- numeric(p); out[sample(1:p, z)] <- 1; out}, FUN.VALUE = numeric(p)
  )
  t(out)
}

# This is the part that takes most time: the prediction data is being generated
modify_and_stack <- function(X, bg, Z) {
  data_list <- vector("list", nrow(Z))
  for (j in seq_len(nrow(Z))) {
    X_mod <- X
    r <- !Z[j, ]
    X_mod[, r] <- bg[, r, drop = FALSE]
    data_list[[j]] <- X_mod
  }
  dplyr::bind_rows(data_list) # much faster than rbind()
}

# Convenience wrapper around mean and weighted.mean
weighted_mean <- function(x, w = NULL, ...) {
  if (is.null(w)) {
    return(mean(x, ...))
  }
  stats::weighted.mean(x, w = w, ...)
}

# Fast way to calculate grouped weighted means
rowmean <- function(x, group, w = NULL) {
  if (is.null(w)) {
    w <- rep(1.0, length(x))
  }
  rowsum(x * w, group) / rowsum(w, group)
}
