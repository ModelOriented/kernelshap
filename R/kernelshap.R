#' Kernel SHAP
#'
#' This function implements the model-agnostic Kernel SHAP algorithm explained in detail in 
#' Covert and Lee (2021). It is an iterative refinement of the original Kernel SHAP algorithm 
#' of Lundberg and Lee (2017).
#' The algorithm is applied to each row in \code{X}. Due to its iterative nature, approximate
#' standard errors of the resulting SHAP values are provided, and convergence is monitored. 
#' The data rows \code{X} to be explained and the
#' background data \code{bg_X} should only contain feature columns required by the
#' prediction function \code{pred_fun}. The latter is a function taking
#' a data structure like \code{X} and \code{bg_X} and providing one numeric 
#' prediction per row.
#' During each iteration, \code{m} subsets are evaluated until the worst standard error of the SHAP 
#' values is small enough relative to the range of the SHAP values. This exactly follows the logic
#' used by Covert and Lee (2021).
#' @param X Matrix or data.frame containing the observations to be explained.
#' Should only contain model features.
#' @param pred_fun A function taking objects like \code{X} or \code{bg_X} as input and providing numeric 
#' predictions. Example: If "fit" denotes a logistic regression fitted via \code{stats::glm}, 
#' and SHAP values should be on the probability scale, then this argument is
#' \code{function(X) predict(fit, X, type = "response")}.
#' @param bg_X Matrix or data.frame used as background data to calculate marginal 
#' expectations. Its column structure must be similar to \code{X}.
#' This data should neither be too small nor too large (50-200 rows). A large background
#' data slows down the calculations, while a small data set leads to imprecise SHAP values.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param paired_sampling Logical flag indicating whether to use paired sampling.
#' The default is \code{TRUE}. This means that with every feature subset S,
#' also its complement is evaluated, which leads to faster convergence.
#' @param m Number of feature subsets S to be evaluated during one iteration. 
#' The default, "auto", equals \code{trunc(20 * sqrt(ncol(X)))}. 
#' For the paired sampling strategy, 2m evaluations are done per iteration.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n) / diff(range(beta_n)) < tol, where sigma_n are the standard errors 
#' and beta_n are the SHAP values of a given observation.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, then the algorithm stops.
#' @param verbose Set to \code{FALSE} to suppress messages, warnings, and the progress bar.
#' @param ... Currently unused.
#' @return An object of class "kernelshap" with the following components:
#' \itemize{
#'   \item \code{S}: Matrix with SHAP values.
#'   \item \code{X}: Same as parameter \code{X}.
#'   \item \code{baseline}: The average prediction on the background data.
#'   \item \code{SE}: Standard errors corresponding to \code{S}.
#'   \item \code{n_iter}: Number of iterations until convergence per row.
#'   \item \code{converged}: Logical vector indicating convergence per row.
#' }
#' @export
#' @references
#' \enumerate{
#'   \item Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
#'   \item Scott M. Lundberg and Su-In Lee. A Unified Approach to Interpreting Model Predictions. Advances in Neural Information Processing Systems 30, 2017.
#'}
#' @examples
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' pred_fun <- function(X) stats::predict(fit, X)
#' s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[-1])
#' s
#' 
#' # Matrix input works as well, and pred_fun may contain preprocessing steps.
#' fit <- stats::lm(Sepal.Length ~ ., data = iris[1:4])
#' pred_fun <- function(X) stats::predict(fit, as.data.frame(X))
#' X <- data.matrix(iris[2:4])
#' s <- kernelshap(X[1:3, ], pred_fun = pred_fun, X)
#' s
kernelshap <- function(X, pred_fun, bg_X, bg_w = NULL, 
                       paired_sampling = TRUE, m = "auto", 
                       tol = 0.01, max_iter = 250, verbose = TRUE, ...) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.function(pred_fun),
    dim(X) >= 1L,
    ncol(bg_X) >= 1L,
    nrow(bg_X) >= 2L,
    !is.null(colnames(X)),
    !is.null(colnames(bg_X)),
    colnames(X) == colnames(bg_X)
  )
  if (!is.null(bg_w)) {
    stopifnot(
      length(bg_w) == nrow(bg_X),
      all(bg_w >= 0),
      !all(bg_w == 0)
    )
  }
  if (verbose && nrow(bg_X) > 1000) {
    warning("Your background data 'bg_X' is large, which will slow down the process. Consider using 50-200 rows.")
  }
  if (verbose && nrow(bg_X) < 10) {
    warning("Your background data 'bg_X' is small, which might lead to imprecise SHAP values. Consider using 50-200 rows.")
  }
  preds <- pred_fun(bg_X)
  stopifnot(
    "Predictions of multiple observations should be a vector (and not, e.g. a matrix)" = 
      is.vector(preds),
    "Predictions should be numeric" = is.numeric(preds), 
    length(preds) == nrow(bg_X)
  )
  
  # Initialization
  v0 <- weighted_mean(preds, bg_w)
  p <- ncol(X)
  n <- nrow(X)
  S <- SE <- matrix(0, nrow = n, ncol = p, dimnames = list(NULL, colnames(X)))
  n_iter <- integer(n)
  converged <- logical(n)
  if (m == "auto") {
    m <- trunc(20 * sqrt(p))
  }

  # Handle simple exact case
  if (p == 1L) {
    S <- matrix(pred_fun(X) - v0, ncol = 1, dimnames = list(NULL, colnames(X)))
    out <- list(
      S = S, X = X, baseline = v0, SE = SE, n_iter = n_iter, converged = rep(TRUE, n)
    )
    class(out) <- "kernelshap"
    return(out)
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
      paired = paired_sampling,
      m = m,
      tol = tol,
      max_iter = max_iter
    )
    
    S[i, ] <- res$beta
    SE[i, ] <- res$sigma
    n_iter[i] <- res$n_iter
    converged[i] <- res$converged
    
    if (verbose && n >= 2L) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (verbose && !all(converged)) {
    warning("\nNon-convergence for ", sum(!converged), " rows.")
  }
  out <- list(
    S = S, X = X, baseline = v0, SE = SE, n_iter = n_iter, converged = converged
  )
  class(out) <- "kernelshap"
  out
}

# Little helpers

# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(x, pred_fun, bg_X, bg_w, v0, paired, m, tol, max_iter) {
  v1 <- pred_fun(x)
  p <- ncol(x)
  X <- x[rep(1L, nrow(bg_X)), ]
  
  # Outer loop init
  est_m = list()
  converged <- FALSE
  n_iter <- 0L
  n <- 0L
  A <- matrix(0, nrow = p, ncol = p)
  b <- numeric(p)

  while(!isTRUE(converged) && n_iter < max_iter) {
    Z <- make_Z(m, p)
    
    # Calling get_vZ() is expensive
    vz <- get_vz(
      X,  bg = bg_X, Z = if (paired) rbind(Z, 1-Z) else Z, pred_fun = pred_fun, w = bg_w
    )
    
    # The inner loop can be replaced by vectorized operation, but is negligible
    counter <- 0L
    Atemp <- matrix(0, nrow = p, ncol = p)
    btemp <- numeric(p)
    
    for (i in 1:m) { # i <- 1
      z <- Z[i, ]
      if (paired) {
        Asample <- (tcrossprod(z) + tcrossprod(1 - z)) / 2
        bsample <- (z * vz[i] + (1 - z) * vz[i + m] - v0) / 2
      } else {
        Asample <- tcrossprod(z)
        bsample <- z * (vz[i] - v0)
      }

      # Welford's algorithm to iteratively calculate covariances etc
      n <- n + 1L
      A <- A + (Asample - A) / n
      b <- b + (bsample - b) / n
      counter <- counter + 1L
      Atemp <- Atemp + (Asample - Atemp) / counter
      btemp <- btemp + (bsample - btemp) / counter
    }

    est_m[[length(est_m) + 1L]] <- solver(Atemp, btemp, v1, v0)
    n_iter <- n_iter + 1L

    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(A, b, v1, v0)
      sigma_beta <- m * stats::cov(do.call(rbind, est_m))
      sigma_n <- sqrt(diag(sigma_beta) / n)
      converged <- max(sigma_n) / diff(range(beta_n)) < tol
    }
  }
  list(beta = beta_n, sigma = sigma_n, n_iter = n_iter, converged = converged)
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
  out <- vapply(len_S, function(z) {
    out <- numeric(p); out[sample(1:p, z)] <- 1; out}, FUN.VALUE = numeric(p)
  )
  t(out)
}

# This step takes time: marginal expectations for multiple z are calculated
get_vz <- function(X, bg, Z, pred_fun, w) {
  modify_X <- function(not_z) {
    X[, not_z] <- bg[, not_z, drop = FALSE]
    X
  }
  data_list <- apply(!Z, 1L, modify_X, simplify = FALSE)
  pred_data <- do.call(rbind, data_list)
  # pred_data <- dplyr::bind_rows(data_list) # Faster, but does not work for matrices
  preds <- pred_fun(pred_data)
  rowmean(
    preds, 
    group = rep(1:nrow(Z), each = nrow(bg)), 
    w = if (!is.null(w)) rep(w, times = nrow(Z))
  )
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
