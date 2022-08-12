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
#' a data structure like \code{X} or \code{bg_X} and providing one numeric 
#' prediction per row.
#' During each iteration, \code{m} feature subsets are evaluated until the worst 
#' standard error of the SHAP values is small enough relative to the range of the SHAP values. 
#' This stopping criterion was suggested in Covert and Lee (2021).
#' @param X Feature matrix or data.frame of dimension (n x p) containing the observations 
#' to be explained. Should only contain model features.
#' @param pred_fun A function taking a matrix/data.frame like \code{X} 
#' as input and providing numeric predictions. 
#' Example: If "fit" denotes a logistic regression fitted via \code{stats::glm}, 
#' and SHAP values should be on the probability scale, then this argument is
#' \code{function(X) predict(fit, X, type = "response")}.
#' @param bg_X Matrix or data.frame used as background data to calculate marginal 
#' expectations. Its column structure must be similar to \code{X}.
#' It should neither be too small nor too large (50-200 rows). A large background
#' data slows down the calculations, while a small data set leads to imprecise SHAP values.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param paired_sampling Logical flag indicating whether to use paired sampling.
#' The default is \code{TRUE}. This means that with every feature subset S,
#' also its complement is evaluated, which leads to considerably faster convergence.
#' @param m Number of feature subsets S to be evaluated during one iteration. 
#' The default, "auto", equals \code{trunc(20 * sqrt(ncol(X)))}. 
#' For the paired sampling strategy, 2m evaluations are done per iteration.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n) / diff(range(beta_n)) < tol, where sigma_n are the standard errors 
#' and beta_n are the SHAP values of a given observation.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, then the algorithm stops.
#' @param use_dt Logical flag indicating whether to use the "data.table" package 
#' for expensive data reshapings (default is \code{TRUE}). The flag is silently
#' ignored if \code{X} is a matrix or if the "data.table" package is not available. 
#' @param verbose Set to \code{FALSE} to suppress messages, warnings, and progress bar.
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
#' s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[, -1])
#' s
#' 
#' # Matrix input works as well, and pred_fun may contain preprocessing steps.
#' fit <- stats::lm(Sepal.Length ~ ., data = iris[1:4])
#' pred_fun <- function(X) stats::predict(fit, as.data.frame(X))
#' X <- data.matrix(iris[2:4])
#' s <- kernelshap(X[1:3, ], pred_fun = pred_fun, X)
#' s
kernelshap <- function(X, pred_fun, bg_X, bg_w = NULL, 
                       paired_sampling = TRUE, m = "auto", tol = 0.01, max_iter = 250, 
                       use_dt = TRUE, verbose = TRUE, ...) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.matrix(X) == is.matrix(bg_X),
    is.function(pred_fun),
    dim(X) >= 1L,
    ncol(bg_X) >= 1L,
    nrow(bg_X) >= 2L,
    !is.null(colnames(X)),
    !is.null(colnames(bg_X)),
    colnames(X) == colnames(bg_X)
  )
  bg_n <- nrow(bg_X)
  if (!is.null(bg_w)) {
    stopifnot(
      length(bg_w) == bg_n,
      all(bg_w >= 0),
      !all(bg_w == 0)
    )
  }
  if (verbose && bg_n > 1000L) {
    warning("Your background data 'bg_X' is large, which will slow down the process. Consider using 50-200 rows.")
  }
  if (verbose && bg_n < 10L) {
    warning("Your background data 'bg_X' is small, which might lead to imprecise SHAP values. Consider using 50-200 rows.")
  }
  bg_preds <- pred_fun(bg_X)
  stopifnot(
    "Predictions of multiple observations should be a vector (and not, e.g. a matrix)" = 
      is.vector(bg_preds),
    "Predictions should be numeric" = is.numeric(bg_preds), 
    length(bg_preds) == bg_n
  )
  use_dt <- use_dt &&!is.matrix(X) && requireNamespace("data.table", quietly = TRUE)

  # Initialization
  v0 <- weighted_mean(bg_preds, bg_w)  # Average prediction of background data
  v1 <- pred_fun(X)                    # Vector of predictions of X
  n <- nrow(X)
  nms <- colnames(X)
  if (m == "auto") {
    m <- trunc(20 * sqrt(ncol(X)))
  }

  # Handle simple exact case
  if (ncol(X) == 1L) {
    out <- list(
      S = matrix(v1 - v0, ncol = 1L, dimnames = list(NULL, nms)), 
      X = X, 
      baseline = v0, 
      SE = matrix(0, nrow = n, ncol = 1L, dimnames = list(NULL, nms)), 
      n_iter = integer(n), 
      converged = rep(TRUE, n)
    )
    class(out) <- "kernelshap"
    return(out)
  }

  # Real work: apply Kernel SHAP to each row of X
  if (verbose && n >= 2L) {
    pb <- utils::txtProgressBar(1L, n, style = 3)  
  }
  res <- vector("list", n)
  for (i in seq_len(n)) {
    res[[i]] <- kernelshap_one(
      x = X[i, , drop = FALSE], 
      pred_fun = pred_fun, 
      bg_X = bg_X, 
      bg_w = bg_w, 
      v0 = v0,
      v1 = v1[i],
      paired = paired_sampling,
      m = m,
      tol = tol,
      max_iter = max_iter,
      use_dt = use_dt
    )
    if (verbose && n >= 2L) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Organize output
  S <- do.call(rbind, lapply(res, `[[`, "beta"))
  colnames(S) <- nms
  SE <- do.call(rbind, lapply(res, `[[`, "sigma"))
  colnames(SE) <- nms
  n_iter <- vapply(res, `[[`, "n_iter", FUN.VALUE = integer(1L))
  converged <- vapply(res, `[[`, "converged", FUN.VALUE = logical(1L))
  if (verbose && !all(converged)) {
    warning("\nNon-convergence for ", sum(!converged), " rows.")
  }
  out <- list(
    S = S, X = X, baseline = v0, SE = SE, n_iter = n_iter, converged = converged
  )
  class(out) <- "kernelshap"
  out
}

# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(x, pred_fun, bg_X, bg_w, v0, v1, 
                           paired, m, tol, max_iter, use_dt) {
  X <- x[rep(1L, nrow(bg_X)), , drop = FALSE]
  p <- ncol(X)
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
      X,  
      bg = bg_X, 
      Z = if (paired) rbind(Z, 1 - Z) else Z, 
      pred_fun = pred_fun, 
      w = bg_w,
      use_dt = use_dt
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

# This function is by far the most expensive: data manipulation and prediction
# We further try to work with data.tables
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
