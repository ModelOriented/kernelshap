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
#' a data structure like \code{X} or \code{bg_X} and providing K >= 1 numeric 
#' predictions per row.
#' During each iteration, \code{m} feature subsets are evaluated until the worst 
#' standard error of the SHAP values is small enough relative to the range of the SHAP values. 
#' This stopping criterion was suggested in Covert and Lee (2021). We apply it to each
#' of the K dimensions of the predictions and stop if the stopping criterion fires for all of them.
#' @param X A (n x p) matrix or data.frame of rows to be explained. 
#' Important: The columns should only represent model features, not the response.
#' @param pred_fun A function that takes a data structure like \code{X} or \code{bg_X} 
#' and provides K >= 1 numeric prediction per row.
#' Example: If "fit" denotes a logistic regression fitted via \code{stats::glm}, 
#' and SHAP values should be on the probability scale, then this argument is
#' \code{function(X) predict(fit, X, type = "response")}.
#' @param bg_X The background data used to integrate out "switched off" features. 
#' It should have the same column structure as \code{X}. A good size is around $50-200$ rows.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param paired_sampling Logical flag indicating whether to use paired sampling.
#' The default is \code{TRUE}. This means that with every feature subset S,
#' also its complement is evaluated, which leads to considerably faster convergence.
#' @param m Number of feature subsets S to be evaluated during one iteration. 
#' The default, "auto", equals \code{trunc(20 * sqrt(ncol(X)))}. 
#' For the paired sampling strategy, 2m evaluations are done per iteration.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n) / diff(range(beta_n)) < tol, where the beta_n are the SHAP values 
#' of a given observation and sigma_n their standard errors. For multidimensional
#' predictions, the criterion must be satisfied for each dimension separately.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, then the algorithm stops.
#' @param use_dt Logical flag indicating whether to use the "data.table" package 
#' for expensive data reshapings (default is \code{TRUE}). The flag is silently
#' ignored if \code{X} is a matrix or if the "data.table" package is not available. 
#' @param verbose Set to \code{FALSE} to suppress messages, warnings, and progress bar.
#' @param ... Currently unused.
#' @return An object of class "kernelshap" with the following components:
#' \itemize{
#'   \item \code{S}: (n x p) matrix with SHAP values or, if the predictions are multidimensional,
#'   a list of K such matrices.
#'   \item \code{X}: Same as input argument \code{X}.
#'   \item \code{baseline}: A vector of length K representing the average prediction on the background data.
#'   \item \code{SE}: Standard errors corresponding to \code{S} (and organized like \code{S}).
#'   \item \code{n_iter}: Integer vector of length n providing the number of iterations per row of \code{X}.
#'   \item \code{converged}: Logical vector of length n indicating convergence per row of \code{X}.
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
#' # Multi-output regression or (probabilistic) classification
#' fit <- stats::lm(
#'   as.matrix(iris[1:2]) ~ Petal.Length + Petal.Width + Species, data = iris
#' )
#' fit
#' s <- kernelshap(iris[1:4, 3:5], pred_fun = pred_fun, iris[, 3:5])
#' s
#' 
#' # Matrix input works as well, and pred_fun may contain preprocessing steps
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
  bg_preds <- check_pred(pred_fun(bg_X), n = bg_n)
  use_dt <- use_dt &&!is.matrix(X) && requireNamespace("data.table", quietly = TRUE)

  # Initialization
  n <- nrow(X)
  v0 <- weighted_colMeans(bg_preds, bg_w)  # Average pred of background data: 1 x K
  v1 <- check_pred(pred_fun(X), n = n)     # Vector of predictions of X:      n x K
  nms <- colnames(X)
  if (m == "auto") {
    m <- trunc(20 * sqrt(ncol(X)))
  }

  # Handle simple exact case
  if (ncol(X) == 1L) {
    S <- matrix(v1 - v0, ncol = 1L, dimnames = list(NULL, nms))
    SE <- matrix(0, nrow = n, ncol = 1L, dimnames = list(NULL, nms))
    
    out <- list(
      S = S, 
      X = X, 
      baseline = v0, 
      SE = SE, 
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
      v1 = v1[i, , drop = FALSE],
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
  S <- reorganize_list(lapply(res, `[[`, "beta"), nms = nms)
  SE <- reorganize_list(lapply(res, `[[`, "sigma"), nms = nms)
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
  Asum <- matrix(0, nrow = p, ncol = p)
  bsum <- numeric(p)
  n_Z <- m * (1L + paired)
  v0_ext <- v0[rep(1L, times = n_Z), , drop = FALSE]             #  (n_Z x K)
  
  while(!isTRUE(converged) && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    Z <- make_Z(m, p)
    if (paired) {
      Z <- rbind(Z, 1 - Z)
    }
    
    # Calling get_vz() is expensive                                  (n_Z x K)
    vz <- get_vz(X, bg = bg_X, Z, pred_fun = pred_fun, w = bg_w, use_dt = use_dt)
    
    Atemp <- crossprod(Z) / n_Z                                   #  (p x p)
    btemp <- crossprod(Z, (vz - v0_ext)) / n_Z                    #  (p x K)
    Asum <- Asum + Atemp                                          #  (p x p)
    bsum <- bsum + btemp                                          #  (p x K)
    est_m[[n_iter]] <- solver(Atemp, btemp, v1, v0)               #  (p x K)
    
    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(Asum / n_iter, bsum / n_iter, v1, v0)      #  (p x K)
      sigma_n <- get_sigma(est_m, iter = n_iter)                  #  (p x K)
      converged <- all(conv_crit(sigma_n, beta_n) < tol)
    }
  }
  list(beta = beta_n, sigma = sigma_n, n_iter = n_iter, converged = converged)
}

# Calculates regression coefficients
solver <- function(A, b, v1, v0) {
  p <- ncol(A)
  Ainv <- solve(A)
  s <- (matrix(colSums(crossprod(Ainv, b)), nrow = 1L) - v1 + v0) / sum(Ainv)  # (1 x K)
  s <- s[rep(1L, times = p), , drop = FALSE]                                   # (p x K) 
  Ainv %*% (b - s)                                                             # (p x K)
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
# We try to work with matrices or data.tables
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

# Average per Z (optimized for speed)
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
  apply(sig, 2L, FUN = max) / apply(bet, 2L, FUN = function(z) diff(range(z)))
}

# Turn list of n (p x K) matrices into K lists of (n x p) matrices. Reduce if K = 1.
reorganize_list <- function(alist, nms) {
  out <- abind1(alist)
  dimnames(out) <- list(NULL, nms, NULL)
  out <- asplit(out, MARGIN = 3L)
  if (length(out) == 1L) {
    return(out[[1L]])
  }
  out
}

# Checks and reshapes predictions to matrix with n rows
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

