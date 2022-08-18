#' Kernel SHAP
#'
#' Implements a multidimensional version of the Kernel SHAP algorithm explained in detail in 
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
#' This stopping criterion was suggested in Covert and Lee (2021). In the multioutput case,
#' the criterion must be fulfilled for each dimension separately until iteration stops.
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
#' @param verbose Set to \code{FALSE} to suppress messages, warnings, and progress bar.
#' @param ... Currently unused.
#' @return An object of class "kernelshap" with the following components:
#' \itemize{
#'   \item \code{S}: (n x p) matrix with SHAP values or, if the model output has dimension K > 1,
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
#' # Multioutput regression (or probabilistic classification)
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
kernelshap <- function(X, pred_fun, bg_X, bg_w = NULL, paired_sampling = TRUE, 
                       m = "auto", tol = 0.01, max_iter = 250, verbose = TRUE, ...) {
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
    stopifnot(length(bg_w) == bg_n, all(bg_w >= 0), !all(bg_w == 0))
  }
  if (verbose && bg_n > 1000L) {
    warning("Your background data 'bg_X' is large, which will slow down the process. Consider using 50-200 rows.")
  }
  if (verbose && bg_n < 10L) {
    warning("Your background data 'bg_X' is small, which might lead to imprecise SHAP values. Consider using 50-200 rows.")
  }
  bg_preds <- check_pred(pred_fun(bg_X), n = bg_n)

  # Initialization
  n <- nrow(X)
  v0 <- weighted_colMeans(bg_preds, bg_w)  # Average pred of background data: 1 x K
  v1 <- check_pred(pred_fun(X), n = n)     # Vector of predictions of X:      n x K
  nms <- colnames(X)
  if (m == "auto") {
    m <- trunc(20 * sqrt(ncol(X)))
  }
  bg_Xm <- bg_X[rep(seq_len(bg_n), times = m * (1L + paired_sampling)), , drop = FALSE]
  
  # Handle simple exact case
  if (ncol(X) == 1L) {
    S <- v1 - v0[rep(1L, n), , drop = FALSE]
    SE <- matrix(numeric(n), dimnames = list(NULL, nms))
    if (ncol(v1) > 1L) {
      SE <- replicate(ncol(v1), SE, simplify = FALSE)
      S <- lapply(
        asplit(S, MARGIN = 2L), function(M) as.matrix(M, dimnames = list(NULL, nms))
      )
    } else {
      colnames(S) <- nms      
    }
    out <- list(
      S = S, 
      X = X, 
      baseline = as.vector(v0), 
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
      X = X[rep(i, times = nrow(bg_Xm)), , drop = FALSE], 
      pred_fun = pred_fun, 
      bg_X = bg_Xm, 
      bg_w = bg_w, 
      v0 = v0,
      v1 = v1[i, , drop = FALSE],
      paired = paired_sampling,
      m = m,
      tol = tol,
      max_iter = max_iter
    )
    if (verbose && n >= 2L) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Organize output
  converged <- vapply(res, `[[`, "converged", FUN.VALUE = logical(1L))
  if (verbose && !all(converged)) {
    warning("\nNon-convergence for ", sum(!converged), " rows.")
  }
  out <- list(
    S = reorganize_list(lapply(res, `[[`, "beta"), nms = nms), 
    X = X, 
    baseline = as.vector(v0), 
    SE = reorganize_list(lapply(res, `[[`, "sigma"), nms = nms), 
    n_iter = vapply(res, `[[`, "n_iter", FUN.VALUE = integer(1L)),
    converged = converged
  )
  class(out) <- "kernelshap"
  out
}

# Kernel SHAP algorithm for a single row x with paired sampling
kernelshap_one <- function(X, pred_fun, bg_X, bg_w, v0, v1, 
                           paired, m, tol, max_iter) {
  p <- ncol(X)
  est_m = list()
  converged <- FALSE
  n_iter <- 0L
  Asum <- matrix(0, nrow = p, ncol = p)
  bsum <- numeric(p)
  n_Z <- m * (1L + paired)
  v0_ext <- v0[rep(1L, n_Z), , drop = FALSE]                      #  (n_Z x K)
  
  while(!isTRUE(converged) && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    Z <- make_Z(m, p)                                             #  (p x m)
    if (paired) {
      Z <- cbind(Z, 1 - Z)                                        #  (p x n_Z)
    }
    
    # Calling get_vz() is expensive                               
    vz <- get_vz(X, bg = bg_X, Z, pred_fun = pred_fun, w = bg_w)  # (n_Z x K)
    
    Atemp <- tcrossprod(Z) / n_Z                                  #  (p x p)
    btemp <- Z %*% (vz - v0_ext) / n_Z                            #  (p x K)
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
