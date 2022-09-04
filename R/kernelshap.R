#' Kernel SHAP
#'
#' Implements a multidimensional version of the Kernel SHAP algorithm explained in detail in 
#' Covert and Lee (2021). It is an iterative refinement of the original Kernel SHAP algorithm 
#' of Lundberg and Lee (2017).
#' The algorithm is applied to each row in \code{X}. Its behaviour depends on the number of features p:
#' \itemize{
#'   \item 2 <= p <= 5: Exact Kernel SHAP values are returned. (Exact regarding the given background data.)
#'   \item p > 5: Sampling version of Kernel SHAP. The algorithm iterates until Kernel SHAP values are sufficiently accurate. 
#'   Approximate standard errors of the SHAP values are returned. 
#'   \item p = 1: Exact Shapley values are returned.
#' }
#' 
#' \code{X} should only contain feature columns required by the
#' prediction function \code{pred_fun}. The latter is a function taking
#' a data structure like \code{X} or \code{bg_X} and providing K >= 1 numeric 
#' predictions per row. The background data \code{bg_X} must contain the same column names
#' as \code{X} (additional columns are silently dropped).
#' During each iteration, \code{m} feature subsets are evaluated until the worst 
#' standard error of the SHAP values is small enough relative to the range of the SHAP values. 
#' This stopping criterion was suggested in Covert and Lee (2021). In the multioutput case,
#' the criterion must be fulfilled for each dimension separately until iteration stops.
#' @param X A (n x p) matrix, data.frame, tibble or data.table of rows to be explained. 
#' Important: The columns should only represent model features, not the response.
#' @param pred_fun A function that takes a data structure like \code{X} 
#' and provides K >= 1 numeric predictions per row.
#' Example: If "fit" denotes a logistic regression fitted via \code{stats::glm}, 
#' and SHAP values should be on the probability scale, then this argument is
#' \code{function(X) predict(fit, X, type = "response")}.
#' @param bg_X The background data used to integrate out "switched off" features. 
#' It should contain the same columns as \code{X}. A good size is around 50 to 200 rows.
#' Columns not in \code{X} are silently dropped and the columns are arranged into
#' the order as they appear in \code{X}.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param paired_sampling Logical flag indicating whether to use paired sampling.
#' The default is \code{TRUE}. This means that with every feature subset S,
#' also its complement is evaluated, which leads to considerably faster convergence.
#' @param m Number of feature subsets S to be evaluated during one iteration. 
#' The default, "auto", equals \code{max(trunc(20*sqrt(p)), 5*p)}, where p is the
#' number of features. 
#' For the paired sampling strategy, 2m evaluations are done per iteration.
#' @param exact If \code{TRUE} (default) and the number of features p is at most 5,
#' the algorithm will produce exact Kernel SHAP values. In this case, the arguments
#' \code{m}, \code{paired_sampling}, \code{tol}, and \code{max_iter} are ignored.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n) / diff(range(beta_n)) < tol, where the beta_n are the SHAP values 
#' of a given observation and sigma_n their standard errors. For multidimensional
#' predictions, the criterion must be satisfied for each dimension separately.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, then the algorithm stops.
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to loop over rows
#' to be explained. Must register backend beforehand, e.g. via \code{doMC}. See
#' example below. Parallelization automatically disables the progress bar.
#' @param verbose Set to \code{FALSE} to suppress messages, warnings, and the progress bar.
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
#' @import doRNG
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
#' # In parallel
#' require(doMC)
#' registerDoMC(cores = 4)
#' system.time(kernelshap(iris[, -1], pred_fun = pred_fun, iris[, -1]))
#' system.time(kernelshap(iris[, -1], pred_fun = pred_fun, iris[, -1], parallel = TRUE))
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
kernelshap <- function(X, pred_fun, bg_X, bg_w = NULL, 
                       paired_sampling = TRUE, m = "auto", exact = TRUE, 
                       tol = 0.01, max_iter = 250, parallel = FALSE, 
                       verbose = TRUE, seed = NULL, ...) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.matrix(X) == is.matrix(bg_X),
    (n <- nrow(X)) >= 1L,
    (bg_n <- nrow(bg_X)) >= 2L,
    (p <- ncol(X)) >= 1L,
    ncol(bg_X) >= 1L,
    !is.null(nms <- colnames(X)),
    !is.null(colnames(bg_X)),
    all(nms %in% colnames(bg_X)),
    is.function(pred_fun)
  )
  if (!is.null(bg_w)) {
    stopifnot(length(bg_w) == bg_n, all(bg_w >= 0), !all(bg_w == 0))
  }
  if (verbose) {
    check_bg_size(bg_n)
  }
  
  # Calculate v0 and v1
  bg_X <- bg_X[, colnames(X), drop = FALSE]
  bg_preds <- check_pred(pred_fun(bg_X), n = bg_n)
  v0 <- weighted_colMeans(bg_preds, bg_w)  # Average pred of background data: 1 x K
  v1 <- check_pred(pred_fun(X), n = n)     # Predictions on X:                n x K
  
  # For p = 1, exact Shapley values are returned
  if (p == 1L) {
    if (verbose) {
      message("Calculating exact Shapley values")
    }
    return(case_p1(n = n, nms = nms, v0 = v0, v1 = v1, X = X))
  }
  
  # Calculate m, S, probs
  if (exact && p <= length(Z_exact)) {
    if (verbose) {
      message("Calculating exact Kernel SHAP values")
    }
    m <- nrow(Z_exact[[p]])
    paired_sampling <- FALSE  
  } else {
    if (verbose) {
      message("Calculating Kernel SHAP values by iterative sampling")
    }
    exact <- FALSE
  }
  if (m == "auto") {
    m <- max(trunc(20 * sqrt(p)), 5L * p)
  }
  S <- 1:(p - 1)
  probs <- (p - 1) / (choose(p, S) * S * (p - S))
  
  # Allocate replicated version of the background data
  bg_Xm <- bg_X[rep(seq_len(bg_n), times = m * (1L + paired_sampling)), , drop = FALSE]
  
  # Real work: apply Kernel SHAP to each row of X
  if (isTRUE(parallel)) {
    res <- foreach(i = seq_len(n)) %dorng% kernelshap_one(
      X = X[rep(i, times = nrow(bg_Xm)), , drop = FALSE],
      pred_fun = pred_fun, 
      bg_X = bg_Xm, 
      bg_w = bg_w, 
      v0 = v0,
      v1 = v1[i, , drop = FALSE],
      paired = paired_sampling,
      m = m,
      S = S,
      probs = probs,
      exact = exact,
      tol = tol,
      max_iter = max_iter
    )
  } else {
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
        S = S,
        probs = probs,
        exact = exact,
        tol = tol,
        max_iter = max_iter
      )
      if (verbose && n >= 2L) {
        utils::setTxtProgressBar(pb, i)
      }
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
                           paired, m, S, probs, exact, tol, max_iter) {
  p <- ncol(X)
  est_m = list()
  converged <- FALSE
  n_iter <- 0L
  Asum <- matrix(0, nrow = p, ncol = p)                           #  (p x p)
  bsum <- matrix(0, nrow = p, ncol = ncol(v0))                    #  (p x K)
  n_Z <- m * (1L + paired)
  v0_ext <- v0[rep(1L, n_Z), , drop = FALSE]                      #  (n_Z x K)
  
  while(!isTRUE(converged) && n_iter < max_iter) {
    n_iter <- n_iter + 1L
    
    # Create Z matrix of dimension ->                             #  (n_Z x p)
    if (exact) {
      Z <- Z_exact[[p]]
    } else {
      Z <- sample_Z(S, m, probs)
      if (paired) {
        Z <- rbind(Z, 1 - Z)
      }
    }
    
    # Calling get_vz() is expensive                               #  (n_Z x K)
    vz <- get_vz(X = X, bg = bg_X, Z = Z, pred_fun = pred_fun, w = bg_w)
    
    # Least-squares with constraint that beta_1 + ... + beta_p = v_1 - v_0. 
    # The additional constraint beta_0 = v_0 is dealt via offset
    Atemp <- crossprod(Z) / n_Z                                   #  (p x p)
    btemp <- crossprod(Z, (vz - v0_ext)) / n_Z                    #  (p x K)
    Asum <- Asum + Atemp                                          #  (p x p)
    bsum <- bsum + btemp                                          #  (p x K)
    est_m[[n_iter]] <- R <- solver(Atemp, btemp, v1, v0)          #  (p x K)
    
    if (exact) {
      return(list(beta = R, sigma = 0 * R, n_iter = 1L, converged = TRUE))
    }
    
    # Covariance calculation would fail in the first iteration
    if (n_iter >= 2L) {
      beta_n <- solver(Asum / n_iter, bsum / n_iter, v1, v0)      #  (p x K)
      sigma_n <- get_sigma(est_m, iter = n_iter)                  #  (p x K)
      converged <- all(conv_crit(sigma_n, beta_n) < tol)
    }
  }
  list(beta = beta_n, sigma = sigma_n, n_iter = n_iter, converged = converged)
}
