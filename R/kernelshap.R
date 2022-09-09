#' Kernel SHAP
#'
#' Implements a multidimensional version of the Kernel SHAP algorithm explained in 
#' detail in Covert and Lee (2021). It is an iterative refinement of the original 
#' Kernel SHAP algorithm of Lundberg and Lee (2017). The algorithm is applied to each 
#' row in \code{X}. Its default behaviour depends on the number of features p:
#' \itemize{
#'   \item 2 <= p <= 8: Exact Kernel SHAP values are returned. 
#'   (Exact regarding the given background data.)
#'   \item p > 8: Sampling version of Kernel SHAP. 
#'   The algorithm iterates until Kernel SHAP values are sufficiently accurate. 
#'   Approximate standard errors of the SHAP values are returned. 
#'   To improve convergence, a hybrid (partly-exact) sampling scheme to produce "on-off" 
#'   vectors z is used.
#'   \item p = 1: Exact Shapley values are returned.
#' }
#' 
#' During each iteration, \code{m} on-off vectors (feature subsets) are evaluated until 
#' the worst standard error of the SHAP values is small enough relative to the range of 
#' the SHAP values. This stopping criterion was suggested in Covert and Lee (2021) and 
#' uses the fact that SHAP values and their standard errors are all on the scale of the 
#' predictions. In the multi-output case, the criterion must be fulfilled for each 
#' dimension separately until iteration stops.
#' 
#' @importFrom doRNG %dorng%
#' 
#' @param object Fitted model object.
#' @param X A (n x p) matrix, data.frame, tibble or data.table of rows to be explained. 
#' Important: The columns should only represent model features, not the response.
#' @param bg_X Background data used to integrate out "switched off" features, 
#' often a subset of the training data (around 100 to 200 rows)
#' It should contain the same columns as \code{X}. Columns not in \code{X} are silently 
#' dropped and the columns are arranged into the order as they appear in \code{X}.
#' @param pred_fun Prediction function of the form \code{function(object, X, ...)},
#' providing K >= 1 numeric predictions per row. Its first argument represents the
#' model \code{object}, its second argument a data structure like \code{X}. 
#' (The names of the first two arguments do not matter.) Additional (named)
#' arguments are passed via \code{...}. The default, \code{stats::predict}, will
#' work in most cases. Some exceptions (classes "ranger" and mlr3 "Learner")
#' are handled separately. In other cases, the function must be specified manually.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param sampling_strategy Sampling strategy to draw m binary "on-off" vectors z of 
#' length p used to calculate Kernel SHAP. Currently, one of:
#' \itemize{
#'   \item \code{"auto"} (default): "exact" for up to 8 features, and "hybrid" otherwise.
#'   \item \code{"simple"}: Each z is drawn randomly according to Kernel SHAP weights.
#'   This is done in two steps: First, the value s of sum(z) drawn from the Kernel weight
#'   distribution (normalized to the range from 1 to p-1), see Covert and Lee (2021). 
#'   Then, s random positions of z are set to 1.
#'   Note that this strategy is strictly worse than paired sampling, so there is 
#'   virtually no reason to use it except for studying properties of Kernel SHAP.
#'   \item \code{"paired"}: Here, m/2 vectors z are drawn as with the "simple" strategy. 
#'   Then, all m/2 complements 1-z are added to the pool of vectors. 
#'   This paired sampling strategy converges faster than the simple method, 
#'   see Covert and Lee (2021).
#'   \item \code{"exact"}: All possible binary vectors z are enumerated and weighted
#'   according to the SHAP kernel weight distribution. This produces exact Kernel SHAP
#'   values with respect to the given background data. The algorithm works with large
#'   prediction data having (2^p-2) * nrow(bg_X) rows. Thus, we recommend this option
#'   up to p=10.
#'   \item \code{"hybrid"}: Sampling is done partly exact and partly via sampling.
#'   Kernel SHAP weighting puts >=75% mass to "on-off" vectors z with sum(z) either 1 
#'   or p-1. (Similarly, it puts >=92.7% mass to vectors with sum(z) <= 2 or >= p-2.)
#'   For large p, the probabilities are even much closer to 1.
#'   Thus, random sampling from the Kernel weight distribution yields many identical z 
#'   with sum(z) = 1 or sum(z) = p-1. This leaves only a few other z. 
#'   The hybrid strategy improves this inefficient behaviour: it first creates all 2p 
#'   vectors z with sum(z) equals to 1 or p-1 ("hybrid degree 1"). 
#'   Then, if m is large enough, it adds all p(1-p) vectors z with sum(z) = 2 or 
#'   sum(z) = p-2 ("degree 2"). Finally, the remaining 
#'   vectors are sampled according to the paired strategy. Convergence of this strategy 
#'   is expected to be much faster than of paired sampling.
#' }
#' @param paired_sampling Deprecated, set \code{sampling_strategy = "paired"}.
#' @param exact Deprecated, set \code{sampling_strategy = "exact"}.
#' @param m Number of on-off vectors to be evaluated during one iteration. 
#' The default, \code{NULL}, equals \code{max(30, pp*(pp+5), 4*p)}, where p is the
#' number of features and pp is min(p, 14). Ignored if exact calculations are done.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n) / diff(range(beta_n)) < tol, where the beta_n are the SHAP values 
#' of a given observation and sigma_n their standard errors. For multidimensional
#' predictions, the criterion must be satisfied for each dimension separately.
#' The stopping criterion uses the fact that standard errors and SHAP values are all
#' on the same scale. Ignored if exact calculations are done.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, the algorithm stops. Ignored if exact calculations are done.
#' @param parallel If \code{TRUE}, use parallel \code{foreach::foreach()} to loop over rows
#' to be explained. Must register backend beforehand, e.g. via "doFuture" package, 
#' see Readme for an example. Parallelization automatically disables the progress bar.
#' @param parallel_args A named list of arguments passed to \code{foreach::foreach()}, see
#' \code{?foreach::foreach}. Ideally, this is \code{NULL} (default). Only relevant
#' if \code{parallel = TRUE}. Example on Windows: if \code{object} is a generalized
#' additive model fitted with package "mgcv", then one might need to set
#' \code{parallel_args = list(.packages = "mgcv")}.
#' @param verbose Set to \code{FALSE} to suppress messages and the progress bar.
#' @param ... Additional arguments passed to \code{pred_fun(object, X, ...)}.
#' @return An object of class "kernelshap" with the following components:
#' \itemize{
#'   \item \code{S}: (n x p) matrix with SHAP values or, if the model output has dimension K > 1,
#'   a list of K such matrices.
#'   \item \code{X}: Same as input argument \code{X}.
#'   \item \code{baseline}: A vector of length K representing the average prediction on the background data.
#'   \item \code{SE}: Standard errors corresponding to \code{S} (and organized like \code{S}).
#'   \item \code{n_iter}: Integer vector of length n providing the number of iterations per row of \code{X}.
#'   \item \code{converged}: Logical vector of length n indicating convergence per row of \code{X}.
#'   \item \code{m}: Integer providing the effective number of on-off vectors used per iteration.
#' }
#' @references
#' \enumerate{
#'   \item Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
#'   \item Scott M. Lundberg and Su-In Lee. A Unified Approach to Interpreting Model Predictions. Advances in Neural Information Processing Systems 30, 2017.
#'}
#' @export
#' @examples
#' # Linear regression
#' fit <- stats::lm(Sepal.Length ~ ., data = iris)
#' s <- kernelshap(fit, iris[1:2, -1], bg_X = iris)
#' s
#' 
#' # Multivariate model
#' fit <- stats::lm(
#'   as.matrix(iris[1:2]) ~ Petal.Length + Petal.Width + Species, data = iris
#' )
#' s <- kernelshap(fit, iris[1:4, 3:5], bg_X = iris)
#' s
#'
#' # Matrix input works as well, and pred_fun can be overwritten
#' fit <- stats::lm(Sepal.Length ~ ., data = iris[1:4])
#' pred_fun <- function(fit, X) stats::predict(fit, as.data.frame(X))
#' X <- data.matrix(iris[2:4])
#' s <- kernelshap(fit, X[1:3, ], bg_X = X, pred_fun = pred_fun)
#' s
#'
#' # Logistic regression
#' fit <- stats::glm(
#'   I(Species == "virginica") ~ Sepal.Length + Sepal.Width, 
#'   data = iris, 
#'   family = binomial
#' )
#' 
#' # On scale of linear predictor
#' s <- kernelshap(fit, iris[1:2], bg_X = iris)
#' s
#' 
#' # On scale of response (probability)
#' s <- kernelshap(fit, iris[1:2], bg_X = iris, type = "response")
#' s
#' 
kernelshap <- function(object, ...){
  UseMethod("kernelshap")
}

#' @describeIn kernelshap Default Kernel SHAP method.
#' @export
kernelshap.default <- function(object, X, bg_X, pred_fun = stats::predict, bg_w = NULL, 
                               sampling_strategy = c("auto", "hybrid", "exact", "paired", "simple"),
                               paired_sampling = NULL, exact = NULL,
                               m = NULL, tol = 0.01, max_iter = 250, parallel = FALSE, 
                               parallel_args = NULL, verbose = TRUE, ...) {
  sampling_strategy <- match.arg(sampling_strategy)
  if (!is.null(paired_sampling)) {
    warning("The argument 'paired_sampling' is deprecated. It is not used anymore 
    and will be removed in version 0.4.0. Use sampling_strategy = 'paired' instead.")
  }
  if (!is.null(exact)) {
    warning("The argument 'exact' is deprecated. It is not used anymore 
    and will be removed in version 0.4.0. Use sampling_strategy = 'exact' instead.")
  }
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
  bg_preds <- check_pred(pred_fun(object, bg_X, ...), n = bg_n)
  v0 <- weighted_colMeans(bg_preds, bg_w)            # Average pred of bg data: 1 x K
  v1 <- check_pred(pred_fun(object, X, ...), n = n)  # Predictions on X:        n x K
  
  # For p = 1, exact Shapley values are returned
  if (p == 1L) {
    if (verbose) {
      message("Calculating exact Shapley values")
    }
    return(case_p1(n = n, nms = nms, v0 = v0, v1 = v1, X = X))
  }
  
  # Select strategy "auto"
  if (sampling_strategy == "auto") {
    sampling_strategy <- if (p <= 8) "exact" else "hybrid"
  }
  
  m <- fix_m(sampling_strategy, p = p, m = m)
  
  if (sampling_strategy == "exact" && p > 10L) {
    warning("Exact calculations with more than ten features might take long because 
    prediction data sets have about 2^p * nrow(bg_X) rows. Consider setting exact = FALSE")  
  }
  
  if (sampling_strategy == "hybrid" && hybrid_degree(p = p, m = m) == 0L) {
    stop("m is too small for the hybrid strategy, please make it larger")
  }
  
  if (verbose) {
    message(summarize_strategy(sampling_strategy, p = p, m = m))
  }
  
  # Precalculate certain objects that are identical for each row of X
  bg_Xm <- bg_X[rep(seq_len(bg_n), times = m), , drop = FALSE]
  precalc <- switch(
    sampling_strategy,
    exact = input_exact(p),
    hybrid = all_pairs(p),
    NULL
  )
  
  # Real work: apply Kernel SHAP to each row of X
  if (isTRUE(parallel)) {
    parallel_args <- c(list(i = seq_len(n)), parallel_args)
    res <- do.call(foreach::foreach, parallel_args) %dorng% kernelshap_one(
      object = object,
      X = X[rep(i, times = nrow(bg_Xm)), , drop = FALSE],
      bg_X = bg_Xm,
      pred_fun = pred_fun,
      bg_w = bg_w, 
      v0 = v0,
      v1 = v1[i, , drop = FALSE],
      sampling_strategy = sampling_strategy,
      m = m,
      precalc = precalc,
      tol = tol,
      max_iter = max_iter,
      ...
    )
  } else {
    if (verbose && n >= 2L) {
      pb <- utils::txtProgressBar(1L, n, style = 3)  
    }
    res <- vector("list", n)
    for (i in seq_len(n)) {
      res[[i]] <- kernelshap_one(
        object = object,
        X = X[rep(i, times = nrow(bg_Xm)), , drop = FALSE],
        bg_X = bg_Xm,
        pred_fun = pred_fun,
        bg_w = bg_w, 
        v0 = v0,
        v1 = v1[i, , drop = FALSE],
        sampling_strategy = sampling_strategy,
        m = m,
        precalc = precalc,
        tol = tol,
        max_iter = max_iter,
        ...
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
    converged = converged,
    m = m
  )
  class(out) <- "kernelshap"
  out
}

#' @describeIn kernelshap Kernel SHAP method for "ranger" models, see Readme for an example.
#' @export
kernelshap.ranger <- function(object, X, bg_X,
                              pred_fun = function(m, X, ...) stats::predict(m, X, ...)$predictions, 
                              bg_w = NULL, 
                              sampling_strategy = c("auto", "hybrid", "exact", "paired", "simple"),
                              paired_sampling = NULL, exact = NULL,
                              m = NULL, tol = 0.01, max_iter = 250, parallel = FALSE, 
                              parallel_args = NULL, verbose = TRUE, ...) {
  kernelshap.default(
    object = object, 
    X = X, 
    bg_X = bg_X, 
    pred_fun = pred_fun, 
    bg_w = bg_w, 
    sampling_strategy = sampling_strategy,
    paired_sampling = paired_sampling,
    exact = exact,
    m = m, 
    tol = tol, 
    max_iter = max_iter,
    parallel = parallel,
    parallel_args = parallel_args, 
    verbose = verbose, 
    ...
  )
}

#' @describeIn kernelshap Kernel SHAP method for "mlr3" models, see Readme for an example.
#' @export
kernelshap.Learner <- function(object, X, bg_X,
                               pred_fun = function(m, X) m$predict_newdata(X)$response, 
                               bg_w = NULL, 
                               sampling_strategy = c("auto", "hybrid", "exact", "paired", "simple"),
                               paired_sampling = NULL, exact = NULL,
                               m = NULL, tol = 0.01, max_iter = 250, parallel = FALSE, 
                               parallel_args = NULL, verbose = TRUE, ...) {
  kernelshap.default(
    object = object, 
    X = X, 
    bg_X = bg_X, 
    pred_fun = pred_fun, 
    bg_w = bg_w, 
    sampling_strategy = sampling_strategy,
    paired_sampling = paired_sampling,
    exact = exact,
    m = m, 
    tol = tol, 
    max_iter = max_iter,
    parallel = parallel,
    parallel_args = parallel_args, 
    verbose = verbose, 
    ...
  )
}

