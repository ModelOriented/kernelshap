#' Kernel SHAP
#'
#' Multidimensional refinement of the Kernel SHAP Algorithm described in Covert and Lee (2021), 
#' in the following abbreviated by "CL21". 
#' The function allows to calculate Kernel SHAP values in an exact way, by iterative sampling 
#' as in CL21, or by a hybrid of these two options. As soon as sampling is involved, 
#' the algorithm iterates until convergence, and standard errors are provided.
#' The default behaviour depends on the number of features p, see also Details below:
#' \itemize{
#'   \item 2 <= p <= 8: Exact Kernel SHAP values are returned (for the given background data). 
#'   \item p > 8: Hybrid (partly exact) iterative version of Kernel SHAP
#'   (degree 2 up to p = 16, degree 1 for p > 16, see Details).
#'   \item p = 1: Exact Shapley values are returned (independent of the background data).
#' }
#' 
#' The iterative Kernel SHAP sampling algorithm (CL21) works by randomly sample 
#' m on-off vectors z so that their sum follows the SHAP Kernel weight distribution 
#' (renormalized to the range from 1 to p-1). Based on these vectors, many predictions 
#' are formed. Then, Kernel SHAP values are derived as the solution of a constrained 
#' linear regression. This is done multiple times until convergence, see CL21 for details.
#' 
#' A drawback of this strategy is that many (at least 75%) of the z vectors will have 
#' sum(z) equal to 1 or p-1, producing many duplicates. Similarly, at least 92% of 
#' the mass will be used for the p(p+1) possible vectors with sum(z) in 1, 2, p-1, p-2. 
#' This inefficiency can be fixed by a hybrid strategy, combining exact calculations with sampling.
#' 
#' The hybrid algorithm has two steps:
#' \enumerate{
#'   \item Step 1 (exact part): There are 2p different on-off vectors z with sum(z) equals to 
#'   1 or p-1, covering a large proportion of the Kernel SHAP distribution. 
#'   The degree 1 hybrid will list those vectors and use them according to their weights 
#'   in the upcoming calculations. Depending on p, we can also go a step further to 
#'   a degree 2 hybrid by adding all p(p-1) vectors with sum(z) equals to 2 or p-2
#'   to the process etc. The necessary predictions are obtained along with other 
#'   calculations similar to those described in CL21.
#'   \item Step 2 (sampling part): The remaining weight is filled by sampling vectors z
#'   according to Kernel SHAP weights renormalized to the values not yet covered by Step 1. 
#'   Together with the results from Step 1 - correctly weighted - this now forms a
#'   complete iteration as in CL21. The difference is that most mass is covered by exact calculations. 
#'   Afterwards, the algorithm iterates until convergence. The output of Step 1 is reused
#'   in every iteration, leading to an extremely efficient strategy.
#' }
#' 
#' If p is sufficiently small, all possible 2^p-2 on-off vectors z can be evaluated.
#' In this case, no sampling is required and the algorithm returns exact Kernel SHAP values 
#' with respect to the given background data. Since \code{kernelshap()} calculates predictions 
#' on data with MN rows (N is the background data size and M the number of z vectors),
#' p should not be much higher than 10 for exact calculations. 
#' For similar reasons, degree 2 hybrids are limited to p up to 30-40.
#' 
#' @importFrom doRNG %dorng%
#' 
#' @param object Fitted model object.
#' @param X A (n x p) matrix, data.frame, tibble or data.table of rows to be explained. 
#' Important: The columns should only represent model features, not the response.
#' @param bg_X Background data used to integrate out "switched off" features, 
#' often a subset of the training data (typically 50 to 500 rows)
#' It should contain the same columns as \code{X}.
#' In cases with a natural "off" value (like MNIST digits), 
#' this can also be a single row with all values set to the off value.
#' @param pred_fun Prediction function of the form \code{function(object, X, ...)},
#' providing K >= 1 numeric predictions per row. Its first argument represents the
#' model \code{object}, its second argument a data structure like \code{X} and \code{bg_X}. 
#' (The names of the first two arguments do not matter.) Additional (named)
#' arguments are passed via \code{...}. The default, \code{stats::predict}, will
#' work in most cases. Some exceptions (classes "ranger" and mlr3 "Learner")
#' are handled separately. In other cases, the function must be specified manually.
#' @param feature_names Optional vector of column names in \code{X} used to calculate 
#' SHAP values. By default, this equals \code{colnames(X)}. Not supported for matrix
#' \code{X}.
#' @param bg_w Optional vector of case weights for each row of \code{bg_X}.
#' @param exact If \code{TRUE}, the algorithm will produce exact Kernel SHAP values
#' with respect to the background data. In this case, the arguments \code{hybrid_degree}, 
#' \code{m}, \code{paired_sampling}, \code{tol}, and \code{max_iter} are ignored.
#' The default is \code{TRUE} up to eight features, and \code{FALSE} otherwise. 
#' @param hybrid_degree Integer controlling the exactness of the hybrid strategy. For
#' 4 <= p <= 16, the default is 2, otherwise it is 1. Ignored if \code{exact = TRUE}.
#' \itemize{
#'   \item \code{0}: Pure sampling strategy not involving any exact part. It is strictly
#'   worse than the hybrid strategy and should therefore only be used for 
#'   studying properties of the Kernel SHAP algorithm.
#'   \item \code{1}: Uses all 2p on-off vectors z with sum(z) equal to 1 or p-1 for the exact 
#'   part, which covers at least 75% of the mass of the Kernel weight distribution. 
#'   The remaining mass is covered by sampling.
#'   \item \code{2}: Uses all p(p+1) on-off vectors z with sum(z) equal to 1, p-1, 2, or p-2. 
#'   This covers at least 92% of the mass of the Kernel weight distribution. 
#'   The remaining mass is covered by sampling. Convergence usually happens in the 
#'   minimal possible number of iterations of two.
#'   \item \code{k>2}: Uses all on-off vectors with sum(z) in 1,...,k and p-1,...,p-k.
#' }
#' @param paired_sampling Logical flag indicating whether to do the sampling in a paired
#' manner. This means that with every on-off vector z, also 1-z is considered.
#' CL21 shows its superiority compared to standard sampling, therefore the 
#' default (\code{TRUE}) should usually not be changed except for studying properties
#' of Kernel SHAP algorithms. Ignored if \code{exact = TRUE}.
#' @param m Even number of on-off vectors sampled during one iteration. 
#' The default is 2p, except when \code{hybrid_degree == 0}. Then it is set to 8p. 
#' Ignored if \code{exact = TRUE}.
#' @param tol Tolerance determining when to stop. The algorithm keeps iterating until
#' max(sigma_n)/diff(range(beta_n)) < tol, where the beta_n are the SHAP values 
#' of a given observation and sigma_n their standard errors, see CL21. For multidimensional
#' predictions, the criterion must be satisfied for each dimension separately.
#' The stopping criterion uses the fact that standard errors and SHAP values are all
#' on the same scale. Ignored if \code{exact = TRUE}.
#' @param max_iter If the stopping criterion (see \code{tol}) is not reached after 
#' \code{max_iter} iterations, the algorithm stops. Ignored if \code{exact = TRUE}.
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
#'   \item \code{m}: Integer providing the effective number of sampled on-off vectors used per iteration.
#'   \item \code{m_exact}: Integer providing the effective number of exact on-off vectors used per iteration.
#'   \item \code{prop_exact}: Proportion of the Kernel SHAP weight distribution covered by exact calculations.
#'   \item \code{exact}: Logical flag indicating whether calculations are exact or not.
#'   \item \code{txt}: Summary text.
#' }
#' @references
#' \enumerate{
#'   \item Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
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
#' summary(s)
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
#' # Non-feature columns can be dropped via 'feature_names'
#' fit <- stats::lm(Sepal.Length ~ . - Species, data = iris)
#' s <- kernelshap(
#'   fit, 
#'   iris[1:2, ], 
#'   bg_X = iris, 
#'   feature_names = c("Sepal.Width", "Petal.Length", "Petal.Width")
#' )
#' s
kernelshap <- function(object, ...){
  UseMethod("kernelshap")
}

#' @describeIn kernelshap Default Kernel SHAP method.
#' @export
kernelshap.default <- function(object, X, bg_X, pred_fun = stats::predict, 
                               feature_names = colnames(X), bg_w = NULL, 
                               exact = length(feature_names) <= 8L, 
                               hybrid_degree = 1L + length(feature_names) %in% 4:16, 
                               paired_sampling = TRUE, 
                               m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)), 
                               tol = 0.005, max_iter = 100L, parallel = FALSE, 
                               parallel_args = NULL, verbose = TRUE, ...) {
  stopifnot(
    is.matrix(X) || is.data.frame(X),
    is.matrix(bg_X) || is.data.frame(bg_X),
    is.matrix(X) == is.matrix(bg_X),
    dim(X) >= 1L,
    dim(bg_X) >= 1L,
    !is.null(colnames(X)),
    !is.null(colnames(bg_X)),
    (p <- length(feature_names)) >= 1L,
    all(feature_names %in% colnames(X)),
    all(feature_names %in% colnames(bg_X)),
    is.function(pred_fun),
    exact %in% c(TRUE, FALSE),
    p == 1L || exact || hybrid_degree %in% 0:(p / 2),
    paired_sampling %in% c(TRUE, FALSE),
    "m must be even" = trunc(m / 2) == m / 2
  )
  n <- nrow(X)
  bg_n <- nrow(bg_X)
  if (!is.null(bg_w)) {
    stopifnot(length(bg_w) == bg_n, all(bg_w >= 0), !all(bg_w == 0))
  }
  if (is.matrix(X) && !identical(colnames(X), feature_names)) {
    stop("If X is a matrix, feature_names must equal colnames(X)")  
  }
  
  # Calculate v0 and v1
  bg_preds <- check_pred(pred_fun(object, bg_X, ...), n = bg_n)
  v0 <- weighted_colMeans(bg_preds, bg_w)            # Average pred of bg data: 1 x K
  v1 <- check_pred(pred_fun(object, X, ...), n = n)  # Predictions on X:        n x K
  
  # For p = 1, exact Shapley values are returned
  if (p == 1L) {
    return(
      case_p1(n = n, nms = feature_names, v0 = v0, v1 = v1, X = X, verbose = verbose)
    )
  }
  
  # Drop unnecessary columns in bg_X. If X is matrix, also column order is relevant
  if (!identical(colnames(bg_X), feature_names)) {
    bg_X <- bg_X[, feature_names, drop = FALSE]
  }
  
  # Precalculations for the real Kernel SHAP
  if (exact || hybrid_degree >= 1L) {
    precalc <- if (exact) input_exact(p) else input_partly_exact(p, hybrid_degree)
    m_exact <- nrow(precalc[["Z"]])
    prop_exact <- sum(precalc[["w"]])
    precalc[["bg_X_exact"]] <- bg_X[rep(seq_len(bg_n), times = m_exact), , drop = FALSE]
  } else {
    precalc <- list()
    m_exact <- 0L
    prop_exact <- 0
  }
  if (!exact) {
    precalc[["bg_X_m"]] <- bg_X[rep(seq_len(bg_n), times = m), , drop = FALSE]  
  }
  
  # Some infos
  txt <- summarize_strategy(p, exact = exact, deg = hybrid_degree)
  if (verbose) {
    message(txt)
  }
  if (verbose && max(m, m_exact) * bg_n > 2e5) {
    warning("\nPredictions on large data sets with ", max(m, m_exact), "x", bg_n,
            " observations are being done. Consider reducing the computational burden ",
            "(e.g. exact = FALSE, low hybrid_degree, smaller background data, smaller m)")
  }
  
  # Apply Kernel SHAP to each row of X
  if (isTRUE(parallel)) {
    parallel_args <- c(list(i = seq_len(n)), parallel_args)
    res <- do.call(foreach::foreach, parallel_args) %dorng% kernelshap_one(
      x = X[i, , drop = FALSE], 
      v1 = v1[i, , drop = FALSE], 
      object = object,
      pred_fun = pred_fun,
      feature_names = feature_names,
      bg_w = bg_w, 
      exact = exact,
      deg = hybrid_degree,
      paired = paired_sampling,
      m = m,
      tol = tol,
      max_iter = max_iter,
      v0 = v0,
      precalc = precalc,
      ...
    )
  } else {
    if (verbose && n >= 2L) {
      pb <- utils::txtProgressBar(1L, n, style = 3)  
    }
    res <- vector("list", n)
    for (i in seq_len(n)) {
      res[[i]] <- kernelshap_one(
        x = X[i, , drop = FALSE], 
        v1 = v1[i, , drop = FALSE], 
        object = object,
        pred_fun = pred_fun,
        feature_names = feature_names,
        bg_w = bg_w, 
        exact = exact,
        deg = hybrid_degree,
        paired = paired_sampling,
        m = m,
        tol = tol,
        max_iter = max_iter,
        v0 = v0,
        precalc = precalc,
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
    S = reorganize_list(lapply(res, `[[`, "beta"), nms = feature_names), 
    X = X, 
    baseline = as.vector(v0), 
    SE = reorganize_list(lapply(res, `[[`, "sigma"), nms = feature_names), 
    n_iter = vapply(res, `[[`, "n_iter", FUN.VALUE = integer(1L)),
    converged = converged,
    m = m,
    m_exact = m_exact,
    prop_exact = prop_exact,
    exact = exact || trunc(p / 2) == hybrid_degree,
    txt = txt
  )
  class(out) <- "kernelshap"
  out
}

#' @describeIn kernelshap Kernel SHAP method for "ranger" models, see Readme for an example.
#' @export
kernelshap.ranger <- function(object, X, bg_X,
                              pred_fun = function(m, X, ...) stats::predict(m, X, ...)$predictions,
                              feature_names = colnames(X), 
                              bg_w = NULL, exact = length(feature_names) <= 8L, 
                              hybrid_degree = 1L + length(feature_names) %in% 4:16, 
                              paired_sampling = TRUE, 
                              m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)), 
                              tol = 0.005, max_iter = 100L, parallel = FALSE, 
                              parallel_args = NULL, verbose = TRUE, ...) {
  kernelshap.default(
    object = object, 
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    exact = exact,
    hybrid_degree = hybrid_degree,
    paired_sampling = paired_sampling,
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
                               feature_names = colnames(X),
                               bg_w = NULL, exact = length(feature_names) <= 8L,
                               hybrid_degree = 1L + length(feature_names) %in% 4:16,
                               paired_sampling = TRUE,
                               m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)),
                               tol = 0.005, max_iter = 100L, parallel = FALSE,
                               parallel_args = NULL, verbose = TRUE, ...) {
  kernelshap.default(
    object = object,
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    exact = exact,
    hybrid_degree = hybrid_degree,
    paired_sampling = paired_sampling,
    m = m,
    tol = tol,
    max_iter = max_iter,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    ...
  )
}

