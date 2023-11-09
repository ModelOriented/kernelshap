#' Permutation SHAP
#'
#' Exact permutation SHAP values with respect to a background dataset.
#'
#' @inheritParams kernelshap
#' @returns
#'   An object of class "permshap" with the following components:
#'   - `S`: \eqn{(n \times p)} matrix with SHAP values or, if the model output has
#'     dimension \eqn{K > 1}, a list of \eqn{K} such matrices.
#'   - `X`: Same as input argument `X`.
#'   - `baseline`: Vector of length K representing the average prediction on the
#'     background data.
#' @export
#' @examples
#' # MODEL ONE: Linear regression
#' fit <- lm(Sepal.Length ~ ., data = iris)
#'
#' # Select rows to explain (only feature columns)
#' X_explain <- iris[1:2, -1]
#'
#' # Select small background dataset (could use all rows here because iris is small)
#' set.seed(1)
#' bg_X <- iris[sample(nrow(iris), 100), ]
#'
#' # Calculate SHAP values
#' s <- permshap(fit, X_explain, bg_X = bg_X)
#' s
#'
#' # MODEL TWO: Multi-response linear regression
#' fit <- lm(as.matrix(iris[1:2]) ~ Petal.Length + Petal.Width + Species, data = iris)
#' s <- permshap(fit, iris[1:4, 3:5], bg_X = bg_X)
#' s
#'
#' # Non-feature columns can be dropped via 'feature_names'
#' s <- permshap(
#'   fit,
#'   iris[1:4, ],
#'   bg_X = bg_X,
#'   feature_names = c("Petal.Length", "Petal.Width", "Species")
#' )
#' s
permshap <- function(object, ...) {
  UseMethod("permshap")
}

#' @describeIn permshap Default permutation SHAP method.
#' @export
permshap.default <- function(object, X, bg_X, pred_fun = stats::predict,
                             feature_names = colnames(X), bg_w = NULL,
                             parallel = FALSE, parallel_args = NULL,
                             verbose = TRUE, ...) {
  basic_checks(X = X, bg_X = bg_X, feature_names = feature_names, pred_fun = pred_fun)
  p <- length(feature_names)
  stopifnot(p <= 14L)
  n <- nrow(X)
  bg_n <- nrow(bg_X)
  if (!is.null(bg_w)) {
    bg_w <- prep_w(bg_w, bg_n = bg_n)
  }
  if (is.matrix(X) && !identical(colnames(X), feature_names)) {
    stop("If X is a matrix, feature_names must equal colnames(X)")
  }
  
  if (verbose) {
    message("Exact permutation SHAP values")
  }
  
  # Baseline
  bg_preds <- align_pred(pred_fun(object, bg_X[, colnames(X), drop = FALSE], ...))
  v0 <- wcolMeans(bg_preds, bg_w)            # Average pred of bg data: 1 x K
  
  # Drop unnecessary columns in bg_X. If X is matrix, also column order is relevant
  # Predictions will never be applied directly to bg_X anymore
  if (!identical(colnames(bg_X), feature_names)) {
    bg_X <- bg_X[, feature_names, drop = FALSE]
  }
  
  # Precalculations that are identical for each row to be explained
  Z <- exact_Z(p, feature_names = feature_names, keep_extremes = TRUE)
  m_exact <- nrow(Z)
  precalc <- list(
    Z = Z,
    Z_code = rowpaste(Z),
    bg_X_rep = bg_X[rep(seq_len(bg_n), times = m_exact), , drop = FALSE]
  )
  
  if (m_exact * bg_n > 2e5) {
    warning("\nPredictions on large data sets with ", m_exact, "x", bg_n,
            " observations are being done.\n",
            "Consider reducing the computational burden (e.g. use smaller X_bg)")
  }
  
  # Apply permutation SHAP to each row of X
  if (isTRUE(parallel)) {
    parallel_args <- c(list(i = seq_len(n)), parallel_args)
    res <- do.call(foreach::foreach, parallel_args) %dopar% permshap_one(
      x = X[i, , drop = FALSE],
      object = object,
      pred_fun = pred_fun,
      bg_w = bg_w,
      precalc = precalc,
      ...
    )
  } else {
    if (verbose && n >= 2L) {
      pb <- utils::txtProgressBar(max = n, style = 3)
    }
    res <- vector("list", n)
    for (i in seq_len(n)) {
      res[[i]] <- permshap_one(
        x = X[i, , drop = FALSE],
        object = object,
        pred_fun = pred_fun,
        bg_w = bg_w,
        precalc = precalc,
        ...
      )
      if (verbose && n >= 2L) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }
  out <- list(S = reorganize_list(res), X = X, baseline = as.vector(v0))
  class(out) <- "permshap"
  out
}

#' @describeIn permshap Permutation SHAP method for "ranger" models, see Readme for an example.
#' @export
permshap.ranger <- function(object, X, bg_X,
                            pred_fun = function(m, X, ...) stats::predict(m, X, ...)$predictions,
                            feature_names = colnames(X),
                            bg_w = NULL, parallel = FALSE, parallel_args = NULL,
                            verbose = TRUE, ...) {
  permshap.default(
    object = object,
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    ...
  )
}

#' @describeIn permshap Permutation SHAP method for "mlr3" models, see Readme for an example.
#' @export
permshap.Learner <- function(object, X, bg_X,
                             pred_fun = NULL,
                             feature_names = colnames(X),
                             bg_w = NULL, parallel = FALSE, parallel_args = NULL,
                             verbose = TRUE, ...) {
  if (is.null(pred_fun)) {
    pred_fun <- mlr3_pred_fun(object, X = X)
  }
  permshap.default(
    object = object,
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    ...
  )
}