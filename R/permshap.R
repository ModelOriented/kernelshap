#' Permutation SHAP
#'
#' Exact permutation SHAP algorithm with respect to a background dataset,
#' see Strumbelj and Kononenko. The function works for up to 14 features.
#' For more than eight features, we recommend [kernelshap()] due to its higher speed.
#'
#' @inheritParams kernelshap
#' @returns
#'   An object of class "kernelshap" with the following components:
#'   - `S`: \eqn{(n \times p)} matrix with SHAP values or, if the model output has
#'     dimension \eqn{K > 1}, a list of \eqn{K} such matrices.
#'   - `X`: Same as input argument `X`.
#'   - `baseline`: Vector of length K representing the average prediction on the
#'     background data.
#'   - `bg_X`: The background data.
#'   - `bg_w`: The background case weights.
#'   - `m_exact`: Integer providing the effective number of exact on-off vectors used.
#'   - `exact`: Logical flag indicating whether calculations are exact or not
#'     (currently always `TRUE`).
#'   - `txt`: Summary text.
#'   - `predictions`: \eqn{(n \times K)} matrix with predictions of `X`.
#'   - `algorithm`: "permshap".
#' @references
#'   1. Erik Strumbelj and Igor Kononenko. Explaining prediction models and individual
#'     predictions with feature contributions. Knowledge and Information Systems 41, 2014.
#' @export
#' @examples
#' # MODEL ONE: Linear regression
#' fit <- lm(Sepal.Length ~ ., data = iris)
#'
#' # Select rows to explain (only feature columns)
#' X_explain <- iris[-1]
#'
#' # Calculate SHAP values
#' s <- permshap(fit, X_explain)
#' s
#'
#' # MODEL TWO: Multi-response linear regression
#' fit <- lm(as.matrix(iris[, 1:2]) ~ Petal.Length + Petal.Width + Species, data = iris)
#' s <- permshap(fit, iris[3:5])
#' s
#'
#' # Note 1: Feature columns can also be selected 'feature_names'
#' # Note 2: Especially when X is small, pass a sufficiently large background data bg_X
#' s <- permshap(
#'   fit,
#'   iris[1:4, ],
#'   bg_X = iris,
#'   feature_names = c("Petal.Length", "Petal.Width", "Species")
#' )
#' s
permshap <- function(object, ...) {
  UseMethod("permshap")
}

#' @describeIn permshap Default permutation SHAP method.
#' @export
permshap.default <- function(
    object,
    X,
    bg_X = NULL,
    pred_fun = stats::predict,
    feature_names = colnames(X),
    bg_w = NULL,
    bg_n = 200L,
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    ...
  ) {
  p <- length(feature_names)
  if (p <= 1L) {
    stop("Case p = 1 not implemented. Use kernelshap() instead.")
  }
  if (p > 14L) {
    stop("Permutation SHAP only supported for up to 14 features")
  }

  txt <- "Exact permutation SHAP"
  if (verbose) {
    message(txt)
  }

  basic_checks(X = X, feature_names = feature_names, pred_fun = pred_fun)
  prep_bg <- prepare_bg(X = X, bg_X = bg_X, bg_n = bg_n, bg_w = bg_w, verbose = verbose)
  bg_X <- prep_bg$bg_X
  bg_w <- prep_bg$bg_w
  bg_n <- nrow(bg_X)
  n <- nrow(X)

  # Baseline and predictions on explanation data
  bg_preds <- align_pred(pred_fun(object, bg_X, ...))
  v0 <- wcolMeans(bg_preds, w = bg_w)         # Average pred of bg data: 1 x K
  v1 <- align_pred(pred_fun(object, X, ...))  # Predictions on X:        n x K

  # Drop unnecessary columns in bg_X. If X is matrix, also column order is relevant
  # Predictions will never be applied directly to bg_X anymore
  if (!identical(colnames(bg_X), feature_names)) {
    bg_X <- bg_X[, feature_names, drop = FALSE]
  }

  # Precalculations that are identical for each row to be explained
  Z <- exact_Z(p, feature_names = feature_names, keep_extremes = TRUE)
  m_exact <- nrow(Z) - 2L  # We won't evaluate vz for first and last row
  precalc <- list(
    Z = Z,
    Z_code = rowpaste(Z),
    bg_X_rep = rep_rows(bg_X, rep.int(seq_len(bg_n), m_exact))
  )

  if (m_exact * bg_n > 2e5) {
    warning_burden(m_exact, bg_n = bg_n)
  }

  # Apply permutation SHAP to each row of X
  if (isTRUE(parallel)) {
    parallel_args <- c(list(i = seq_len(n)), parallel_args)
    res <- do.call(foreach::foreach, parallel_args) %dopar% permshap_one(
      x = X[i, , drop = FALSE],
      v1 = v1[i, , drop = FALSE],
      object = object,
      pred_fun = pred_fun,
      bg_w = bg_w,
      v0 = v0,
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
        v1 = v1[i, , drop = FALSE],
        object = object,
        pred_fun = pred_fun,
        bg_w = bg_w,
        v0 = v0,
        precalc = precalc,
        ...
      )
      if (verbose && n >= 2L) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }
  if (verbose) {
    cat("\n")
  }
  out <- list(
    S = reorganize_list(res),
    X = X,
    baseline = as.vector(v0),
    bg_X = bg_X,
    bg_w = bg_w,
    m_exact = m_exact,
    exact = TRUE,
    txt = txt,
    predictions = v1,
    algorithm = "permshap"
  )
  class(out) <- "kernelshap"
  out
}

#' @describeIn permshap Permutation SHAP method for "ranger" models, see Readme for an example.
#' @export
permshap.ranger <- function(
    object,
    X,
    bg_X = NULL,
    pred_fun = NULL,
    feature_names = colnames(X),
    bg_w = NULL,
    bg_n = 200L,
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    survival = c("chf", "prob"),
    ...
  ) {

  if (is.null(pred_fun)) {
    pred_fun <- create_ranger_pred_fun(object$treetype, survival = match.arg(survival))
  }

  permshap.default(
    object = object,
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    bg_n = bg_n,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    ...
  )
}

