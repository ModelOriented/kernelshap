#' Permutation SHAP
#'
#' @description
#' Permutation SHAP algorithm with respect to a background dataset,
#' see Strumbelj and Kononenko (2014) for the basic idea.
#'
#' By default, for up to p=8 features, exact SHAP values are returned
#' (exact with respect to the selected background data).
#'
#' Otherwise, the sampling process iterates until the resulting values
#' are sufficiently precise, and standard errors are provided.
#'
#' During each iteration, the algorithm cycles twice through a random permutation:
#' It starts with all feature components "turned on" (i.e., taking them
#' from the observation to be explained), then gradually turning off components
#' according to the permutation (i.e., marginalizing them over the background data).
#' When all components are turned off, the algorithm - one by one - turns the components
#' back on, until all components are turned on again. This antithetic scheme allows to
#' evaluate Shapley's formula 2p times with each permutation, using a total of
#' 2p + 1 evaluations of marginal means.
#'
#' For models with interactions up to order two, one can show that
#' even a single iteration provides exact SHAP values (with respect to the
#' given background dataset).
#'
#' The Python implementation "shap" uses a similar approach, but without
#' providing standard errors, and without early stopping. To mimic its behavior,
#' we would need to set `max_iter = p` in R, and `max_eval = (2*p+1)*p` in Python.
#'
#' For faster convergence, we use balanced permutations in the sense that
#' p subsequent permutations each start with a different feature.
#' Furthermore, the 2p on-off vectors with sum <=1 or >=p-1 are evaluated only once,
#' similar to the degree 1 hybrid in [kernelshap()] (but covering less weight).
#'
#' @param exact If `TRUE`, the algorithm will produce exact SHAP values
#'   with respect to the background data.
#'   The default is `TRUE` for up to eight features, and `FALSE` otherwise.
#' @param low_memory If `FALSE` (default up to p = 15), the algorithm evaluates p
#'   predictions together, reducing the number of calls to `predict()`.
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
#'   - `m_exact`: Number of on-off vectors evaluated once per row of `X`.
#'   - `exact`: Logical flag indicating whether calculations are exact or not.
#'   - `txt`: Summary text.
#'   - `predictions`: \eqn{(n \times K)} matrix with predictions of `X`.
#'   - `algorithm`: "permshap".
#'   - `m`: Number of sampled on-off vectors evaluated per iteration (if not exact).
#'   - `SE`: Standard errors corresponding to `S` (if not exact).
#'   - `n_iter`: Integer vector of length n providing the number of iterations
#'     per row of `X` (if not exact).
#'   - `converged`: Logical vector of length n indicating convergence per row of `X`
#'     (if not exact).
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
    exact = length(feature_names) <= 8L,
    low_memory = length(feature_names) > 15L,
    tol = 0.01,
    max_iter = 10L * length(feature_names),
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    ...) {
  p <- length(feature_names)
  if (p <= 1L) {
    stop("Case p = 1 not implemented. Use kernelshap() instead.")
  }
  if (exact && p > 14L) {
    stop("Exact permutation SHAP only supported for up to 14 features")
  }
  if (!exact && p < 4L) {
    stop("Sampling version of permutation SHAP only supported for p >= 4 features")
  }

  txt <- paste(if (exact) "Exact" else "Sampling version of", "permutation SHAP")
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
  v0 <- wcolMeans(bg_preds, w = bg_w) # Average pred of bg data: 1 x K
  v1 <- align_pred(pred_fun(object, X, ...)) # Predictions on X:        n x K

  # Pre-calculations that are identical for each row to be explained
  if (exact) {
    Z <- exact_Z(p, feature_names = feature_names)
    m_exact <- nrow(Z) - 2L # We won't evaluate vz for first and last row
    m_eval <- 0L # for consistency with sampling case
    precalc <- list(
      Z = Z,
      bg_X_rep = rep_rows(bg_X, rep.int(seq_len(bg_n), m_exact)),
      positions = positions_for_exact(Z),
      shapley_w = shapley_weights(p, ell = rowSums(Z) - 1) # how many other players?
    )
  } else {
    max_iter <- as.integer(ceiling(max_iter / p) * p) # should be multiple of p
    m_exact <- 2L * p
    m <- 2L * (p - 3L) # inner loop
    m_eval <- if (low_memory) m else m * p # outer loop
    precalc <- list(
      bg_X_rep = rep_rows(bg_X, rep.int(seq_len(bg_n), m_eval)),
      bg_X_balanced = rep_rows(bg_X, rep.int(seq_len(bg_n), m_exact)),
      Z_balanced = exact_Z_balanced(p, feature_names)
    )
  }

  if (max(m_eval, m_exact) * bg_n > 2e5) {
    warning_burden(max(m_eval, m_exact), bg_n = bg_n)
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
      feature_names = feature_names,
      exact = exact,
      low_memory = low_memory,
      tol = tol,
      max_iter = max_iter,
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
        feature_names = feature_names,
        exact = exact,
        low_memory = low_memory,
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
  out <- list(
    S = reorganize_list(lapply(res, `[[`, "beta")),
    X = X,
    baseline = as.vector(v0),
    bg_X = bg_X,
    bg_w = bg_w,
    m_exact = m_exact,
    exact = exact,
    txt = txt,
    predictions = v1,
    algorithm = "permshap"
  )
  if (!exact) {
    out$m <- m
    out$SE <- reorganize_list(lapply(res, `[[`, "sigma"))
    out$n_iter <- vapply(res, `[[`, "n_iter", FUN.VALUE = integer(1L))
    out$converged <- vapply(res, `[[`, "converged", FUN.VALUE = logical(1L))

    if (verbose && !all(out$converged)) {
      warning("\nNon-convergence for ", sum(!out$converged), " rows.")
    }
  }
  if (verbose) {
    cat("\n")
  }
  class(out) <- "kernelshap"
  return(out)
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
    exact = length(feature_names) <= 8L,
    low_memory = length(feature_names) > 15L,
    tol = 0.01,
    max_iter = 10L * length(feature_names),
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    survival = c("chf", "prob"),
    ...) {
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
    exact = exact,
    low_memory = low_memory,
    tol = tol,
    max_iter = max_iter,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    ...
  )
}
