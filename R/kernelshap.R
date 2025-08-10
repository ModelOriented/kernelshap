#' Kernel SHAP
#'
#' @description
#' Efficient implementation of Kernel SHAP, see Lundberg and Lee (2017), and
#' Covert and Lee (2021), abbreviated by CL21.
#' By default, for up to p=8 features, exact SHAP values are returned
#' (with respect to the selected background data).
#' Otherwise, a partly exact hybrid algorithm combining exact calculations and
#' iterative paired sampling is used, see Details.
#'
#' @details
#' The pure iterative Kernel SHAP sampling as in Covert and Lee (2021) works like this:
#'
#' 1. A binary "on-off" vector \eqn{z} is drawn from \eqn{\{0, 1\}^p} according to
#'   a special weighting logic.
#' 2. For each \eqn{j} with \eqn{z_j = 1}, the \eqn{j}-th column of the
#'   original background data is replaced by the corresponding feature value \eqn{x_j}
#'   of the observation to be explained.
#' 3. The average prediction \eqn{v_z} on the data of Step 2 is calculated, and the
#'   average prediction \eqn{v_0} on the background data is subtracted.
#' 4. Steps 1 to 3 are repeated \eqn{m} times. This produces a binary \eqn{m \times p}
#'   matrix \eqn{Z} (each row equals one of the \eqn{z}) and a vector \eqn{v} of
#'   shifted predictions.
#' 5. \eqn{v} is regressed onto \eqn{Z} under the constraint that the sum of the
#'   coefficients equals \eqn{v_1 - v_0}, where \eqn{v_1} is the prediction of the
#'   observation to be explained. The resulting coefficients are the Kernel SHAP values.
#'
#' This is repeated multiple times until convergence, see CL21 for details.
#'
#' To avoid the re-evaluation of identical coalition vectors, we have implemented
#' a hybrid strategy, combining exact calculations with sampling.
#'
#' The hybrid algorithm has two steps:
#' 1. Step 1 (exact part): There are \eqn{2p} different on-off vectors \eqn{z} with
#'   \eqn{\sum z \in \{1, p-1\}}.
#'   The degree 1 hybrid will list those vectors and use them according
#'   to their weights in the upcoming calculations. Depending on \eqn{p}, we can also go
#'   a step further to a degree 2 hybrid by adding all \eqn{p(p-1)} vectors with
#'   \eqn{\sum z \in \{2, p-2\}} to the process etc. The necessary predictions are
#'   obtained along with other calculations similar to those described in CL21.
#' 2. Step 2 (sampling part): The remaining weight is filled by sampling vectors z
#'   according to Kernel SHAP weights normalized to the values not yet covered by Step 1.
#'   Together with the results from Step 1 - correctly weighted - this now forms a
#'   complete iteration as in CL21. The difference is that a significant part of the mass
#'   is covered by exact calculations. Afterwards, the algorithm iterates until
#'   convergence. The output of Step 1 is reused in every iteration.
#'
#' If \eqn{p} is sufficiently small, all possible \eqn{2^p-2} on-off vectors \eqn{z} can be
#' evaluated. In this case, no sampling is required and the algorithm returns exact
#' Kernel SHAP values with respect to the given background data.
#' Since [kernelshap()] calculates predictions on data with \eqn{MN} rows
#' (\eqn{N} is the background data size and \eqn{M} the number of \eqn{z} vectors), \eqn{p}
#' should not be higher than 10 for exact calculations.
#' For similar reasons, degree 2 hybrids should not use \eqn{p} larger than 40.
#'
#' @importFrom doFuture %dofuture%
#'
#' @param object Fitted model object.
#' @param X \eqn{(n \times p)} matrix or `data.frame` with rows to be explained.
#'   The columns should only represent model features, not the response
#'   (but see `feature_names` on how to overrule this).
#' @param bg_X Background data used to integrate out "switched off" features,
#'   often a subset of the training data (typically 50 to 500 rows).
#'   In cases with a natural "off" value (like MNIST digits),
#'   this can also be a single row with all values set to the off value.
#'   If no `bg_X` is passed (the default) and if `X` is sufficiently large,
#'   a random sample of `bg_n` rows from `X` serves as background data.
#' @param pred_fun Prediction function of the form `function(object, X, ...)`,
#'   providing \eqn{K \ge 1} predictions per row. Its first argument
#'   represents the model `object`, its second argument a data structure like `X`.
#'   Additional (named) arguments are passed via `...`.
#'   The default, [stats::predict()], will work in most cases.
#' @param feature_names Optional vector of column names in `X` used to calculate
#'   SHAP values. By default, this equals `colnames(X)`.
#' @param bg_w Optional vector of case weights for each row of `bg_X`.
#'   If `bg_X = NULL`, must be of same length as `X`. Set to `NULL` for no weights.
#' @param bg_n If `bg_X = NULL`: Size of background data to be sampled from `X`.
#' @param exact If `TRUE`, the algorithm will produce exact SHAP values
#'   with respect to the background data.
#'   The default is `TRUE` for up to eight features, and `FALSE` otherwise.
#' @param hybrid_degree Integer controlling the exactness of the hybrid strategy. For
#'   \eqn{4 \le p \le 16}, the default is 2, otherwise it is 1.
#'   Ignored if `exact = TRUE`.
#'   - `0`: Pure sampling strategy not involving any exact part. It is strictly
#'     worse than the hybrid strategy and should therefore only be used for
#'     studying properties of the Kernel SHAP algorithm.
#'   - `1`: Uses all \eqn{2p} on-off vectors \eqn{z} with \eqn{\sum z \in \{1, p-1\}}
#'     for the exact part. The remaining mass is covered by random sampling.
#'   - `2`: Uses all \eqn{p(p+1)} on-off vectors \eqn{z} with
#'     \eqn{\sum z \in \{1, 2, p-2, p-1\}}. The remaining mass is covered by sampling.
#'     Usually converges fast.
#'   - `k>2`: Uses all on-off vectors with
#'     \eqn{\sum z \in \{1, \dots, k, p-k, \dots, p-1\}}.
#' @param m Even number of on-off vectors sampled during one iteration.
#'   The default is \eqn{2p}, except when `hybrid_degree == 0`.
#'   Then it is set to \eqn{8p}. Ignored if `exact = TRUE`.
#' @param tol Tolerance determining when to stop. As in CL21, the algorithm keeps
#'   iterating until \eqn{\textrm{max}(\sigma_n)/(\textrm{max}(\beta_n) - \textrm{min}(\beta_n)) < \textrm{tol}},
#'   where the \eqn{\beta_n} are the SHAP values of a given observation,
#'   and \eqn{\sigma_n} their standard errors.
#'   For multidimensional predictions, the criterion must be satisfied for each
#'   dimension separately. The stopping criterion uses the fact that standard errors
#'   and SHAP values are all on the same scale. Ignored if `exact = TRUE`.
#'   For `permshap()`, the default is 0.01, while for `kernelshap()` it is set to 0.005.
#' @param max_iter If the stopping criterion (see `tol`) is not reached after
#'   `max_iter` iterations, the algorithm stops. Ignored if `exact = TRUE`.
#' @param parallel If `TRUE`, use [foreach::foreach()] and `%dofuture%` to loop over rows
#'   to be explained. Must register backend beforehand, e.g., `plan(multisession)`,
#'   see README for an example. Currently disables the progress bar.
#' @param parallel_args Named list of arguments passed to
#'   `foreach::foreach(.options.future = ...)`, ideally `NULL` (default).
#'   Only relevant if `parallel = TRUE`.
#'   Example on Windows: if `object` is a GAM fitted with package 'mgcv',
#'   then one might need to set `parallel_args = list(packages = "mgcv")`.
#'   Similarly, if the model has been fitted with `ranger()`, then it might be necessary
#'   to pass `parallel_args = list(packages = "ranger")`.
#' @param verbose Set to `FALSE` to suppress messages and the progress bar.
#' @param seed Optional integer random seed. Note that it changes the global seed.
#' @param survival Should cumulative hazards ("chf", default) or survival
#'   probabilities ("prob") per time be predicted? Only in `ranger()` survival models.
#' @param ... Additional arguments passed to `pred_fun(object, X, ...)`.
#' @returns
#'   An object of class "kernelshap" with the following components:
#'   - `S`: \eqn{(n \times p)} matrix with SHAP values or, if the model output has
#'     dimension \eqn{K > 1}, a list of \eqn{K} such matrices.
#'   - `X`: Same as input argument `X`.
#'   - `baseline`: Vector of length K representing the average prediction on the
#'     background data.
#'   - `bg_X`: The background data.
#'   - `bg_w`: The background case weights.
#'   - `m_exact`: Number of on-off vectors evaluated for exact calculations.
#'   - `prop_exact`: Proportion of the Kernel SHAP weight distribution covered by
#'     exact calculations.
#'   - `exact`: Logical flag indicating whether calculations are exact or not.
#'   - `txt`: Summary text.
#'   - `predictions`: \eqn{(n \times K)} matrix with predictions of `X`.
#'   - `algorithm`: "kernelshap".
#'   - `m`: Number of sampled on-off vectors evaluated per iteration (if not exact).
#'   - `SE`: Standard errors corresponding to `S` (if not exact).
#'   - `n_iter`: Integer vector of length n providing the number of iterations
#'     per row of `X` (if not exact).
#'   - `converged`: Logical vector of length n indicating convergence per row of `X`
#'     (if not exact).
#' @references
#'   1. Scott M. Lundberg and Su-In Lee. A unified approach to interpreting model
#'     predictions. Proceedings of the 31st International Conference on Neural
#'     Information Processing Systems, 2017.
#'   2. Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value
#'     Estimation Using Linear Regression. Proceedings of The 24th International
#'     Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
#' @export
#' @examples
#' # MODEL ONE: Linear regression
#' fit <- lm(Sepal.Length ~ ., data = iris)
#'
#' # Select rows to explain (only feature columns)
#' X_explain <- iris[-1]
#'
#' # Calculate SHAP values
#' s <- kernelshap(fit, X_explain)
#' s
#'
#' # MODEL TWO: Multi-response linear regression
#' fit <- lm(as.matrix(iris[, 1:2]) ~ Petal.Length + Petal.Width + Species, data = iris)
#' s <- kernelshap(fit, iris[3:5])
#' s
#'
#' # Note 1: Feature columns can also be selected 'feature_names'
#' # Note 2: Especially when X is small, pass a sufficiently large background data bg_X
#' s <- kernelshap(
#'   fit,
#'   iris[1:4, ],
#'   bg_X = iris,
#'   feature_names = c("Petal.Length", "Petal.Width", "Species")
#' )
#' s
kernelshap <- function(object, ...) {
  UseMethod("kernelshap")
}

#' @describeIn kernelshap Default Kernel SHAP method.
#' @export
kernelshap.default <- function(
    object,
    X,
    bg_X = NULL,
    pred_fun = stats::predict,
    feature_names = colnames(X),
    bg_w = NULL,
    bg_n = 200L,
    exact = length(feature_names) <= 8L,
    hybrid_degree = 1L + length(feature_names) %in% 4:16,
    m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)),
    tol = 0.005,
    max_iter = 100L,
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    seed = NULL,
    ...) {
  p <- length(feature_names)
  basic_checks(X = X, feature_names = feature_names, pred_fun = pred_fun)
  stopifnot(
    exact %in% c(TRUE, FALSE),
    p == 1L || exact || hybrid_degree %in% 0:(p / 2),
    "m must be even" = trunc(m / 2) == m / 2
  )
  prep_bg <- prepare_bg(X = X, bg_X = bg_X, bg_n = bg_n, bg_w = bg_w, verbose = verbose)
  bg_X <- prep_bg$bg_X
  bg_w <- prep_bg$bg_w
  bg_n <- nrow(bg_X)
  n <- nrow(X)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Calculate v1 and v0
  bg_preds <- align_pred(pred_fun(object, bg_X, ...))
  v0 <- wcolMeans(bg_preds, bg_w) # Average pred of bg data: 1 x K
  v1 <- align_pred(pred_fun(object, X, ...)) # Predictions on X:        n x K

  # For p = 1, exact Shapley values are returned
  if (p == 1L) {
    out <- case_p1(
      n = n, feature_names = feature_names, v0 = v0, v1 = v1, X = X, verbose = verbose
    )
    return(out)
  }

  txt <- summarize_strategy(p, exact = exact, deg = hybrid_degree)
  if (verbose) {
    message(txt)
  }

  # Pre-calculations that are identical for each row to be explained
  if (exact || hybrid_degree >= 1L) {
    if (exact) {
      precalc <- input_exact(p, feature_names = feature_names)
    } else {
      precalc <- input_partly_exact(
        p = p, deg = hybrid_degree, feature_names = feature_names
      )
    }
    Z <- precalc[["Z"]]
    m_exact <- nrow(Z)
    prop_exact <- sum(precalc[["w"]])
    precalc[["bg_exact_rep"]] <- rep_rows(bg_X, rep.int(seq_len(bg_n), m_exact))
    g <- rep_each(m_exact, each = bg_n)
    precalc[["Z_exact_rep"]] <- Z[g, , drop = FALSE]
  } else {
    precalc <- list()
    m_exact <- 0L
    prop_exact <- 0
  }
  if (!exact) {
    precalc[["bg_sampling_rep"]] <- rep_rows(bg_X, rep.int(seq_len(bg_n), m))
  }

  if (max(m, m_exact) * bg_n > 2e5) {
    warning_burden(max(m, m_exact), bg_n = bg_n)
  }

  # Apply Kernel SHAP to each row of X
  if (isTRUE(parallel)) {
    future_args <- c(list(seed = TRUE), parallel_args)
    parallel_args <- c(list(i = seq_len(n)), list(.options.future = future_args))
    res <- do.call(foreach::foreach, parallel_args) %dofuture% kernelshap_one(
      x = X[i, , drop = FALSE],
      v1 = v1[i, , drop = FALSE],
      object = object,
      pred_fun = pred_fun,
      feature_names = feature_names,
      bg_w = bg_w,
      exact = exact,
      deg = hybrid_degree,
      m = m,
      tol = tol,
      max_iter = max_iter,
      v0 = v0,
      precalc = precalc,
      bg_n = bg_n,
      ...
    )
  } else {
    if (verbose && n >= 2L) {
      pb <- utils::txtProgressBar(max = n, style = 3)
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
        m = m,
        tol = tol,
        max_iter = max_iter,
        v0 = v0,
        precalc = precalc,
        bg_n = bg_n,
        ...
      )
      if (verbose && n >= 2L) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }

  # Organize output
  exact <- exact || trunc(p / 2) == hybrid_degree

  out <- list(
    S = reorganize_list(lapply(res, `[[`, "beta")),
    X = X,
    baseline = as.vector(v0),
    bg_X = bg_X,
    bg_w = bg_w,
    m_exact = m_exact,
    prop_exact = prop_exact,
    exact = exact,
    txt = txt,
    predictions = v1,
    algorithm = "kernelshap"
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

#' @describeIn kernelshap Kernel SHAP method for "ranger" models, see Readme for an example.
#' @export
kernelshap.ranger <- function(
    object,
    X,
    bg_X = NULL,
    pred_fun = NULL,
    feature_names = colnames(X),
    bg_w = NULL,
    bg_n = 200L,
    exact = length(feature_names) <= 8L,
    hybrid_degree = 1L + length(feature_names) %in% 4:16,
    m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)),
    tol = 0.005,
    max_iter = 100L,
    parallel = FALSE,
    parallel_args = NULL,
    verbose = TRUE,
    seed = NULL,
    survival = c("chf", "prob"),
    ...) {
  if (is.null(pred_fun)) {
    pred_fun <- create_ranger_pred_fun(object$treetype, survival = match.arg(survival))
  }

  kernelshap.default(
    object = object,
    X = X,
    bg_X = bg_X,
    pred_fun = pred_fun,
    feature_names = feature_names,
    bg_w = bg_w,
    bg_n = bg_n,
    exact = exact,
    hybrid_degree = hybrid_degree,
    m = m,
    tol = tol,
    max_iter = max_iter,
    parallel = parallel,
    parallel_args = parallel_args,
    verbose = verbose,
    seed = seed,
    ...
  )
}
