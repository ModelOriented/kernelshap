#' Kernel SHAP
#'
#' Efficient implementation of Kernel SHAP, see Lundberg and Lee (2017), and
#' Covert and Lee (2021), abbreviated by CL21.
#' For up to \eqn{p=8} features, the resulting Kernel SHAP values are exact regarding
#' the selected background data. For larger \eqn{p}, an almost exact
#' hybrid algorithm involving iterative sampling is used, see Details.
#'
#' Pure iterative Kernel SHAP sampling as in Covert and Lee (2021) works like this:
#'
#' 1. A binary "on-off" vector \eqn{z} is drawn from \eqn{\{0, 1\}^p}
#'   such that its sum follows the SHAP Kernel weight distribution
#'   (normalized to the range \eqn{\{1, \dots, p-1\}}).
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
#' A drawback of this strategy is that many (at least 75%) of the \eqn{z} vectors will
#' have \eqn{\sum z \in \{1, p-1\}}, producing many duplicates. Similarly, at least 92%
#' of the mass will be used for the \eqn{p(p+1)} possible vectors with
#' \eqn{\sum z \in \{1, 2, p-2, p-1\}}.
#' This inefficiency can be fixed by a hybrid strategy, combining exact calculations
#' with sampling.
#'
#' The hybrid algorithm has two steps:
#' 1. Step 1 (exact part): There are \eqn{2p} different on-off vectors \eqn{z} with
#'   \eqn{\sum z \in \{1, p-1\}}, covering a large proportion of the Kernel SHAP
#'   distribution. The degree 1 hybrid will list those vectors and use them according
#'   to their weights in the upcoming calculations. Depending on \eqn{p}, we can also go
#'   a step further to a degree 2 hybrid by adding all \eqn{p(p-1)} vectors with
#'   \eqn{\sum z \in \{2, p-2\}} to the process etc. The necessary predictions are
#'   obtained along with other calculations similar to those described in CL21.
#' 2. Step 2 (sampling part): The remaining weight is filled by sampling vectors z
#'   according to Kernel SHAP weights renormalized to the values not yet covered by Step 1.
#'   Together with the results from Step 1 - correctly weighted - this now forms a
#'   complete iteration as in CL21. The difference is that most mass is covered by exact
#'   calculations. Afterwards, the algorithm iterates until convergence.
#'   The output of Step 1 is reused in every iteration, leading to an extremely
#'   efficient strategy.
#'
#' If \eqn{p} is sufficiently small, all possible \eqn{2^p-2} on-off vectors \eqn{z} can be
#' evaluated. In this case, no sampling is required and the algorithm returns exact
#' Kernel SHAP values with respect to the given background data.
#' Since [permshap()] calculates predictions on data with \eqn{MN} rows
#' (\eqn{N} is the background data size and \eqn{M} the number of \eqn{z} vectors), \eqn{p}
#' should not be much higher than 10 for exact calculations.
#' For similar reasons, degree 2 hybrids should not use \eqn{p} much larger than 40.
#'
#' @importFrom foreach %dopar%
#'
#' @param object Fitted model object.
#' @param X \eqn{(n \times p)} matrix or `data.frame` with rows to be explained.
#'   The columns should only represent model features, not the response
#'   (but see `feature_names` on how to overrule this).
#' @param bg_X Background data used to integrate out "switched off" features,
#'   often a subset of the training data (typically 50 to 500 rows)
#'   It should contain the same columns as `X`.
#'   In cases with a natural "off" value (like MNIST digits),
#'   this can also be a single row with all values set to the off value.
#' @param pred_fun Prediction function of the form `function(object, X, ...)`,
#'   providing \eqn{K \ge 1} numeric predictions per row. Its first argument
#'   represents the model `object`, its second argument a data structure like `X`.
#'   Additional (named) arguments are passed via `...`.
#'   The default, [stats::predict()], will work in most cases.
#' @param feature_names Optional vector of column names in `X` used to calculate
#'   SHAP values. By default, this equals `colnames(X)`. Not supported if `X`
#'   is a matrix.
#' @param bg_w Optional vector of case weights for each row of `bg_X`.
#' @param exact If `TRUE`, the algorithm will produce exact Kernel SHAP values
#'   with respect to the background data. In this case, the arguments `hybrid_degree`,
#'   `m`, `paired_sampling`, `tol`, and `max_iter` are ignored.
#'   The default is `TRUE` up to eight features, and `FALSE` otherwise.
#' @param hybrid_degree Integer controlling the exactness of the hybrid strategy. For
#'   \eqn{4 \le p \le 16}, the default is 2, otherwise it is 1.
#'   Ignored if `exact = TRUE`.
#'   - `0`: Pure sampling strategy not involving any exact part. It is strictly
#'     worse than the hybrid strategy and should therefore only be used for
#'     studying properties of the Kernel SHAP algorithm.
#'   - `1`: Uses all \eqn{2p} on-off vectors \eqn{z} with \eqn{\sum z \in \{1, p-1\}}
#'     for the exact part, which covers at least 75% of the mass of the Kernel weight
#'     distribution. The remaining mass is covered by random sampling.
#'   - `2`: Uses all \eqn{p(p+1)} on-off vectors \eqn{z} with
#'     \eqn{\sum z \in \{1, 2, p-2, p-1\}}. This covers at least 92% of the mass of the
#'     Kernel weight distribution. The remaining mass is covered by sampling.
#'     Convergence usually happens in the minimal possible number of iterations of two.
#'   - `k>2`: Uses all on-off vectors with
#'     \eqn{\sum z \in \{1, \dots, k, p-k, \dots, p-1\}}.
#' @param paired_sampling Logical flag indicating whether to do the sampling in a paired
#'   manner. This means that with every on-off vector \eqn{z}, also \eqn{1-z} is
#'   considered. CL21 shows its superiority compared to standard sampling, therefore the
#'   default (`TRUE`) should usually not be changed except for studying properties
#'   of Kernel SHAP algorithms. Ignored if `exact = TRUE`.
#' @param m Even number of on-off vectors sampled during one iteration.
#'   The default is \eqn{2p}, except when `hybrid_degree == 0`.
#'   Then it is set to \eqn{8p}. Ignored if `exact = TRUE`.
#' @param tol Tolerance determining when to stop. Following CL21, the algorithm keeps
#'   iterating until \eqn{\textrm{max}(\sigma_n)/(\textrm{max}(\beta_n) - \textrm{min}(\beta_n)) < \textrm{tol}},
#'   where the \eqn{\beta_n} are the SHAP values of a given observation,
#'   and \eqn{\sigma_n} their standard errors.
#'   For multidimensional predictions, the criterion must be satisfied for each
#'   dimension separately. The stopping criterion uses the fact that standard errors
#'   and SHAP values are all on the same scale. Ignored if `exact = TRUE`.
#' @param max_iter If the stopping criterion (see `tol`) is not reached after
#'   `max_iter` iterations, the algorithm stops. Ignored if `exact = TRUE`.
#' @param parallel If `TRUE`, use parallel [foreach::foreach()] to loop over rows
#'   to be explained. Must register backend beforehand, e.g., via {doFuture} package,
#'   see README for an example. Parallelization automatically disables the progress bar.
#' @param parallel_args Named list of arguments passed to [foreach::foreach()].
#'   Ideally, this is `NULL` (default). Only relevant if `parallel = TRUE`.
#'   Example on Windows: if `object` is a GAM fitted with package {mgcv},
#'   then one might need to set `parallel_args = list(.packages = "mgcv")`.
#' @param verbose Set to `FALSE` to suppress messages and the progress bar.
#' @param ... Additional arguments passed to `pred_fun(object, X, ...)`.
#' @returns
#'   An object of class "permshap" with the following components:
#'   - `S`: \eqn{(n \times p)} matrix with SHAP values or, if the model output has
#'     dimension \eqn{K > 1}, a list of \eqn{K} such matrices.
#'   - `X`: Same as input argument `X`.
#'   - `baseline`: Vector of length K representing the average prediction on the
#'     background data.
#'   - `SE`: Standard errors corresponding to `S` (and organized like `S`).
#'   - `n_iter`: Integer vector of length n providing the number of iterations
#'     per row of `X`.
#'   - `converged`: Logical vector of length n indicating convergence per row of `X`.
#'   - `m`: Integer providing the effective number of sampled on-off vectors used
#'     per iteration.
#'   - `m_exact`: Integer providing the effective number of exact on-off vectors used
#'     per iteration.
#'   - `prop_exact`: Proportion of the Kernel SHAP weight distribution covered by
#'     exact calculations.
#'   - `exact`: Logical flag indicating whether calculations are exact or not.
#'   - `txt`: Summary text.
#'   - `predictions`: \eqn{(n \times K)} matrix with predictions of `X`.
#' @references
#'   1. Scott M. Lundberg and Su-In Lee. A unified approach to interpreting model
#'     predictions. Proceedings of the 31st International Conference on Neural
#'     Information Processing Systems, 2017.
#'   2. Ian Covert and Su-In Lee. Improving permshap: Practical Shapley Value
#'     Estimation Using Linear Regression. Proceedings of The 24th International
#'     Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.
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
#' summary(s)
#'
#' # Non-feature columns can be dropped via 'feature_names'
#' s <- permshap(
#'   fit,
#'   iris[1:4, ],
#'   bg_X = bg_X,
#'   feature_names = c("Petal.Length", "Petal.Width", "Species")
#' )
#' s
permshap <- function(object, ...){
  UseMethod("permshap")
}

#' @describeIn permshap Default Kernel SHAP method.
#' @export
permshap.default <- function(object, X, bg_X, pred_fun = stats::predict,
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
    all(feature_names %in% colnames(bg_X)),  # not necessary, but clearer
    all(colnames(X) %in% colnames(bg_X)),
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

  # Calculate v1 and v0
  v1 <- check_pred(pred_fun(object, X, ...), n = n)  # Predictions on X:        n x K
  bg_preds <- check_pred(
    pred_fun(object, bg_X[, colnames(X), drop = FALSE], ...), n = bg_n
  )
  v0 <- weighted_colMeans(bg_preds, bg_w)            # Average pred of bg data: 1 x K

  # For p = 1, exact Shapley values are returned
  if (p == 1L) {
    return(
      case_p1(n = n, nms = feature_names, v0 = v0, v1 = v1, X = X, verbose = verbose)
    )
  }

  # Drop unnecessary columns in bg_X. If X is matrix, also column order is relevant
  # In what follows, predictions will never be applied directly to bg_X anymore
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
  if (max(m, m_exact) * bg_n > 2e5) {
    warning("\nPredictions on large data sets with ", max(m, m_exact), "x", bg_n,
            " observations are being done.\n",
            "Consider reducing the computational burden (e.g. use smaller X_bg)")
  }

  # Apply Kernel SHAP to each row of X
  if (isTRUE(parallel)) {
    parallel_args <- c(list(i = seq_len(n)), parallel_args)
    res <- do.call(foreach::foreach, parallel_args) %dopar% permshap_one(
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
      res[[i]] <- permshap_one(
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
    txt = txt,
    predictions = v1
  )
  class(out) <- "permshap"
  out
}

#' @describeIn permshap Kernel SHAP method for "ranger" models, see Readme for an example.
#' @export
permshap.ranger <- function(object, X, bg_X,
                              pred_fun = function(m, X, ...) stats::predict(m, X, ...)$predictions,
                              feature_names = colnames(X),
                              bg_w = NULL, exact = length(feature_names) <= 8L,
                              hybrid_degree = 1L + length(feature_names) %in% 4:16,
                              paired_sampling = TRUE,
                              m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)),
                              tol = 0.005, max_iter = 100L, parallel = FALSE,
                              parallel_args = NULL, verbose = TRUE, ...) {
  permshap.default(
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

#' @describeIn permshap Kernel SHAP method for "mlr3" models, see Readme for an example.
#' @export
permshap.Learner <- function(object, X, bg_X,
                               pred_fun = function(m, X) m$predict_newdata(X)$response,
                               feature_names = colnames(X),
                               bg_w = NULL, exact = length(feature_names) <= 8L,
                               hybrid_degree = 1L + length(feature_names) %in% 4:16,
                               paired_sampling = TRUE,
                               m = 2L * length(feature_names) * (1L + 3L * (hybrid_degree == 0L)),
                               tol = 0.005, max_iter = 100L, parallel = FALSE,
                               parallel_args = NULL, verbose = TRUE, ...) {
  permshap.default(
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

