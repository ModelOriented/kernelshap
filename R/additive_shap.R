#' Additive SHAP
#'
#' Exact additive SHAP assuming feature independence. The implementation
#' works for models fitted via
#' - [lm()],
#' - [glm()],
#' - [mgcv::gam()],
#' - [mgcv::bam()],
#' - `gam::gam()`,
#' - [survival::coxph()], and
#' - [survival::survreg()].
#' 
#' The SHAP values are extracted via `predict(object, newdata = X, type = "terms")`,
#' a logic heavily inspired by `fastshap:::explain.lm(..., exact = TRUE)`.
#' Models with interactions (specified via `:` or `*`), or with terms of
#' multiple features like `log(x1/x2)` are not supported.
#' 
#' Note that the SHAP values obtained by [additive_shap()] are expected to
#' match those of [permshap()] and [kernelshap()] as long as their background
#' data equals the full training data (which is typically not feasible).
#'
#' @inheritParams kernelshap
#' @param X Dataframe with rows to be explained. Will be used like
#'   `predict(object, newdata = X, type = "terms")`.
#' @param ... Currently unused.
#' @returns
#'   An object of class "kernelshap" with the following components:
#'   - `S`: \eqn{(n \times p)} matrix with SHAP values.
#'   - `X`: Same as input argument `X`.
#'   - `baseline`: The baseline.
#'   - `exact`: `TRUE`.
#'   - `txt`: Summary text.
#'   - `predictions`: Vector with predictions of `X` on the scale of "terms".
#'   - `algorithm`: "additive_shap".
#' @export
#' @examples
#' # MODEL ONE: Linear regression
#' fit <- lm(Sepal.Length ~ ., data = iris)
#' s <- additive_shap(fit, head(iris))
#' s
#' 
#' # MODEL TWO: More complicated (but not very clever) formula
#' fit <- lm(
#'   Sepal.Length ~ poly(Sepal.Width, 2) + log(Petal.Length) + log(Sepal.Width),
#'   data = iris
#' )
#' s_add <- additive_shap(fit, head(iris))
#' s_add
#' 
#' # Equals kernelshap()/permshap() when background data is full training data
#' s_kernel <- kernelshap(
#'  fit, head(iris[c("Sepal.Width", "Petal.Length")]), bg_X = iris
#' )
#' all.equal(s_add$S, s_kernel$S)
additive_shap <- function(object, X, verbose = TRUE, ...) {
  stopifnot(
    inherits(object, c("lm", "glm", "gam", "bam", "Gam", "coxph", "survreg"))
  )
  if (any(attr(stats::terms(object), "order") > 1)) {
    stop("Additive SHAP not appropriate for models with interactions.")
  }
  
  txt <- "Exact additive SHAP via predict(..., type = 'terms')"
  if (verbose) {
    message(txt)
  }
  
  S <- stats::predict(object, newdata = X, type = "terms")
  rownames(S) <- NULL
  
  # Baseline value
  b <- as.vector(attr(S, "constant"))
  if (is.null(b)) {
    b <- 0
  }
  
  # Which columns of X are used in each column of S?
  s_names <- colnames(S)
  cols_used <- lapply(s_names, function(z) all.vars(stats::reformulate(z)))
  if (any(lengths(cols_used) > 1L)) {
    stop("The formula contains terms with multiple features (not supported).")
  }
  
  # Collapse all columns in S using the same column in X and rename accordingly
  mapping <- split(
    s_names, factor(unlist(cols_used), levels = colnames(X)), drop = TRUE
  )
  S <- do.call(
    cbind,
    lapply(mapping, function(z) rowSums(S[, z, drop = FALSE], na.rm = TRUE))
  )
  
  structure(
    list(
      S = S,
      X = X,
      baseline = b,
      exact = TRUE,
      txt = txt,
      predictions = b + rowSums(S),
      algorithm = "additive_shap"
    ),
    class = "kernelshap"
  )
}
