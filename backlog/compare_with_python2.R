library(kernelshap)

n <- 1000

X <- data.frame(
  x1 = seq(1:n) / 100,
  x2 = log(1:n),
  x3 = sqrt(1:n),
  x4 = sin(1:n)
)
head(X)

pf <- function(model, newdata) {
  x <- newdata
  x[, 1] * x[, 2] * x[, 3] + x[, 4]
}
ks <- kernelshap(pf, X, bg_X = X, pred_fun = pf)
ks
es <- permshap(pf, X, bg_X = X, pred_fun = pf)
es

library(shapr)

specs <- function(x) {
  feature_specs <- list()
  feature_specs$labels <- c("x1", "x2", "x3", "x4")
  feature_specs$classes <- setNames(rep("numeric", 4), feature_specs$labels)
  feature_specs$factor_levels <- list(x1 = NULL, x2 = NULL, x3 = NULL, x4 = NULL)
  feature_specs
}


sh <- explain(pf, X, X,
  approach = "independence",
  predict_model = pf, phi0 = es$baseline,
  get_model_specs = specs,
  max_n_coalitions = 100
)
