fit <- stats::lm(Sepal.Length ~ ., data = iris)
pred_fun <- function(X) stats::predict(fit, X)
s <- kernelshap(iris[1:2, -1], pred_fun = pred_fun, iris[-1])

test_that("get_* functions get the right element", {
  expect_equal(s$S, ks_shap_values(s))
  expect_equal(s$X, ks_feature_values(s))
  expect_equal(s$SE, ks_standard_errors(s))
  expect_equal(s$baseline, ks_baseline(s))
  expect_equal(s$converged, ks_converged(s))
  expect_equal(s$n_iter, ks_n_iter(s))
})

test_that("get_* functions fail for wrong objects", {
  expect_error(ks_shap_values(1))
  expect_error(ks_feature_values(1))
  expect_error(ks_standard_errors(1))
  expect_error(ks_baseline(1))
  expect_error(ks_converged(1))
  expect_error(ks_n_iter(1))
})

test_that("is_kernelshap() works", {
  expect_true(is.kernelshap(s))
  expect_false(is.kernelshap(1))
})
