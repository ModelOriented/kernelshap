fit <- lm(Sepal.Length ~ ., data = iris)

set.seed(1)

shap <- list(
  kernelshap(
    fit, iris[1:2, -1L], bg_X = iris, verbose = FALSE, exact = FALSE, hybrid_degree = 1
  ),
  permshap(fit, iris[1:2, -1L], bg_X = iris, verbose = FALSE),
  additive_shap(fit, iris, verbose = FALSE)
)

test_that("is.kernelshap() works", {
  for (s in shap) {
    expect_true(is.kernelshap(s))
    expect_false(is.kernelshap(1))
  }
})

test_that("print() and summary() do not give an error", {
  for (s in shap) {
    capture_output(expect_no_error(print(s)))
    capture_output(expect_no_error(summary(s)))
    capture_output(expect_no_error(summary(s, compact = TRUE)))
  }
})

