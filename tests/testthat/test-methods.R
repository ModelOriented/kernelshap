fit <- stats::lm(Sepal.Length ~ ., data = iris)
s <- kernelshap(fit, iris[1:2, -1L], bg_X = iris[-1L], verbose = FALSE)

test_that("is_kernelshap() works", {
  expect_true(is.kernelshap(s))
  expect_false(is.kernelshap(1))
})

