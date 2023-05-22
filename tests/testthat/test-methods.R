fit <- lm(Sepal.Length ~ ., data = iris)
s <- kernelshap(fit, iris[1:2, -1L], bg_X = iris[-1L], verbose = FALSE)

test_that("is_kernelshap() works", {
  expect_true(is.kernelshap(s))
  expect_false(is.kernelshap(1))
})

test_that("print() and summary() do not give an error", {
  capture_output(expect_no_error(print(s)))
  capture_output(expect_no_error(summary(s)))
})

