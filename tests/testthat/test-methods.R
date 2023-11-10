fit <- lm(Sepal.Length ~ ., data = iris)
set.seed(1)
s <- kernelshap(
  fit, iris[1:2, -1L], bg_X = iris[-1L], verbose = FALSE, exact = FALSE, hybrid_degree = 1
)

test_that("is.kernelshap() works", {
  expect_true(is.kernelshap(s))
  expect_false(is.kernelshap(1))
})

test_that("print() and summary() do not give an error", {
  capture_output(expect_no_error(print(s)))
  capture_output(expect_no_error(summary(s)))
  capture_output(expect_no_error(summary(s, compact = TRUE)))
})

# permshap
s <- permshap(fit, iris[1:2, -1L], bg_X = iris[-1L], verbose = FALSE)
test_that("is.permshap() works", {
  expect_true(is.permshap(s))
  expect_false(is.permshap(1))
})

test_that("print() and summary() also works for permshap", {
  capture_output(expect_no_error(print(s)))
  capture_output(expect_no_error(summary(s)))
})
