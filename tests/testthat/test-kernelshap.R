# Model with non-linearities and interactions
fit <- stats::lm(Sepal.Length ~ poly(Petal.Width, 2) * Species, data = iris)
x <- c("Petal.Width", "Species")
preds <- unname(predict(fit, iris))
s <- kernelshap(fit, iris[c(1, 51, 101), x], bg_X = iris)

test_that("Baseline equals average prediction on background data", {
  expect_equal(s$baseline, mean(iris$Sepal.Length))
})

test_that("SHAP + baseline = prediction", {
  expect_equal(rowSums(s$S) + s$baseline, preds[c(1, 51, 101)])
})

test_that("Non-exact calculation is similar to exact", {
  s1 <- kernelshap(fit, iris[c(1, 51, 101), x], bg_X = iris, hybrid_degree = 1)
  expect_equal(s$S, s1$S)
})

test_that("Decomposing a single row works", {
  s <- kernelshap(fit, iris[1L, x], bg_X = iris)
  
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1])
})

test_that("Background data can contain additional columns", {
  ks4 <- kernelshap(fit, iris[1L, x], bg_X = cbind(d = 1, iris))
  expect_true(is.kernelshap(ks4))
})

fit <- stats::lm(Sepal.Length ~ stats::poly(Petal.Width, 2), data = iris)
x <- "Petal.Width"
preds <- unname(stats::predict(fit, iris))

test_that("Special case p = 1 works", {
  s <- kernelshap(fit, iris[1:5, x, drop = FALSE], bg_X = iris)
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
  expect_equal(s$SE[1L], 0)
})

fit <- stats::lm(Sepal.Length ~ ., data = iris[1:4])
X <- data.matrix(iris[2:4])
pred_fun <- function(m, X) stats::predict(m, as.data.frame(X))
preds <- unname(pred_fun(fit, X))
s <- kernelshap(fit, X[1:3, ], pred_fun = pred_fun, bg_X = X)

test_that("Matrix input is fine", {
  expect_true(is.kernelshap(s))
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:3])
})

test_that("Matrix input works if bg data containts extra columns", {
  ks5 <- kernelshap(fit, X[1:3, ], pred_fun = pred_fun, bg_X = cbind(d = 1, X))
  expect_true(is.kernelshap(ks5))
})

## Now with case weights
fit <- stats::lm(
  Sepal.Length ~ poly(Petal.Width, 2) * Species, data = iris, weights = Petal.Length
)
x <- c("Petal.Width", "Species")
preds <- unname(stats::predict(fit, iris))
s <- kernelshap(fit, iris[1:5, x], bg_X = iris, bg_w = iris$Petal.Length)

test_that("Baseline equals weighted average prediction on background data", {
  expect_equal(s$baseline, stats::weighted.mean(iris$Sepal.Length, iris$Petal.Length))
})

test_that("SHAP + baseline = prediction works with case weights", {
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
})

test_that("Decomposing a single row works with case weights", {
  s <- kernelshap(fit, iris[1, x], bg_X = iris, bg_w = iris$Petal.Length)
  expect_equal(s$baseline, stats::weighted.mean(iris$Sepal.Length, iris$Petal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1])
})

fit <- stats::lm(
  Sepal.Length ~ stats::poly(Petal.Width, 2), data = iris, weights = Petal.Length
)
x <- "Petal.Width"
preds <- unname(stats::predict(fit, iris))

test_that("Special case p = 1 works with case weights", {
  s <- kernelshap(
    fit, iris[1:5, x, drop = FALSE], bg_X = iris, bg_w = iris$Petal.Length
  )
  
  expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, iris$Petal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
})

fit <- stats::lm(
  Sepal.Length ~ . , data = iris[c(1, 3, 4)], weights = iris$Sepal.Width
)
X <- data.matrix(iris[3:4])
preds <- unname(pred_fun(fit, X))

test_that("Matrix input is fine with case weights", {
  s <- kernelshap(fit, X[1:3, ], pred_fun = pred_fun, bg_X = X, bg_w = iris$Sepal.Width)
  
  expect_true(is.kernelshap(s))
  expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, iris$Sepal.Width))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:3])
})

test_that("kernelshap works for large p", {
  set.seed(9L)
  X <- data.frame(matrix(rnorm(20000L), ncol = 100L))
  y <- rnorm(200L)
  fit <- lm(y ~ ., data = cbind(y = y, X))
  s <- kernelshap(fit, X[1L, ], bg_X = X)

  expect_equal(s$baseline, mean(y))
  expect_equal(rowSums(s$S) + s$baseline, unname(stats::predict(fit, X[1L, ])))
})

