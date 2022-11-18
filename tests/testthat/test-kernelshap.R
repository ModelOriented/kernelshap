# Model with non-linearities and interactions
fit <- stats::lm(
  Sepal.Length ~ poly(Petal.Width, 2) * Species + Petal.Length, data = iris
)
x <- c("Petal.Width", "Species", "Petal.Length")
preds <- unname(predict(fit, iris))
s <- kernelshap(fit, iris[c(1, 51, 101), x], bg_X = iris)

test_that("Baseline equals average prediction on background data", {
  expect_equal(s$baseline, mean(iris$Sepal.Length))
})

test_that("SHAP + baseline = prediction", {
  expect_equal(rowSums(s$S) + s$baseline, preds[c(1, 51, 101)])
})

test_that("Exact hybrid calculation is similar to exact (non-hybrid)", {
  s1 <- kernelshap(
    fit, iris[c(1, 51, 101), x], bg_X = iris, exact = FALSE, hybrid_degree = 1
  )
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

test_that("Background data can contain only one single row", {
  expect_true(is.kernelshap(kernelshap(fit, iris[1L, x], bg_X = iris[150L, ])))
  expect_true(is.kernelshap(kernelshap(fit, iris[1:10, x], bg_X = iris[150L, ])))
})

test_that("feature_names can drop columns from SHAP calculations", {
  s_f <- kernelshap(fit, iris[c(1, 51, 101), ], bg_X = iris, feature_names = x)
  expect_equal(within(unclass(s), rm(X)), within(unclass(s_f), rm(X)))
})

test_that("feature_names can rearrange column names in result", {
  s_f2 <- kernelshap(fit, iris[c(1, 51, 101), ], bg_X = iris, feature_names = rev(x))
  expect_equal(s$S, s_f2$S[, x])
})

test_that("feature_names must be in colnames(X) and colnames(bg_X)", {
  expect_error(kernelshap(fit, iris, bg_X = cbind(iris, a = 1), feature_names = "a"))
  expect_error(kernelshap(fit, cbind(iris, a = 1), bg_X = iris, feature_names = "a"))
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

set.seed(9L)
X <- data.frame(matrix(rnorm(20000L), ncol = 100L))
y <- X[, 1L] * X[, 2L] * X[, 3L]
fit <- lm(y ~ X1:X2:X3 + ., data = cbind(y = y, X))
s <- kernelshap(fit, X[1L, ], bg_X = X)

test_that("kernelshap works for large p (hybrid case)", {
  expect_equal(s$baseline, mean(y))
  expect_equal(rowSums(s$S) + s$baseline, unname(stats::predict(fit, X[1L, ])))
})


