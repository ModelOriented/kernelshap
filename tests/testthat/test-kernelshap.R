# Model with non-linearities and interactions
fit <- stats::lm(Sepal.Length ~ poly(Petal.Width, 2) * Species, data = iris)
pred_fun <- function(X) stats::predict(fit, X)
x <- c("Petal.Width", "Species")
preds <- unname(pred_fun(iris))
s <- kernelshap(iris[1:5, x], pred_fun = pred_fun, bg_X = iris[, x])

test_that("Baseline equals average prediction on background data", {
  expect_equal(s$baseline, mean(iris$Sepal.Length))
})

test_that("SHAP + baseline = prediction", {
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
})

test_that("Decomposing a single row works", {
  s <- kernelshap(iris[1L, x], pred_fun = pred_fun, bg_X = iris[, x])
  
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1])
})

test_that("Background data can contain additional columns", {
  expect_true(
    is.kernelshap(
      kernelshap(iris[1L, x], pred_fun = pred_fun, bg_X = cbind(d = 1, iris[, x]))
    )
  )
})

fit <- stats::lm(Sepal.Length ~ stats::poly(Petal.Width, 2), data = iris)
x <- "Petal.Width"
preds <- unname(pred_fun(iris[x]))

test_that("Special case p = 1 works", {
  s <- kernelshap(iris[1:5, x, drop = FALSE], pred_fun = pred_fun, bg_X = iris[x])
  
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
  expect_equal(s$SE[1L], 0)
})

fit <- stats::lm(Sepal.Length ~ ., data = iris[1:4])
X <- data.matrix(iris[2:4])
pred_fun2 <- function(X) stats::predict(fit, as.data.frame(X))
preds <- unname(pred_fun2(X))
s <- kernelshap(X[1:3, ], pred_fun = pred_fun2, X)

test_that("Matrix input is fine", {
  expect_true(is.kernelshap(s))
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:3])
})

test_that("Matrix input works if bg data containts extra columns", {
  expect_true(is.kernelshap(kernelshap(X[1:3, ], pred_fun = pred_fun2, cbind(d = 1, X))))
})

## Now with case weights
fit <- stats::lm(
  Sepal.Length ~ poly(Petal.Width, 2) * Species, 
  data = iris, 
  weights = Petal.Length
)
x <- c("Petal.Width", "Species")
preds <- unname(pred_fun(iris))
s <- kernelshap(
  iris[1:5, x], pred_fun = pred_fun, bg_X = iris[, x], bg_w = iris$Petal.Length
)

test_that("Baseline equals weighted average prediction on background data", {
  expect_equal(s$baseline, stats::weighted.mean(iris$Sepal.Length, iris$Petal.Length))
})

test_that("SHAP + baseline = prediction works with case weights", {
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
})

test_that("Decomposing a single row works with case weights", {
  s <- kernelshap(
    iris[1, x], pred_fun = pred_fun, bg_X = iris[, x], bg_w = iris$Petal.Length
  )
  
  expect_equal(s$baseline, stats::weighted.mean(iris$Sepal.Length, iris$Petal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1])
})

fit <- stats::lm(
  Sepal.Length ~ stats::poly(Petal.Width, 2), data = iris, weights = Petal.Length
)
x <- "Petal.Width"
preds <- unname(pred_fun(iris[x]))

test_that("Special case p = 1 works with case weights", {
  s <- kernelshap(
    iris[1:5, x, drop = FALSE], 
    pred_fun = pred_fun, 
    bg_X = iris[x], 
    bg_w = iris$Petal.Length
  )
  
  expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, iris$Petal.Length))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:5])
})

fit <- stats::lm(
  Sepal.Length ~ . , data = iris[c(1, 3, 4)], weights = iris$Sepal.Width
)
X <- data.matrix(iris[3:4])
preds <- unname(pred_fun2(X))

test_that("Matrix input is fine with case weights", {
  s <- kernelshap(X[1:3, ], pred_fun = pred_fun2, X, bg_w = iris$Sepal.Width)
  
  expect_true(is.kernelshap(s))
  expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, iris$Sepal.Width))
  expect_equal(rowSums(s$S) + s$baseline, preds[1:3])
})

test_that("kernelshap works for large p", {
  set.seed(9L)
  X <- data.frame(matrix(rnorm(20000L), ncol = 100L))
  y <- rnorm(200L)
  fit <- lm(y ~ ., data = cbind(y = y, X))
  s <- kernelshap(X[1L, ], function(X) predict(fit, X), X)

  expect_equal(s$baseline, mean(y))
  expect_equal(rowSums(s$S) + s$baseline, unname(predict(fit, X[1L, ])))
})

