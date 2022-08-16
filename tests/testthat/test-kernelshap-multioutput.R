#===========================================================
# Tests for multi-output model
#===========================================================

# Model with non-linearities and interactions
y <- iris$Sepal.Length
Y <- as.matrix(iris[, c("Sepal.Length", "Sepal.Width")])

fity <- stats::lm(y ~ poly(Petal.Width, 2) * Species, data = iris)
fitY <- stats::lm(Y ~ poly(Petal.Width, 2) * Species, data = iris)

pred_funy <- function(X) stats::predict(fity, X)
pred_funY <- function(X) stats::predict(fitY, X)

x <- c("Petal.Width", "Species")

predsy <- unname(pred_funy(iris))
predsY <- unname(pred_funY(iris))

sy <- kernelshap(iris[1:5, x], pred_fun = pred_funy, bg_X = iris[, x])
sY <- kernelshap(iris[1:5, x], pred_fun = pred_funY, bg_X = iris[, x])

test_that("Baseline equals average prediction on background data", {
  expect_equal(sY$baseline, unname(colMeans(Y)))
})

test_that("SHAP + baseline = prediction", {
  expect_equal(rowSums(sY$S[[1L]]) + sY$baseline[1L], predsY[1:5, 1L])
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[1:5, 2L])
})

test_that("First dimension of multioutput model equals single output (approx)", {
  expect_equal(sY$baseline[1L], sy$baseline)
  expect_equal(sY$S[[1L]], sy$S)
})

test_that("Decomposing a single row works", {
  sY <- kernelshap(iris[1L, x], pred_fun = pred_funY, bg_X = iris[, x])
  
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(rowSums(sY$S[[1L]]) + sY$baseline[1L], predsY[1L, 1L])
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[1L, 2L])
})

fitY <- stats::lm(Y ~ stats::poly(Petal.Width, 2), data = iris)
x <- "Petal.Width"
predsY <- unname(pred_funY(iris[x]))

test_that("Special case p = 1 works", {
  sY <- kernelshap(iris[1:5, x, drop = FALSE], pred_fun = pred_funY, bg_X = iris[x])
  
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(unname(rowSums(sY$S[[2L]]) + sY$baseline[2L]), predsY[1:5, 2L])
  expect_equal(sY$SE[[1L]][1L], 0)
})

fitY <- stats::lm(Y ~ Petal.Length + Petal.Width, data = iris[1:4])
X <- data.matrix(iris[2:4])
pred_fun2 <- function(X) stats::predict(fitY, as.data.frame(X))
predsY <- unname(pred_fun2(X))
sY <- kernelshap(X[1:3, ], pred_fun = pred_fun2, X)

test_that("Matrix input is fine", {
  expect_true(is.kernelshap(sY))
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[1:3, 2L])
})

## Now with case weights
fitY <- stats::lm(
  Y ~ poly(Petal.Width, 2) * Species, data = iris, weights = Petal.Length
)
x <- c("Petal.Width", "Species")
predsY <- unname(pred_funY(iris))
sY <- kernelshap(
  iris[5:10, x], pred_fun = pred_funY, bg_X = iris[, x], bg_w = iris$Petal.Length
)

test_that("Baseline equals weighted average prediction on background data", {
  expect_equal(sY$baseline[1L], stats::weighted.mean(Y[, 1L], iris$Petal.Length))
  expect_equal(sY$baseline[2L], stats::weighted.mean(Y[, 2L], iris$Petal.Length))
})

test_that("SHAP + baseline = prediction works with case weights", {
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[5:10, 2L])
})

