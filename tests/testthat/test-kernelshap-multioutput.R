#===========================================================
# Tests for multi-output model
#===========================================================

# Model with non-linearities and interactions
y <- iris$Sepal.Length
Y <- as.matrix(iris[, c("Sepal.Length", "Sepal.Width")])

fity <- lm(y ~ poly(Petal.Width, degree = 2L) * Species, data = iris)
fitY <- lm(Y ~ poly(Petal.Width, degree = 2L) * Species, data = iris)

x <- c("Petal.Width", "Species")

predsy <- unname(predict(fity, iris))
predsY <- unname(predict(fitY, iris))

sy <- kernelshap(fity, iris[1:5, x], bg_X = iris, verbose = FALSE)
sY <- kernelshap(fitY, iris[1:5, x], bg_X = iris, verbose = FALSE)

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
  sY <- kernelshap(fitY, iris[1L, x], bg_X = iris, verbose = FALSE)
  
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(rowSums(sY$S[[1L]]) + sY$baseline[1L], predsY[1L, 1L])
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[1L, 2L])
})

fitY <- lm(Y ~ poly(Petal.Width, degree = 2L), data = iris)
x <- "Petal.Width"
predsY <- unname(predict(fitY, iris))

test_that("Special case p = 1 works", {
  sY <- kernelshap(fitY, iris[1:5, x, drop = FALSE], bg_X = iris, verbose = FALSE)
  
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(unname(rowSums(sY$S[[2L]]) + sY$baseline[2L]), predsY[1:5, 2L])
  expect_equal(sY$SE[[1L]][1L], 0)
})

fitY <- lm(Y ~ Petal.Length + Petal.Width, data = iris[1:4])
X <- data.matrix(iris[2:4])
pred_fun <- function(fit, X) predict(fit, as.data.frame(X))
predsY <- unname(pred_fun(fitY, X))
sY <- kernelshap(fitY, X[1:3, ], pred_fun = pred_fun, bg_X = X, verbose = FALSE)

test_that("Matrix input is fine", {
  expect_true(is.kernelshap(sY))
  expect_equal(sY$baseline, unname(colMeans(Y)))
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[1:3, 2L])
})

## Now with case weights
fitY <- lm(
  Y ~ poly(Petal.Width, degree = 2L) * Species, data = iris, weights = Petal.Length
)
x <- c("Petal.Width", "Species")
predsY <- unname(predict(fitY, iris))
sY <- kernelshap(
  fitY, 
  iris[5:10, x], 
  pred_fun = predict, 
  bg_X = iris, 
  bg_w = iris$Petal.Length, 
  verbose = FALSE
)

test_that("Baseline equals weighted average prediction on background data", {
  expect_equal(sY$baseline[1L], weighted.mean(Y[, 1L], iris$Petal.Length))
  expect_equal(sY$baseline[2L], weighted.mean(Y[, 2L], iris$Petal.Length))
})

test_that("SHAP + baseline = prediction works with case weights", {
  expect_equal(rowSums(sY$S[[2L]]) + sY$baseline[2L], predsY[5:10, 2L])
})

