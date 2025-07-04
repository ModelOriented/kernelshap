# Model with non-linearities and interactions
y <- iris$Sepal.Length
Y <- as.matrix(iris[, c("Sepal.Length", "Sepal.Width")])

fit_y <- lm(y ~ poly(Petal.Width, degree = 2L) * Species, data = iris)
fit_Y <- lm(Y ~ poly(Petal.Width, degree = 2L) * Species, data = iris)

x <- c("Petal.Width", "Species")
J <- c(1L, 51L, 101L)

preds_y <- unname(predict(fit_y, iris))
preds_Y <- unname(predict(fit_Y, iris))

shap_y <- list(
  kernelshap(fit_y, iris[J, x], bg_X = iris, verbose = FALSE),
  permshap(fit_y, iris[J, x], bg_X = iris, verbose = FALSE),
  permshap(fit_y, iris[J, x], bg_X = iris, exact = FALSE, verbose = FALSE)
)

shap_Y <- list(
  kernelshap(fit_Y, iris[J, x], bg_X = iris, verbose = FALSE),
  permshap(fit_Y, iris[J, x], bg_X = iris, verbose = FALSE),
  permshap(fit_Y, iris[J, x], bg_X = iris, exact = FALSE, verbose = FALSE)
)

test_that("Baseline equals average prediction on background data", {
  for (i in 1:3) {
    expect_equal(shap_Y[[i]]$baseline, unname(colMeans(Y)))
  }
})

test_that("SHAP + baseline = prediction", {
  for (i in 1:3) {
    s <- shap_Y[[i]]
    expect_equal(rowSums(s$S[[1L]]) + s$baseline[1L], preds_Y[J, 1L])
    expect_equal(rowSums(s$S[[2L]]) + s$baseline[2L], preds_Y[J, 2L])
  }
})

test_that("First dimension of multioutput model equals single output", {
  for (i in 1:3) {
    expect_equal(shap_Y[[i]]$baseline[1L], shap_y[[i]]$baseline)
    expect_equal(shap_Y[[i]]$S[[1L]], shap_y[[i]]$S)
  }
})
