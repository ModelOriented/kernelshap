# Model with non-linearities and interactions
fit <- lm(
  Sepal.Length ~ poly(Petal.Width, degree = 2L) * Species,
  data = iris,
  weights = Petal.Length
)
x <- c("Petal.Width", "Species")
preds <- unname(predict(fit, iris))
J <- c(1L, 51L, 101L)
w <- iris$Petal.Length

shap <- list(
  kernelshap(fit, iris[x], bg_X = iris, bg_w = w, verbose = FALSE),
  permshap(fit, iris[x], bg_X = iris, bg_w = w, verbose = FALSE)
)

test_that("constant weights gives same as no weights", {
  shap_unweighted <- list(
    kernelshap(fit, iris[x], bg_X = iris, verbose = FALSE),
    permshap(fit, iris[x], bg_X = iris, verbose = FALSE)
  )
  
  w2 <- rep(3, nrow(iris))
  shap2 <- list(
    kernelshap(fit, iris[x], bg_X = iris, bg_w = w2, verbose = FALSE),
    permshap(fit, iris[x], bg_X = iris, bg_w = w2, verbose = FALSE)
  )
  
  for (i in seq_along(shap))
    expect_equal(shap2[[i]]$S, shap_unweighted[[i]]$S)
})

test_that("baseline equals average prediction on background data", {
  for (s in shap)
    expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, w))
})

test_that("SHAP + baseline = prediction for exact mode", {
  for (s in shap)
    expect_equal(rowSums(s$S) + s$baseline, preds)
})

test_that("Decomposing a single row works", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[1L, x], bg_X = iris, bg_w = w, verbose = FALSE)
    expect_equal(s$baseline, weighted.mean(iris$Sepal.Length, w))
    expect_equal(rowSums(s$S) + s$baseline, preds[1])
  }
})

test_that("auto-selection of background data works", {
  # Here, the background data equals the full X
  shap2 <- list(
    kernelshap(fit, iris[x], bg_w = w, verbose = FALSE),
    permshap(fit, iris[x], bg_w = w, verbose = FALSE)
  )
  
  for (i in 1:2) {
    expect_equal(shap$S, shap2$S)
  }
})

test_that("selection of bg_X can be controlled via bg_n", {
  n <- 20L
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris, bg_w = w, verbose = FALSE, bg_n = n)
    expect_equal(nrow(s$bg_X), n)
  }
})

test_that("weights must have correct length", {
  for (algo in c(kernelshap, permshap)) {
    expect_error(algo(fit, iris[J, ], bg_X = iris, bg_w = 1:3, verbose = FALSE))
  }
})

test_that("weights can't be all 0", {
  for (algo in c(kernelshap, permshap)) {
    expect_error(
      algo(fit, iris[J, ], bg_X = iris, bg_w = rep(0, nrow(iris)), verbose = FALSE)
    )
  }
})

test_that("weights can't be negative", {
  for (algo in c(kernelshap, permshap)) {
    expect_error(
      algo(fit, iris[J, ], bg_X = iris, bg_w = rep(-1, nrow(iris)), verbose = FALSE)
    )
  }
})

