# Model with non-linearities and interactions
fit <- lm(
  Sepal.Length ~ poly(Petal.Width, degree = 2L) * Species + Petal.Length + Sepal.Width,
  data = iris
)
x <- c("Petal.Width", "Species", "Petal.Length", "Sepal.Width")
J <- c(1L, 51L, 101L)
preds <- unname(predict(fit, iris[J, ]))

shap <- list(
  # Exact
  kernelshap(fit, iris[J, x], bg_X = iris, verbose = FALSE),
  permshap(fit, iris[J, x], bg_X = iris, verbose = FALSE),
  # Sampling
  kernelshap(
    fit, iris[J, x],
    bg_X = iris, exact = FALSE, hybrid_degree = 0, verbose = FALSE
  ),
  permshap(fit, iris[J, x], bg_X = iris, exact = FALSE, verbose = FALSE),
  permshap(fit, iris[J, x], bg_X = iris, exact = FALSE, low_memory = TRUE, verbose = FALSE)
)

test_that("baseline equals average prediction on background data", {
  for (s in shap) {
    expect_equal(s$baseline, mean(iris$Sepal.Length))
  }
})

test_that("SHAP + baseline = prediction", {
  for (s in shap) {
    expect_equal(rowSums(s$S) + s$baseline, preds)
  }
})

test_that("Exact and sampling modes agree with interactions of order 2", {
  expect_equal(shap[[1L]]$S, shap[[2L]]$S) # exact ks vs exact ps
  expect_equal(shap[[1L]]$S, shap[[3L]]$S) # exact ks vs sampling ks
  expect_equal(shap[[2L]]$S, shap[[4L]]$S) # exact ps vs sampling ps
  expect_equal(shap[[4L]], shap[[5L]]) # low/high-memory
  expect_true(all(shap[[3L]]$n_iter == 2L)) # ks stops after second iteration
  expect_true(all(shap[[4L]]$n_iter == length(x))) # ps stops after p iteration
})

test_that("permshap() in sampling mode requires at least 4 features", {
  expect_error(
    permshap(
      fit, iris[1:3, x],
      bg_X = iris, exact = FALSE, feature_names = x[1:2], verbose = FALSE
    )
  )
})

test_that("kernelshap() with max_iter = 1 works", {
  ks <- kernelshap(
    fit, iris[J, x],
    bg_X = iris, exact = FALSE,
    hybrid_degree = 0L, max_iter = 1L, verbose = FALSE
  )
  expect_equal(ks$S, shap[[1L]]$S) # should be the same as exact ks for simple model
  expect_equal(ks$n_iter, rep(1L, 3L))
  expect_true(all(is.na(ks$SE)))
  expect_true(all(!ks$converged))
})

test_that("auto-selection of background data works", {
  # Here, the background data equals the full X
  shap2 <- list(
    kernelshap(fit, iris[x], verbose = FALSE),
    permshap(fit, iris[x], verbose = FALSE)
  )

  for (i in 1:2) {
    expect_equal(shap$S, shap2$S)
  }
})

test_that("missing bg_X gives error if X is very small", {
  for (algo in c(kernelshap, permshap)) {
    expect_error(algo(fit, iris[1:10, x], verbose = FALSE))
  }
})

test_that("missing bg_X gives warning if X is quite small", {
  for (algo in c(kernelshap, permshap)) {
    expect_warning(algo(fit, iris[1:30, x], verbose = FALSE))
  }
})

test_that("selection of bg_X can be controlled via bg_n", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[x], verbose = FALSE, bg_n = 20L)
    expect_equal(nrow(s$bg_X), 20L)
  }
})

test_that("verbose is chatty", {
  for (algo in c(kernelshap, permshap)) {
    capture_output(expect_message(algo(fit, iris[J, x], bg_X = iris, verbose = TRUE)))
  }
})

test_that("large background data cause warning", {
  # Takes a bit of time, thus only for one algo
  large_bg <- iris[rep(1:150, 230), ]
  expect_warning(
    kernelshap(fit, iris[1L, x], bg_X = large_bg, verbose = FALSE)
  )
})

test_that("Decomposing a single row works", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[1L, x], bg_X = iris, verbose = FALSE)
    expect_equal(s$baseline, mean(iris$Sepal.Length))
    expect_equal(rowSums(s$S) + s$baseline, preds[1])
  }
})

test_that("Background data can contain additional columns", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[1L, x], bg_X = cbind(d = 1, iris), verbose = FALSE)
    expect_true(is.kernelshap(s))
  }
})

test_that("Background data can contain only one single row", {
  for (algo in c(kernelshap, permshap)) {
    expect_no_error(algo(fit, iris[1L, x], bg_X = iris[150L, ], verbose = FALSE))
  }
})

test_that("feature_names can drop columns from SHAP calculations", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[J, ], bg_X = iris, feature_names = x, verbose = FALSE)
    expect_equal(colnames(s$S), x)
  }
})

test_that("feature_names can rearrange column names in result", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[J, ], bg_X = iris, feature_names = rev(x), verbose = FALSE)
    expect_equal(colnames(s$S), rev(x))
  }
})

test_that("feature_names must be in colnames(X) and colnames(bg_X)", {
  for (algo in c(kernelshap, permshap)) {
    expect_error(algo(fit, iris, bg_X = cbind(iris, a = 1), feature_names = "a"))
    expect_error(algo(fit, cbind(iris, a = 1), bg_X = iris, feature_names = "a"))
  }
})

test_that("Matrix input is fine", {
  X <- data.matrix(iris)
  pred_fun <- function(m, X) {
    data <- as.data.frame(X) |>
      transform(Species = factor(Species, labels = levels(iris$Species)))
    predict(m, data)
  }

  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, X[J, x], pred_fun = pred_fun, bg_X = X, verbose = FALSE)

    expect_equal(s$baseline, mean(iris$Sepal.Length)) # baseline is mean of bg
    expect_equal(rowSums(s$S) + s$baseline, preds) # sum shap = centered preds
    expect_no_error( # additional cols in bg are ok
      algo(fit, X[J, x], pred_fun = pred_fun, bg_X = cbind(d = 1, X), verbose = FALSE)
    )
  }
})

test_that("Special case p = 1 works only for kernelshap()", {
  capture_output(
    expect_message(
      s <- kernelshap(fit, X = iris[J, ], bg_X = iris, feature_names = "Petal.Width")
    )
  )
  expect_equal(s$baseline, mean(iris$Sepal.Length))
  expect_equal(unname(rowSums(s$S)) + s$baseline, preds)

  expect_error( # Not implemented
    permshap(
      fit, iris[J, ],
      bg_X = iris, verbose = FALSE, feature_names = "Petal.Width"
    )
  )
})

test_that("exact hybrid kernelshap() is similar to exact (non-hybrid)", {
  s1 <- kernelshap(
    fit, iris[J, x],
    bg_X = iris, exact = FALSE, hybrid_degree = 1L, verbose = FALSE
  )
  expect_equal(s1$S, shap[[1L]]$S)
})

# test_that("kernelshap works for large p (hybrid case)", {
#   set.seed(9L)
#   X <- data.frame(matrix(rnorm(20000L), ncol = 100L))
#   y <- X[, 1L] * X[, 2L] * X[, 3L]
#   fit <- lm(y ~ X1:X2:X3 + ., data = cbind(y = y, X))
#   s <- kernelshap(fit, X[1L, ], bg_X = X, verbose = FALSE)
#
#   expect_equal(s$baseline, mean(y))
#   expect_equal(rowSums(s$S) + s$baseline, unname(predict(fit, X[1L, ])))
# })

test_that("kernelshap and permshap work for models with high-order interactions", {
  # Expected: Python output
  # import numpy as np
  # import shap # 0.47.2
  #
  # X = np.array(
  #   [
  #     np.arange(1, 101) / 100,
  #     np.log(np.arange(1, 101)),
  #     np.sqrt(np.arange(1, 101)),
  #     np.sin(np.arange(1, 101)),
  #     (np.arange(1, 101) / 100) ** 2,
  #     np.cos(np.arange(1, 101)),
  #   ]
  # ).T
  #
  #
  # def predict(X):
  #   return X[:, 0] * X[:, 1] * X[:, 2] * X[:, 3] + X[:, 4] + X[:, 5]
  #
  #
  # ks = shap.explainers.Kernel(predict, X, nsamples=10000)
  # es = shap.explainers.Exact(predict, X)
  #
  # print("Exact Kernel SHAP:\n", ks(X[0:2]).values)
  # print("Exact (Permutation) SHAP:\n", es(X[0:2]).values)
  #
  # # Exact Kernel SHAP:
  # #  [[-1.19621609 -1.24184808 -0.9567848   3.87942037 -0.33825     0.54562519]
  # #  [-1.64922699 -1.20770105 -1.18388581  4.54321217 -0.33795    -0.41082395]]
  # # Exact (Permutation) SHAP:
  # #  [[-1.19621609 -1.24184808 -0.9567848   3.87942037 -0.33825     0.54562519]
  # #  [-1.64922699 -1.20770105 -1.18388581  4.54321217 -0.33795    -0.41082395]]

  expected <- rbind(
    c(-1.19621609, -1.24184808, -0.9567848, 3.87942037, -0.33825, 0.54562519),
    c(-1.64922699, -1.20770105, -1.18388581, 4.54321217, -0.33795, -0.41082395)
  )

  n <- 100

  X <- data.frame(
    x1 = seq(1:n) / 100,
    x2 = log(1:n),
    x3 = sqrt(1:n),
    x4 = sin(1:n),
    x5 = (seq(1:n) / 100)^2,
    x6 = cos(1:n)
  )

  pf <- function(model, newdata) {
    x <- newdata
    x[, 1] * x[, 2] * x[, 3] * x[, 4] + x[, 5] + x[, 6]
  }
  ks <- kernelshap(pf, head(X, 2), bg_X = X, pred_fun = pf, verbose = FALSE)
  expect_equal(unname(ks$S), expected)

  ps <- permshap(pf, head(X, 2), bg_X = X, pred_fun = pf, verbose = FALSE)
  expect_equal(unname(ps$S), expected)

  # Sampling versions of KernelSHAP is quite close
  set.seed(1)
  ksh2 <- kernelshap(
    pf,
    head(X, 1),
    bg_X = X,
    pred_fun = pf,
    hybrid_degree = 2,
    exact = FALSE,
    m = 1000,
    max_iter = 100,
    tol = 0.001,
    verbose = FALSE
  )
  expect_equal(
    c(ksh2$S),
    c(-1.202621, -1.245655, -0.9510531, 3.880391, -0.3384063, 0.5492909),
    tolerance = 1e-4
  )

  set.seed(1)
  ksh1 <- kernelshap(
    pf,
    head(X, 1),
    bg_X = X,
    pred_fun = pf,
    hybrid_degree = 1,
    exact = FALSE,
    m = 1000,
    max_iter = 1000,
    tol = 0.002,
    verbose = FALSE
  )
  expect_equal(
    c(ksh1$S),
    c(-1.203028, -1.240056, -0.9574796, 3.89151, -0.3360694, 0.5370693),
    tolerance = 1e-3
  )

  set.seed(1)
  ksh0 <- suppressWarnings(
    kernelshap(
      pf,
      head(X, 1),
      bg_X = X,
      pred_fun = pf,
      hybrid_degree = 0,
      exact = FALSE,
      m = 10000,
      max_iter = 10000,
      tol = 0.005,
      verbose = FALSE
    )
  )
  expect_equal(
    c(ksh0$S),
    c(-1.212142, -1.238348, -0.9795299, 3.937577, -0.3726428, 0.5570322),
    tolerance = 1e-3
  )
})


test_that("Random seed works", {
  n <- 100

  X <- data.frame(
    x1 = seq(1:n) / 100,
    x2 = log(1:n),
    x3 = sqrt(1:n),
    x4 = sin(1:n),
    x5 = (seq(1:n) / 100)^2,
    x6 = cos(1:n)
  )

  pf <- function(model, newdata, ...) {
    x <- newdata
    x[, 1] * x[, 2] * x[, 3] * x[, 4] + x[, 5] + x[, 6]
  }

  for (algo in c(permshap, kernelshap)) {
    s1a <- algo(pf, head(X, 2), bg_X = X, pred_fun = pf, verbose = FALSE, seed = 1, exact = FALSE, hybrid_degree = 0)
    s1b <- algo(pf, head(X, 2), bg_X = X, pred_fun = pf, verbose = FALSE, seed = 1, exact = FALSE, hybrid_degree = 0)
    s2 <- algo(pf, head(X, 2), bg_X = X, pred_fun = pf, verbose = FALSE, seed = 2, exact = FALSE, hybrid_degree = 0)
    expect_equal(s1a, s1b)
    expect_false(identical(s1a$S, s2$S))
  }
})
