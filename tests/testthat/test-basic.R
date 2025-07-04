# Model with non-linearities and interactions
fit <- lm(
  Sepal.Length ~ poly(Petal.Width, degree = 2L) * Species + Petal.Length,
  data = iris
)
x <- c("Petal.Width", "Species", "Petal.Length")
preds <- unname(predict(fit, iris))
J <- c(1L, 51L, 101L)

shap <- list(
  # Exact
  kernelshap(fit, iris[x], bg_X = iris, verbose = FALSE),
  permshap(fit, iris[x], bg_X = iris, verbose = FALSE),
  # Sampling
  kernelshap(
    fit, iris[x],
    bg_X = iris, exact = FALSE, hybrid_degree = 0, verbose = FALSE
  ),
  permshap(fit, iris[x], bg_X = iris, exact = FALSE, verbose = FALSE)
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
  expect_true(all(shap[[3L]]$n_iter == 2L)) # ks stops after second iteraction
  expect_true(all(shap[[4L]]$n_iter == 1L)) # ps stops after first iteration
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

test_that("using foreach (non-parallel) gives the same as normal mode", {
  for (algo in c(kernelshap, permshap)) {
    s <- algo(fit, iris[J, x], bg_X = iris, verbose = FALSE)
    s2 <- suppressWarnings(
      algo(fit, iris[J, x], bg_X = iris, verbose = FALSE, parallel = TRUE)
    )
    expect_equal(s, s2)
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
    expect_equal(rowSums(s$S) + s$baseline, preds[J]) # sum shap = centered preds
    expect_no_error( # additional cols in bg are ok
      algo(fit, X[J, x], pred_fun = pred_fun, bg_X = cbind(d = 1, X), verbose = FALSE)
    )
    expect_error( # feature_names are less flexible
      algo(fit, X[J, ],
        pred_fun = pred_fun, bg_X = X,
        verbose = FALSE, feature_names = "Sepal.Width"
      )
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
  expect_equal(unname(rowSums(s$S)) + s$baseline, preds[J])

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
  expect_equal(s1$S, shap[[1L]]$S[J, ])
})

test_that("kernelshap works for large p (hybrid case)", {
  set.seed(9L)
  X <- data.frame(matrix(rnorm(20000L), ncol = 100L))
  y <- X[, 1L] * X[, 2L] * X[, 3L]
  fit <- lm(y ~ X1:X2:X3 + ., data = cbind(y = y, X))
  s <- kernelshap(fit, X[1L, ], bg_X = X, verbose = FALSE)

  expect_equal(s$baseline, mean(y))
  expect_equal(rowSums(s$S) + s$baseline, unname(predict(fit, X[1L, ])))
})
