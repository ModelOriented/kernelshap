# Model with non-linearities and interactions
y <- mtcars$mpg
Y <- as.matrix(mtcars[, c("mpg", "cyl")])

fit_y <- lm(y ~ poly(disp, degree = 2L) * vs + drat + wt, data = mtcars)
fit_Y <- lm(Y ~ poly(disp, degree = 2L) * vs + drat + wt, data = mtcars)

x <- c("disp", "vs", "drat", "wt")
J <- c(1L, 10L, 20L)

preds_y <- unname(predict(fit_y, mtcars[J, ]))
preds_Y <- unname(predict(fit_Y, mtcars[J, ]))

shap_y <- list(
  kernelshap(fit_y, mtcars[J, x], bg_X = mtcars, verbose = FALSE),
  kernelshap(
    fit_y, mtcars[J, x],
    bg_X = mtcars, exact = FALSE, hybrid_degree = 1L, verbose = FALSE
  ),
  permshap(fit_y, mtcars[J, x], bg_X = mtcars, verbose = FALSE),
  permshap(fit_y, mtcars[J, x], bg_X = mtcars, exact = FALSE, verbose = FALSE),
  permshap(
    fit_y, mtcars[J, x],
    bg_X = mtcars, exact = FALSE, low_memory = TRUE, verbose = FALSE
  )
)

shap_Y <- list(
  kernelshap(fit_Y, mtcars[J, x], bg_X = mtcars, verbose = FALSE),
  kernelshap(
    fit_Y, mtcars[J, x],
    bg_X = mtcars,
    exact = FALSE, hybrid_degree = 1L, verbose = FALSE
  ),
  permshap(fit_Y, mtcars[J, x], bg_X = mtcars, verbose = FALSE),
  permshap(fit_Y, mtcars[J, x], bg_X = mtcars, exact = FALSE, verbose = FALSE),
  permshap(
    fit_Y, mtcars[J, x],
    bg_X = mtcars, exact = FALSE, low_memory = TRUE, verbose = FALSE
  )
)

test_that("Exact and sampling modes agree with interactions of order 2", {
  expect_equal(shap_Y[[1L]]$S, shap_Y[[2L]]$S) # exact ks vs sampling ks
  expect_equal(shap_Y[[1L]]$S, shap_Y[[3L]]$S) # exact ks vs exact ps
  expect_equal(shap_Y[[3L]]$S, shap_Y[[4L]]$S) # exact ps vs sampling ps
  expect_equal(shap_Y[[4L]]$S, shap_Y[[5L]]$S) # low vs high memory ps
})

test_that("Baseline equals average prediction on background data", {
  for (i in 1:5) {
    expect_equal(shap_Y[[i]]$baseline, unname(colMeans(Y)))
  }
})

test_that("SHAP + baseline = prediction", {
  for (i in 1:5) {
    s <- shap_Y[[i]]
    expect_equal(rowSums(s$S[[1L]]) + s$baseline[1L], preds_Y[, 1L])
    expect_equal(rowSums(s$S[[2L]]) + s$baseline[2L], preds_Y[, 2L])
  }
})

test_that("First dimension of multioutput model equals single output", {
  for (i in 1:5) {
    expect_equal(shap_Y[[i]]$baseline[1L], shap_y[[i]]$baseline)
    expect_equal(shap_Y[[i]]$S[[1L]], shap_y[[i]]$S)
  }
})
