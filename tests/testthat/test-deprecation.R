test_that("Deprecation warnings work", {
  fit <- lm(Sepal.Length ~ ., data = iris)

  for (algo in c(kernelshap, permshap)) {
    expect_warning(
      algo(fit, X = iris[1:2, -1], bg_X = iris, parallel = TRUE, verbose = FALSE)
    )
    expect_warning(
      algo(
        fit,
        X = iris[1:2, -1],
        bg_X = iris,
        parallel_args = list(a = 1),
        verbose = FALSE
      )
    )
  }
})
