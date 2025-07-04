test_that("sum of kernel weights is 1", {
  for (p in 2:10) {
    expect_equal(sum(kernel_weights(p)), 1.0)
  }
})

test_that("Sum of kernel weights is 1, even for subset of domain", {
  expect_equal(sum(kernel_weights(10L, S = 2:5)), 1.0)
})

p <- 10L
m <- 100L

test_that("Random z have right output dim and the sums are between 1 and p-1", {
  Z <- sample_Z(p, m = m, feature_names = LETTERS[1:p])

  expect_equal(dim(Z), c(m, p))
  expect_true(all(rowSums(Z) %in% 1:(p - 1L)))
})

test_that("Random z have right output dim and the sums are in subset S", {
  S <- 2:3
  Z <- sample_Z(p, m = m, feature_names = LETTERS[1:p], S = S)

  expect_equal(dim(Z), c(m, p))
  expect_true(all(rowSums(Z) %in% S))
})

test_that("Sampling input structure is ok (deg = 0)", {
  input <- input_sampling(
    p,
    m = m, deg = 0L, feature_names = LETTERS[1:p]
  )

  expect_equal(dim(input$Z), c(m, p))
  expect_equal(sum(input$w), 1.0)
  expect_equal(dim(input$A), c(p, p))
  expect_equal(unname(diag(input$A)), rep(0.5, p))
})

test_that("Sampling input structure is ok (deg = 1)", {
  input <- input_sampling(
    p,
    m = m, deg = 1L, feature_names = LETTERS[1:p]
  )

  expect_equal(dim(input$Z), c(m, p))
  expect_true(sum(input$w) < 1.0)
  expect_equal(dim(input$A), c(p, p))
  expect_true(all(diag(input$A) < 0.5))
})

test_that("Sampling input input structure ok (deg = 2)", {
  input <- input_sampling(
    p,
    m = m, deg = 2L, feature_names = LETTERS[1:p]
  )

  expect_equal(dim(input$Z), c(m, p))
  expect_true(sum(input$w) < 1.0)
  expect_equal(dim(input$A), c(p, p))
  expect_true(all(diag(input$A) < 0.5))
})

test_that("Partly exact A, w, Z equal exact for sufficiently large deg", {
  for (p in 2:10) {
    pa <- input_partly_exact(p, deg = trunc(p / 2), feature_names = LETTERS[1:p])
    ex <- input_exact(p, feature_names = LETTERS[1:p])
    pa_rs <- rowSums(pa$Z)
    ex_rs <- rowSums(ex$Z)

    expect_equal(pa$A, ex$A)
    expect_equal(pa$w[order(pa_rs)], ex$w[order(ex_rs)])
    expect_equal(tabulate(pa_rs), tabulate(ex_rs))
  }
})

test_that("hybrid weights sum to 1 for different p and degree 1", {
  deg <- 1L
  expect_error(input_sampling(2L, deg = deg, feature_names = LETTERS[1:p]))
  expect_error(input_sampling(3L, deg = deg, feature_names = LETTERS[1:p]))

  for (p in 4:20) {
    pa <- input_partly_exact(p, deg = deg, feature_names = LETTERS[1:p])
    sa <- input_sampling(
      p,
      m = 100L, deg = deg, feature_names = LETTERS[1:p]
    )
    expect_equal(sum(pa$w) + sum(sa$w), 1.0)
  }
})

test_that("hybrid weights sum to 1 for different p and degree 2", {
  deg <- 2L
  expect_error(input_sampling(4L, deg = deg, feature_names = LETTERS[1:p]))
  expect_error(input_sampling(5L, deg = deg, feature_names = LETTERS[1:p]))

  for (p in 6:20) {
    pa <- input_partly_exact(p, deg = deg, feature_names = LETTERS[1:p])
    sa <- input_sampling(
      p,
      m = 100L, deg = deg, feature_names = LETTERS[1:p]
    )
    expect_equal(sum(pa$w) + sum(sa$w), 1L)
  }
})

test_that("partly_exact_Z(p, k) fails for bad p or k", {
  expect_error(partly_exact_Z(0L, k = 1L, feature_names = LETTERS[1:p]))
  expect_error(partly_exact_Z(5L, k = 3L, feature_names = LETTERS[1:p]))
  expect_error(partly_exact_Z(5L, k = 0L, feature_names = LETTERS[1:p]))
})

test_that("input_partly_exact(p, deg) fails for bad p or deg", {
  expect_error(input_partly_exact(2L, deg = 0L, feature_names = LETTERS[1:p]))
  expect_error(input_partly_exact(5L, deg = 3L, feature_names = LETTERS[1:p]))
})
