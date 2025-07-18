test_that("check_convergence() works", {
  beta <- cbind(c(-1, 2)) # range 3
  sigma <- cbind(c(0.1, 0.3)) # max 0.3

  expect_true(check_convergence(beta, sigma, tol = 0.1))
  expect_false(check_convergence(beta, sigma, tol = 0.01))
})

test_that("get_sigma() works", {
  K <- 1L
  p <- 2L
  n_iter <- 3L
  est <- array(dim = c(n_iter, p, K), dimnames = list(NULL, c("x1", "x2"), "K"))
  est[1, , ] <- c(1, 1)
  est[2, , ] <- c(1, 2)
  est[3, , ] <- c(1, 3)
  # Standard error without Bessel correction in the variance
  expected <- cbind(K = c(x1 = 0, x2 = sd(1:3) * sqrt(2) / 3))
  expect_equal(get_sigma(est), expected)
  expect_error(get_sigma(est[1L, , , drop = FALSE]))
})

test_that("head_list(x) = head(x) for matrix x", {
  x <- cbind(1:10, 2:11)
  expect_equal(head_list(x), utils::head(x))
})

test_that("head_list(x)[[1L]] = head(x[[1L]]) for list of matries x", {
  x1 <- cbind(1:10, 2:11)
  x2 <- cbind(1:7, 2:8)
  x <- list(x1, x2)
  expect_equal(head_list(x)[[1L]], utils::head(x[[1L]]))
})

test_that("reorganize_list() fails for non-list inputs", {
  expect_error(reorganize_list(alist = 1:10))
})

test_that("wcolMeans() gives the same as colMeans() in unweighted case", {
  X <- cbind(1:3, 2:4)
  expect_equal(c(wcolMeans(X)), colMeans(X))
  expect_equal(c(wcolMeans(X, w = c(1, 1, 1))), colMeans(X))
  expect_equal(c(wcolMeans(X, w = c(2, 2, 2))), colMeans(X))
})

test_that("exact_Z() works", {
  for (p in 2:8) {
    nm <- LETTERS[1:p]
    r <- exact_Z(p, feature_names = nm)
    row_sum <- 0
    for (j in p:1) {
      row_sum <- row_sum + r[, j] * 2^(p - j)
    }
    expect_equal(colnames(r), nm)
    expect_equal(dim(r), c(2^p, p))
    expect_equal(row_sum, 0:(2^p - 1))
  }
  expect_error(exact_Z(1, "1"))
})

test_that("rep_rows() gives the same as usual subsetting (except rownames)", {
  setrn <- function(x) {
    rownames(x) <- 1:nrow(x)
    x
  }

  expect_equal(rep_rows(iris, 1), iris[1, ])
  expect_equal(rep_rows(iris, 2:1), setrn(iris[2:1, ]))
  expect_equal(rep_rows(iris, c(1, 1, 1)), setrn(iris[c(1, 1, 1), ]))

  ir <- iris[1, ]
  ir$y <- list(list(a = 1, b = 2))
  expect_equal(rep_rows(ir, c(1, 1)), setrn(ir[c(1, 1), ]))
})

test_that("rep_rows() gives the same as usual subsetting for matrices", {
  ir <- data.matrix(iris[1:4])

  expect_equal(rep_rows(ir, c(1, 1, 2)), ir[c(1, 1, 2), ])
  expect_equal(rep_rows(ir, 1), ir[1, , drop = FALSE])
})


# Unit tests copied from {hstats}

test_that("rep_each() works", {
  expect_equal(rep_each(3, 10), rep_each(3L, 10L))
  expect_equal(rep_each(3, 10), rep(1:3, each = 10))
  expect_true(is.integer(rep_each(100, 100)))
})

test_that("wrowmean_vector() works for 1D matrices", {
  x2 <- cbind(a = 6:1)
  out2 <- wrowmean_vector(x2, ngroups = 2L)
  expec <- rowsum(x2, group = rep(1:2, each = 3)) / 3
  rownames(expec) <- NULL

  expect_error(wrowmean_vector(matrix(1:4, ncol = 2L)))
  expect_equal(out2, expec)

  expect_equal(wrowmean_vector(x2, ngroups = 3L), cbind(a = c(5.5, 3.5, 1.5)))

  # Constant weights have no effect
  expect_equal(wrowmean_vector(x2, ngroups = 2L, w = c(1, 1, 1)), out2)
  expect_equal(wrowmean_vector(x2, ngroups = 2L, w = c(4, 4, 4)), out2)

  # Non-constant weights
  a <- weighted.mean(6:4, 1:3)
  b <- weighted.mean(3:1, 1:3)
  out2 <- wrowmean_vector(x2, ngroups = 2L, w = 1:3)
  expect_equal(out2, cbind(a = c(a, b)))
})

test_that("align_pred() works", {
  expect_error(align_pred(c("A", "B")))
  expect_error(align_pred(factor(c("A", "B"))))
  expect_equal(align_pred(1:4), as.matrix(1:4))
})
