# Helper functions
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

test_that("exact_Z() works for both kernel- and permshap", {
  p <- 5
  nm <- LETTERS[1:p]
  r1 <- exact_Z(p, feature_names = nm, keep_extremes = TRUE)
  r2 <- exact_Z(p, feature_names = nm, keep_extremes = FALSE)
  expect_equal(r1[2:(nrow(r1) - 1L), ], r2)
  expect_equal(colnames(r1), nm)
  expect_equal(dim(r1), c(2^p, p))
})

# Unit tests copied from {hstats}

test_that("rep_each() works", {
  expect_equal(rep_each(3, 10), rep_each(3L, 10L))
  expect_equal(rep_each(3, 10), rep(1:3, each = 10))
  expect_true(is.integer(rep_each(100, 100)))
})

test_that("fdummy() works", {
  x <- c("A", "A", "C", "D")
  mm <- matrix(model.matrix(~ x + 0), ncol = 3, dimnames = list(NULL, c("A", "C", "D")))
  expect_equal(fdummy(x), mm)
})

test_that("fdummy() works for singletons", {
  x <- c("A")
  expect_equal(fdummy(x), cbind(A = 1))
  expect_true(is.matrix(fdummy(x)))
})

test_that("fdummy() respects factor level order", {
  x1 <- factor(c("A", "A", "C", "D"))
  x2 <- factor(x1, levels = rev(levels(x1)))
  d1 <- fdummy(x1)
  d2 <- fdummy(x2)
  expect_equal(d1, d2[, colnames(d1)])
  expect_equal(colnames(d1), rev(colnames(d2)))
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

test_that("rowmean_factor() works for factor input", {
  x <- factor(c("C", "A", "C", "C", "A", "A"))
  out <- rowmean_factor(x, ngroups = 2L)
  
  expect_true(is.matrix(out))
  expect_equal(out, cbind(A = c(1/3, 2/3), C = c(2/3, 1/3)))
})
