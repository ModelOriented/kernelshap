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

test_that("regoranize_list() fails for non-list inputs", {
  expect_error(reorganize_list(alist = 1:10))
})

test_that("weighted_colMeans() gives the same as colMeans() in unweighted case", {
  X <- cbind(1:3, 2:4)
  expect_equal(c(weighted_colMeans(X)), colMeans(X))
  expect_equal(c(weighted_colMeans(X, w = c(1, 1, 1))), colMeans(X))
  expect_equal(c(weighted_colMeans(X, w = c(2, 2, 2))), colMeans(X))
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