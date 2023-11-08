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
