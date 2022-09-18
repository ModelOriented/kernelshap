fit <- stats::lm(Sepal.Length ~ ., data = iris)
s <- kernelshap(fit, iris[1:2, -1], bg_X = iris[-1])

test_that("get_* functions get the right element", {
  expect_equal(s$S, ks_extract(s, "S"))
  expect_equal(s$X, ks_extract(s, "X"))
  expect_equal(s$SE, ks_extract(s, "SE"))
  expect_equal(s$baseline, ks_extract(s, "baseline"))
  expect_equal(s$n_iter, ks_extract(s, "n_iter"))
  expect_equal(s$converged, ks_extract(s, "converged"))
  expect_equal(s$m, ks_extract(s, "m"))
  expect_equal(s$m_exact, ks_extract(s, "m_exact"))
  expect_equal(s$prop_exact, ks_extract(s, "prop_exact"))
  expect_equal(s$txt, ks_extract(s, "txt"))
})

test_that("ks_extract() fails for wrong objects", {
  expect_error(ks_extract(1))
})

test_that("ks_extract() fails for wrong elements", {
  expect_error(ks_extract(s, "x"))
})

test_that("is_kernelshap() works", {
  expect_true(is.kernelshap(s))
  expect_false(is.kernelshap(1))
})
